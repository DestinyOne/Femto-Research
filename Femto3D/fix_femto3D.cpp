/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
    Contributing authors: Paul Crozier (Original fix ttm author)
                          Carolyn Phillips (Original fix ttm author)
                          Sergey Starikov (fix ttm mod author)
                          Vasily Pisarev (fix ttm mod author)
                          Weirong Yuan (fix femto3D author, Center for Materials Under eXtreme Environments(CMUXE))
------------------------------------------------------------------------- */

#include "potential_file_reader.h"

#include "lmptype.h"
#include <mpi.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "fix_femto3D.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "comm.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"
#include "citeme.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define MAXLINE 1024

static const char cite_fix_femto3D[] =
"fix femto3D command:\n\n"
"@article{Pisarev2014,\n"
"author = {Pisarev, V. V. and Starikov, S. V.},\n"
"title = {{Atomistic simulation of ion track formation in UO2.}},\n"
"journal = {J.~Phys.:~Condens.~Matter},\n"
"volume = {26},\n"
"number = {47},\n"
"pages = {475401},\n"
"year = {2014}\n"
"}\n\n";

/* ---------------------------------------------------------------------- */

FixFEMTO3D::FixFEMTO3D(LAMMPS* lmp, int narg, char** arg) :
  Fix(lmp, narg, arg)
{
  // MPI_Comm world = MPI_COMM_WORLD;
  MPI_Comm_rank(world, &pid);
  MPI_Comm_size(world, &numP);
  if (lmp->citeme) lmp->citeme->add(cite_fix_femto3D);

  if (narg != 7 && narg != 8) error->all(FLERR, "Illegal fix femto3D command");

  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;
  extvector = 1;
  nevery = 1;
  restart_peratom = 1;
  restart_global = 1;

  seed = utils::inumeric(FLERR, arg[3], false, lmp);
  if (seed <= 0) error->all(FLERR, "Invalid random number seed in fix femto3D command");

  // Output every this many timesteps, 0 = no dump
  nfileevery = utils::inumeric(FLERR, arg[4], false, lmp);
  if (nfileevery < 0) error->all(FLERR, "Invalid output parameter in fix femto3D command");

  // fp_parameter
  read_parameter(arg[5]);

  // fp_tablelist
  read_tablelist(arg[6]);

  if (premode == 0) {
    // fp_outlist
    fp_outlist = arg[7];
  }

  // t_surface is determined by electronic temperature (not constant)
  // t_surface_l = nxnodes; //surface_l;
  // t_surface_r = 0; //surface_r;
  duration = 0.0;
  // if (gamma_s < 0.0) error->all(FLERR,"Fix femto3D gamma_s must be >= 0.0");
  // if (E_stop < 0.0) error->all(FLERR,"Fix femto3D E_stop must be >= 0.0");
  if (nxnodes <= 0 || nynodes <= 0 || nznodes <= 0) error->all(FLERR, "Fix femto3D number of nodes must be > 0");
  if (ionic_density <= 0.0) error->all(FLERR, "Fix femto3D ionic_density must be > 0.0");
  if (surface_l < 0) error->all(FLERR, "Surface coordinates must be >= 0");
  if (surface_l >= surface_r) error->all(FLERR, "Left surface coordinate must be less than right surface coordinate");
  if ((bulk_ttm != 0) && (bulk_ttm != 1)) error->all(FLERR, "bulk_ttm should be 0 or 1");
  if ((premode != 0) && (premode != 1)) error->all(FLERR, "premode should be 0 or 1");
  // initialize Marsaglia RNG with processor-unique seed
  random = new RanMars(lmp, seed + comm->me);
  // allocate 3d grid variables
  total_nnodes = nxnodes * nynodes * nznodes;
  memory->create(nsum, nxnodes, nynodes, nznodes, "femto3D:nsum");
  memory->create(sum_mass_vsq, nxnodes, nynodes, nznodes, "femto3D:sum_mass_vsq");
  memory->create(T_a, nxnodes, nynodes, nznodes, "femto3D:T_a");
  if (premode == 0) {
    memory->create(Activated, nxnodes, nynodes, nznodes, "femto3D:Activated");

    memory->create(sum_mass_v, nxnodes, nynodes, nznodes, "femto3D:sum_mass_v");
    memory->create(sum_mass, nxnodes, nynodes, nznodes, "femto3D:sum_mass");
    memory->create(average_v, nxnodes, nynodes, nznodes, "femto3D:average_v");
    memory->create(T_electron_old, nxnodes, nynodes, nznodes, "femto3D:T_electron_old");
    memory->create(T_electron_first, nxnodes, nynodes, nznodes, "femto3D:T_electron_first");
    memory->create(T_electron, nxnodes, nynodes, nznodes, "femto3D:T_electron");

    memory->create(net_energy_transfer, nxnodes, nynodes, nznodes, "femto3D:net_energy_transfer");
    memory->create(energy_conduction, nxnodes, nynodes, nznodes, "femto3D:energy_conduction");
    memory->create(mult_factor, nxnodes, nynodes, nznodes, "femto3D:mult_factor");
    memory->create(skin_layer, nxnodes, nynodes, nznodes, "femto3D:skin_layer");
    memory->create(ke_real, nxnodes, nynodes, nznodes, "femto3D:ke_real");
    memory->create(CeT, nxnodes, nynodes, nznodes, "femto3D:CeT");
    memory->create(GT, nxnodes, nynodes, nznodes, "femto3D:GT");

    memory->create(Genergy, nxnodes, "femto3D:Genergy");
    memory->create(Kenergy, nxnodes, "femto3D:Kenergy");

    if (bulk_ttm == 1) {
      memory->create(Ta_bulk, bxsize, nynodes, nznodes, "femto3D:Ta_bulk");
      memory->create(Te_bulk, bxsize, nynodes, nznodes, "femto3D:Te_bulk");
      memory->create(CeT_bulk, bxsize, nynodes, nznodes, "femto3D:CeT_bulk");
      memory->create(ke_real_bulk, bxsize, nynodes, nznodes, "femto3D:ke_real_bulk");
      memory->create(E_melt_buffer, bxsize, nynodes, nznodes, "femto3D:E_melt_buffer");
      memory->create(GT_bulk, bxsize, nynodes, nznodes, "femto3D:GT_bulk");
      memory->create(ki_bulk, bxsize, nynodes, nznodes, "femto3D:ki_bulk");
      memory->create(x_max, nynodes, nznodes, "femto3D:x_max");
      memory->create(Tes, nynodes, nznodes, "femto3D:Tes");
      memory->create(mult_factor_bulk, bxsize, nynodes, nznodes, "femto3D:mult_factor_bulk");
      memory->create(skin_layer_bulk, bxsize, nynodes, nznodes, "femto3D:skin_layer_bulk");
    }
  }

  if (premode == 1) {
    if (bulk_ttm == 1) {
      memory->create(x_max, nynodes, nznodes, "femto3D:x_max");
    }
  }

  grow_arrays(atom->nmax); // allocation/reallocation of Langevin force array

  for (int i = 0; i < atom->nmax; i++) {
    flangevin[i][0] = 0.0;
    flangevin[i][1] = 0.0;
    flangevin[i][2] = 0.0;
  }

  atom->add_callback(Atom::GROW);
  atom->add_callback(Atom::RESTART);

  // epsilon0 = 8.854187817e-21; // C^2*ps^2/gram/Angstrom^3
  epsilon0 = 1.418597153e-40; // C^2/eV/Angstrom
  speed_light = 2.99792458e+6; // Angstroms/picosecond
  Pi = MY_PI;
  e_charge = 1.60217662e-19; // Electron charge in coulombs (C)
  // m_e = 9.10938e-31; // rest mass of an electron in kg
  m_e = 9.10938e-31 / e_charge * 1.0e+4; // rest mass of an electron in eV*ps^2/A^2
  k_b = force->boltz; // Boltzmann constant in eV/K
  // h_bar = 1.0545718e-34; // J*s Plank's constant
  h_bar = force->hplanck / 2.0 / Pi;
  Na = 6.02214086e+23; // 1/mol

  crit_num_f = 0.1;
  dx = domain->xprd / nxnodes;
  dy = domain->yprd / nynodes;
  dz = domain->zprd / nznodes;
  del_vol = dx * dy * dz;
  if (bulk_ttm == 1) {
    Latent_melt *= ionic_density / Na; // in eV/A^3
    dx_bulk = bulk_thick / bxsize;
  }

  if (premode == 0) {
    // set initial electron temperatures from user input file
    set_initial_temperatures();
  }

}

/* ---------------------------------------------------------------------- */

FixFEMTO3D::~FixFEMTO3D()
{
  if (premode == 0) {

    delete random;
    memory->destroy(nsum);
    memory->destroy(Activated);
    memory->destroy(sum_mass_vsq);
    memory->destroy(sum_mass_v);
    memory->destroy(sum_mass);
    memory->destroy(average_v);
    memory->destroy(T_electron_first);
    memory->destroy(T_electron_old);
    memory->destroy(T_electron);
    if (bulk_ttm == 1) {
      memory->destroy(Te_bulk);
      memory->destroy(Ta_bulk);
      memory->destroy(CeT_bulk);
      memory->destroy(ke_real_bulk);
      memory->destroy(GT_bulk);
      memory->destroy(ki_bulk);
      memory->destroy(x_max);
      memory->destroy(Tes);
      memory->destroy(mult_factor_bulk);
      memory->destroy(skin_layer_bulk);
      memory->destroy(E_melt_buffer);
    }
    memory->destroy(T_a);
    memory->destroy(flangevin);
    memory->destroy(net_energy_transfer);
    memory->destroy(energy_conduction);
    memory->destroy(ZTe);
    memory->destroy(CeTe);
    memory->destroy(GTe);
    memory->destroy(KeTe);
    memory->destroy(ReflecTe);
    memory->destroy(PenTe);
    memory->destroy(ke_real);
    memory->destroy(CeT);
    memory->destroy(GT);
    memory->destroy(Genergy);
    memory->destroy(Kenergy);
    memory->destroy(mult_factor);
    memory->destroy(skin_layer);
    atom->delete_callback(id, 0);
    atom->delete_callback(id, 1);
  }

  if (premode == 1) {

    delete random;
    memory->destroy(nsum);
    memory->destroy(sum_mass_vsq);
    memory->destroy(ITemp);
    memory->destroy(ITempBulk);
    memory->destroy(T_a);
    if (bulk_ttm == 1) {
      memory->destroy(x_max);
    }
    memory->destroy(flangevin);
    atom->delete_callback(id, 0);
    atom->delete_callback(id, 1);
  }

}

/* ---------------------------------------------------------------------- */

int FixFEMTO3D::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixFEMTO3D::init()
{
  if (strcmp(update->unit_style, "metal") != 0)
    error->all(FLERR, "fix femto3D should use unit style 'metal'");
  if (domain->dimension == 2)
    error->all(FLERR, "Cannot use fix femto3D with 2d simulation");
  if (domain->boundary[0][0] != 1 || domain->boundary[0][1] != 1 ||
    domain->boundary[1][0] != 0 || domain->boundary[1][1] != 0 ||
    domain->boundary[2][0] != 0 || domain->boundary[2][1] != 0)
    error->all(FLERR, "Use wrong boundaries with fix femto3D, should be f p p");
  if (domain->triclinic)
    error->all(FLERR, "Cannot use fix femto3D with triclinic box");
  if (utils::strmatch(update->integrate_style,"^respa"))
    nlevels_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels;

  if (premode == 0) {
    read_outlist(fp_outlist, static_cast<int>(1000.0 * duration));
    update_Ta();
    // Initialize parameters with the initial temperature
    update_parameters();
    if (bulk_ttm == 1) {
      update_parameters_bulk();
    }
    for (int ixnode = 0; ixnode < nxnodes; ixnode++) {
      Genergy[ixnode] = 0.0;
      Kenergy[ixnode] = 0.0;
    }

    hasrun = true;
    if (NewPulse) duration = 0.0;
  }

}

/* ---------------------------------------------------------------------- */

void FixFEMTO3D::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style,"^verlet")) {
    post_force_setup(vflag);
  } else {
    (dynamic_cast<Respa *>(update->integrate))->copy_flevel_f(nlevels_respa-1);
    post_force_respa_setup(vflag,nlevels_respa-1,0);
    (dynamic_cast<Respa *>(update->integrate))->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixFEMTO3D::post_force(int vflag)
{
  // int pid;
  // MPI_Comm_rank(world,&pid);
  //for (int ixnode = 0; ixnode < bxsize; ixnode++)
    //for (int iynode = 0; iynode < nynodes; iynode++)
      //for (int iznode = 0; iznode < nznodes; iznode++)
      //if (pid == 0 && update->ntimestep <= 203) {
        //printf("Te is %f in x = %d \n", Te_bulk[ixnode][iynode][iznode], ixnode);
      //}
        // if (pid == 0 && Activated[ixnode][iynode][iznode] == 1 && update->ntimestep <= 203) {
         // printf("Te is %f in x = %d \n", T_electron[ixnode][iynode][iznode], ixnode);
        // }
  double** x = atom->x;
  double** v = atom->v;
  double** f = atom->f;
  double* mass = atom->mass;
  double* rmass = atom->rmass;
  int* type = atom->type;
  int* mask = atom->mask;
  int nlocal = atom->nlocal;
  double gamma1, gamma2;

  if (premode == 0) {
    double xmax = 0.0;
    xmin = domain->boxhi[0];

    int nmark;
    double bound_thick_temp;
    double v_b = 0.0;

    if (bulk_ttm == 1) {
      for (int iynode = 0; iynode < nynodes; iynode++)
        for (int iznode = 0; iznode < nznodes; iznode++) {
          x_max[iynode][iznode] = domain->boxlo[0];
          Tes[iynode][iznode] = 0;
        }
    }

    update_Ta();
    ChangeType();

    // Calculate the rear boundary velocity for TTM boundary force
    if (bulk_ttm == 1) {
      for (int iynode = 0; iynode < nynodes; iynode++)
        for (int iznode = 0; iznode < nznodes; iznode++) {
          xmax += x_max[iynode][iznode];
        }
      xmax /= nynodes * nznodes; // average rear surface position

      double crit = ionic_density * domain->yprd * domain->zprd * bound_thick;// have to account for initial thermal expansion
      bound_thick_temp = bound_thick;
      do {
        nmark = 0;
        v_b = 0;
        for (int i = 0; i < nlocal; i++) {
          if (mask[i] & groupbit) {
            if (rmass) massone = rmass[i];
            else massone = mass[type[i]];
            double xscale = (x[i][0] - domain->boxlo[0]) / domain->xprd;
            double yscale = (x[i][1] - domain->boxlo[1]) / domain->yprd;
            double zscale = (x[i][2] - domain->boxlo[2]) / domain->zprd;
            int ixnode = static_cast<int>(xscale * nxnodes);
            int iynode = static_cast<int>(yscale * nynodes);
            int iznode = static_cast<int>(zscale * nznodes);
            // Eliminate numerical errors
            if (ixnode > nxnodes - 1) ixnode = nxnodes - 1;
            if (iynode > nynodes - 1) iynode = nynodes - 1;
            if (iznode > nznodes - 1) iznode = nznodes - 1;
            if (ixnode < 0) ixnode = 0;
            if (iynode < 0) iynode = 0;
            if (iznode < 0) iznode = 0;

            if (x[i][0] > (xmax - bound_thick_temp)) {
              v_b += v[i][0];
              nmark++;
            }
          }
        }
        MPI_Allreduce(MPI_IN_PLACE, &nmark, 1, MPI_INT, MPI_SUM, world);
        if (nmark < crit) bound_thick_temp *= 1.01;
      } while (nmark < crit);

      MPI_Allreduce(MPI_IN_PLACE, &v_b, 1, MPI_DOUBLE, MPI_SUM, world);
      v_b /= nmark;
    }

    update_parameters();

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        // printf("Atom ID is %d, local ID is %d, core ID is %d\n", atom->tag[i], i, pid);
        if (rmass) massone = rmass[i];
        else massone = mass[type[i]];
        double xscale = (x[i][0] - domain->boxlo[0]) / domain->xprd;
        double yscale = (x[i][1] - domain->boxlo[1]) / domain->yprd;
        double zscale = (x[i][2] - domain->boxlo[2]) / domain->zprd;
        int ixnode = static_cast<int>(xscale * nxnodes);
        int iynode = static_cast<int>(yscale * nynodes);
        int iznode = static_cast<int>(zscale * nznodes);
        // Eliminate numerical errors
        if (ixnode > nxnodes - 1) ixnode = nxnodes - 1;
        if (iynode > nynodes - 1) iynode = nynodes - 1;
        if (iznode > nznodes - 1) iznode = nznodes - 1;
        if (ixnode < 0) ixnode = 0;
        if (iynode < 0) iynode = 0;
        if (iznode < 0) iznode = 0;
        if (T_electron[ixnode][iynode][iznode] < 0)
          error->all(FLERR, "Electronic temperature dropped below zero");
        if (Activated[ixnode][iynode][iznode] == 0) {
          flangevin[i][0] = 0.0;
          flangevin[i][1] = 0.0;
          flangevin[i][2] = 0.0;
        }
        else {
          double G = GT[ixnode][iynode][iznode];
          if (T_electron[ixnode][iynode][iznode] < T_a[ixnode][iynode][iznode]) {
            double ci = 3.0 * nsum[ixnode][iynode][iznode] / del_vol * force->boltz;
            double upper = ci * CeT[ixnode][iynode][iznode] / (ci + CeT[ixnode][iynode][iznode]) / update->dt; // upper limit
            if (G > upper) {
              G = upper;
            }
          }

          double gamma_p = G * massone / 3.0 / force->boltz / ionic_density;
          gamma1 = -gamma_p / force->ftm2v;
          gamma2 = sqrt(2.0 * force->boltz * gamma_p / update->dt / force->mvv2e) / force->ftm2v * sqrt(T_electron[ixnode][iynode][iznode]);

          flangevin[i][0] = gamma1 * (v[i][0] - average_v[ixnode][iynode][iznode]) + gamma2 * random->gaussian(); // random->gaussian(); //
          flangevin[i][1] = gamma1 * v[i][1] + gamma2 * random->gaussian(); //sqrt(12.0)*(random->uniform()-0.5);
          flangevin[i][2] = gamma1 * v[i][2] + gamma2 * random->gaussian(); //sqrt(12.0)*(random->uniform()-0.5);

          double x_surf = dx * double(t_surface_l) + domain->boxlo[0];
          if (x_surf < xmin) x_surf = xmin;
          double x_at = x[i][0];
          int right_xnode = ixnode + 1;
          int right_ynode = iynode + 1;
          int right_znode = iznode + 1;
          if (right_xnode == nxnodes) right_xnode = nxnodes - 1;
          if (right_ynode == nynodes) right_ynode = 0;
          if (right_znode == nznodes) right_znode = 0;
          int left_xnode = ixnode - 1;
          int left_ynode = iynode - 1;
          int left_znode = iznode - 1;
          if (left_xnode == -1) left_xnode = 0;
          if (left_ynode == -1) left_ynode = nynodes - 1;
          if (left_znode == -1) left_znode = nznodes - 1;
          double T_i = T_electron[ixnode][iynode][iznode];
          double T_ir = T_electron[right_xnode][iynode][iznode];
          // if(Tes[iynode][iznode])
          double T_iu = T_electron[ixnode][right_ynode][iznode];
          double T_if = T_electron[ixnode][iynode][right_znode];
          double T_il = T_electron[left_xnode][iynode][iznode];
          double T_id = T_electron[ixnode][left_ynode][iznode];
          double T_ib = T_electron[ixnode][iynode][left_znode];
          double Ni = (double)nsum[ixnode][iynode][iznode] / del_vol;
          double C_i = CeT[ixnode][iynode][iznode] / Ni * ionic_density;
          double C_ir = CeT[right_xnode][iynode][iznode] / Ni * ionic_density;
          double C_iu = CeT[ixnode][right_ynode][iznode] / Ni * ionic_density;
          double C_if = CeT[ixnode][iynode][right_znode] / Ni * ionic_density;
          double C_il = CeT[left_xnode][iynode][iznode] / Ni * ionic_density;
          double C_id = CeT[ixnode][left_ynode][iznode] / Ni * ionic_density;
          double C_ib = CeT[ixnode][iynode][left_znode] / Ni * ionic_density;
          double factor0, factor1, factor2;

          // Elliminate the error when Te = 0
          if (Activated[right_xnode][iynode][iznode] == 0 && Activated[left_xnode][iynode][iznode] == 0) {
            factor0 = 0.0;
          }
          else if (Activated[right_xnode][iynode][iznode] == 1 && Activated[left_xnode][iynode][iznode] == 0) {
            factor0 = (C_ir * T_ir - C_i * T_i) / dx;
          }
          else if (Activated[right_xnode][iynode][iznode] == 0 && Activated[left_xnode][iynode][iznode] == 1) {
            factor0 = (C_i * T_i - C_il * T_il) / dx;
          }
          else {
            factor0 = (C_ir * T_ir - C_il * T_il) / 2.0 / dx;
          }

          if (factor0 > 0.0) factor0 = 0.0;

          if (Activated[ixnode][right_ynode][iznode] == 0 && Activated[ixnode][left_ynode][iznode] == 0) {
            factor1 = 0.0;
          }
          else if (Activated[ixnode][right_ynode][iznode] == 1 && Activated[ixnode][left_ynode][iznode] == 0) {
            factor1 = (C_iu * T_iu - C_i * T_i) / dy;
          }
          else if (Activated[ixnode][right_ynode][iznode] == 0 && Activated[ixnode][left_ynode][iznode] == 1) {
            factor1 = (C_i * T_i - C_id * T_id) / dy;
          }
          else {
            factor1 = (C_iu * T_iu - C_id * T_id) / 2.0 / dy;
          }

          if (Activated[ixnode][iynode][right_znode] == 0 && Activated[ixnode][iynode][left_znode] == 0) {
            factor2 = 0.0;
          }
          else if (Activated[ixnode][iynode][right_znode] == 1 && Activated[ixnode][iynode][left_znode] == 0) {
            factor2 = (C_if * T_if - C_i * T_i) / dz;
          }
          else if (Activated[ixnode][iynode][right_znode] == 0 && Activated[ixnode][iynode][left_znode] == 1) {
            factor2 = (C_i * T_i - C_ib * T_ib) / dz;
          }
          else {
            factor2 = (C_if * T_if - C_ib * T_ib) / 2.0 / dz;
          }

          // ------------ Only for Gold (Mean Free Path Mod)------------------
          // double FE = 5.5115;
          // free_path = sqrt(2.0*sqrt(force->boltz*force->boltz*T_i*T_i+FE*FE)*e_charge/9.10938e-31)/0.129/1.e5;
          // printf("MFP is %f \n", free_path);
          // -----------------------------------------------------------------

          // if (update->ntimestep >= 1000 && ixnode <= t_surface_l + 1 ) factor0=0;

          double diff_x = (x_at - x_surf);
          if (x_at >= x_surf) {
            // free_path = 20.0; // Starikov 2014
            // if (diff_x < free_path) diff_x = free_path;
            // flangevin[i][0] -= pres_factor/ionic_density*(C_i*T_i*free_path/(diff_x+free_path)/(diff_x+free_path) + diff_x/(diff_x+free_path)*factor0);
            flangevin[i][0] -= pres_factor / ionic_density * factor0;
            flangevin[i][1] -= pres_factor / ionic_density * factor1;
            flangevin[i][2] -= pres_factor / ionic_density * factor2;
          }
          else {
            flangevin[i][0] -= 0;
            flangevin[i][1] -= 0;
            flangevin[i][2] -= 0;
          }

          // if (update->ntimestep >= 1000) printf("650 in core %d at step %d\n",pid,update->ntimestep);

          double forcess = sqrt(flangevin[i][0] * flangevin[i][0] + flangevin[i][1] * flangevin[i][1] + flangevin[i][2] * flangevin[i][2]);
          if (forcess > 100) printf("The force of atom %d is %f in x = %d; Te = %f %f %f; x_at, x_surf = %f, %f\n", atom->tag[i], forcess, ixnode, T_i, T_ir, T_il, x_at, x_surf);
          // if (abs(factor0) > 0.1) printf("The factor0 of atom %d is %f in x = %d; Te = %f %f %f; x_at, x_surf = %f, %f\n", atom->tag[i], factor0, ixnode, T_i, T_ir, T_il, x_at, x_surf);
        }
        /*
        if (T_electron[ixnode][iynode][iznode] < T_a[ixnode][iynode][iznode]) {
          flangevin[i][0] = 0;
          flangevin[i][1] = 0;
          flangevin[i][2] = 0;
        }
        */

        f[i][0] += flangevin[i][0];
        f[i][1] += flangevin[i][1];
        f[i][2] += flangevin[i][2];

        // pressure from the bulk
        if (bulk_ttm == 1 && x[i][0] > (xmax - bound_thick_temp)) {
          f[i][0] += F_0 * domain->yprd * domain->zprd / nmark;
          f[i][0] -= massone * ionic_density * A_cross * v_s * v_b / force->ftm2v;
          // f[i][0] -= massone*ionic_density*A_cross*v_s*v[i][0]/force->ftm2v;
        }

      }
    }
  }

  if (premode == 1) {

    double v_b = 0.0;
    double xmax = 0.0;
    int nmark;
    double bound_thick_temp;

    // Calculate the rear boundary velocity for TTM boundary force
    if (bulk_ttm == 1) {

      for (int iynode = 0; iynode < nynodes; iynode++)
        for (int iznode = 0; iznode < nznodes; iznode++)
          x_max[iynode][iznode] = domain->boxlo[0];

      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
          if (rmass) massone = rmass[i];
          else massone = mass[type[i]];
          double xscale = (x[i][0] - domain->boxlo[0]) / domain->xprd;
          double yscale = (x[i][1] - domain->boxlo[1]) / domain->yprd;
          double zscale = (x[i][2] - domain->boxlo[2]) / domain->zprd;
          int ixnode = static_cast<int>(xscale * nxnodes);
          int iynode = static_cast<int>(yscale * nynodes);
          int iznode = static_cast<int>(zscale * nznodes);
          // Eliminate numerical errors
          if (ixnode > nxnodes - 1) ixnode = nxnodes - 1;
          if (iynode > nynodes - 1) iynode = nynodes - 1;
          if (iznode > nznodes - 1) iznode = nznodes - 1;
          if (ixnode < 0) ixnode = 0;
          if (iynode < 0) iynode = 0;
          if (iznode < 0) iznode = 0;
          if (x_max[iynode][iznode] < x[i][0]) x_max[iynode][iznode] = x[i][0];
        }
      }

      MPI_Allreduce(MPI_IN_PLACE, &x_max[0][0], nynodes * nznodes, MPI_DOUBLE, MPI_MAX, world);

      for (int iynode = 0; iynode < nynodes; iynode++)
        for (int iznode = 0; iznode < nznodes; iznode++) {
          xmax += x_max[iynode][iznode];
        }
      xmax /= nynodes * nznodes; // average rear surface position

      double crit = ionic_density * domain->yprd * domain->zprd * bound_thick;// have to account for initial thermal expansion
      bound_thick_temp = bound_thick;

      do {
        nmark = 0;
        v_b = 0.0;
        for (int i = 0; i < nlocal; i++) {
          if (mask[i] & groupbit) {
            if (rmass) massone = rmass[i];
            else massone = mass[type[i]];
            double xscale = (x[i][0] - domain->boxlo[0]) / domain->xprd;
            double yscale = (x[i][1] - domain->boxlo[1]) / domain->yprd;
            double zscale = (x[i][2] - domain->boxlo[2]) / domain->zprd;
            int ixnode = static_cast<int>(xscale * nxnodes);
            int iynode = static_cast<int>(yscale * nynodes);
            int iznode = static_cast<int>(zscale * nznodes);
            // Eliminate numerical errors
            if (ixnode > nxnodes - 1) ixnode = nxnodes - 1;
            if (iynode > nynodes - 1) iynode = nynodes - 1;
            if (iznode > nznodes - 1) iznode = nznodes - 1;
            if (ixnode < 0) ixnode = 0;
            if (iynode < 0) iynode = 0;
            if (iznode < 0) iznode = 0;

            if (x[i][0] > (xmax - bound_thick_temp)) {
              v_b += v[i][0];
              nmark++;
            }
          }
        }
        MPI_Allreduce(MPI_IN_PLACE, &nmark, 1, MPI_INT, MPI_SUM, world);
        if (nmark < crit) bound_thick_temp *= 1.01;
      } while (nmark < crit);

      MPI_Allreduce(MPI_IN_PLACE, &v_b, 1, MPI_DOUBLE, MPI_SUM, world);
      v_b /= nmark;
    }


    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        // printf("Atom ID is %d, local ID is %d, core ID is %d\n", atom->tag[i], i, pid);
        if (rmass) massone = rmass[i];
        else massone = mass[type[i]];
        double xscale = (x[i][0] - domain->boxlo[0]) / domain->xprd;
        double yscale = (x[i][1] - domain->boxlo[1]) / domain->yprd;
        double zscale = (x[i][2] - domain->boxlo[2]) / domain->zprd;
        int ixnode = static_cast<int>(xscale * nxnodes);
        int iynode = static_cast<int>(yscale * nynodes);
        int iznode = static_cast<int>(zscale * nznodes);
        // Eliminate numerical errors
        if (ixnode > nxnodes - 1) ixnode = nxnodes - 1;
        if (iynode > nynodes - 1) iynode = nynodes - 1;
        if (iznode > nznodes - 1) iznode = nznodes - 1;
        if (ixnode < 0) ixnode = 0;
        if (iynode < 0) iynode = 0;
        if (iznode < 0) iznode = 0;
        int TempN;
        if (ixnode < surface_l)
          TempN = 0;
        else if (ixnode >= surface_r)
          TempN = surface_r - surface_l - 1;
        else
          TempN = ixnode - surface_l;

        double G = 1.e-7 * GStrength;
        double gamma_p = G * massone / 3.0 / force->boltz / ionic_density;
        gamma1 = -gamma_p / force->ftm2v;
        gamma2 = sqrt(2.0 * force->boltz * gamma_p / update->dt / force->mvv2e) / force->ftm2v * sqrt(ITemp[TempN]);

        flangevin[i][0] = gamma1 * v[i][0] + gamma2 * random->gaussian();
        flangevin[i][1] = gamma1 * v[i][1] + gamma2 * random->gaussian();
        flangevin[i][2] = gamma1 * v[i][2] + gamma2 * random->gaussian();

        f[i][0] += flangevin[i][0];
        f[i][1] += flangevin[i][1];
        f[i][2] += flangevin[i][2];


        // pressure from the bulk
        if (bulk_ttm == 1 && x[i][0] > (xmax - bound_thick_temp)) {
          f[i][0] += F_0 * domain->yprd * domain->zprd / nmark;
          f[i][0] -= massone * ionic_density * A_cross * v_s * v_b / force->ftm2v;
        }

      }
    }

  }

}

/* ---------------------------------------------------------------------- */
void FixFEMTO3D::end_of_step()
{
  //int numP, pid;
  //MPI_Comm_size(world, &numP);
  //MPI_Comm_rank(world, &pid);
  double** x = atom->x;
  double** v = atom->v;
  double* mass = atom->mass;
  double* rmass = atom->rmass;
  int* type = atom->type;
  int* mask = atom->mask;
  int nlocal = atom->nlocal;

  if (premode == 0) {
    update_Ta();
    /*
      for (int ixnode = 0; ixnode < nxnodes; ixnode++)
        for (int iynode = 0; iynode < nynodes; iynode++)
          for (int iznode = 0; iznode < nznodes; iznode++) {
            if (Activated[ixnode][iynode][iznode] == 1) {
              T_electron_first[ixnode][iynode][iznode] = T_electron[ixnode][iynode][iznode] - update->dt/CeT[ixnode][iynode][iznode]* net_energy_transfer[ixnode][iynode][iznode]/del_vol;
            } else {
              T_a[ixnode][iynode][iznode] = 0.0;
              T_electron_first[ixnode][iynode][iznode] = 0.0;
            }
          }
    */
    // Tempout();

    double el_specific_heat = 1.0; // MAX
    double sh_min = MyCe(T_init, ionic_density, T_init);
    double el_ke = 0.0; // MIN
    update_parameters(); // Update for Ti changes

    for (int ixnode = 0; ixnode < nxnodes; ixnode++)
      for (int iynode = 0; iynode < nynodes; iynode++)
        for (int iznode = 0; iznode < nznodes; iznode++) {
          // Validation of the inactivated cells
          if (Activated[ixnode][iynode][iznode] == 0 && net_energy_transfer[ixnode][iynode][iznode] >= 1.e-5) printf("Bugs here! (Activated cells), energy = %f\n", net_energy_transfer[ixnode][iynode][iznode]);
          // Get the highest ke
          if (Activated[ixnode][iynode][iznode] == 1 && ke_real[ixnode][iynode][iznode] > el_ke)
            el_ke = ke_real[ixnode][iynode][iznode];
          // Get the lowest Ce
          if (Activated[ixnode][iynode][iznode] == 1 && CeT[ixnode][iynode][iznode] < el_specific_heat)
            el_specific_heat = CeT[ixnode][iynode][iznode];
        }
    // el_specific_heat = MAX(sh_min,el_specific_heat);
    // num_inner_timesteps = # of inner steps (thermal solves)
    // required this MD step to maintain a stable explicit solve
    double duration_temp;
    int num_inner_timesteps = 1;
    double inner_dt = update->dt;
    double stability_criterion = 0.0;
    int mxnodei, mxnodef, msize, mdisp;
    int rsize[numP], disp[numP];
    mxnodei = (pid * nxnodes) / numP;
    mxnodef = ((pid + 1) * nxnodes) / numP - 1;
    mdisp = mxnodei * nynodes * nznodes;
    msize = (mxnodef - mxnodei + 1) * nynodes * nznodes;
    MPI_Allgather(&msize, 1, MPI_INT, rsize, 1, MPI_INT, world);
    MPI_Allgather(&mdisp, 1, MPI_INT, disp, 1, MPI_INT, world);

    if (bulk_ttm != 1) bxsize = 1;
    double Te_bulk_first[bxsize][nynodes][nznodes];
    double Ta_bulk_first[bxsize][nynodes][nznodes];

    for (int ixnode = 0; ixnode < nxnodes; ixnode++)
      for (int iynode = 0; iynode < nynodes; iynode++)
        for (int iznode = 0; iznode < nznodes; iznode++) {
          T_electron_first[ixnode][iynode][iznode] = T_electron[ixnode][iynode][iznode];
        }

    if (bulk_ttm == 1) {
      for (int ixnode = 0; ixnode < bxsize; ixnode++)
        for (int iynode = 0; iynode < nynodes; iynode++)
          for (int iznode = 0; iznode < nznodes; iznode++) {
            Te_bulk_first[ixnode][iynode][iznode] = Te_bulk[ixnode][iynode][iznode];
            Ta_bulk_first[ixnode][iynode][iznode] = Ta_bulk[ixnode][iynode][iznode];
          }
    }

    do {

      duration_temp = duration;

      for (int ixnode = 0; ixnode < nxnodes; ixnode++)
        for (int iynode = 0; iynode < nynodes; iynode++)
          for (int iznode = 0; iznode < nznodes; iznode++) {
            T_electron[ixnode][iynode][iznode] = T_electron_first[ixnode][iynode][iznode];
          }

      if (bulk_ttm == 1) {
        for (int ixnode = 0; ixnode < bxsize; ixnode++)
          for (int iynode = 0; iynode < nynodes; iynode++)
            for (int iznode = 0; iznode < nznodes; iznode++) {
              Te_bulk[ixnode][iynode][iznode] = Te_bulk_first[ixnode][iynode][iznode];
              Ta_bulk[ixnode][iynode][iznode] = Ta_bulk_first[ixnode][iynode][iznode];
            }
      }

      update_parameters(); // Update for Te changes
      if (bulk_ttm == 1) {
        update_parameters_bulk(); // Update for Te/a_bulk changes
      }

      stability_criterion = 1.0 - 2.0 * inner_dt / el_specific_heat * (el_ke * (1.0 / dx / dx + 1.0 / dy / dy + 1.0 / dz / dz));

      if (stability_criterion < 0.0) {
        inner_dt = 1. / 4 * el_specific_heat / (el_ke * (1.0 / dx / dx + 1.0 / dy / dy + 1.0 / dz / dz));
      }
      num_inner_timesteps = static_cast<unsigned int>(update->dt / inner_dt) + 1;
      // num_inner_timesteps = 1;
      inner_dt = update->dt / double(num_inner_timesteps);

      if (pid == 0 && num_inner_timesteps > 500) {// 1000000
        char str[128];
        sprintf(str, "Too many inner timesteps: %d", num_inner_timesteps);
        error->warning(FLERR, str);
      }
      for (int ixnode = 0; ixnode < nxnodes; ixnode++)
        for (int iynode = 0; iynode < nynodes; iynode++)
          for (int iznode = 0; iznode < nznodes; iznode++)
            T_electron_old[ixnode][iynode][iznode] = T_electron[ixnode][iynode][iznode];

      for (int ith_inner_timestep = 0; ith_inner_timestep < num_inner_timesteps; ith_inner_timestep++) {
        // compute new electron T profile
        duration_temp += inner_dt;
        if (duration_temp <= 4.0 * width) laser(duration_temp);
        // Divide work into several parts and send them to different cores
        for (int ixnode = mxnodei; ixnode <= mxnodef; ixnode++)
          for (int iynode = 0; iynode < nynodes; iynode++)
            for (int iznode = 0; iznode < nznodes; iznode++) {
              int right_xnode = ixnode + 1;
              int right_ynode = iynode + 1;
              int right_znode = iznode + 1;
              if (right_xnode == nxnodes) right_xnode = nxnodes - 1;
              if (right_ynode == nynodes) right_ynode = 0;
              if (right_znode == nznodes) right_znode = 0;
              int left_xnode = ixnode - 1;
              int left_ynode = iynode - 1;
              int left_znode = iznode - 1;
              if (left_xnode == -1) left_xnode = 0;
              if (left_ynode == -1) left_ynode = nynodes - 1;
              if (left_znode == -1) left_znode = nznodes - 1;

              int cr_vac = 1;
              if (Activated[ixnode][iynode][iznode] == 0) cr_vac = 0;
              int cr_v_l_x = 1;
              if (Activated[left_xnode][iynode][iznode] == 0) cr_v_l_x = 0;
              int cr_v_r_x = 1;
              if (Activated[right_xnode][iynode][iznode] == 0) cr_v_r_x = 0;
              int cr_v_l_y = 1;
              if (Activated[ixnode][left_ynode][iznode] == 0) cr_v_l_y = 0;
              int cr_v_r_y = 1;
              if (Activated[ixnode][right_ynode][iznode] == 0) cr_v_r_y = 0;
              int cr_v_l_z = 1;
              if (Activated[ixnode][iynode][left_znode] == 0) cr_v_l_z = 0;
              int cr_v_r_z = 1;
              if (Activated[ixnode][iynode][right_znode] == 0) cr_v_r_z = 0;
              // double free_path = 377;
              if (cr_vac != 0) {

                int indic = 1;
                // -----------  Considering the great temperature gradient ------------------------------------
                double flux_l = 0.0, flux_r = 0.0;
                double f_L = 0.1; //flux limiter, ranges from 0.1 to 0.01
                if (cr_v_r_x == 1) {
                  flux_r = (ke_real[ixnode][iynode][iznode] + ke_real[right_xnode][iynode][iznode]) / 2.0 * (T_electron_old[right_xnode][iynode][iznode] - T_electron_old[ixnode][iynode][iznode]) / dx;
                  double over_LT = abs(T_electron_old[right_xnode][iynode][iznode] - T_electron_old[ixnode][iynode][iznode]) / dx / ((T_electron_old[right_xnode][iynode][iznode] + T_electron_old[ixnode][iynode][iznode]) / 2.0);

                  if (indic && free_path * over_LT > 0.1) {
                    double LT_r = 1.0 / over_LT; // electron temperature gradient length
                    double T_e = (T_electron_old[right_xnode][iynode][iznode] + T_electron_old[ixnode][iynode][iznode]) / 2.0;
                    double Zcharge = interpolation(ZTe, Zsize, T_e, 0.0);
                    double Ni = (nsum[right_xnode][iynode][iznode] + nsum[ixnode][iynode][iznode]) / 2.0 / del_vol;
                    double q_L = f_L * Zcharge * Ni * force->boltz * T_e * sqrt(force->boltz * T_e / m_e);
                    flux_r /= 1 + abs(flux_r) / q_L;
                  }

                }

                if (cr_v_l_x == 1) {
                  flux_l = (ke_real[ixnode][iynode][iznode] + ke_real[left_xnode][iynode][iznode]) / 2.0 * (T_electron_old[left_xnode][iynode][iznode] - T_electron_old[ixnode][iynode][iznode]) / dx;
                  double over_LT = abs(T_electron_old[left_xnode][iynode][iznode] - T_electron_old[ixnode][iynode][iznode]) / dx / ((T_electron_old[left_xnode][iynode][iznode] + T_electron_old[ixnode][iynode][iznode]) / 2.0);

                  if (indic && free_path * over_LT > 0.1) {
                    double LT_l = 1.0 / over_LT; // electron temperature gradient length
                    double T_e = (T_electron_old[left_xnode][iynode][iznode] + T_electron_old[ixnode][iynode][iznode]) / 2.0;
                    double Zcharge = interpolation(ZTe, Zsize, T_e, 0.0);
                    double Ni = (nsum[left_xnode][iynode][iznode] + nsum[ixnode][iynode][iznode]) / 2.0 / del_vol;
                    double q_L = f_L * Zcharge * Ni * force->boltz * T_e * sqrt(force->boltz * T_e / m_e);
                    flux_l /= 1 + abs(flux_l) / q_L;
                  }

                }
                // ------------------------------------------------------------------------------
                              // Energy conduction within Electron subsystem
                energy_conduction[ixnode][iynode][iznode] = ((cr_v_r_x * flux_r +
                  cr_v_l_x * flux_l) / dx +
                  (cr_v_r_y * (ke_real[ixnode][iynode][iznode] + ke_real[ixnode][right_ynode][iznode]) / 2.0 * (T_electron_old[ixnode][right_ynode][iznode] - T_electron_old[ixnode][iynode][iznode]) / dy +
                    cr_v_l_y * (ke_real[ixnode][iynode][iznode] + ke_real[ixnode][left_ynode][iznode]) / 2.0 * (T_electron_old[ixnode][left_ynode][iznode] - T_electron_old[ixnode][iynode][iznode]) / dy) / dy +
                  (cr_v_r_z * (ke_real[ixnode][iynode][iznode] + ke_real[ixnode][iynode][right_znode]) / 2.0 * (T_electron_old[ixnode][iynode][right_znode] - T_electron_old[ixnode][iynode][iznode]) / dz +
                    cr_v_l_z * (ke_real[ixnode][iynode][iznode] + ke_real[ixnode][iynode][left_znode]) / 2.0 * (T_electron_old[ixnode][iynode][left_znode] - T_electron_old[ixnode][iynode][iznode]) / dz) / dz);
                T_electron[ixnode][iynode][iznode] = T_electron_old[ixnode][iynode][iznode] + inner_dt / CeT[ixnode][iynode][iznode] * energy_conduction[ixnode][iynode][iznode]
                  ;

                if (duration_temp <= 4.0 * width && skin_layer[ixnode][iynode][iznode] == 0) printf("Bug!!! skin_layer is zero! ixnode = %d\n Te = %f\n", ixnode, T_electron_old[ixnode][iynode][iznode]);

                // Energy trasfer from laser and to lattice
                if (duration_temp <= 4.0 * width) {
                  T_electron[ixnode][iynode][iznode] += inner_dt / CeT[ixnode][iynode][iznode] * (mult_factor[ixnode][iynode][iznode] / skin_layer[ixnode][iynode][iznode] - net_energy_transfer[ixnode][iynode][iznode] / del_vol);
                }
                else {
                  T_electron[ixnode][iynode][iznode] += inner_dt / CeT[ixnode][iynode][iznode] * (0 - net_energy_transfer[ixnode][iynode][iznode] / del_vol);
                }

                // Heat transfer from/to the bulk
                if (bulk_ttm == 1 && ixnode == Tes[iynode][iznode]) {

                  double ketr = (ke_real[ixnode][iynode][iznode] * dx_bulk + ke_real_bulk[0][iynode][iznode] * dx) / (dx_bulk + dx);
                  double EC_bulk = ketr * (T_electron_old[ixnode][iynode][iznode] - Te_bulk[0][iynode][iznode]) / (dx_bulk + dx) / dx;
                  energy_conduction[ixnode][iynode][iznode] -= EC_bulk;
                  T_electron[ixnode][iynode][iznode] -= inner_dt / CeT[ixnode][iynode][iznode] * EC_bulk;
                }

              }
              else {
                if (Activated[ixnode][iynode][iznode] == 1) printf("Bug here abc\n");
                T_electron[ixnode][iynode][iznode] = T_electron_old[ixnode][iynode][iznode];
              }

            }

        if (bulk_ttm == 1) {
          int mxnodei_bulk, mxnodef_bulk, msize_bulk, mdisp_bulk;
          int rsize_bulk[numP], disp_bulk[numP];
          double ci = 3.0 * ionic_density * force->boltz;
          mxnodei_bulk = (pid * bxsize) / numP;
          mxnodef_bulk = ((pid + 1) * bxsize) / numP - 1;
          mdisp_bulk = mxnodei_bulk * nynodes * nznodes;
          msize_bulk = (mxnodef_bulk - mxnodei_bulk + 1) * nynodes * nznodes;
          MPI_Allgather(&msize_bulk, 1, MPI_INT, rsize_bulk, 1, MPI_INT, world);
          MPI_Allgather(&mdisp_bulk, 1, MPI_INT, disp_bulk, 1, MPI_INT, world);
          double Te_bulk_buffer[mxnodef_bulk - mxnodei_bulk + 1][nynodes][nznodes];
          double Ta_bulk_buffer[mxnodef_bulk - mxnodei_bulk + 1][nynodes][nznodes];
          // Bulk Temperature calculation
          if (duration_temp <= 4.0 * width) laser_bulk();
          for (int i = 0; i < mxnodef_bulk - mxnodei_bulk + 1; i++)
            for (int iynode = 0; iynode < nynodes; iynode++)
              for (int iznode = 0; iznode < nznodes; iznode++) {
                int ixnode = i + mxnodei_bulk;
                int right_xnode = ixnode + 1;
                int right_ynode = iynode + 1;
                int right_znode = iznode + 1;
                // if (right_xnode == bxsize) right_xnode = 0;
                if (right_ynode == nynodes) right_ynode = 0;
                if (right_znode == nznodes) right_znode = 0;
                int left_xnode = ixnode - 1;
                int left_ynode = iynode - 1;
                int left_znode = iznode - 1;
                // if (left_xnode == -1) left_xnode = bxsize - 1;
                if (left_ynode == -1) left_ynode = nynodes - 1;
                if (left_znode == -1) left_znode = nznodes - 1;

                if (ixnode == 0) {
                  double ketr = (ke_real[Tes[iynode][iznode]][iynode][iznode] * dx_bulk + ke_real_bulk[0][iynode][iznode] * dx) / (dx_bulk + dx);
                  double ketr2 = (ke_real_bulk[1][iynode][iznode] + ke_real_bulk[0][iynode][iznode]) / 2.0;

                  Te_bulk_buffer[i][iynode][iznode] = Te_bulk[0][iynode][iznode] + inner_dt / CeT_bulk[0][iynode][iznode] * (ketr * (T_electron_old[Tes[iynode][iznode]][iynode][iznode] - Te_bulk[0][iynode][iznode]) / (dx_bulk + dx) + ketr2 * (Te_bulk[1][iynode][iznode] - Te_bulk[0][iynode][iznode]) / dx_bulk) / dx_bulk;

                  double kitr = 0.5 * (ki_bulk[0][iynode][iznode] + ki_bulk[1][iynode][iznode]);
                  Ta_bulk_buffer[i][iynode][iznode] = Ta_bulk[0][iynode][iznode] + inner_dt / ci * (kitr * (Ta_bulk[1][iynode][iznode] - Ta_bulk[0][iynode][iznode]) / dx_bulk) / dx_bulk;

                }
                else if (ixnode == bxsize - 1) {
                  double ketr = (ke_real_bulk[bxsize - 1][iynode][iznode] + ke_real_bulk[bxsize - 2][iynode][iznode]) / 2.0;

                  Te_bulk_buffer[i][iynode][iznode] = Te_bulk[bxsize - 1][iynode][iznode] + inner_dt / CeT_bulk[bxsize - 1][iynode][iznode] * (ketr * (Te_bulk[bxsize - 2][iynode][iznode] - Te_bulk[bxsize - 1][iynode][iznode]) / dx_bulk) / dx_bulk;

                  double kitr = 0.5 * (ki_bulk[bxsize - 2][iynode][iznode] + ki_bulk[bxsize - 1][iynode][iznode]);
                  Ta_bulk_buffer[i][iynode][iznode] = Ta_bulk[bxsize - 1][iynode][iznode] + inner_dt / ci * (kitr * (Ta_bulk[bxsize - 2][iynode][iznode] - Ta_bulk[bxsize - 1][iynode][iznode]) / dx_bulk) / dx_bulk;

                }
                else {
                  // printf("ixnode = %d, bxsize = %d\n",ixnode,bxsize);
                  double ketr1 = (ke_real_bulk[left_xnode][iynode][iznode] + ke_real_bulk[ixnode][iynode][iznode]) / 2.0;

                  double ketr2 = (ke_real_bulk[right_xnode][iynode][iznode] + ke_real_bulk[ixnode][iynode][iznode]) / 2.0;

                  Te_bulk_buffer[i][iynode][iznode] = Te_bulk[ixnode][iynode][iznode] + inner_dt / CeT_bulk[ixnode][iynode][iznode] * (ketr1 * (Te_bulk[left_xnode][iynode][iznode] - Te_bulk[ixnode][iynode][iznode]) / dx_bulk + ketr2 * (Te_bulk[right_xnode][iynode][iznode] - Te_bulk[ixnode][iynode][iznode]) / dx_bulk) / dx_bulk;

                  double kitr1 = 0.5 * (ki_bulk[left_xnode][iynode][iznode] + ki_bulk[ixnode][iynode][iznode]);
                  double kitr2 = 0.5 * (ki_bulk[right_xnode][iynode][iznode] + ki_bulk[ixnode][iynode][iznode]);
                  Ta_bulk_buffer[i][iynode][iznode] = Ta_bulk[ixnode][iynode][iznode] + inner_dt / ci * (kitr1 * (Ta_bulk[left_xnode][iynode][iznode] - Ta_bulk[ixnode][iynode][iznode]) / dx_bulk + kitr2 * (Ta_bulk[right_xnode][iynode][iznode] - Ta_bulk[ixnode][iynode][iznode]) / dx_bulk) / dx_bulk;
                }

                double ketr1 = (ke_real_bulk[ixnode][left_ynode][iznode] + ke_real_bulk[ixnode][iynode][iznode]) / 2.0;

                double ketr2 = (ke_real_bulk[ixnode][right_ynode][iznode] + ke_real_bulk[ixnode][iynode][iznode]) / 2.0;

                double ketr3 = (ke_real_bulk[ixnode][iynode][left_znode] + ke_real_bulk[ixnode][iynode][iznode]) / 2.0;

                double ketr4 = (ke_real_bulk[ixnode][iynode][right_znode] + ke_real_bulk[ixnode][iynode][iznode]) / 2.0;

                Te_bulk_buffer[i][iynode][iznode] += inner_dt / CeT_bulk[ixnode][iynode][iznode] * ((ketr1 * (Te_bulk[ixnode][left_ynode][iznode] - Te_bulk[ixnode][iynode][iznode]) / dy + ketr2 * (Te_bulk[ixnode][right_ynode][iznode] - Te_bulk[ixnode][iynode][iznode]) / dy) / dy +
                  (ketr3 * (Te_bulk[ixnode][iynode][left_znode] - Te_bulk[ixnode][iynode][iznode]) / dz + ketr4 * (Te_bulk[ixnode][iynode][right_znode] - Te_bulk[ixnode][iynode][iznode]) / dz) / dz);

                double kitr1 = 0.5 * (ki_bulk[ixnode][left_ynode][iznode] + ki_bulk[ixnode][iynode][iznode]);
                double kitr2 = 0.5 * (ki_bulk[ixnode][right_ynode][iznode] + ki_bulk[ixnode][iynode][iznode]);
                double kitr3 = 0.5 * (ki_bulk[ixnode][iynode][left_znode] + ki_bulk[ixnode][iynode][iznode]);
                double kitr4 = 0.5 * (ki_bulk[ixnode][iynode][right_znode] + ki_bulk[ixnode][iynode][iznode]);

                Ta_bulk_buffer[i][iynode][iznode] += inner_dt / ci * ((kitr1 * (Ta_bulk[ixnode][left_ynode][iznode] - Ta_bulk[ixnode][iynode][iznode]) / dy + kitr2 * (Ta_bulk[ixnode][right_ynode][iznode] - Ta_bulk[ixnode][iynode][iznode]) / dy) / dy +
                  (kitr3 * (Ta_bulk[ixnode][iynode][left_znode] - Ta_bulk[ixnode][iynode][iznode]) / dz + kitr4 * (Ta_bulk[ixnode][iynode][right_znode] - Ta_bulk[ixnode][iynode][iznode]) / dz) / dz);

                double Ee2Ea = GT_bulk[ixnode][iynode][iznode] * (Te_bulk[ixnode][iynode][iznode] - Ta_bulk[ixnode][iynode][iznode]);
                Te_bulk_buffer[i][iynode][iznode] -= inner_dt / CeT_bulk[ixnode][iynode][iznode] * Ee2Ea;
                Ta_bulk_buffer[i][iynode][iznode] += inner_dt / ci * Ee2Ea;
                // laser input
                if (duration_temp <= 4.0 * width) {
                  Te_bulk_buffer[i][iynode][iznode] += inner_dt / CeT_bulk[ixnode][iynode][iznode] * mult_factor_bulk[ixnode][iynode][iznode] / skin_layer_bulk[ixnode][iynode][iznode];
                }
                bool flag_melt = ((Ta_bulk_buffer[i][iynode][iznode] - T_melt) * (Ta_bulk[ixnode][iynode][iznode] - T_melt) <= 0.0);
                if (flag_melt) {
                  E_melt_buffer[ixnode][iynode][iznode] += (Ta_bulk_buffer[i][iynode][iznode] - T_melt) * ci;
                  if (E_melt_buffer[ixnode][iynode][iznode] < 0.0) {
                    Ta_bulk_buffer[i][iynode][iznode] = T_melt - (0.0 - E_melt_buffer[ixnode][iynode][iznode]) / ci;
                    E_melt_buffer[ixnode][iynode][iznode] = 0.0;
                  }
                  else if (E_melt_buffer[ixnode][iynode][iznode] > Latent_melt) {
                    Ta_bulk_buffer[i][iynode][iznode] = T_melt + (E_melt_buffer[ixnode][iynode][iznode] - Latent_melt) / ci;
                    E_melt_buffer[ixnode][iynode][iznode] = Latent_melt;
                  }
                  else {
                    Ta_bulk_buffer[i][iynode][iznode] = T_melt;
                  }
                }

              }
          MPI_Allgatherv(&Te_bulk_buffer[0][0][0], msize_bulk, MPI_DOUBLE, &Te_bulk[0][0][0], rsize_bulk, disp_bulk, MPI_DOUBLE, world);
          MPI_Allgatherv(&Ta_bulk_buffer[0][0][0], msize_bulk, MPI_DOUBLE, &Ta_bulk[0][0][0], rsize_bulk, disp_bulk, MPI_DOUBLE, world);
          MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &E_melt_buffer[0][0][0], rsize_bulk, disp_bulk, MPI_DOUBLE, world);
        }

        //  Update all new temperatures for all cores
        MPI_Allgatherv(&T_electron[mxnodei][0][0], msize, MPI_DOUBLE, &T_electron_old[0][0][0], rsize, disp, MPI_DOUBLE, world);

        for (int ixnode = 0; ixnode < nxnodes; ixnode++)
          for (int iynode = 0; iynode < nynodes; iynode++)
            for (int iznode = 0; iznode < nznodes; iznode++)
              T_electron[ixnode][iynode][iznode] = T_electron_old[ixnode][iynode][iznode]; // unify Te bewteen cores

        update_parameters(); // Update for Te changes
        if (bulk_ttm == 1) {
          update_parameters_bulk(); // Update for Te/a_bulk changes
        }

        for (int ixnode = mxnodei; ixnode <= mxnodef; ixnode++)
          for (int iynode = 0; iynode < nynodes; iynode++)
            for (int iznode = 0; iznode < nznodes; iznode++) {
              if (Activated[ixnode][iynode][iznode] == 0)
                continue;
              if (CeT[ixnode][iynode][iznode] < el_specific_heat)
                el_specific_heat = CeT[ixnode][iynode][iznode];
              if (ke_real[ixnode][iynode][iznode] > el_ke)
                el_ke = ke_real[ixnode][iynode][iznode];
            }

      }

      MPI_Allreduce(MPI_IN_PLACE, &el_ke, 1, MPI_DOUBLE, MPI_MAX, world);
      MPI_Allreduce(MPI_IN_PLACE, &el_specific_heat, 1, MPI_DOUBLE, MPI_MIN, world);
      //el_specific_heat = MAX(sh_min,el_specific_heat);

      stability_criterion = 1.0 - 2.0 * inner_dt / el_specific_heat * (el_ke * (1.0 / dx / dx + 1.0 / dy / dy + 1.0 / dz / dz));
      // stability_criterion = 1.0;
      //printf("Here is debugging time!\n");
      // if (pid == 0) printf("Ce = %f;   ke = %f;   s_c = %f;   i_steps = %d\n",el_specific_heat,el_ke, stability_criterion, num_inner_timesteps);
    } while (stability_criterion < 0.0);

    // if (pid == 0) printf("inner_step is %d in step %d\n", num_inner_timesteps, update->ntimestep);

    duration = duration_temp;
    // if (pid == 0) printf("Duration is %f \n", duration);

    // output nodal temperatures for current timestep
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &energy_conduction[0][0][0], rsize, disp, MPI_DOUBLE, world);
    Tempout();
    Otherout();
  }

  if (premode == 1) {
    if ((nfileevery) && !(update->ntimestep % nfileevery)) {
      for (int ixnode = 0; ixnode < nxnodes; ixnode++)
        for (int iynode = 0; iynode < nynodes; iynode++)
          for (int iznode = 0; iznode < nznodes; iznode++) {
            nsum[ixnode][iynode][iznode] = 0;
            sum_mass_vsq[ixnode][iynode][iznode] = 0.0;
          }

      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
          if (rmass) massone = rmass[i];
          else massone = mass[type[i]];
          double xscale = (x[i][0] - domain->boxlo[0]) / domain->xprd;
          double yscale = (x[i][1] - domain->boxlo[1]) / domain->yprd;
          double zscale = (x[i][2] - domain->boxlo[2]) / domain->zprd;
          int ixnode = static_cast<int>(xscale * nxnodes);
          int iynode = static_cast<int>(yscale * nynodes);
          int iznode = static_cast<int>(zscale * nznodes);
          // Eliminate numerical errors
          if (ixnode > nxnodes - 1) ixnode = nxnodes - 1;
          if (iynode > nynodes - 1) iynode = nynodes - 1;
          if (iznode > nznodes - 1) iznode = nznodes - 1;
          if (ixnode < 0) ixnode = 0;
          if (iynode < 0) iynode = 0;
          if (iznode < 0) iznode = 0;

          nsum[ixnode][iynode][iznode] += 1;
          double vsq = v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2];
          sum_mass_vsq[ixnode][iynode][iznode] += massone * vsq;
        }
      }

      MPI_Allreduce(MPI_IN_PLACE, &nsum[0][0][0], total_nnodes, MPI_INT, MPI_SUM, world);
      MPI_Allreduce(MPI_IN_PLACE, &sum_mass_vsq[0][0][0], total_nnodes, MPI_DOUBLE, MPI_SUM, world);

      double crit_num = crit_num_f * ionic_density * del_vol;

      for (int ixnode = 0; ixnode < nxnodes; ixnode++)
        for (int iynode = 0; iynode < nynodes; iynode++)
          for (int iznode = 0; iznode < nznodes; iznode++) {
            if ((double)nsum[ixnode][iynode][iznode] > crit_num) {
              T_a[ixnode][iynode][iznode] = sum_mass_vsq[ixnode][iynode][iznode] / (3.0 * force->boltz * nsum[ixnode][iynode][iznode] / force->mvv2e);
            }
          }

      if (comm->me == 0) {
        FILE* fp = fopen(fp_Ta_out.c_str(), "a");
        fprintf(fp, BIGINT_FORMAT, update->ntimestep);
        fprintf(fp, "\n--------------------------------------------------------------------");
        double crit_num = crit_num_f * ionic_density * del_vol;
        for (int ixnode = 0; ixnode < nxnodes; ixnode++) {
          int num = 0;
          double ptemp = 0.0;
          for (int iynode = 0; iynode < nynodes; iynode++)
            for (int iznode = 0; iznode < nznodes; iznode++) {
              if ((double)nsum[ixnode][iynode][iznode] > crit_num) {
                num++;
                ptemp += T_a[ixnode][iynode][iznode];
              }
            }
          if (num != 0) {
            fprintf(fp, "\n%d\t%f", ixnode, ptemp);
          }
        }
        fprintf(fp, "\n\n");
        fclose(fp);
      }

    }
  }

}


/* ---------------------------------------------------------------------- */

void FixFEMTO3D::post_force_setup(int vflag)
{
  // int pid;
  // MPI_Comm_rank(world,&pid);
  double** f = atom->f;
  int* mask = atom->mask;
  int nlocal = atom->nlocal;
  // apply langevin forces that have been stored from previous run
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      f[i][0] += flangevin[i][0];
      f[i][1] += flangevin[i][1];
      f[i][2] += flangevin[i][2];
    }
  }
  // if (pid == 0) printf("Here is part post_force_setup!\n");
}

/* ---------------------------------------------------------------------- */

void FixFEMTO3D::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa - 1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixFEMTO3D::post_force_respa_setup(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa - 1) post_force_setup(vflag);
}

/* ---------------------------------------------------------------------- */

void FixFEMTO3D::reset_dt()
{

}

void FixFEMTO3D::set_initial_temperatures()
{
  // Initialize all cell temperatures
  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++) {
        Activated[ixnode][iynode][iznode] = 0;
        T_electron[ixnode][iynode][iznode] = T_init;
      }

  if (bulk_ttm == 1) {
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++) {
        Tes[iynode][iznode] = 0;
        x_max[iynode][iznode] = domain->boxlo[0];
        for (int ixnode = 0; ixnode < bxsize; ixnode++) {
          Te_bulk[ixnode][iynode][iznode] = T_init;
          Ta_bulk[ixnode][iynode][iznode] = T_init;
          E_melt_buffer[ixnode][iynode][iznode] = 0.0;
        }
      }
  }

}

/* ---------------------------------------------------------------------- */

void FixFEMTO3D::laser(double duration_temp)
{
  double T_e, T_i, Ni;
  reflectivity = 0.0;
  for (int ixnode = t_surface_l; ixnode < t_surface_r; ixnode++) // during laser, surface doesn't move
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++) {
        T_e = T_electron[ixnode][iynode][iznode];
        T_i = T_a[ixnode][iynode][iznode];
        //Ni = nsum[ixnode][iynode][iznode]/del_vol;
        Ni = ionic_density; // During laser, Ni is not changing
        skin_layer[ixnode][iynode][iznode] = interpolation(PenTe, Pensize, T_e, 0.0);

        if (ixnode == t_surface_l) {
          double refl = interpolation(ReflecTe, Reflecsize, T_e, 0.0, 1.0);
          mult_factor[ixnode][iynode][iznode] = 1.0 - refl; // Absorption
          reflectivity += 1 - mult_factor[ixnode][iynode][iznode];
          // if (pid == 0) printf("multi = %f;   skin = %f\n",mult_factor[ixnode][iynode][iznode],skin_layer[ixnode][iynode][iznode]);
          mult_factor[ixnode][iynode][iznode] *= 0.9395 * intensity * exp(-2.77259 * (duration_temp - 2.0 * width) * (duration_temp - 2.0 * width) / width / width); // guassian in time;
        }
        else {
          mult_factor[ixnode][iynode][iznode] = mult_factor[ixnode - 1][iynode][iznode] * exp((-1.0) * dx / skin_layer[ixnode - 1][iynode][iznode]);
        }

      }
  reflectivity /= nynodes * nznodes;
  // skin_layer = 180; // Angstrom
}

void FixFEMTO3D::laser_bulk()
{
  // int numP, pid;
  // MPI_Comm_size(world, &numP);
  // MPI_Comm_rank(world, &pid);
  double T_e, T_i, Ni;
  for (int ixnode = 0; ixnode < bxsize; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++) {
        T_e = Te_bulk[ixnode][iynode][iznode];
        T_i = Ta_bulk[ixnode][iynode][iznode];
        Ni = ionic_density;
        skin_layer_bulk[ixnode][iynode][iznode] = interpolation(PenTe, Pensize, T_e, 0.0);
        if (skin_layer_bulk[ixnode][iynode][iznode] < 0.0)
          skin_layer_bulk[ixnode][iynode][iznode] = 0.0;

        if (ixnode == 0) {
          if (pid == 0 && t_surface_r - 1 != Tes[iynode][iznode]) printf("t_surface_r = %d\tTes = %d\n", t_surface_r, Tes[iynode][iznode]);
          mult_factor_bulk[0][iynode][iznode] = mult_factor[Tes[iynode][iznode]][iynode][iznode] * exp((-1.0) * dx / skin_layer[Tes[iynode][iznode]][iynode][iznode]);
        }
        else {
          mult_factor_bulk[ixnode][iynode][iznode] = mult_factor_bulk[ixnode - 1][iynode][iznode] * exp((-1.0) * dx / skin_layer_bulk[ixnode - 1][iynode][iznode]);
        }
      }
}


double FixFEMTO3D::MyCe(double T_e, double Ni, double T_i) {
  int isFromTable = 1;
  double Ce;
  if (isFromTable) {
    Ce = interpolation(CeTe, Cesize, T_e, 0.0) / ionic_density * Ni;
  }
  else {
    Ce = 0;
    // printf("Cep = %lg\n", Ce);
  }
  return Ce;
}

double FixFEMTO3D::MyKet(double T_e, double Ni, double T_i) {

  int isFromTable = 0;
  double ket;
  double pp = 1.0;
  if (T_i > T_e) {
    //ket = 0;
    //return ket;
    pp = 0.2;
  }
  if (isFromTable) {
    ket = interpolation(KeTe, Kesize, T_e, 0.0);
  }
  else { //Numerical Calculation

    int Rethfeld = 1;
    int Petrov = 2;
    int KeMode = Rethfeld;

    if (KeMode == Rethfeld) {
      // Rethfeld
      double thetae, thetal;
      thetae = T_e / 64200;
      thetal = T_i / 64200;
      ket = 353.0 / e_charge / 1.e22 * thetae * pow((thetae * thetae + 0.16), 1.25) * (thetae * thetae + 0.44)
        / sqrt(thetae * thetae + 0.092) / (thetae * thetae + 0.16 * thetal);
    }
    else if (KeMode == Petrov) {
      // Petrov ----------------------------
      double kee, t, x, k0, kes;
      x = Ni / ionic_density;
      t = 5.62e-5 * T_e / x;
      kee = 1.e4 * (1 + 0.03 * sqrt(t) - 0.2688 * t + 0.9722 * t * t) * pow(x, 4 / 3) / 9.294 / t;
      k0 = 131 * t * (1 + 3.07 * t * t) / (1 + 1.08 * pow(t, 2.07));
      if (x > 0.88) {
        double xrt = 19.3 / 19.5;
        double ybar = 1.6678 * pow(x, 8.84) / (1 + 0.6678 * pow(x, 4.92));
        kes = 4.5876e4 * pow(x / xrt, 4 / 3) * ybar / T_i * k0;
      }
      else {
        double xl = 0.887179 - 0.0328321 * (T_i / 1000 - 1.337)
          - 0.0030982 * pow(T_i / 1000 - 1.337, 2) - 0.000164884 * pow(T_i / 1000 - 1.337, 3);
        double rl = 148.5 + 119.3 * T_i / 1000 * (15337 / (14000 + T_i));
        kes = k0 * 3254 / rl * x * pow(x / xl, 2);
      }
      ket = 1.0 / (1.0 / kee + 1.0 / kes) / (e_charge * 1.e22); // Convert to LAMMPS units
    }
    else {
      error->all(FLERR, "Wrong Ke Mode");
    }

  }
  ket *= pp;
  return ket;
}

double FixFEMTO3D::MyG(double T_e, double Ni, double T_i) {
  int isFromTable = 1;
  double G;
  double pp = 1.0;
  if (T_e < T_i) {
    pp = 5.0;
  }

  if (isFromTable) {
    // Lin's data
    G = interpolation(GTe, Gsize, T_e, 0.0);
  }
  else {
    // Constant
    G = 3.e16 / e_charge / 1.e42;

    /*
    // Chen's data
    if (Ni > 0.88*ionic_density) {
      G = 2.2e16 /e_charge/1.e42;
    } else {
      G = 2.6e16 /e_charge/1.e42;
    }
    G *= 9.756e-5*(T_e+T_i) + 1;
    */

    /*
    // Ashitkov
    double Tev = T_e * force->boltz;
    double ka = 4;
    G = (0.2 + 4.3/ka * pow(Tev,3.6) / (1+pow(Tev,3.5)+0.9*pow(Tev,4.1) ) ) * pow((Ni/ionic_density),5./3) * 1.e17/e_charge/1.e42;
    */

  }
  // printf("G = %lg\n", G);
  G *= pp;
  return G;
}

void FixFEMTO3D::ChangeType()
{
  double** x = atom->x;
  double** v = atom->v;
  double** f = atom->f;
  double* mass = atom->mass;
  double* rmass = atom->rmass;
  int* type = atom->type;
  int* mask = atom->mask;
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      double xscale = (x[i][0] - domain->boxlo[0]) / domain->xprd;
      double yscale = (x[i][1] - domain->boxlo[1]) / domain->yprd;
      double zscale = (x[i][2] - domain->boxlo[2]) / domain->zprd;
      int ixnode = static_cast<int>(xscale * nxnodes);
      int iynode = static_cast<int>(yscale * nynodes);
      int iznode = static_cast<int>(zscale * nznodes);
      // Eliminate numerical errors
      if (ixnode > nxnodes - 1) ixnode = nxnodes - 1;
      if (iynode > nynodes - 1) iynode = nynodes - 1;
      if (iznode > nznodes - 1) iznode = nznodes - 1;
      if (ixnode < 0) ixnode = 0;
      if (iynode < 0) iynode = 0;
      if (iznode < 0) iznode = 0;
      double TeV = T_electron[ixnode][iynode][iznode] * force->boltz;
      // printf("%f %f\n",T_electron[ixnode][iynode][iznode], force->boltz);

      if (TeV < 1.0)
        type[i] = 1; // Au_0.1
      else if (TeV >= 1.0 && TeV < 2.0)
        type[i] = 2; // Au_1.5
      else if (TeV >= 2.0 && TeV < 4.0)
        type[i] = 3; // Au_3.0
      else if (TeV >= 4.0 && TeV < 5.0)
        type[i] = 4; // Au_4.5
      else
        type[i] = 5; // Au_6.0

      /*
      if(TeV < 1.0)
        type[i] = MAX(type[i],1); // Au_0.1
      else if(TeV >= 1.0 && TeV < 2.0)
        type[i] = MAX(type[i],2); // Au_1.5
      else if(TeV >= 2.0 && TeV < 4.0)
        type[i] = MAX(type[i],3); // Au_3.0
      else if(TeV >= 4.0 && TeV < 5.0)
        type[i] = MAX(type[i],4); // Au_4.5
      else
        type[i] = 5; // Au_6.0
      */
      // type[i] = 1;
      // printf("TeV = %f, type = %d\n", TeV, type[i]);

    }
  }
}

void FixFEMTO3D::update_Ta()
{
  double** x = atom->x;
  double** v = atom->v;
  double** f = atom->f;
  double* mass = atom->mass;
  double* rmass = atom->rmass;
  int* type = atom->type;
  int* mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++) {
        nsum[ixnode][iynode][iznode] = 0;
        sum_mass_vsq[ixnode][iynode][iznode] = 0.0;
        sum_mass_v[ixnode][iynode][iznode] = 0.0;
        sum_mass[ixnode][iynode][iznode] = 0.0;
        net_energy_transfer[ixnode][iynode][iznode] = 0.0;
        average_v[ixnode][iynode][iznode] = 0.0;
        T_electron_old[ixnode][iynode][iznode] = T_electron[ixnode][iynode][iznode];
      }

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      double xscale = (x[i][0] - domain->boxlo[0]) / domain->xprd;
      double yscale = (x[i][1] - domain->boxlo[1]) / domain->yprd;
      double zscale = (x[i][2] - domain->boxlo[2]) / domain->zprd;
      int ixnode = static_cast<int>(xscale * nxnodes);
      int iynode = static_cast<int>(yscale * nynodes);
      int iznode = static_cast<int>(zscale * nznodes);
      // Eliminate numerical errors
      if (ixnode > nxnodes - 1) ixnode = nxnodes - 1;
      if (iynode > nynodes - 1) iynode = nynodes - 1;
      if (iznode > nznodes - 1) iznode = nznodes - 1;
      if (ixnode < 0) ixnode = 0;
      if (iynode < 0) iynode = 0;
      if (iznode < 0) iznode = 0;

      nsum[ixnode][iynode][iznode] += 1;
      double vsq = v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2];
      sum_mass[ixnode][iynode][iznode] += massone;
      sum_mass_v[ixnode][iynode][iznode] += massone * v[i][0]; // momentum in laser direction
      sum_mass_vsq[ixnode][iynode][iznode] += massone * vsq;
      net_energy_transfer[ixnode][iynode][iznode] += flangevin[i][0] * v[i][0] + flangevin[i][1] * v[i][1] + flangevin[i][2] * v[i][2];
      if (x[i][0] < xmin) xmin = x[i][0];
      if (bulk_ttm == 1 && x_max[iynode][iznode] < x[i][0]) x_max[iynode][iznode] = x[i][0];
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &xmin, 1, MPI_DOUBLE, MPI_MIN, world);
  MPI_Allreduce(MPI_IN_PLACE, &nsum[0][0][0], total_nnodes, MPI_INT, MPI_SUM, world);
  MPI_Allreduce(MPI_IN_PLACE, &sum_mass_vsq[0][0][0], total_nnodes, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(MPI_IN_PLACE, &sum_mass_v[0][0][0], total_nnodes, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(MPI_IN_PLACE, &sum_mass[0][0][0], total_nnodes, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(MPI_IN_PLACE, &net_energy_transfer[0][0][0], total_nnodes, MPI_DOUBLE, MPI_SUM, world);

  if (bulk_ttm == 1) {
    MPI_Allreduce(MPI_IN_PLACE, &x_max[0][0], nynodes * nznodes, MPI_DOUBLE, MPI_MAX, world);
  }

  // Activation
  if (duration <= 4.0 * width) {
    crit_num_f = 0.5;
  }
  else {
    crit_num_f = 0.1;
  }

  double crit_num = crit_num_f * ionic_density * del_vol;
  // if (duration <= 5.0) crit_num = 0.05*ionic_density*del_vol; // for the surface atoms in the early time 4.0*width

  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++) {
        if ((double)nsum[ixnode][iynode][iznode] > crit_num) {
          // Activated[ixnode][iynode][iznode] = 1;
          average_v[ixnode][iynode][iznode] = sum_mass_v[ixnode][iynode][iznode] / sum_mass[ixnode][iynode][iznode];
          T_a[ixnode][iynode][iznode] = (sum_mass_vsq[ixnode][iynode][iznode] - sum_mass[ixnode][iynode][iznode] * average_v[ixnode][iynode][iznode] * average_v[ixnode][iynode][iznode]) / (3.0 * force->boltz * nsum[ixnode][iynode][iznode] / force->mvv2e);
          // T_a[ixnode][iynode][iznode] = sum_mass_vsq[ixnode][iynode][iznode]/(3.0*force->boltz*nsum[ixnode][iynode][iznode]/force->mvv2e);
          // if(pid == 0) printf("Ta = %f at ( %d, %d, %d )\n",T_a[ixnode][iynode][iznode],ixnode,iynode,iznode);
          // T_a[ixnode][iynode][iznode] = 300.0;
          if (Activated[ixnode][iynode][iznode] == 0) {
            int neigh = 0;
            double temp = 0.0;
            int right_xnode = ixnode + 1;
            int right_ynode = iynode + 1;
            int right_znode = iznode + 1;
            if (right_xnode == nxnodes) right_xnode = nxnodes - 1;
            if (right_ynode == nynodes) right_ynode = 0;
            if (right_znode == nznodes) right_znode = 0;
            int left_xnode = ixnode - 1;
            int left_ynode = iynode - 1;
            int left_znode = iznode - 1;
            if (left_xnode == -1) left_xnode = 0;
            if (left_ynode == -1) left_ynode = nynodes - 1;
            if (left_znode == -1) left_znode = nznodes - 1;
            if (Activated[right_xnode][iynode][iznode] == 1) {
              neigh += 1;
              temp += T_electron_old[right_xnode][iynode][iznode];
            }
            if (Activated[ixnode][right_ynode][iznode] == 1) {
              neigh += 1;
              temp += T_electron_old[ixnode][right_ynode][iznode];
            }
            if (Activated[ixnode][iynode][right_znode] == 1) {
              neigh += 1;
              temp += T_electron_old[ixnode][iynode][right_znode];
            }
            if (Activated[left_xnode][iynode][iznode] == 1) {
              neigh += 1;
              temp += T_electron_old[left_xnode][iynode][iznode];
            }
            if (Activated[ixnode][left_ynode][iznode] == 1) {
              neigh += 1;
              temp += T_electron_old[ixnode][left_ynode][iznode];
            }
            if (Activated[ixnode][iynode][left_znode] == 1) {
              neigh += 1;
              temp += T_electron_old[ixnode][iynode][left_znode];
            }
            if (neigh >= 1) {
              T_electron[ixnode][iynode][iznode] = temp / neigh;
            }
            else {
              T_electron[ixnode][iynode][iznode] = T_a[ixnode][iynode][iznode];
            }
          }
        }
      }

  if (duration <= 4.0 * width) {
    t_surface_l = nxnodes;
    t_surface_r = 0;
  }

  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++) {
        if ((double)nsum[ixnode][iynode][iznode] > crit_num) {
          Activated[ixnode][iynode][iznode] = 1;
          // Surface Movement
          if (duration <= 4.0 * width && ixnode < t_surface_l) t_surface_l = ixnode;
          if (duration <= 4.0 * width && ixnode >= t_surface_r) t_surface_r = ixnode + 1;
          if (bulk_ttm == 1 && Tes[iynode][iznode] < ixnode) Tes[iynode][iznode] = ixnode;
        }
        else {
          Activated[ixnode][iynode][iznode] = 0;
        }
      }

  // Assume smooth surface during laser irradiation
  if (duration <= 4.0 * width) {
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++) {
        if (bulk_ttm == 1) Tes[iynode][iznode] = t_surface_r - 1;
        for (int ixnode = t_surface_l; ixnode < t_surface_r; ixnode++)
          Activated[ixnode][iynode][iznode] = 1;
      }
  }

}

void FixFEMTO3D::update_parameters()
{
  int mxnodei, mxnodef, msize, mdisp, ixnode;
  int rsize[numP], disp[numP];
  mxnodei = (pid * nxnodes) / numP;
  mxnodef = ((pid + 1) * nxnodes) / numP - 1;
  mdisp = mxnodei * nynodes * nznodes;
  msize = (mxnodef - mxnodei + 1) * nynodes * nznodes;
  MPI_Allgather(&msize, 1, MPI_INT, rsize, 1, MPI_INT, world);
  MPI_Allgather(&mdisp, 1, MPI_INT, disp, 1, MPI_INT, world);
  double T_e, T_i, Ni;
  for (int i = 0; i < mxnodef - mxnodei + 1; i++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++) {
        ixnode = i + mxnodei;
        if (Activated[ixnode][iynode][iznode] == 1) {
          T_e = T_electron[ixnode][iynode][iznode];
          T_i = T_a[ixnode][iynode][iznode];
          Ni = (double)nsum[ixnode][iynode][iznode] / del_vol;
          CeT[ixnode][iynode][iznode] = MyCe(T_e, Ni, T_i);
          ke_real[ixnode][iynode][iznode] = MyKet(T_e, Ni, T_i);
          GT[ixnode][iynode][iznode] = MyG(T_e, Ni, T_i);
        }
        else {
          CeT[ixnode][iynode][iznode] = 0.0;
          GT[ixnode][iynode][iznode] = 0.0;
          ke_real[ixnode][iynode][iznode] = 0.0;
        }
      }
  //MPI_Allgatherv(&ke_real_buffer[0][0][0],msize,MPI_DOUBLE,&ke_real[0][0][0],rsize,disp,MPI_DOUBLE,world);
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &ke_real[0][0][0], rsize, disp, MPI_DOUBLE, world);
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &CeT[0][0][0], rsize, disp, MPI_DOUBLE, world);
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &GT[0][0][0], rsize, disp, MPI_DOUBLE, world);
}

void FixFEMTO3D::update_parameters_bulk()
{
  int mxnodei, mxnodef, msize, mdisp, ixnode;
  int rsize[numP], disp[numP];
  mxnodei = (pid * bxsize) / numP;
  mxnodef = ((pid + 1) * bxsize) / numP - 1;
  mdisp = mxnodei * nynodes * nznodes; // nxnodes
  msize = (mxnodef - mxnodei + 1) * nynodes * nznodes;
  MPI_Allgather(&msize, 1, MPI_INT, rsize, 1, MPI_INT, world);
  MPI_Allgather(&mdisp, 1, MPI_INT, disp, 1, MPI_INT, world);
  double T_e, T_i, Ni;
  for (int i = 0; i < mxnodef - mxnodei + 1; i++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++) {
        ixnode = i + mxnodei;
        T_e = Te_bulk[ixnode][iynode][iznode];
        T_i = Ta_bulk[ixnode][iynode][iznode];
        Ni = ionic_density;
        CeT_bulk[ixnode][iynode][iznode] = MyCe(T_e, Ni, T_i);
        ke_real_bulk[ixnode][iynode][iznode] = MyKet(T_e, Ni, T_i);
        GT_bulk[ixnode][iynode][iznode] = MyG(T_e, Ni, T_i);
        ki_bulk[ixnode][iynode][iznode] = ki_function(T_i);
      }
  //MPI_Allgatherv(&ke_real_buffer[0][0][0],msize,MPI_DOUBLE,&ke_real[0][0][0],rsize,disp,MPI_DOUBLE,world);
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &ke_real_bulk[0][0][0], rsize, disp, MPI_DOUBLE, world);
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &CeT_bulk[0][0][0], rsize, disp, MPI_DOUBLE, world);
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &GT_bulk[0][0][0], rsize, disp, MPI_DOUBLE, world);
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &ki_bulk[0][0][0], rsize, disp, MPI_DOUBLE, world);
}

double FixFEMTO3D::ki_function(double T_i)
{
  // dummy output, the ion conductivity is neglegible compared to electrons'
  double ki = 0.0 * T_i;
  return ki;
}


/* ----------------------------------------------------------------------
   memory usage of 3d grid
------------------------------------------------------------------------- */

double FixFEMTO3D::memory_usage()
{
  double bytes = 0.0;

  if (premode == 0) {
    bytes += 2 * total_nnodes * sizeof(int); // Activated, nsum
    bytes += 15 * total_nnodes * sizeof(double); // sum_mass_vsq, sum_mass_v, sum_mass, average_v, T_electron_first, T_electron_old, T_electron, T_a, net_energy_transfer, energy_conduction, ke_real, CeT, GT, mult_factor, skin_layer
    bytes += 6 * 2 * 1000 * sizeof(double); // parameter tables, size = 1000, could be different, ZTe, CeTe, GTe, KeTe, ReflecTe, PenTe
    bytes += 2 * nxnodes * sizeof(double); // Kenergy, Genergy for the calculation of energy transfer
    if (bulk_ttm == 1) {
      bytes += 9 * bxsize * nynodes * nznodes * sizeof(double); // Te_bulk, Ta_bulk, CeT_bulk, ke_real_bulk, GT_bulk, ki_bulk, mul_factor_bulk, skin_layer_bulk, E_melt_buffer
      bytes += nynodes * nznodes * sizeof(double); // xmax
      bytes += nynodes * nznodes * sizeof(int); // Tes
    }
  }

  if (premode == 1) {
    bytes += 1 * total_nnodes * sizeof(int); // nsum
    bytes += 2 * total_nnodes * sizeof(double); // sum_mass_vsq, T_a
  }

  return bytes;
}

/* ---------------------------------------------------------------------- */

void FixFEMTO3D::grow_arrays(int ngrow)
{
  flangevin = NULL;
  memory->grow(flangevin, ngrow, 3, "femto3D:flangevin");
  // zero out the flangevin array
  for (int i = 0; i < atom->nmax; i++) {
    flangevin[i][0] = 0.0;
    flangevin[i][1] = 0.0;
    flangevin[i][2] = 0.0;
  }
}


/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixFEMTO3D::write_restart(FILE* fpr)
{
  double* rlist;
  if (bulk_ttm == 0) {
    memory->create(rlist, 2 + 2 * nxnodes * nynodes * nznodes, "femto3D:rlist");
  }
  else {
    memory->create(rlist, 2 + 2 * nxnodes * nynodes * nznodes + 2 * bxsize * nynodes * nznodes, "femto3D:rlist");
  }
  int n = 0;
  rlist[n++] = static_cast<double> (seed);
  rlist[n++] = duration;

  if (premode == 0) {
    for (int ixnode = 0; ixnode < nxnodes; ixnode++)
      for (int iynode = 0; iynode < nynodes; iynode++)
        for (int iznode = 0; iznode < nznodes; iznode++)
          rlist[n++] = T_electron[ixnode][iynode][iznode];
    for (int ixnode = 0; ixnode < nxnodes; ixnode++)
      for (int iynode = 0; iynode < nynodes; iynode++)
        for (int iznode = 0; iznode < nznodes; iznode++)
          rlist[n++] = static_cast<double>(Activated[ixnode][iynode][iznode]);
    if (bulk_ttm == 1) {
      for (int ixnode = 0; ixnode < bxsize; ixnode++)
        for (int iynode = 0; iynode < nynodes; iynode++)
          for (int iznode = 0; iznode < nznodes; iznode++)
            rlist[n++] = Te_bulk[ixnode][iynode][iznode];
      for (int ixnode = 0; ixnode < bxsize; ixnode++)
        for (int iynode = 0; iynode < nynodes; iynode++)
          for (int iznode = 0; iznode < nznodes; iznode++)
            rlist[n++] = Ta_bulk[ixnode][iynode][iznode];
    }
  }

  if (premode == 1) {
    for (int ixnode = 0; ixnode < nxnodes; ixnode++) {
      int TempN;
      if (ixnode < surface_l)
        TempN = 0;
      else if (ixnode >= surface_r)
        TempN = surface_r - surface_l - 1;
      else
        TempN = ixnode - surface_l;
      for (int iynode = 0; iynode < nynodes; iynode++)
        for (int iznode = 0; iznode < nznodes; iznode++)
          rlist[n++] = ITemp[TempN];
    }
    for (int ixnode = 0; ixnode < nxnodes; ixnode++)
      for (int iynode = 0; iynode < nynodes; iynode++)
        for (int iznode = 0; iznode < nznodes; iznode++) {
          double crit_num = crit_num_f * ionic_density * del_vol;
          if ((double)nsum[ixnode][iynode][iznode] > crit_num)
            rlist[n++] = static_cast<double>(1);
          else
            rlist[n++] = static_cast<double>(0);
        }
    if (bulk_ttm == 1) {
      for (int ixnode = 0; ixnode < bxsize; ixnode++)
        for (int iynode = 0; iynode < nynodes; iynode++)
          for (int iznode = 0; iznode < nznodes; iznode++)
            rlist[n++] = ITempBulk[ixnode];
      for (int ixnode = 0; ixnode < bxsize; ixnode++)
        for (int iynode = 0; iynode < nynodes; iynode++)
          for (int iznode = 0; iznode < nznodes; iznode++)
            rlist[n++] = ITempBulk[ixnode];
    }
  }

  if (comm->me == 0) {

    int size = n * sizeof(double);
    fwrite(&size, sizeof(int), 1, fpr);
    fwrite(rlist, sizeof(double), n, fpr);
  }
  memory->destroy(rlist);
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixFEMTO3D::restart(char* buf)
{
  int n = 0;
  double* rlist = (double*)buf;
  // the seed must be changed from the initial seed
  seed = static_cast<int> (0.5 * rlist[n++]);
  duration = rlist[n++];

  if (premode == 0) {
    for (int ixnode = 0; ixnode < nxnodes; ixnode++)
      for (int iynode = 0; iynode < nynodes; iynode++)
        for (int iznode = 0; iznode < nznodes; iznode++)
          T_electron[ixnode][iynode][iznode] = rlist[n++];
    for (int ixnode = 0; ixnode < nxnodes; ixnode++)
      for (int iynode = 0; iynode < nynodes; iynode++)
        for (int iznode = 0; iznode < nznodes; iznode++)
          Activated[ixnode][iynode][iznode] = static_cast<int>(rlist[n++]);
    if (bulk_ttm == 1) {
      for (int ixnode = 0; ixnode < bxsize; ixnode++)
        for (int iynode = 0; iynode < nynodes; iynode++)
          for (int iznode = 0; iznode < nznodes; iznode++)
            Te_bulk[ixnode][iynode][iznode] = rlist[n++];
      for (int ixnode = 0; ixnode < bxsize; ixnode++)
        for (int iynode = 0; iynode < nynodes; iynode++)
          for (int iznode = 0; iznode < nznodes; iznode++)
            Ta_bulk[ixnode][iynode][iznode] = rlist[n++];
    }
  }
  else if (premode == 1) {
    // It is empty here
  }
  delete random;
  random = new RanMars(lmp, seed + comm->me);
}
/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixFEMTO3D::pack_restart(int i, double* buf)
{
  buf[0] = 4;
  buf[1] = flangevin[i][0];
  buf[2] = flangevin[i][1];
  buf[3] = flangevin[i][2];
  return 4;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixFEMTO3D::unpack_restart(int nlocal, int nth)
{
  double** extra = atom->extra;
  // skip to Nth set of extra values
  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int>(extra[nlocal][m]);
  m++;
  flangevin[nlocal][0] = extra[nlocal][m++];
  flangevin[nlocal][1] = extra[nlocal][m++];
  flangevin[nlocal][2] = extra[nlocal][m++];
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixFEMTO3D::maxsize_restart()
{
  return 4;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixFEMTO3D::size_restart(int nlocal)
{
  return 4;
}

/* ----------------------------------------------------------------------
   read Z, G, C, K from a user-specified file called by all procs
------------------------------------------------------------------------- */

int FixFEMTO3D::read_data_table(const std::string& filename, double*** iTe)
{
  PotentialFileReader reader(lmp, filename, "ttm/femto3D data table");

  int size = reader.next_values(1).next_int();
  int count = 0;
  double** xTe;
  memory->create(xTe, size, 2, "femto3D:read_data_table");
  while (count < size) {
    auto line = reader.next_values(2);
    if (!line.has_next())
      break;
    xTe[count][0] = line.next_double();
    xTe[count][1] = line.next_double();
    if (xTe[count][0] < 0.0 || xTe[count][1] < 0.0) error->one(FLERR, "Fix femto3D electron temperatures or read_data_table must be >= 0.0");
    count++;
  }
  if (size != count) error->one(FLERR, "Fix femto3D read_data_table has size problems, register size is greater than actual size");
  *iTe = xTe;
  return count;
}

int FixFEMTO3D::read_temperature_table(const std::string& filename, double** iTe)
{
  PotentialFileReader reader(lmp, filename, "ttm/femto3D temperature table");

  int size = reader.next_values(1).next_int();

  int count = 0;
  double* ITemp;
  memory->create(ITemp, size, "femto3D:read_temperature_table");
  while (count < size) {
    ITemp[count] = reader.next_values(1).next_double();
    if (ITemp[count] < 0.0) error->one(FLERR, "Fix femto3D input temperatures must be >= 0.0");
    count++;
  }
  if (size != count) error->one(FLERR, "Fix femto3D read_temperature_table has size problems, register size is greater than actual size");
  *iTe = ITemp;
  return count;
}

/* ----------------------------------------------------------------------
   read tables from a 2D array
------------------------------------------------------------------------- */

double FixFEMTO3D::interpolation(double** table, int size, double A, double min, double max)
{
  double B;
  int i;
  for (int j = 0; j < size; j++) {
    i = j;
    if (table[j][0] >= A) break;
  }

  if (i == 0)
    B = table[i][1];
  else
    B = A * (table[i][1] - table[i - 1][1]) / (table[i][0] - table[i - 1][0]) +
    (table[i][0] * table[i - 1][1] - table[i][1] * table[i - 1][0]) / (table[i][0] - table[i - 1][0]);

  if (B < min) B = min;
  if (B > max) B = max;

  return B;
}

/* ----------------------------------------------------------------------
   read parameters from a user-specified file
   only called by all procs
------------------------------------------------------------------------- */

void FixFEMTO3D::read_parameter(const std::string& filename)
{

  if (comm->me == 0) {

    try {
      PotentialFileReader reader(lmp, filename, "ttm/femto3D parameter");

      // Preset temperature mode (Activated = 1; Inactivated = 0)
      premode = reader.next_values(1).next_int();
      // printf("premode: %d\n", premode);

      if (premode == 0) {
        // Number of thermal solve grid points in the x, y, z directions
        auto values = reader.next_values(3);
        nxnodes = values.next_int();
        nynodes = values.next_int();
        nznodes = values.next_int();

        // coordinate of 1st and 2nd surface in x-direction (in box units) - constant
        auto values2 = reader.next_values(2);
        surface_l = values2.next_int();
        surface_r = values2.next_int();

        // Initial material temperature
        T_init = reader.next_values(1).next_double();

        // average intensity of pulse (source of energy) (metal units)
        intensity = reader.next_values(1).next_double();

        // width of pulse (picoseconds)
        width = reader.next_values(1).next_double();

        // laser wavelength (A)
        wavelength_laser = reader.next_values(1).next_double();

        // factor of electronic pressure (PF) Pe = PF*Ce*Te
        pres_factor = reader.next_values(1).next_double();

        // effective free path of electrons (angstrom)
        free_path = reader.next_values(1).next_double();

        // ionic density (ions*angstrom^{-3})
        ionic_density = reader.next_values(1).next_double();

        // NewPulse, 1 for activate, 0 for inactivate
        NewPulse = reader.next_values(1).next_int();

        // bulk_ttm, 1 for activate, 0 for inactivate
        bulk_ttm = reader.next_values(1).next_int();


        if (bulk_ttm == 1) {
          // Atomic mass of the bulk target(grams/mole)
          massT = reader.next_values(1).next_double();

          // Target lattice structure (1 for bcc, 2 for fcc)
          lstr = reader.next_values(1).next_int();

          // Target melting point
          T_melt = reader.next_values(1).next_double();

          // Target latent heat of fusion, in eV/mol
          Latent_melt = reader.next_values(1).next_double();

          // Target Boundary thickness
          bound_thick = reader.next_values(1).next_double();

          // Static Stress at the boundary
          F_0 = reader.next_values(1).next_double();

          // Atomic cross-section
          A_cross = reader.next_values(1).next_double();

          // Speed of sound in the target
          v_s = reader.next_values(1).next_double();

          // bulk thickness
          bulk_thick = reader.next_values(1).next_double();

          // In the bulk, number of grid points
          bxsize = reader.next_values(1).next_int();
        }
      }
      else if (premode == 1) {
        // Number of thermal solve grid points in the x, y, z directions
        auto values3 = reader.next_values(3);
        nxnodes = values3.next_int();
        nynodes = values3.next_int();
        nznodes = values3.next_int();

        // printf("nodes: %d %d %d\n", nxnodes, nynodes, nznodes);

        // coordinate of 1st and 2nd surface in x-direction (in box units) - constant
        auto values4 = reader.next_values(2);
        surface_l = values4.next_int();
        surface_r = values4.next_int();

        // ionic density (ions*angstrom^{-3})
        ionic_density = reader.next_values(1).next_double();

        // GStrength, how strong the electron lattice coupling is
        GStrength = reader.next_values(1).next_double();

        // bulk_ttm, 1 for activate, 0 for inactivate
        bulk_ttm = reader.next_values(1).next_int();

        if (bulk_ttm == 1) {
          // Target Boundary thickness
          bound_thick = reader.next_values(1).next_double();

          // In the bulk, number of grid points
          bxsize = reader.next_values(1).next_int();
        }
      }
    }
    catch (std::exception& e) {
      error->one(FLERR, e.what());
    }
  }

  MPI_Bcast(&premode, 1, MPI_INT, 0, world);

  if (premode == 0) {
    MPI_Bcast(&nxnodes, 1, MPI_INT, 0, world);
    MPI_Bcast(&nynodes, 1, MPI_INT, 0, world);
    MPI_Bcast(&nznodes, 1, MPI_INT, 0, world);
    MPI_Bcast(&surface_l, 1, MPI_INT, 0, world);
    MPI_Bcast(&surface_r, 1, MPI_INT, 0, world);

    MPI_Bcast(&T_init, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&intensity, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&width, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&wavelength_laser, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&pres_factor, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&free_path, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&ionic_density, 1, MPI_DOUBLE, 0, world);

    MPI_Bcast(&NewPulse, 1, MPI_INT, 0, world);
    MPI_Bcast(&bulk_ttm, 1, MPI_INT, 0, world);

    if (bulk_ttm == 1) {
      MPI_Bcast(&massT, 1, MPI_DOUBLE, 0, world);
      MPI_Bcast(&lstr, 1, MPI_INT, 0, world);
      MPI_Bcast(&T_melt, 1, MPI_DOUBLE, 0, world);
      MPI_Bcast(&Latent_melt, 1, MPI_DOUBLE, 0, world);
      MPI_Bcast(&bound_thick, 1, MPI_DOUBLE, 0, world);
      MPI_Bcast(&F_0, 1, MPI_DOUBLE, 0, world);
      MPI_Bcast(&A_cross, 1, MPI_DOUBLE, 0, world);
      MPI_Bcast(&v_s, 1, MPI_DOUBLE, 0, world);
      MPI_Bcast(&bulk_thick, 1, MPI_DOUBLE, 0, world);
      MPI_Bcast(&bxsize, 1, MPI_INT, 0, world);
    }
  }
  else if (premode == 1) {

    MPI_Bcast(&nxnodes, 1, MPI_INT, 0, world);
    MPI_Bcast(&nynodes, 1, MPI_INT, 0, world);
    MPI_Bcast(&nznodes, 1, MPI_INT, 0, world);

    MPI_Bcast(&surface_l, 1, MPI_INT, 0, world);
    MPI_Bcast(&surface_r, 1, MPI_INT, 0, world);

    MPI_Bcast(&ionic_density, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&GStrength, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&bulk_ttm, 1, MPI_INT, 0, world);

    if (bulk_ttm == 1) {
      MPI_Bcast(&bound_thick, 1, MPI_DOUBLE, 0, world);
      MPI_Bcast(&bxsize, 1, MPI_INT, 0, world);
    }
  }
}

void FixFEMTO3D::read_tablelist(const std::string& f_tablelist)
{
  if (comm->me == 0) {

    try {
      PotentialFileReader reader(lmp, f_tablelist, "ttm/femto3D table list");

      if (premode == 0) {

        fp_Z = reader.next_values(1).next_string();
        fp_Ce = reader.next_values(1).next_string();
        fp_G = reader.next_values(1).next_string();
        fp_Ke = reader.next_values(1).next_string();
        fp_Reflec = reader.next_values(1).next_string();
        fp_Pen = reader.next_values(1).next_string();

        Zsize = read_data_table(fp_Z, &ZTe);
        Cesize = read_data_table(fp_Ce, &CeTe);
        Gsize = read_data_table(fp_G, &GTe);
        Kesize = read_data_table(fp_Ke, &KeTe);
        Reflecsize = read_data_table(fp_Reflec, &ReflecTe);
        Pensize = read_data_table(fp_Pen, &PenTe);
      }
      else if (premode == 1) {

        fp_Temp = reader.next_string();
        fp_Ta_out = reader.next_string();
        Tsize = read_temperature_table(fp_Temp, &ITemp);
        if (Tsize != surface_r - surface_l) error->all(FLERR, "Tsize does not match");
        // printf("%s\n", fp_Temp.c_str());

        if (bulk_ttm == 1) {
          fp_TempBulk = reader.next_string();
          Tsize_bulk = read_temperature_table(fp_TempBulk, &ITempBulk);
          if (Tsize_bulk != bxsize) error->all(FLERR, "Tsize_bulk does not match");
        }
      }
    }
    catch (std::exception& e) {
      error->one(FLERR, e.what());
    }
  }

  if (premode == 0) {
    MPI_Bcast(&Zsize, 1, MPI_INT, 0, world);
    MPI_Bcast(&Cesize, 1, MPI_INT, 0, world);
    MPI_Bcast(&Gsize, 1, MPI_INT, 0, world);
    MPI_Bcast(&Kesize, 1, MPI_INT, 0, world);
    MPI_Bcast(&Reflecsize, 1, MPI_INT, 0, world);
    MPI_Bcast(&Pensize, 1, MPI_INT, 0, world);

    MPI_Bcast(&ZTe[0][0], 2 * Zsize, MPI_DOUBLE, 0, world);
    MPI_Bcast(&CeTe[0][0], 2 * Cesize, MPI_DOUBLE, 0, world);
    MPI_Bcast(&GTe[0][0], 2 * Gsize, MPI_DOUBLE, 0, world);
    MPI_Bcast(&KeTe[0][0], 2 * Kesize, MPI_DOUBLE, 0, world);
    MPI_Bcast(&ReflecTe[0][0], 2 * Reflecsize, MPI_DOUBLE, 0, world);
    MPI_Bcast(&PenTe[0][0], 2 * Pensize, MPI_DOUBLE, 0, world);
  }
  else if (premode == 1) {
    MPI_Bcast(&Tsize, 1, MPI_INT, 0, world);
    MPI_Bcast(ITemp, Tsize, MPI_DOUBLE, 0, world);
    if (bulk_ttm == 1) {
      MPI_Bcast(&Tsize_bulk, 1, MPI_INT, 0, world);
      MPI_Bcast(ITempBulk, Tsize_bulk, MPI_DOUBLE, 0, world);
    }
  }
}

void FixFEMTO3D::read_outlist(const std::string& fp_outlist, int nowdur)
{
  if (comm->me == 0) {
    try {
      PotentialFileReader reader(lmp, fp_outlist, "ttm/femto3D output list");
      FILE * fp;
      bool flag = (nowdur <= 100); // create new output files if restart before 0.1 ps
      if (nfileevery > 0) {
        fp_Ta_out = reader.next_string();
        fp_Te_out = reader.next_string();
        fp_Et_out = reader.next_string();
        fp_laser_out = reader.next_string();
        fp_Ta_early_out = reader.next_string();
        fp_Te_early_out = reader.next_string();
        if (bulk_ttm == 1) {
          fp_Ta_bulk = reader.next_string();
          fp_Te_bulk = reader.next_string();
        }

        if (flag) {
          fp = fopen(fp_Ta_out.c_str(), "w");
          fclose(fp);
          fp = fopen(fp_Te_out.c_str(), "w");
          fclose(fp);
          fp = fopen(fp_Et_out.c_str(), "w");
          fclose(fp);
          fp = fopen(fp_laser_out.c_str(), "w");
          fclose(fp);
          fp = fopen(fp_Ta_early_out.c_str(), "w");
          fclose(fp);
          fp = fopen(fp_Te_early_out.c_str(), "w");

          if (bulk_ttm == 1) {
            fp = fopen(fp_Ta_bulk.c_str(), "w");
            fclose(fp);
            fp = fopen(fp_Te_bulk.c_str(), "w");
            fclose(fp);
          }

        }
      }
    } catch (std::exception &e) {
      error->one(FLERR,e.what());
    }
  }
}

void FixFEMTO3D::Tempout() {
  if (comm->me == 0) {
    FILE * fp;
    // output nodal temperatures for current timestep
    if ((nfileevery) && !(update->ntimestep % nfileevery)) {
      fp = fopen(fp_Ta_out.c_str(), "a");
      fprintf(fp, BIGINT_FORMAT, update->ntimestep);
      fprintf(fp, "\n--------------------------------------------------------------------");
      for (int ixnode = 0; ixnode < nxnodes; ixnode++) {
        int num = 0;
        double ptemp = 0.0;
        for (int iynode = 0; iynode < nynodes; iynode++)
          for (int iznode = 0; iznode < nznodes; iznode++) {
            if (Activated[ixnode][iynode][iznode] == 1) {
              num++;
              ptemp += T_a[ixnode][iynode][iznode];
            }
          }
        if (num != 0) {
          ptemp /= num;
          fprintf(fp, "\n%d\t%f", ixnode, ptemp);
        }
      }
      fprintf(fp, "\n\n");
      fclose(fp);

      fp = fopen(fp_Te_out.c_str(), "a");
      fprintf(fp, BIGINT_FORMAT, update->ntimestep);
      fprintf(fp, "\n--------------------------------------------------------------------");
      for (int ixnode = 0; ixnode < nxnodes; ixnode++) {
        int num = 0;
        double ptemp = 0.0;
        for (int iynode = 0; iynode < nynodes; iynode++)
          for (int iznode = 0; iznode < nznodes; iznode++) {
            if (Activated[ixnode][iynode][iznode] == 1) {
              num++;
              // ptemp += T_electron_first[ixnode][iynode][iznode];
              ptemp += T_electron[ixnode][iynode][iznode];
            }
          }
        if (num != 0) {
          ptemp /= num;
          fprintf(fp, "\n%d\t%f", ixnode, ptemp);
        }
      }
      fprintf(fp, "\n\n");
      fclose(fp);

      if (bulk_ttm == 1) {

        double Ee2Ea, dTe_bulk = 0.0, dTa_bulk = 0.0;
        double ci = 3 * ionic_density * force->boltz;

        fp = fopen(fp_Ta_bulk.c_str(), "a");
        fprintf(fp, BIGINT_FORMAT, update->ntimestep);
        fprintf(fp, "\n--------------------------------------------------------------------");
        for (int ixnode = 0; ixnode < bxsize; ixnode++) {
          double ptemp = 0.0;
          for (int iynode = 0; iynode < nynodes; iynode++)
            for (int iznode = 0; iznode < nznodes; iznode++) {
              // Ee2Ea = GT_bulk[ixnode][iynode][iznode]*(Te_bulk[ixnode][iynode][iznode] - Ta_bulk[ixnode][iynode][iznode]);
              // dTa_bulk = +update->dt/ci * Ee2Ea;
              ptemp += Ta_bulk[ixnode][iynode][iznode] + dTa_bulk;
            }
          ptemp /= nynodes * nznodes;
          fprintf(fp, "\n%d\t%f", ixnode, ptemp);
        }
        fprintf(fp, "\n\n");
        fclose(fp);

        fp = fopen(fp_Te_bulk.c_str(), "a");
        fprintf(fp, BIGINT_FORMAT, update->ntimestep);
        fprintf(fp, "\n--------------------------------------------------------------------");
        for (int ixnode = 0; ixnode < bxsize; ixnode++) {
          double ptemp = 0.0;
          for (int iynode = 0; iynode < nynodes; iynode++)
            for (int iznode = 0; iznode < nznodes; iznode++) {
              // Ee2Ea = GT_bulk[ixnode][iynode][iznode]*(Te_bulk[ixnode][iynode][iznode] - Ta_bulk[ixnode][iynode][iznode]);
              // dTe_bulk = -update->dt/CeT_bulk[ixnode][iynode][iznode] * Ee2Ea;
              ptemp += Te_bulk[ixnode][iynode][iznode] + dTe_bulk;
            }
          ptemp /= nynodes * nznodes;
          fprintf(fp, "\n%d\t%f", ixnode, ptemp);
        }
        fprintf(fp, "\n\n");
        fclose(fp);

      }
    }
  }
}

void FixFEMTO3D::Otherout() {

  if (comm->me == 0) {
    FILE * fp;

    for (int ixnode = 0; ixnode < nxnodes; ixnode++) {
      double Gtemp = 0.0;
      double Ktemp = 0.0;
      for (int iynode = 0; iynode < nynodes; iynode++)
        for (int iznode = 0; iznode < nznodes; iznode++) {
          if (Activated[ixnode][iynode][iznode] == 1) {
            Gtemp -= net_energy_transfer[ixnode][iynode][iznode] / del_vol;
            Ktemp += energy_conduction[ixnode][iynode][iznode];
          }
        }
      Genergy[ixnode] -= Gtemp / nynodes / nznodes;
      Kenergy[ixnode] += Ktemp / nynodes / nznodes;
    }

    if ((nfileevery) && !(update->ntimestep % nfileevery)) {
      fp = fopen(fp_Et_out.c_str(), "a");
      fprintf(fp, BIGINT_FORMAT, update->ntimestep);
      fprintf(fp, "\tG\tK\n--------------------------------------------------------------------");
      double GTotal = 0.0;
      double KTotal = 0.0;
      for (int ixnode = 0; ixnode < nxnodes; ixnode++) {
        if ((abs(Genergy[ixnode]) > 1.0e-5) || (abs(Kenergy[ixnode]) > 1.0e-5)) {
          fprintf(fp, "\n%d\t%f\t%f", ixnode, Genergy[ixnode] / nfileevery, Kenergy[ixnode] / nfileevery); // Engergy per volume per time
        }
        GTotal += Genergy[ixnode] * del_vol * nynodes * nznodes * update->dt;
        KTotal += Kenergy[ixnode] * del_vol * nynodes * nznodes * update->dt;
      }
      fprintf(fp, "\nTotal:\t%f\t%f", GTotal, KTotal);
      for (int ixnode = 0; ixnode < nxnodes; ixnode++) {
        // Reset Genergy and Kenergy every nfileevery
        Genergy[ixnode] = 0.0;
        Kenergy[ixnode] = 0.0;
      }
      fprintf(fp, "\n\n");
      fclose(fp);
    }

    if (duration <= 4.0 * width + 0.0005) {
      fp = fopen(fp_laser_out.c_str(), "a");
      fprintf(fp, BIGINT_FORMAT, update->ntimestep); //BIGINT_FORMAT, update->ntimestep, static_cast<int> (1000.0*duration))
      fprintf(fp, "\t%f", reflectivity);
      fprintf(fp, "\n--------------------------------------------------------------------");
      for (int ixnode = 0; ixnode < nxnodes; ixnode++) {
        int num = 0;
        double dtemp = 0.0;
        double Qtemp = 0.0;
        for (int iynode = 0; iynode < nynodes; iynode++)
          for (int iznode = 0; iznode < nznodes; iznode++) {
            if (Activated[ixnode][iynode][iznode] == 1) {
              num++;
              dtemp += skin_layer[ixnode][iynode][iznode];
              Qtemp += mult_factor[ixnode][iynode][iznode];
            }
          }
        if (num != 0) {
          dtemp /= num;
          Qtemp /= num;
          fprintf(fp, "\n%d\t%f\t%f", ixnode, dtemp, Qtemp);
        }
      }
      fprintf(fp, "\n\n");
      fclose(fp);
    }

    if (duration <= 1.0) { // 1 picosecond
      fp = fopen(fp_Te_early_out.c_str(), "a");
      fprintf(fp, BIGINT_FORMAT, update->ntimestep);
      fprintf(fp, "\n--------------------------------------------------------------------");
      for (int ixnode = 0; ixnode < nxnodes; ixnode++) {
        int num = 0;
        double ptemp = 0.0;
        for (int iynode = 0; iynode < nynodes; iynode++)
          for (int iznode = 0; iznode < nznodes; iznode++) {
            if (Activated[ixnode][iynode][iznode] == 1) {
              num++;
              ptemp += T_electron[ixnode][iynode][iznode];
            }
          }
        if (num != 0) {
          ptemp /= num;
          fprintf(fp, "\n%d\t%f", ixnode, ptemp);
        }
      }
      fprintf(fp, "\n\n");
      fclose(fp);

      fp = fopen(fp_Ta_early_out.c_str(), "a");
      fprintf(fp, BIGINT_FORMAT, update->ntimestep);
      fprintf(fp, "\n--------------------------------------------------------------------");
      for (int ixnode = 0; ixnode < nxnodes; ixnode++) {
        int num = 0;
        double ptemp = 0.0;
        for (int iynode = 0; iynode < nynodes; iynode++)
          for (int iznode = 0; iznode < nznodes; iznode++) {
            if (Activated[ixnode][iynode][iznode] == 1) {
              num++;
              ptemp += T_a[ixnode][iynode][iznode];
            }
          }
        if (num != 0) {
          ptemp /= num;
          fprintf(fp, "\n%d\t%f", ixnode, ptemp);
        }
      }
      fprintf(fp, "\n\n");
      fclose(fp);
    }
  }
}
