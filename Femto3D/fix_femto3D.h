/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(femto3D,FixFEMTO3D)

#else

#ifndef LMP_FIX_FEMTO3D_H
#define LMP_FIX_FEMTO3D_H

#include "fix.h"

namespace LAMMPS_NS {

class FixFEMTO3D : public Fix {
  public:
    FixFEMTO3D(class LAMMPS *, int, char **);
    ~FixFEMTO3D();
    int setmask();
    void init();
    void setup(int);
    void post_force(int);
    void post_force_respa(int, int, int);
    void post_force_setup(int);
    void post_force_respa_setup(int, int, int);
    void end_of_step();
    void reset_dt();
    void write_restart(FILE *);
    void restart(char *);
    int pack_restart(int, double *);
    void unpack_restart(int, int);
    int size_restart(int);
    int maxsize_restart();
    double memory_usage();
    void grow_arrays(int);

  private:

    // Constant
    double epsilon0, speed_light, Pi, e_charge, m_e, k_b, h_bar, Na;
    // Parameters Calculation
    double ***CeT, ***ke_real, ***GT;
    double ***CeT_bulk, ***ke_real_bulk, ***GT_bulk, ***ki_bulk;
    double dx, dy, dz, dx_bulk, del_vol;

    bool hasrun=false;
    int pid, numP;
    int nfileevery;
    int nlevels_respa;
    int seed, writenn, premode, NewPulse;
    int Cesize, Gsize, Kesize, Reflecsize, Pensize, Zsize;
    class RanMars *random;
    FILE *fp_Te_out, *fp_Ta_out, *fp_Et_out, *fp_Tidiff_out, *fp_laser_out, *fp_Te_bulk, *fp_Ta_bulk, *fp_Z, *fp_Ce, *fp_G, *fp_Ke, *fp_Pen, *fp_Reflec, *fp_Te_in, *fp_Te_early_out, *fp_Ta_early_out, *fp_parameter, *fp_tablelist, *fp_outlist;
    int nxnodes,nynodes,nznodes,total_nnodes;
    int ***Activated, ***nsum;
    double *Kenergy, *Genergy;
    double **CeTe, **GTe, **KeTe, **ReflecTe, **PenTe, **ZTe;
    double *ratio;
    double **flangevin;
    double **x_max;
    int **Tes;
    double ***T_electron,***T_electron_old,***T_electron_first,***T_a,***Ta_bulk,***Te_bulk;
    double ***net_energy_transfer,***energy_conduction,***skin_layer,***mult_factor,***skin_layer_bulk,***mult_factor_bulk,***E_melt_buffer,***sum_mass_vsq, ***sum_mass_v, ***sum_mass, ***average_v;
    int surface_l,surface_r,t_surface_l,t_surface_r;
    double electronic_density, reflectivity;
    double intensity,width,duration,wavelength_laser,surface_double,massT,massone;
    double ttm_dt;
    double pres_factor,free_path,ionic_density,crit_num_f;
    double T_init,T_melt,Latent_melt,xmin;
    double F_0, A_cross, v_s, bound_thick, bulk_thick;
    int bulk_ttm, bxsize, lstr;

    double GStrength, *ITemp, *ITempBulk;
    FILE *fp_Temp, *fp_TempBulk;
    int Tsize, Tsize_bulk;

    void set_initial_temperatures();
    void read_parameter(FILE *);
    void read_tablelist(FILE *);
    void read_outlist(FILE *, int);
    int read_in_table(FILE *, double ***);
    int read_in_table2(FILE *, double **);
    double read_tables(double **, int, double);
    void Tempout();
    void Otherout();
    void laser(double);
    void laser_bulk();
    void update_Ta();
    void ChangeType();
    void update_parameters();
    void update_parameters_bulk();
    double ki_function(double);
    double MyCe(double,double,double);
    double MyKet(double,double,double);
    double MyG(double,double,double);
};

}

#endif
#endif
