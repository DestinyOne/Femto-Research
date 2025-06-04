# Femto-Research
Codes for femtosecond laser interaction research. Note Femto3D works in LAMMPS with version (29 Aug 2024).

If you have any questions, please email yuanweirong78@gmail.com.

## `Femto3D`
This is the modified LAMMPS function based on its two-temperature function. Put `fix_femto3D.cpp` and `fix_femto3D.h` in the same folder as `fix_ttm_mod.cpp` (`src/EXTRA-FIX`) and compile in the same way. Put `femto3D.cmke` in the `cmake\presets` folder inside LAMMPS folder. `installation help.txt` is a file useful for the installation Reference: https://docs.lammps.org/fix_ttm.html#fix-ttm-mod-command

## `Preheating`
A target generated in MD is preheated to a wanted temperature and then relaxes to equilibrium. The output file for further use is a LAMMPS restart file for further simulation. The default output from the sample input code is `restart.equil`. Example: `mpiexec -np 60 PATH/TO/YOUR/lmp_mpi -in Au_TTM_pre.in` with all necessary parameter files in the same folder.

## `Common Runs`
The LAMMPS input file requires a restart file which would be generated from the Preheating LAMMPS codes. The restart file by default should be renamed from `restart.equil` to `TTM_restart.0`. Example: `mpiexec -np 60 PATH/TO/YOUR/lmp_mpi -in Au_TTM_restart.in` with all necessary parameter files in the same folder.

## `Au Parameters`
Example parameters for Au. They are required for common runs. The units are based on the metal-style of LAMMPS.

Ke: eV / (ps * A * K); Ce: eV / (A^3 * K); G: eV / (A^3 * ps * K).

## `Analysis`
The `read.cpp` and `read2.cpp` are to process the dump files from LAMMPS outputs and they should run in order. The rest MATLAB files are to generate figures using the temperature files from LAMMPS and the output files from `read.cpp` and `read2.cpp`.

## Example commands to run LAMMPS simulations
`mpiexec -np 60 PATH/TO/YOUR/lmp_mpi -in Au_TTM_restart.in` as in `run.sh`.

## Tips
- It is recommended to run `Preheating` and `Common Runs` in seperate folders.
- Use `tmux` in a remote server for long runs.
- Modify `void LAMMPS_NS::FixFEMTO3D::ChangeType()` for flexible atoms with different potentials.

## Reference
- [Paper: Ablation study in gold irradiated by single femtosecond laser pulse with electron temperature dependent interatomic potential and electron–phonon coupling factor](https://iopscience.iop.org/article/10.1088/1555-6611/abdcb8/meta)
