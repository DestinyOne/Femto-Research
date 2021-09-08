# Femto-Research
Codes for femtosecond laser interaction research

## `Femto3D`
This is the modified LAMMPS function based on its two-temperature function.

## `Preheating`
A target generated in MD is preheated to a wanted temperature and then relaxes to equilibrium. The output file for further use is a LAMMPS restart file for further simulation. The default output from the sample input code is `restart.equil.mpiio`.

## `Common Runs`
The LAMMPS input file requires a restart file which would be generated from the Preheating LAMMPS codes. The restart file by default should be renamed from `restart.equil.mpiio` to `TTM_restart_0.mpiio`.

## `Au Parameters`

## `Analysis`
The `read.cpp` and `read2.cpp` are to process the dump files from LAMMPS outputs and they should run in order. The rest MATLAB files are to generate figures using the temperature files from LAMMPS and the output files from `read.cpp` and `read2.cpp`.

## Example commands to run LAMMPS simulations
`mpiexec -np 60 PATH/TO/YOUR/lmp_mpi -in Au_TTM_restart.in` as in `run.sh`.

## Tips
- It is recommended to run `Preheating` and `Common Runs` in seperate folders.
- Use `tmux` in a remote server for long runs.
