# Femto-Research
Codes for femtosecond laser interaction research

## Femto3D
LAMMPS function

## Femto3D LAMMPS function
This is the modified LAMMPS function based on its two-temperature function. The
corresponding files include ’fix_femto3D.cpp’ and ’fix_femto3D.h’ in the ’ ’ folder.

## Preheating LAMMPS codes
These files are in the ’ ’ folder. A target generated in MD is preheated to a
wanted temperature and then relaxes to equilibrium. The output file for further use is a
LAMMPS restart file for further simulation.

## Running LAMMPS codes
These files are in the ’ ’ folder. The LAMMPS input file requires a restart
file which would be generated from the Preheating LAMMPS codes.

## Analysis codes
These files are in the ’ ’ folder. The ’read.cpp’ and ’read2.cpp’ are to process the
dump files from LAMMPS outputs and they should run in order. The rest MATLAB files
are to generate figures using the temperature files from LAMMPS and the output files from
’read.cpp’ and ’read2.cpp’.
