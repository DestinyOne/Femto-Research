echo "Start --- `date`" >> time.txt
# mpiexec -np 16 ~/Research5/mylammps/src/lmp_mpi -sf gpu -pk gpu 1 -in TTM.in
mpiexec -np 60 ~/Research5/mylammps/src/lmp_mpi -in Au_TTM_restart.in
tail -1 log.lammps >> time.txt
echo -e "End   --- `date`\n" >> time.txt
