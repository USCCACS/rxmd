#!/bin/sh

exec='srun -n 1024 --mpi=pmi2 ./rxmd --vprocs 32 32 1 --xyz_pto_tionly'

make -j rxmd 
${exec} --start_from_pto --random_velocity --treq 100 --mdmode 1 --ntime_step 5000 --dt 0.5 | tee log01
${exec} --treq 150 --mdmode 5 --sstep 10 --ntime_step 20000 --dt 0.5 | tee log02
