#!/bin/sh

exec='srun -n 256 --mpi=pmi2 ./rxmd --vprocs 16 16 1 --xyz_pto_tionly'

make -j rxmd 
${exec} --start_from_pto --random_velocity --treq 150 --mdmode 1 --sstep 200 --ntime_step 1000 --dt 0.2
${exec} --treq 10 --mdmode 7 --sstep 10 --ntime_step 10000 --dt 1.0
