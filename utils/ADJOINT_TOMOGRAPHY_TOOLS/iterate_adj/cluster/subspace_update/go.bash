#!/bin/bash -v
#BSUB -o OUTPUT_FILES/%J.o
#BSUB -a mpich_gm
#BSUB -J subspace_update

echo "$LSB_MCPU_HOSTS" > OUTPUT_FILES/lsf_machines
echo "$LSB_JOBID" > OUTPUT_FILES/jobid
remap_lsf_machines.pl OUTPUT_FILES/lsf_machines >OUTPUT_FILES/machines

mpirun.lsf --gm-no-shmem --gm-copy-env src/subspace_update
