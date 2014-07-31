#!/bin/bash -v
#BSUB -o OUTPUT_FILES/%J.o
#BSUB -a openmpi
#BSUB -J sum_kernels

echo "$LSB_MCPU_HOSTS" > OUTPUT_FILES/lsf_machines
echo "$LSB_JOBID" > OUTPUT_FILES/jobid
remap_lsf_machines.pl OUTPUT_FILES/lsf_machines >OUTPUT_FILES/machines

#mpirun src/sum_kernels
mpirun.lsf src/sum_kernels

# Caltech command, used with mpich_gm
#mpirun.lsf --gm-no-shmem --gm-copy-env src/sum_kernels

