#!/bin/bash -v
#BSUB -o OUTPUT_FILES/%J.o
#BSUB -a mpich_gm
#BSUB -J add_model

echo "$LSB_MCPU_HOSTS" > OUTPUT_FILES/lsf_machines
echo "$LSB_JOBID" > OUTPUT_FILES/jobid
remap_lsf_machines.pl OUTPUT_FILES/lsf_machines >OUTPUT_FILES/machines

mpirun.lsf --gm-no-shmem --gm-copy-env src/add_model

