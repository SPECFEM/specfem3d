#!/bin/bash -v
#BSUB -o OUTPUT_FILES/%J.o
#BSUB -a mpich_gm
#BSUB -J go_mesher_solver_lsf

BASEMPIDIR=/scratch/$USER/DATABASES_MPI

echo "$LSB_MCPU_HOSTS" > OUTPUT_FILES/lsf_machines
echo "$LSB_JOBID" > OUTPUT_FILES/jobid

remap_lsf_machines.pl OUTPUT_FILES/lsf_machines >OUTPUT_FILES/machines

shmux -M50 -Sall -c "mkdir -p $BASEMPIDIR" - < OUTPUT_FILES/machines >/dev/null

mpirun.lsf --gm-no-shmem --gm-copy-env $PWD/remap old_machines 150
