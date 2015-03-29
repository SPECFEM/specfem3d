#!/bin/bash
#PBS -S /bin/bash

## job name and output file
#PBS -N sum_kernel
#PBS -j oe
#PBS -o OUTPUT_FILES/job.o

###########################################################
# USER PARAMETERS

## e.g. 144 CPUs ( 18*8  ), walltime 2 hour
#PBS -l nodes=18:ppn=8,walltime=2:00:00

numnodes=144

###########################################################

#(takes about 10 minutes...)

date

cd $PBS_O_WORKDIR

# obtain lsf job information
cat $PBS_NODEFILE > OUTPUT_FILES/compute_nodes
echo "$PBS_JOBID" > OUTPUT_FILES/jobid

mpiexec -np $numnodes $PWD/xsum_kernels

