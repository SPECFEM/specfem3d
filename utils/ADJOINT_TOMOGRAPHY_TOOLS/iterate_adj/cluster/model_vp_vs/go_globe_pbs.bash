#!/bin/bash
#PBS -S /bin/bash

## job name and output file
#PBS -N model_update
#PBS -j oe
#PBS -o OUTPUT_FILES/job.o

###########################################################
# USER PARAMETERS

## 144 CPUs ( 18*8  ), walltime 5 hour
#PBS -l nodes=18:ppn=8,walltime=5:00:00

numnodes=144

# model update percentage 3%
percentage=0.03


###########################################################

# (takes about 10 min...)

date

cd $PBS_O_WORKDIR

# obtain lsf job information
cat $PBS_NODEFILE > OUTPUT_FILES/compute_nodes
echo "$PBS_JOBID" > OUTPUT_FILES/jobid

# steepest descent step
mpiexec -np $numnodes $PWD/xadd_model_iso $percentage


echo "done successfully"
date
