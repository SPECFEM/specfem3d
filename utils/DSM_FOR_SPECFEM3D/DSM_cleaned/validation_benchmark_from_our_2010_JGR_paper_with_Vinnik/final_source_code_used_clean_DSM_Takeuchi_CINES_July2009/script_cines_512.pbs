#!/bin/bash
#PBS -N dsmti_r05
#PBS -m n
#PBS -l select=64:ncpus=8:mpiprocs=8
#PBS -l walltime=02:00:00

cd $PBS_O_WORKDIR

# save information about the job
cat $PBS_NODEFILE | tr '\n' ',' | sed s/,$// > host_list_${PBS_JOBID}
touch job_ID_is_${PBS_JOBID}

/usr/pbs/bin/mpiexec -n 512 dsmti < data/input_IASP91_regular_nocrust_nod410_nod670.inf

