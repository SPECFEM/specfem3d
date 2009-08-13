#!/bin/bash
#PBS -S /bin/bash

## job name and output file
#PBS -N go_mesher
#PBS -j oe
#PBS -o OUTPUT_FILES/job.o

## 25 CPUs ( 1 + 3*8 ), walltime 10 hours
##PBS -l nodes=3:ppn=8+1:ppn=1,walltime=10:00:00

## 4 CPUs ( 4  ), walltime 1 hour
#PBS -l nodes=1:ppn=4,walltime=1:00:00

cd /home/cmorency/SPECFEM3D_SESAME/

if [ -z $USER ]; then
	echo "could not run go_mesher_...bash as no USER env is set"
	exit 2
fi
BASEMPIDIR=/home/cmorency/SPECFEM3D_SESAME/OUTPUT_FILES

# script to run the mesher and the solver
# read DATA/Par_file to get information about the run
# compute total number of nodes needed
NPROC_XI=`grep NPROC_XI DATA/Par_file | cut -d = -f 2 `
NPROC_ETA=`grep NPROC_ETA DATA/Par_file | cut -d = -f 2`

# total number of nodes is the product of the values read
numnodes=$(( $NPROC_XI * $NPROC_ETA ))

cp DATA/Par_file OUTPUT_FILES/

# obtain lsf job information
cat $PBS_NODEFILE > OUTPUT_FILES/compute_nodes
echo "$PBS_JOBID" > OUTPUT_FILES/jobid

echo starting MPI mesher on $numnodes processors
echo " "

sleep 2 
mpiexec -np $numnodes $PWD/xgenerate_databases

echo "done "
