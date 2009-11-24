#!/bin/bash
#PBS -S /bin/bash

## job name and output file
#PBS -N go_solver
#PBS -j oe
#PBS -o OUTPUT_FILES/job.o

## 4 CPUs ( 4  ), walltime 1 hour
#PBS -l nodes=1:ppn=4,walltime=1:00:00

cd $PBS_O_WORKDIR

if [ -z $USER ]; then
	echo "could not run go_solver_...bash as no USER env is set"
	exit 2
fi

# script to run the mesher and the solver
# read DATA/Par_file to get information about the run
# compute total number of nodes needed
NPROC_XI=`grep NPROC_XI DATA/Par_file | cut -d = -f 2 `
NPROC_ETA=`grep NPROC_ETA DATA/Par_file | cut -d = -f 2`

# total number of nodes is the product of the values read
numnodes=$(( $NPROC_XI * $NPROC_ETA ))

rm -r -f OUTPUT_FILES
mkdir OUTPUT_FILES
cp DATA/Par_file OUTPUT_FILES/
cp DATA/CMTSOLUTION OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/

# obtain lsf job information
cat $PBS_NODEFILE > OUTPUT_FILES/compute_nodes
echo "$PBS_JOBID" > OUTPUT_FILES/jobid

echo starting run in current directory $PWD
echo " "

sleep 2 
mpiexec -np $numnodes $PWD/xspecfem3D

echo "finished successfully"

