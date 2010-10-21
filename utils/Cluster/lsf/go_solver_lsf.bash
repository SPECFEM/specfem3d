#!/bin/bash

## job name and output file
#BSUB -J go_solver
#BSUB -o OUTPUT_FILES/%J.o
#BSUB -e OUTPUT_FILES/%J.e

###########################################################
# USER PARAMETERS

## 4 CPUs ( 4  ), walltime 1 hour
#BSUB -n 4
#BSUB -R span[ptile=4]
#BSUB -q short_parallel

#BSUB -a openmpi

###########################################################

if [ -z $USER ]; then
	echo "could not run go_solver_...bash as no USER env is set"
	exit 2
fi

# script to run the mesher and the solver
# read DATA/Par_file to get information about the run
# compute total number of nodes needed
NPROC=`grep NPROC DATA/Par_file | cut -d = -f 2 `

# total number of nodes is the product of the values read
numnodes=$NPROC

rm -r -f OUTPUT_FILES
mkdir OUTPUT_FILES
cp DATA/Par_file OUTPUT_FILES/
cp DATA/CMTSOLUTION OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/

# obtain lsf job information
cat $BSUB_DJOB_NODEFILE > OUTPUT_FILES/compute_nodes
echo "$LSB_JOBID" > OUTPUT_FILES/jobid

echo starting run in current directory $PWD
echo " "

sleep 2 
mpirun.lsf $PWD/xspecfem3D

echo "finished successfully"

