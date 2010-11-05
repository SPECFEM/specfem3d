#!/bin/bash

## job name and output file
#BSUB -J go_solver
#BSUB -o in_out_files/OUTPUT_FILES/%J.o
#BSUB -e in_out_files/OUTPUT_FILES/%J.e

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
# read Par_file to get information about the run
# compute total number of nodes needed
NPROC=`grep NPROC in_data_files/Par_file | cut -d = -f 2 `

# total number of nodes is the product of the values read
numnodes=$NPROC

rm -r -f in_out_files/OUTPUT_FILES
mkdir in_out_files/OUTPUT_FILES
cp in_data_files/Par_file in_out_files/OUTPUT_FILES/
cp in_data_files/CMTSOLUTION in_out_files/OUTPUT_FILES/
cp in_data_files/STATIONS in_out_files/OUTPUT_FILES/

# obtain lsf job information
cat $BSUB_DJOB_NODEFILE > in_out_files/OUTPUT_FILES/compute_nodes
echo "$LSB_JOBID" > in_out_files/OUTPUT_FILES/jobid

echo starting run in current directory $PWD
echo " "

sleep 2
cd bin/
mpirun.lsf ./xspecfem3D

echo "finished successfully"

