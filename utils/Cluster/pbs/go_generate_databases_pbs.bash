#!/bin/bash
#PBS -S /bin/bash

## job name and output file
#PBS -N go_mesher
#PBS -j oe
#PBS -o OUTPUT_FILES/job.o

###########################################################
# USER PARAMETERS

## 4 CPUs ( 4  ), walltime 1 hour
#PBS -l nodes=1:ppn=4,walltime=1:00:00

###########################################################

cd $PBS_O_WORKDIR


# script to run the mesher and the solver
# read DATA/Par_file to get information about the run
# compute total number of nodes needed
NPROC=`grep NPROC DATA/Par_file | cut -d = -f 2 `

# total number of nodes is the product of the values read
numnodes=$NPROC

mkdir -p OUTPUT_FILES

# backup files used for this simulation
cp DATA/Par_file OUTPUT_FILES/

rm -rf OUTPUT_FILES/src
mkdir OUTPUT_FILES/src
cp -p *.f90 OUTPUT_FILES/src/
cp -p *.c OUTPUT_FILES/src/
cp -p *.in OUTPUT_FILES/src/
cp -p *.h OUTPUT_FILES/src/

# obtain lsf job information
cat $PBS_NODEFILE > OUTPUT_FILES/compute_nodes
echo "$PBS_JOBID" > OUTPUT_FILES/jobid

echo starting MPI mesher on $numnodes processors
echo " "

sleep 2 
mpiexec -np $numnodes $PWD/xgenerate_databases

echo "done "
