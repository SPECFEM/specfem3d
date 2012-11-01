#!/bin/bash
#PBS -S /bin/bash

## job name and output file
#PBS -N go_decomposer
#PBS -j oe
#PBS -o OUTPUT_FILES/$PBS_JOBID.o

###########################################################
# USER PARAMETERS

## 1 CPU, walltime 1 hour
#PBS -l nodes=1:ppn=1,walltime=1:00:00
##PBS -q debug

###########################################################

cd $PBS_O_WORKDIR

# script to run the mesher and the solver
# read Par_file to get information about the run
# compute total number of nodes needed
NPROC=`grep NPROC DATA/Par_file | cut -d = -f 2 `

# total number of nodes is the product of the values read
numnodes=$NPROC

echo starting decomposer for $numnodes partitions
echo " "

# save a copy
cp go_decomposer_pbs.bash OUTPUT_FILES/

# USER CHANGE MESH DIRECTORY
MESHDIR=examples/homogeneous_halfspace_HEX8/MESH/

cd bin/
./xdecompose_mesh $numnodes ../$MESHDIR ../OUTPUT_FILES/DATABASES_MPI/

echo "done "
