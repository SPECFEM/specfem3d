#!/bin/bash
#PBS -S /bin/bash

## job name and output file
#PBS -N go_decomposer
#PBS -j oe
#PBS -o OUTPUT_FILES/$PBS_JOBID.o

## group/others read .o file
#PBS -Wumask=0022

###########################################################
# USER PARAMETERS

## 1 CPU, walltime 1 hour
#PBS -l nodes=1:ppn=1,walltime=1:00:00
## queue name will depend on the cluster
#PBS -q debug

###########################################################

cd $PBS_O_WORKDIR

# script to run the mesher and the solver
# read Par_file to get information about the run
# compute total number of nodes needed
NPROC=`grep NPROC DATA/Par_file | cut -d = -f 2`

# path to database files for mesh (relative to bin/)
LOCALPATH=`grep LOCAL_PATH DATA/Par_file | cut -d = -f 2`

echo starting decomposer for $NPROC partitions
echo " "

mkdir -p OUTPUT_FILES/

# save a copy
cp go_decomposer_pbs.bash OUTPUT_FILES/

# USER CHANGE MESH DIRECTORY
MESHDIR=examples/homogeneous_halfspace/MESH

cd bin/
mkdir -p $LOCALPATH
./xdecompose_mesh $NPROC ../$MESHDIR/ $LOCALPATH/

echo "done "
