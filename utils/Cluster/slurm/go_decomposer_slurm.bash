#!/bin/bash

#SBATCH -p debug
#SBATCH --ntasks=1
#SBATCH -t 60

#SBATCH --output=%j.o
#SBATCH --job-name=go_decomposer

umask 0022

cd $SLURM_SUBMIT_DIR

# script to decompose mesh
# read Par_file to get information about the run
# compute total number of partitions needed
NPROC=`cat DATA/Par_file | egrep "^NPROC" | awk '{ print $3 }'`

mkdir -p OUTPUT_FILES
LOCALPATH=`cat DATA/Par_file | egrep "^LOCAL_PATH" | awk '{ print $3 }'`
mkdir -p $LOCALPATH

# USER: CHANGE MESH DIRECTORY
#MESHDIR=./EXAMPLES/homogeneous_halfspace/MESH
MESHDIR=./EXAMPLES/homogeneous_halfspace/MESH-default

# save a copy
cp go_decomposer_slurm.bash OUTPUT_FILES/

# obtain slurm job information
echo "$SLURM_JOBID" > OUTPUT_FILES/jobid

echo starting decomposer for $NPROC partitions
echo " "

./bin/xdecompose_mesh $NPROC $MESHDIR/ $LOCALPATH/

JOBID=$(<OUTPUT_FILES/jobid)
mv $JOBID.o OUTPUT_FILES/

echo "done "
