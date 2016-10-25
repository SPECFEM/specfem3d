#!/bin/bash

#SBATCH -p debug
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH -t 60

#SBATCH --output=OUTPUT_FILES/%j.o
#SBATCH --job-name=go_decomposer

umask 0022

cd $SLURM_SUBMIT_DIR

# script to decompose mesh
# read Par_file to get information about the run
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

# path to database files for mesh
LOCALPATH=`grep ^LOCAL_PATH DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

echo starting decomposer for $NPROC partitions
echo " "

mkdir -p OUTPUT_FILES/

# save a copy
cp go_decomposer_slurm.bash OUTPUT_FILES/

# USER: CHANGE MESH DIRECTORY
#MESHDIR=./EXAMPLES/homogeneous_halfspace/MESH
MESHDIR=./EXAMPLES/homogeneous_halfspace/MESH-default

mkdir -p $LOCALPATH
./bin/xdecompose_mesh $NPROC $MESHDIR/ $LOCALPATH/

echo "done "
