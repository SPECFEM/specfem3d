#!/bin/bash

#SBATCH -p debug
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH -t 60

#SBATCH --output=OUTPUT_FILES/%j.o
#SBATCH --job-name=go_database

umask 0022

cd $SLURM_SUBMIT_DIR

# script to generate databases
# read Par_file to get information about the run
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

mkdir -p OUTPUT_FILES

# backup files used for this simulation
cp go_generate_databases_slurm.bash OUTPUT_FILES/
cp DATA/Par_file OUTPUT_FILES/

# save a complete copy of source files
#rm -rf OUTPUT_FILES/src
#cp -rp ./src OUTPUT_FILES/

# obtain job information
cat $SLURM_JOB_NODELIST > OUTPUT_FILES/compute_nodes
echo "$SLURM_JOBID" > OUTPUT_FILES/jobid

echo starting MPI mesher on $NPROC processors
echo " "

sleep 2
mpiexec -np $NPROC ./bin/xgenerate_databases

echo "done "
