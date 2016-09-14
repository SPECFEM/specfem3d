#!/bin/bash

#SBATCH -p debug
#SBATCH -n 4
#SBATCH -t 60

#SBATCH --output=OUTPUT_FILES/%j.o
#SBATCH --job-name=go_solver

cd $SLURM_SUBMIT_DIR

# script to run the solver
# read Par_file to get information about the run
# compute total number of nodes needed
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

mkdir -p OUTPUT_FILES

# backup files used for this simulation
cp go_solver_slurm.bash OUTPUT_FILES/
cp DATA/Par_file OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/
cp DATA/CMTSOLUTION OUTPUT_FILES/
cp src/shared/constants.h OUTPUT_FILES/

# save a complete copy of source files
#rm -rf OUTPUT_FILES/src
#cp -rp ./src OUTPUT_FILES/

# obtain job information
cat $SLURM_JOB_NODELIST > OUTPUT_FILES/compute_nodes
echo "$SLURM_JOBID" > OUTPUT_FILES/jobid

echo starting run in current directory $PWD
echo " "

sleep 2
mpiexec -np $NPROC ./bin/xspecfem3D

cp DATA/STATIONS_FILTERED OUTPUT_FILES/

echo "finished successfully"
