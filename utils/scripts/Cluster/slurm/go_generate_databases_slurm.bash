#!/bin/bash

#SBATCH -p debug
#SBATCH --ntasks=4
#SBATCH -t 60

#SBATCH --output=%j.o
#SBATCH --job-name=go_database

umask 0022

cd $SLURM_SUBMIT_DIR

# script to generate databases
# read Par_file to get information about the run
# compute total number of nodes needed
NPROC=`cat DATA/Par_file | egrep "^NPROC" | awk '{ print $3 }'`
# check the type of model
MODEL=`cat DATA/Par_file | egrep "^MODEL" | awk '{ print $3 }'`

mkdir -p OUTPUT_FILES

# backup tomography files if any for this simulation
if [[ "${MODEL}" == "tomo" ]]; then
    cp -r DATA/tomo_files OUTPUT_FILES/
fi

# backup files used for database generation
cp go_generate_databases_slurm.bash OUTPUT_FILES/
cp DATA/Par_file OUTPUT_FILES/

# save a complete copy of source files
#rm -rf OUTPUT_FILES/src
#cp -rp ./src OUTPUT_FILES/

# obtain job information
echo "$SLURM_JOB_NODELIST" > OUTPUT_FILES/compute_nodes
echo "$SLURM_JOBID" > OUTPUT_FILES/jobid

echo starting MPI database generation on $NPROC processors
echo " "

sleep 2
mpiexec -np $NPROC ./bin/xgenerate_databases

JOBID=$(<OUTPUT_FILES/jobid)
mv $JOBID.o OUTPUT_FILES/

echo "done "
