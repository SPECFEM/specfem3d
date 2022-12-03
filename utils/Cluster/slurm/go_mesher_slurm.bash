#!/bin/bash

#SBATCH -p debug
#SBATCH --ntasks=4
#SBATCH -t 60

#SBATCH --output=%j.o
#SBATCH --job-name=go_mesher

umask 0022

cd $SLURM_SUBMIT_DIR

# script to make mesh internally
# read Par_file to get information about the run
# compute total number of nodes needed
NPROC=`cat DATA/Par_file | egrep "^NPROC" | awk '{ print $3 }'`

mkdir -p OUTPUT_FILES
LOCALPATH=`cat DATA/Par_file | egrep "^LOCAL_PATH" | awk '{ print $3 }'`
mkdir -p $LOCALPATH

# backup files used for mesh generation
cp go_mesher_slurm.bash OUTPUT_FILES/
cp -r DATA/meshfem3D_files OUTPUT_FILES/

# save a complete copy of source files
#rm -rf OUTPUT_FILES/src
#cp -rp ./src OUTPUT_FILES/

# obtain slurm job information
echo "$SLURM_JOB_NODELIST" > OUTPUT_FILES/compute_nodes
echo "$SLURM_JOBID" > OUTPUT_FILES/jobid

echo starting MPI internal mesher on $NPROC processors
echo " "

sleep 2
mpiexec -np $NPROC ./bin/xmeshfem3D

JOBID=$(<OUTPUT_FILES/jobid)
mv $JOBID.o OUTPUT_FILES/

echo "done "
