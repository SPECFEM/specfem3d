#!/bin/bash

#SBATCH -p debug
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH -t 10

#SBATCH --output=OUTPUT_FILES/%j.o
#SBATCH --job-name=go_combine_vol_data

umask 0022

cd $SLURM_SUBMIT_DIR

# script to run the mesher and the solver
# read Par_file to get information about the run
# compute total number of nodes needed
NPROC=`grep NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

# path to database files for mesh (relative to bin/)
LOCALPATH=`grep LOCAL_PATH DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

nmax=$(($NPROC-1))

make xcombine_vol_data

# model variable is vs; output file will be vs.vtk
#./bin/xcombine_vol_data_vtk 0 $nmax vs $LOCALPATH/ $LOCALPATH 0
srun -l /bin/hostname | sort -n | awk '{print $2}' > ./nodes.$SLURM_JOB_ID
mpirun -np 1 -machinefile ./nodes.$SLURM_JOB_ID ./bin/xcombine_vol_data_vtk 0 $nmax vs $LOCALPATH/ $LOCALPATH 0

cp go_combine_vol_data_slurm.bash OUTPUT_FILES/

echo "done "
