#!/bin/bash

#SBATCH -p debug
#SBATCH --ntasks=1
#SBATCH -t 10

#SBATCH --output=%j.o
#SBATCH --job-name=go_combine_vol_data

umask 0022

cd $SLURM_SUBMIT_DIR

# script to combine vol data
# read Par_file to get information about the run
# compute total number of nodes needed
NPROC=`cat DATA/Par_file | egrep "^NPROC" | awk '{ print $3 }'`

# path to database files for mesh (relative to bin/)
LOCALPATH=`cat DATA/Par_file | egrep "^LOCAL_PATH" | awk '{ print $3 }'`

nmax=$(($NPROC-1))

make xcombine_vol_data_vtk

# model variable is vs; output file will be vs.vtk
#./bin/xcombine_vol_data_vtk 0 $nmax vs $LOCALPATH/ $LOCALPATH 0
srun -l /bin/hostname | sort -n | awk '{print $2}' > ./nodes.$SLURM_JOB_ID
mpirun -np 1 -machinefile ./nodes.$SLURM_JOB_ID ./bin/xcombine_vol_data_vtk 0 $nmax vs $LOCALPATH/ $LOCALPATH 0

cp go_combine_vol_data_slurm.bash OUTPUT_FILES/

# obtain and store job information
echo "$SLURM_JOBID" > OUTPUT_FILES/jobid

JOBID=$(<OUTPUT_FILES/jobid)
mv $JOBID.o OUTPUT_FILES/

echo "done "
