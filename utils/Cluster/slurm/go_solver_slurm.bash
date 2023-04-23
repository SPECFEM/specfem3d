#!/bin/bash

#SBATCH -p debug
#SBATCH -n 4
#SBATCH -t 60

#SBATCH --output=%j.o
#SBATCH --job-name=go_solver

umask 0022

cd $SLURM_SUBMIT_DIR

# script to run the solver
# read Par_file to get information about the run
# compute total number of nodes needed
NPROC=`cat DATA/Par_file | egrep "^NPROC" | awk '{ print $3 }'`

FORCESOLUTION=`cat DATA/Par_file | egrep "^USE_FORCE_POINT_SOURCE" | awk '{ print $3 }'`

EXTERNAL_STF=`cat DATA/Par_file | egrep "^USE_EXTERNAL_SOURCE_FILE" | awk '{ print $3 }'`

mkdir -p OUTPUT_FILES

# backup files used for this simulation
cp go_solver_slurm.bash OUTPUT_FILES/
cp DATA/Par_file OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/
cp setup/constants.h OUTPUT_FILES/

if [[ "${FORCESOLUTION}" == ".true." ]]; then
    cp DATA/FORCESOLUTION OUTPUT_FILES/
else
    cp DATA/CMTSOLUTION OUTPUT_FILES/
fi

if [[ "${EXTERNAL_STF}" == ".true." ]]; then
    cp source_time_function.txt OUTPUT_FILES/
fi

# save a complete copy of source files
#rm -rf OUTPUT_FILES/src
#cp -rp ./src OUTPUT_FILES/

# obtain job information
echo "$SLURM_JOB_NODELIST" > OUTPUT_FILES/compute_nodes
echo "$SLURM_JOBID" > OUTPUT_FILES/jobid

echo starting solver on $NPROC processors
echo " "

sleep 2
mpiexec -np $NPROC ./bin/xspecfem3D

cp DATA/STATIONS_FILTERED OUTPUT_FILES/

JOBID=$(<OUTPUT_FILES/jobid)
mv $JOBID.o OUTPUT_FILES/

echo "done"
