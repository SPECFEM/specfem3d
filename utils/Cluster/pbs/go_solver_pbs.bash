#!/bin/bash
#PBS -S /bin/bash

## job name and output file
#PBS -N go_solver
#PBS -j oe
#PBS -o OUTPUT_FILES/$PBS_JOBID.o

## group/others read .o file
#PBS -Wumask=0022

###########################################################
# USER PARAMETERS

## 4 CPUs, walltime 1 hour
#PBS -l nodes=1:ppn=4,walltime=1:00:00
## queue name will depend on the cluster
#PBS -q debug

###########################################################

cd $PBS_O_WORKDIR

# script to run the mesher and the solver
# read Par_file to get information about the run
# compute total number of nodes needed
NPROC=`grep NPROC DATA/Par_file | cut -d = -f 2`

mkdir -p OUTPUT_FILES

# backup files used for this simulation
cp go_solver_pbs.bash OUTPUT_FILES/
cp DATA/Par_file OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/
cp DATA/CMTSOLUTION OUTPUT_FILES/
cp src/shared/constants.h OUTPUT_FILES/

# save a complete copy of source files
#rm -rf OUTPUT_FILES/src
#cp -rp ./src OUTPUT_FILES/

# obtain job information
cat $PBS_NODEFILE > OUTPUT_FILES/compute_nodes
echo "$PBS_JOBID" > OUTPUT_FILES/jobid

echo starting run in current directory $PWD
echo " "

sleep 2
cd bin/
mpiexec -np $NPROC ./xspecfem3D

cp ../DATA/STATIONS_FILTERED ../OUTPUT_FILES/

echo "finished successfully"

