#!/bin/bash
#PBS -S /bin/bash

## job name and output file
#PBS -N go_database
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

# script to generate databases
# read Par_file to get information about the run
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

mkdir -p OUTPUT_FILES

# backup files used for this simulation
cp go_generate_databases_pbs.bash OUTPUT_FILES/
cp DATA/Par_file OUTPUT_FILES/

# save a complete copy of source files
#rm -rf OUTPUT_FILES/src
#cp -rp ./src OUTPUT_FILES/

# obtain job information
cat $PBS_NODEFILE > OUTPUT_FILES/compute_nodes
echo "$PBS_JOBID" > OUTPUT_FILES/jobid

echo starting MPI mesher on $NPROC processors
echo " "

sleep 2
mpiexec -np $NPROC ./bin/xgenerate_databases

echo "done "
