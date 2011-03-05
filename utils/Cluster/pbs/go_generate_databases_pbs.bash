#!/bin/bash
#PBS -S /bin/bash

## job name and output file
#PBS -N go_mesher
#PBS -j oe
#PBS -o in_out_files/OUTPUT_FILES/$PBS_JOBID.o

###########################################################
# USER PARAMETERS

## 4 CPUs ( 4  ), walltime 1 hour
#PBS -l nodes=1:ppn=4,walltime=1:00:00
##PBS -q debug

###########################################################

cd $PBS_O_WORKDIR

# script to run the mesher and the solver
# read Par_file to get information about the run
# compute total number of nodes needed
NPROC=`grep NPROC in_data_files/Par_file | cut -d = -f 2 `

# total number of nodes is the product of the values read
numnodes=$NPROC

mkdir -p in_out_files/OUTPUT_FILES

# backup files used for this simulation
cp in_data_files/Par_file in_out_files/OUTPUT_FILES/

# save a complete copy of source files
#rm -rf in_out_files/OUTPUT_FILES/src
#cp -rp ./src in_out_files/OUTPUT_FILES/

# obtain lsf job information
cat $PBS_NODEFILE > in_out_files/OUTPUT_FILES/compute_nodes
echo "$PBS_JOBID" > in_out_files/OUTPUT_FILES/jobid

echo starting MPI mesher on $numnodes processors
echo " "

sleep 2
cd bin/
mpiexec -np $numnodes ./xgenerate_databases

echo "done "
