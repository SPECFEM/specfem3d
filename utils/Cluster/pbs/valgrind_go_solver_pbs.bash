#!/bin/bash
#
# Valgrind, a suite of tools for debugging and profiling
# http://valgrind.org/
#

# bash script
#PBS -S /bin/bash

# job name
#PBS -N valgrind_go_solver

# joins output and error information
#PBS -j oe

# job output file
#PBS -o in_out_files/OUTPUT_FILES/job.o

###########################################################
# USER PARAMETERS

# Queue
#PBS -q tromp

# 150 CPUs ( 18*8+6 ), walltime 15 hour
#PBS -l nodes=18:ppn=8+1:ppn=6,walltime=15:00:00

# valgrind mpi library
PRELOAD_LIB=/my_valgrind_path/valgrind/lib/valgrind/libmpiwrap-x86-linux.so

###########################################################

cd $PBS_O_WORKDIR

# script to run the mesher and the solver
# read Par_file to get information about the run
# compute total number of nodes needed
NPROC=`grep NPROC DATA/Par_file | cut -d = -f 2 `

# total number of nodes is the product of the values read
numnodes=$NPROC

mkdir -p in_out_files/OUTPUT_FILES

# backup files used for this simulation
cp DATA/Par_file in_out_files/OUTPUT_FILES/
cp DATA/STATIONS in_out_files/OUTPUT_FILES/
cp DATA/CMTSOLUTION in_out_files/OUTPUT_FILES/

# obtain job information
cat $PBS_NODEFILE > in_out_files/OUTPUT_FILES/compute_nodes
echo "$PBS_JOBID" > in_out_files/OUTPUT_FILES/jobid

echo starting run in current directory $PWD
cd bin/

echo " "
echo "run: memory leaks"
echo " "
sleep 2

# memory leaks
LD_PRELOAD=$PRELOAD_LIB mpiexec -np $numnodes valgrind --leak-check=full ./xspecfem3D >& ../in_out_files/OUTPUT_FILES/output.memory-leaks.log

sleep 2
echo " "
echo "run: cache misses"
echo " "

# cache misses
LD_PRELOAD=$PRELOAD_LIB mpiexec -np $numnodes valgrind --tool=cachegrind ./xspecfem3D >& ../in_out_files/OUTPUT_FILES/output.cache-misses.log

cp cachegrind.out.* ../in_out_files/OUTPUT_FILES/

echo " "
echo "finished successfully"

