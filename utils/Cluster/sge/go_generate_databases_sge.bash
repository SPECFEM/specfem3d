#!/bin/sh 
# sun grid engine cluster @ oxford

# use current working directory
#$ -cwd

# merge error output into standard output stream
#$ -j yes
#$ -o OUTPUT_FILES/job.o

###########################################################
# USER PARAMETERS

# specific environment with 180 cpu total, request to use 4 cpus:
#$ -pe make 4

###########################################################


# script to run the mesher and the solver
# read DATA/Par_file to get information about the run
# compute total number of nodes needed
NPROC=`grep NPROC DATA/Par_file | cut -d = -f 2 `

# total number of nodes is the product of the values read
numnodes=$NPROC

mkdir -p OUTPUT_FILES

# backup for files used for this simulation
cp DATA/Par_file OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/
cp DATA/CMTSOLUTION OUTPUT_FILES/

rm -rf OUTPUT_FILES/src
mkdir OUTPUT_FILES/src
cp -p *.f90 OUTPUT_FILES/src/
cp -p *.c OUTPUT_FILES/src/
cp -p *.in OUTPUT_FILES/src/
cp -p *.h OUTPUT_FILES/src/

echo starting run in current directory $PWD
echo " "

# run in parallel with sge as resource manager
#set echo

/opt/SUNWhpc/bin/mprun -x sge ./xgenerate_databases

# or using mpiexec 
#/opt/mpich/bin/mpiexec -np $numnodes ./xgenerate_databases


