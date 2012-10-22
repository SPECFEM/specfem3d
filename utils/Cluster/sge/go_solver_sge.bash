#!/bin/sh
# sun grid engine cluster @ oxford

# use current working directory
#$ -cwd

# merge error output into standard output stream
#$ -j yes
#$ -o in_out_files/OUTPUT_FILES/job.o

###########################################################
# USER PARAMETERS

# specific environment with 180 cpu total, request to use 4 cpus:
#$ -pe make 4

###########################################################


# script to run the mesher and the solver
# read Par_file to get information about the run
# compute total number of nodes needed
NPROC=`grep NPROC DATA/Par_file | cut -d = -f 2 `

# total number of nodes is the product of the values read
numnodes=$NPROC

mkdir -p in_out_files/OUTPUT_FILES

# backup for files used for this simulation
cp DATA/Par_file in_out_files/OUTPUT_FILES/
cp DATA/STATIONS in_out_files/OUTPUT_FILES/
cp DATA/CMTSOLUTION in_out_files/OUTPUT_FILES/

rm -rf in_out_files/OUTPUT_FILES/src
cp -rp ./src in_out_files/OUTPUT_FILES/

echo starting run in current directory $PWD
echo " "

# run in parallel with sge as resource manager
#set echo
cd bin/
/opt/SUNWhpc/bin/mprun -x sge ./xspecfem3D

# or using mpiexec
#/opt/mpich/bin/mpiexec -np $numnodes ./xspecfem3D
