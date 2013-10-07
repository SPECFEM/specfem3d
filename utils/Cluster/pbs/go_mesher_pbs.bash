#!/bin/bash
#PBS -S /bin/bash

## job name and output file
#PBS -N go_mesher
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

# number of cores for the job
NPROC=`grep NPROC DATA/Par_file | cut -d = -f 2`

mkdir -p OUTPUT_FILES
mkdir -p OUTPUT_FILES/DATABASES_MPI

# backup files used for this simulation
cp go_mesher_pbs.bash OUTPUT_FILES/
cp DATA/Par_file OUTPUT_FILES/
cp DATA/meshfem3D_files/Mesh_Par_file OUTPUT_FILES/

# save a complete copy of source files
#rm -rf OUTPUT_FILES/src
#cp -rp ./src OUTPUT_FILES/

# obtain pbs job information
cat $PBS_NODEFILE > OUTPUT_FILES/compute_nodes
echo "$PBS_JOBID" > OUTPUT_FILES/jobid

echo starting MPI internal mesher on $NPROC processors
echo " "

sleep 2
cd bin/
mpiexec -np $NPROC ./xmeshfem3D

echo "done "
