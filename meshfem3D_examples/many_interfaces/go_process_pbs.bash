#!/bin/bash
#
# usage: submit this script in the current example directory like
#            > qsub go_process_pbs.bash
#
#PBS -S /bin/bash

## job name and output file
#PBS -N go_process
#PBS -j oe
#PBS -o in_out_files/OUTPUT_FILES/$PBS_JOBID.o
## PBS -q test

###########################################################
# USER PARAMETERS

## number of processes
## 81 CPUs ( 10*8 + 1  ), walltime 4 hour
#PBS -l nodes=10:ppn=8+1:ppn=1,walltime=4:00:00

###########################################################

cd $PBS_O_WORKDIR

echo "running example: `date`"
currentdir=`pwd`

echo
echo "(will take about 1 h 15 minutes)"
echo

# sets up directory structure in current example directoy
echo
echo "   setting up example..."
echo

# setup
mkdir -p bin
mkdir -p in_out_files/OUTPUT_FILES
mkdir -p in_out_files/DATABASES_MPI

rm -rf in_out_files/OUTPUT_FILES/*
rm -rf in_out_files/DATABASES_MPI/*

mkdir -p in_data_files
mkdir -p in_data_files/meshfem3D_files/

cp Mesh_Par_file in_data_files/meshfem3D_files/
cp example_*.dat in_data_files/meshfem3D_files/

cp Par_file in_data_files/
cp CMTSOLUTION in_data_files/
cp STATIONS in_data_files/

# compilation of executables must have been done prior on front node
# links executables
cd bin/
ln -s ../../../../bin/xmeshfem3D
ln -s ../../../../bin/xgenerate_databases
ln -s ../../../../bin/xspecfem3D
cd ../

# script to run the mesher and the solver
# read Par_file to get information about the run
# compute total number of nodes needed
NPROC=`grep NPROC in_data_files/Par_file | cut -d = -f 2 `

# total number of nodes is the product of the values read
numnodes=$NPROC

# decomposes mesh
echo
echo "  meshing..."
echo
cd bin/
mpiexec -np $numnodes ./xmeshfem3D
cd ../
mv in_out_files/OUTPUT_FILES/output_mesher.txt in_out_files/OUTPUT_FILES/output_meshfem3D.txt

# stores setup
cp in_data_files/Par_file in_out_files/OUTPUT_FILES/
cp in_data_files/CMTSOLUTION in_out_files/OUTPUT_FILES/
cp in_data_files/STATIONS in_out_files/OUTPUT_FILES/

# runs database generation
echo
echo "  running database generation..."
echo
cd bin/
mpiexec -np $numnodes ./xgenerate_databases
cd ../

# runs simulation
echo
echo "  running solver..."
echo
cd bin/
mpiexec -np $numnodes ./xspecfem3D
cd ../

echo
echo "see results in directory: in_out_files/OUTPUT_FILES/"
echo
echo "done"
echo `date`
