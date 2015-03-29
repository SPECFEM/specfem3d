#!/bin/bash
#
# script runs decomposition,database generation and solver
# using this example setup
#
# prior to running this script, you must create the mesh files
# in directory MESH/
#

###################################################

# number of processes
NPROC=4

##################################################

echo "running example: `date`"
currentdir=`pwd`

echo
echo "(will take about 27 minutes)"
echo


# sets up directory structure in current example directory
echo
echo "   setting up example..."
echo

mkdir -p bin
mkdir -p OUTPUT_FILES/DATABASES_MPI

rm -f OUTPUT_FILES/*
rm -rf OUTPUT_FILES/DATABASES_MPI/*

cd $currentdir

# links executables
cd bin/
rm -f ./x*
ln -s ../../../bin/xdecompose_mesh
ln -s ../../../bin/xgenerate_databases
ln -s ../../../bin/xspecfem3D
cd ../

# stores setup
cp DATA/Par_file OUTPUT_FILES/
cp DATA/CMTSOLUTION OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/

# decomposes mesh (THIS WILL NOT WORK UNLESS IN THE bin/ DIR!)
echo
echo "  decomposing mesh..."
echo
./bin/xdecompose_mesh $NPROC ./MESH/ ./OUTPUT_FILES/DATABASES_MPI/

# runs database generation
echo
echo "  running database generation..."
echo
mpiexec -np $NPROC ./bin/xgenerate_databases

# runs simulation
echo
echo "  running solver..."
echo
mpiexec -np $NPROC ./bin/xspecfem3D

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
date


