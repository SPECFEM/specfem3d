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
echo "(will take about 5 minutes)"
echo

# sets up directory structure in current example directoy
echo
echo "   setting up example..."
echo

mkdir -p bin
mkdir -p OUTPUT_FILES/DATABASES_MPI

rm -f OUTPUT_FILES/*
rm -rf OUTPUT_FILES/DATABASES_MPI/*

# compiles executables in root directory
cd ../../
make clean
make > $currentdir/tmp.log
cd $currentdir

# links executables
cd bin/
rm -f ./x*
cp ../../../bin/xdecompose_mesh ./
cp ../../../bin/xgenerate_databases ./
cp ../../../bin/xspecfem3D ./
cd ../

# decomposes mesh
echo
echo "  decomposing mesh..."
echo
./bin/xdecompose_mesh $NPROC MESH/ OUTPUT_FILES/DATABASES_MPI/

# stores setup
cp DATA/Par_file OUTPUT_FILES/
cp DATA/CMTSOLUTION OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/

# runs database generation
echo
echo "  running database generation..."
echo
cd bin/
mpirun -np $NPROC ./xgenerate_databases
cd ../

# runs simulation
echo
echo "  running solver..."
echo
cd bin/
mpirun -np $NPROC ./xspecfem3D
cd ../

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
echo `date`


