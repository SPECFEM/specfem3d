#!/bin/bash
#
# script runs decomposition,database generation and solver
# using this example setup
#
# prior to running this script, you must create the mesh files
# in directory MESH/
# (see section 3.1 "Meshing with CUBIT" in user guide)
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
cd ../../../

rm -fr DATA/*
cd $currentdir
cp -fr DATA/* ../../../DATA/.

cd ../../../

make clean
./configure
make all > $currentdir/tmp.log
cd $currentdir

# links executables
cd bin/
rm -f *
cp ../../../../bin/xdecompose_mesh .
cp ../../../../bin/xgenerate_databases .
cp ../../../../bin/xspecfem3D .
cd ../

# stores setup
cp DATA/Par_file OUTPUT_FILES/
cp DATA/CMTSOLUTION OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/

# decomposes mesh
echo
echo "  decomposing mesh..."
echo
cd bin/
./xdecompose_mesh $NPROC ../MESH-default ../OUTPUT_FILES/DATABASES_MPI/

# runs database generation
echo
echo "  running database generation..."
echo
mpirun -np $NPROC ./xgenerate_databases

# runs simulation
echo
echo "  running solver..."
echo
mpirun -np $NPROC ./xspecfem3D
cd ../

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
echo `date`


