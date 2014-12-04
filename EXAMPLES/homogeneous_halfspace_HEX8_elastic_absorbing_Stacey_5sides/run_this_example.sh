#!/bin/bash

echo "running example: `date`"
currentdir=`pwd`

echo
echo "(will take about 15 minutes)"
echo

# sets up directory structure in current example directoy
echo
echo "   setting up example..."
echo

rm -f -r OUTPUT_FILES

mkdir -p bin
mkdir -p OUTPUT_FILES
mkdir -p OUTPUT_FILES/DATABASES_MPI

# compiles executables in root directory
cd ../../

rm -fr DATA/*
cd $currentdir
cp -fr DATA/* ../../DATA/.

cd ../../


# echo "(won't work if specfem3d is not already configured)"

#make clean
#./configure CC=icc FC=ifort MPIFC=mpif90 --with-mpi --enable-double-precision --with-scotch-dir=/home/cristini/CODES/scotch_6.0.0
make all # > $currentdir/tmp.log
cd $currentdir

# links executables
cd bin/
rm -f *
cp ../../../bin/xdecompose_mesh .
cp ../../../bin/xgenerate_databases .
cp ../../../bin/xspecfem3D .
cd ../

# stores setup
cp DATA/Par_file OUTPUT_FILES/
cp DATA/CMTSOLUTION OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/

# get the number of processors
NPROC=`grep NPROC DATA/Par_file | cut -d = -f 2`

# decomposes mesh using the pre-saved mesh files in MESH-default
echo
echo "  decomposing mesh..."
echo
./bin/xdecompose_mesh $NPROC ./MESH-default ./OUTPUT_FILES/DATABASES_MPI/

# runs database generation
echo
echo "  running database generation on $NPROC processors..."
echo
mpirun -np $NPROC ./bin/xgenerate_databases

# runs simulation
echo
echo "  running solver on $NPROC processors..."
echo
mpirun -np $NPROC ./bin/xspecfem3D

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
date


