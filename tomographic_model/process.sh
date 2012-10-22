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
echo "(will take about 10 minutes)"
echo

# sets up directory structure in current example directoy
echo
echo "   setting up example..."
echo

mkdir -p bin
mkdir -p in_out_files/OUTPUT_FILES
mkdir -p in_out_files/DATABASES_MPI

rm -rf in_out_files/OUTPUT_FILES/*
rm -rf in_out_files/DATABASES_MPI/*

# sets up tomography model file
./create_tomography_model_file.sh
mv tomography_model.xyz DATA/
echo

# compiles executables in root directory
cd ../../
make > tmp.log
cd $currentdir

# links executables
cd bin/
rm -f ./x*
cp ../../../bin/xdecompose_mesh_SCOTCH ./
cp ../../../bin/xgenerate_databases ./
cp ../../../bin/xspecfem3D ./
cd ../

# decomposes mesh
echo
echo "  decomposing mesh..."
echo
./bin/xdecompose_mesh_SCOTCH $NPROC MESH/ in_out_files/DATABASES_MPI/

# stores setup
cp DATA/Par_file in_out_files/OUTPUT_FILES/
cp DATA/CMTSOLUTION in_out_files/OUTPUT_FILES/
cp DATA/STATIONS in_out_files/OUTPUT_FILES/

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
echo "see results in directory: in_out_files/OUTPUT_FILES/"
echo
echo "done"
echo `date`


