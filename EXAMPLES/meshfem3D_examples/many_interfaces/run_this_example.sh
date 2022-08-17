#!/bin/bash
#
# script runs mesher,database generation and solver
# using this example setup
#

###################################################

# number of processes
NPROC=405

##################################################

echo "running example: `date`"
currentdir=`pwd`

echo
echo "(will take about 1 h 15 minutes)"
echo

# sets up directory structure in current example directory
echo
echo "   setting up example..."
echo

# checks if executables were compiled and available
if [ ! -e ../../../bin/xspecfem3D ]; then
  echo "Please compile first all binaries in the root directory, before running this example..."; echo
  exit 1
fi

mkdir -p bin
mkdir -p OUTPUT_FILES/DATABASES_MPI

rm -f OUTPUT_FILES/*
rm -rf OUTPUT_FILES/DATABASES_MPI/*

cd $currentdir

# links executables
cd bin/
ln -s ../../../../bin/xmeshfem3D
ln -s ../../../../bin/xgenerate_databases
ln -s ../../../../bin/xspecfem3D
cd ../

# stores setup
cp DATA/meshfem3D_files/Mesh_Par_file OUTPUT_FILES/
cp DATA/Par_file OUTPUT_FILES/
cp DATA/CMTSOLUTION OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/

# decomposes mesh
echo
echo "running mesher..."
echo
mpirun -np $NPROC ./bin/xmeshfem3D
mv OUTPUT_FILES/output_mesher.txt OUTPUT_FILES/output_meshfem3D.txt

# runs database generation
echo
echo "  running database generation..."
echo
mpirun -np $NPROC ./bin/xgenerate_databases

# runs simulation
echo
echo "  running solver..."
echo
mpirun -np $NPROC ./bin/xspecfem3D

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
date


