#!/bin/bash
#
# script runs mesher,database generation and solver
# using this example setup
#

###################################################

# number of processes
NPROC=81

##################################################

echo "running example: `date`"
currentdir=`pwd`

echo
echo "(will take about 1 h 15 minutes)"
echo

# sets up directory structure in current example directoy
echo
echo "   setting up example..."
echo

mkdir -p bin
mkdir -p OUTPUT_FILES
mkdir -p OUTPUT_FILES/DATABASES_MPI

rm -rf OUTPUT_FILES/*
rm -rf OUTPUT_FILES/DATABASES_MPI/*

mkdir -p DATA
mkdir -p DATA/meshfem3D_files/

cp Mesh_Par_file DATA/meshfem3D_files/
cp example_*.dat DATA/meshfem3D_files/

cp Par_file DATA/
cp CMTSOLUTION DATA/
cp STATIONS DATA/


# compiles executables in root directory
cd ../../../
make > tmp.log
cd $currentdir

# links executables
cd bin/
ln -s ../../../../bin/xmeshfem3D
ln -s ../../../../bin/xgenerate_databases
ln -s ../../../../bin/xspecfem3D
cd ../

# decomposes mesh
echo
echo "  meshing..."
echo
cd bin/
mpirun -np $NPROC ./xmeshfem3D
cd ../
mv OUTPUT_FILES/output_mesher.txt OUTPUT_FILES/output_meshfem3D.txt

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


