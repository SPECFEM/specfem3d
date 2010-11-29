#!/bin/bash
#
# script runs mesher,database generation and solver
# using this example setup
#

###################################################

# number of processes
NPROC=4

##################################################

echo "running example: `date`"
currentdir=`pwd`

echo
echo "(will take about 15 minutes)"
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

mkdir -p in_data_files
mkdir -p in_data_files/meshfem3D_files/

cp Mesh_Par_file in_data_files/meshfem3D_files/
cp interface*.dat in_data_files/meshfem3D_files/

cp Par_file in_data_files/
cp CMTSOLUTION in_data_files/
cp STATIONS in_data_files/


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

# stores setup
cp in_data_files/Par_file in_out_files/OUTPUT_FILES/
cp in_data_files/CMTSOLUTION in_out_files/OUTPUT_FILES/
cp in_data_files/STATIONS in_out_files/OUTPUT_FILES/

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


