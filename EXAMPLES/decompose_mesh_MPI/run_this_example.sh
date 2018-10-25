#!/bin/bash

echo "running example: `date`"
currentdir=`pwd`

# sets up directory structure in current example directory
echo
echo "   setting up example..."
echo

rm -f -r OUTPUT_FILES

mkdir -p bin
mkdir -p OUTPUT_FILES
mkdir -p OUTPUT_FILES/DATABASES_MPI

# stores setup
cp DATA/Par_file OUTPUT_FILES/
cp DATA/CMTSOLUTION OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/

# links executables
cd bin
rm -f xdecompose_mesh_mpi xgenerate_databases xspecfem3D
ln -s ../../../bin/xdecompose_mesh_mpi
ln -s ../../../bin/xgenerate_databases
ln -s ../../../bin/xspecfem3D
cd ..

# get the number of processors, ignoring comments in the Par_file
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2 | cut -d \# -f 1`
echo "The simulation will run on NPROC = " $NPROC " MPI tasks"

# decomposes mesh using the pre-saved mesh files in the MESH directory
echo
echo "  decomposing mesh..."
echo
 mpirun -np $NPROC ./bin/xdecompose_mesh_mpi 2 2 2

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
echo `date`


