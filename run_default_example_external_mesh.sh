#!/bin/bash
#
# script runs mesher,database generation and solver
# using this example setup
#

###################################################

# number of processes
NPROC=4

##################################################

# setup default example
rm -rf DATA
ln -s EXAMPLES/homogeneous_halfspace/DATA DATA

rm -rf OUTPUT_FILES

mkdir OUTPUT_FILES
mkdir OUTPUT_FILES/DATABASES_MPI

# decompose an existing external mesh
echo
echo "running mesh decomposer..."
echo
./bin/xdecompose_mesh $NPROC ./MESH-default ./OUTPUT_FILES/DATABASES_MPI/

# runs database generation
echo
echo "running database generation..."
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
echo `date`

