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
ln -s EXAMPLES/meshfem3D_examples/simple_model/DATA DATA

rm -rf OUTPUT_FILES

mkdir OUTPUT_FILES
mkdir OUTPUT_FILES/DATABASES_MPI

# creates and decomposes mesh
echo
echo "running mesher..."
echo
mpirun -np $NPROC ./bin/xmeshfem3D

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

