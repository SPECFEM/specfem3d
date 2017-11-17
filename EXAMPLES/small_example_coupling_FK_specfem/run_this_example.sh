#!/usr/bin/env bash


# specfem binaries directory
specfem_bin=/mnt/Data1/vmont/GIT/specfem3d/bin/

# get the number of processors, ignoring comments in the Par_file
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2 | cut -d \# -f 1`

mkdir OUTPUT_FILES DATABASES_MPI

# creating mesh with meshfem
echo
echo "  creating  mesh..."
echo

mpirun -np $NPROC $specfem_bin/xmeshfem3D


# runs database generation
echo
echo "  running database generation on $NPROC processors..."
echo
mpirun -np $NPROC $specfem_bin/xgenerate_databases


# runs simulation
echo
echo "  running solver on $NPROC processors..."
echo
mpirun -np $NPROC $specfem_bin/xspecfem3D

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
echo `date`

