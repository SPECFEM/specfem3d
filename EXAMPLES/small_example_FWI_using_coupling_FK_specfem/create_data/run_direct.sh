#!/bin/bash

bin=../../../specfem3d_Git_devel/bin

echo "running example: `date`"
currentdir=`pwd`

# get the number of processors, ignoring comments in the Par_file
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

mpirun -np $NPROC $bin/xinverse_problem_for_model forward
