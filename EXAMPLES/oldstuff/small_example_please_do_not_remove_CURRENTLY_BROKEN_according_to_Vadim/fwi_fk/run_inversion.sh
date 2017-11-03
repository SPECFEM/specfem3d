#!/bin/bash

bin=../../../specfem3d_Git_devel/bin

echo "running example: `date`"
currentdir=`pwd`

# get the number of processors, ignoring comments in the Par_file
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

# define mesh with initial model
./create_mesh.sh
./copy_mesh.sh

# ------ run the first frequency group : 0.001Hz to 0.5Hz -----
# define inversion_fwi_par file
cp DATA/inverse_problem/inversion_fwi_par_1 DATA/inverse_problem/inversion_fwi_par
mpirun -np $NPROC $bin/xinverse_problem_for_model l-bfgs

#------ define initial model for second frequency group
./define_input_model_for_next_frequency.sh frq1 10

# ------ run the first frequency group : 0.001Hz to max frequency  -----
# define inversion_fwi_par file
cp DATA/inverse_problem/inversion_fwi_par_2 DATA/inverse_problem/inversion_fwi_par
mpirun -np $NPROC $bin/xinverse_problem_for_model l-bfgs
