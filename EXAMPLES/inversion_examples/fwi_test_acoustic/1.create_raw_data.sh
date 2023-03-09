#!/bin/bash
#
#
# script to create "true" data from true model
#
###############################################
# USER PARAMETER

# directory holding "raw" data
DATA_DIR=./RAW_DATA

###############################################

echo "create_data: `date`"
currentdir=`pwd`

echo
echo "creating true event data"
echo

# setup model files
mkdir -p DATABASES_MPI OUTPUT_FILES
rm -rf DATABASES_MPI/* OUTPUT_FILES/*

cp -v MODEL_TRUE/proc* DATABASES_MPI/
cp -v MODEL_TRUE/values_from_mesher.h OUTPUT_FILES/

mkdir -p RAW_DATA

# runs forward simulations
./run_inverse_problem_for_model.sh forward

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "done"
echo




