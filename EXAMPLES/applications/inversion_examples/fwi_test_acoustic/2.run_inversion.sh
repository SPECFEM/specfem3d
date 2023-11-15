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
echo "running inversion"
echo

# setup model files
mkdir -p DATABASES_MPI OUTPUT_FILES
rm -rf DATABASES_MPI/* OUTPUT_FILES/*

cp -v MODEL_INIT/proc* DATABASES_MPI/
cp -v MODEL_INIT/values_from_mesher.h OUTPUT_FILES/

if [ ! -e "RAW_DATA" ]; then echo "could not find RAW_DATA/ directory, please check if data available..."; exit 1; fi
echo

# runs inversion
./run_inverse_problem_for_model.sh l-bfgs

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "done"
echo




