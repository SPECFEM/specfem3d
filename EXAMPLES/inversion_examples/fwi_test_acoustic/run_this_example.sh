#!/bin/bash
#
# inversion example
#

echo "running example: `date`"
currentdir=`pwd`

# checks if executables were compiled and available
if [ ! -e ../../../bin/xspecfem3D ]; then
  echo "Please compile first all binaries in the root directory, before running this example..."; echo
  exit 1
fi

# cleans output files
mkdir -p OUTPUT_FILES
rm -rf OUTPUT_FILES/*

mkdir -p DATABASES_MPI
rm -rf DATABASES_MPI/*

mkdir -p MODEL_INIT
rm -rf MODEL_INIT/*

mkdir -p MODEL_TRUE
rm -rf MODEL_TRUE/*

# model setup
./0.setup_models.sh
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# data setup
./1.create_raw_data.sh
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# FWI model inversion
./2.run_inversion.sh
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "done"
echo `date`


