#!/bin/bash

echo "---------------------------------------------------------------------------------"
echo
echo "                          S P E C F E M 3 D   -   Tests"
echo
echo "---------------------------------------------------------------------------------"
echo "This script runs a set of compilation and unit tests in directory tests/"
echo "It may take a few minutes to execute."
echo
echo "Please consider adding more test units to this directory here ..."
echo "Contributions can be sent to: $(tput bold)https://github.com/SPECFEM/specfem3d$(tput sgr0)"
echo

# directories
dir=`pwd`

# changes to subdirectory tests/ if called in root directory SPECFEM3D/
currentdir=`basename $dir`
echo "current directory: $currentdir"
if [ "$currentdir" == "SPECFEM3D" ]; then
cd tests/
dir=`pwd`
fi

# default sub-directories
tests=( compilations \
        decompose_mesh \
        meshfem3D \
        generate_databases \
        specfem3D \
        auxiliaries \
        tomography \
      )

# running tests
echo "main directory: $dir"
echo
date
if [ "$1" != "" ]; then

# specified test directory
echo "test $1 starting"
# runs all bash scripts in specified test-subdirectory
./run_tests.sh $1

# checks exit code
if [[ $? -ne 0 ]]; then
  dir=`basename $testdir`
  echo "ERROR"
  echo "ERROR test failed, please check file results.log in tests/$dir"
  echo "ERROR"
  exit 1
fi

# all test directories
echo
echo "test completed"

else
echo "all tests starting"
# loops over subdirectories
for testdir in ${tests[@]};
do
  testdir=${testdir%*/}

  if [[ "$testdir" == *buildbot* ]]; then
    # skips this test directory
    :
  else
    # runs all bash scripts in test-subdirectory
    ./run_tests.sh $testdir

    # checks exit code
    if [[ $? -ne 0 ]]; then
      dir=`basename $testdir`
      echo "ERROR"
      echo "ERROR test failed, please check file results.log in tests/$dir"
      echo "ERROR"
      exit 1
    fi

    cd $dir/

  fi

done
echo
echo "all tests completed"
fi

echo
date
echo
