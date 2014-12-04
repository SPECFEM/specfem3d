#!/bin/bash

echo "---------------------------------------------------------------------------------"
echo
echo "                          S P E C F E M 3 D   -   Tests"
echo
echo "---------------------------------------------------------------------------------"
echo "This script is running a set of compilation and unit tests in directory tests/"
echo "It will take a short while to execute, grab a coffee and enjoy ;)"
echo
echo "Please consider adding more test units to this directory here ..."
echo "Contributions can be sent to: $(tput bold)http://github.com/geodynamics/specfem3d$(tput sgr0)"
echo

# directories
dir=`pwd`

# changes to subdirectory tests/ if called in root directory SPECFEM3D/
currentdir=`basename $dir`
#echo "current directory: $currentdir"
if [ "$currentdir" == "SPECFEM3D" ]; then
cd tests/
dir=`pwd`
fi

# running tests
echo "main directory: $dir"
echo "all tests starting: `date`"

# loops over subdirectories
for testdir in ./*/
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
date
echo
