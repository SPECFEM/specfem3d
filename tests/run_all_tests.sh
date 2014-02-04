#!/bin/bash

# directories
dir=`pwd`

# changes to subdirectory tests/ if called in root directory SPECFEM3D/
currentdir=`basename $dir`
echo "current directory: $currentdir"
if [ "$currentdir" == "SPECFEM3D" ]; then
cd tests/
dir=`pwd`
fi

# running tests
echo "all tests"
echo "`date`"
echo  
echo "directory: $dir"
echo 

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
echo "`date`"
echo
