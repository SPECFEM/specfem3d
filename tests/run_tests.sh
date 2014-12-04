#!/bin/bash
#############################################################

# USER parameters

# relative location of SPECFEM3D root directory for tests directories (e.g. SPECFEM3D/tests/compilations)
export ROOT=../../

# test directory
testdir=$1

#############################################################

# helper functions
step() {
  # output
  echo -n "$@"
  STEP_OK=0
}

try() {
  # runs command
  "$@"
  # Check if command failed and update $STEP_OK if so.
  local EXIT_CODE=$?
  if [[ $EXIT_CODE -ne 0 ]]; then
    STEP_OK=$EXIT_CODE
    #echo "Command \`$*' failed with exit code $EXIT_CODE."
  fi
  return $EXIT_CODE
}

next() {
  [[ $STEP_OK -eq 0 ]]  && echo "[  $(tput setaf 2)OK$(tput sgr0)  ]" || echo "[$(tput setaf 1)FAILED$(tput sgr0)]"
  #echo
  return $STEP_OK
}

#############################################################

# checks if argument given
if [ "$testdir" == "" ]; then echo "usage: ./run_tests.sh testdir[e.g.=compilations]"; exit 1; fi

# changes to subdirectory tests/ if called in root directory SPECFEM3D/
dir=`pwd`
currentdir=`basename $dir`
if [ "$currentdir" == "SPECFEM3D" ]; then
cd tests/
fi

# checks if test directory exists
if [ ! -e $testdir  ]; then
  #see if argument is e.g. tests/compilations
  dir=`basename $testdir`
  if [ ! -e $dir  ]; then
    echo "test directory given does not exists: $testdir, please check..."
    exit 1
  else
    testdir=$dir
  fi
fi

# running tests
cd $testdir/

#checks if ROOT valid
if [ ! -e "${ROOT}/src/specfem3D/specfem3D.f90" ]; then echo "please check if ROOT set correctly ..."; exit 1; fi

# user output
echo "-------------------------------------"
echo "tests in directory : $testdir "
echo "-------------------------------------"
date

# test names
title="tests in $testdir"

echo "$title" > results.log
date >> results.log
echo >> results.log

files=`echo ./*.sh`
if [ "$files" == "./*.sh" ]; then
  echo "  directory contains no test-scripts"
else

  for file in ./*.sh
  do
    if [[ "$file" == *run_tests.sh* ]]; then
      # skips this run script
      : 
    else
      step "  processing $file : "
      try ./$file
      next

      # checks exit code
      if [[ $? -ne 0 ]]; then exit 1; fi

      echo >> results.log
    fi
  done

fi

echo >> results.log
echo "tests completed" >> results.log
date >> results.log

cd ../


#return
echo
exit 0
