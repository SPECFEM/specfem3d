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
  # starts timer
  timer_start
  # runs command
  "$@"
  # stops timer
  timer_stop
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

timer_start(){
  t_start=`date +'%s'`
  trap '[ -n "$(jobs -pr)" ] && kill $(jobs -pr)' INT QUIT TERM EXIT
  # progress bar
  while true; do
    echo -n "."
    sleep 30.0
  done &
  timer_PID=$!
}

timer_stop(){
  # saves previous return code
  local exit_code=$?
  # stop background process
  kill $timer_PID
  # get time difference (in s)
  local t_end=`date +'%s'`
  local diff=$(($t_end - $t_start))
  # output nice time format info
  local lapsed_time=$(convertsecs $diff)
  echo -n " ($lapsed_time) "
  # return with previous exit code
  return $exit_code
}

convertsecs() {
  h=$(($1/3600))
  m=$((($1/60)%60))
  s=$(($1%60))
  printf "%02dh %02dm %02ds\n" $h $m $s
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
if [ ! -e "${ROOT}/src/specfem3D/specfem3D.F90" ]; then echo "please check if ROOT set correctly ..."; exit 1; fi

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
    if [[ "$file" == *run_tests.sh* ]] || [[ "$file" == *run_this_example*  ]]; then
      # skips this run script
      :
    else
      step "  processing $file : "
      try ./$file
      next

      # checks exit code
      if [[ $? -ne 0 ]]; then echo "***** results.log ******"; cat results.log; exit 1; fi

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
