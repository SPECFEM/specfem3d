#!/bin/bash
#
# builds all executables
#

# getting updated environment (CUDA_HOME, PATH, ..)
if [ -f $HOME/.tmprc ]; then source $HOME/.tmprc; fi

# checks if anything to do
echo "run checks: $RUN_CHECKS"
if [ "$RUN_CHECKS" == "0" ]; then
  echo "  no run checks required, exiting..."
  exit 0
else
  echo "  run checks required, start building..."
fi
echo

# info
echo $WORKDIR
echo
echo "**********************************************************"
echo
echo "configuration test: TESTID=${TESTID} TESTDIR=${TESTDIR} TESTCOV=${TESTCOV} TESTFLAGS=${TESTFLAGS}"
echo
echo "**********************************************************"
echo

echo "compiler versions:"
echo "gcc --version"
gcc --version
echo "gfortran --version"
gfortran --version
echo "mpif90 --version"
mpif90 --version
if [ "$CUDA" == "true" ]; then
  echo "nvcc --version"
  nvcc --version
fi
echo ""

###########################################################
# configuration & compilation
###########################################################
# configuration
echo 'Configure...' && echo -en 'travis_fold:start:configure\\r'
echo "configuration:"

if [ "$TESTCOV" == "1" ]; then
  echo "configuration: for coverage"
  ./configure FC=${FC} MPIFC=${MPIFC} CC=${CC} ${TESTFLAGS} FLAGS_CHECK="-fprofile-arcs -ftest-coverage -O0" CFLAGS="-coverage -O0"
else
  if [ "$CUDA" == "true" ]; then
    echo "configuration: for cuda"
    ./configure FC=${FC} MPIFC=${MPIFC} CC=${CC} ${TESTFLAGS} CUDA_LIB="${CUDA_HOME}/lib64" CUDA_INC="${CUDA_HOME}/include" CUDA_FLAGS="-Xcompiler -Wall,-Wno-unused-function,-Wno-unused-const-variable,-Wfatal-errors -g -G"
  else
    echo "configuration: default"
    ./configure FC=${FC} MPIFC=${MPIFC} CC=${CC} ${TESTFLAGS}
  fi
fi
if [[ $? -ne 0 ]]; then echo "configuration failed:"; cat config.log; echo ""; echo "exiting..."; exit 1; fi

# we output to console
sed -i "s:IMAIN .*:IMAIN = ISTANDARD_OUTPUT:" setup/constants.h

# inversion example
if [ "$TESTID" == "36" ]; then
  sed -i "s:IMAIN .*:IMAIN = 42:" setup/constants.h
fi

# layered example w/ NGLL = 6
if [ "$TESTID" == "28" ]; then
  sed -i "s:NGLLX =.*:NGLLX = 6:" setup/constants.h
fi

echo -en 'travis_fold:end:configure\\r'

# compilation
echo 'Build...' && echo -en 'travis_fold:start:build\\r'
echo "compilation:"

make clean; make -j2 all

if [[ $? -ne 0 ]]; then exit 1; fi
echo -en 'travis_fold:end:build\\r'

echo
echo "done "
echo `date`
echo
