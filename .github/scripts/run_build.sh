#!/bin/bash
#
# builds all executables
#

# getting updated environment (CUDA_HOME, PATH, ..)
if [ -f $HOME/.tmprc ]; then source $HOME/.tmprc; fi

WORKDIR=`pwd`
TESTCOV=${TESTCOV:-}

# info
echo "work directory: $WORKDIR"
echo `date`
echo
echo "**********************************************************"
echo
echo "configuration test: TESTFLAGS=${TESTFLAGS} TESTNGLL=${TESTNGLL} TESTCOV=${TESTCOV}"
echo
echo "**********************************************************"
echo

# compiler infos
echo "compiler versions:"
echo "gcc --version"
gcc --version
echo "gfortran --version"
gfortran --version
echo "mpif90 --version"
mpif90 --version
echo

## ADIOS2
if [ "${ADIOS2}" == "true" ]; then
  echo
  echo "enabling ADIOS2"
  echo
  ADIOS2_CONFIG="${ADIOS2_DIR}/bin/adios2-config"
  adios=(--with-adios2 ADIOS2_CONFIG="$ADIOS2_CONFIG" )
else
  adios=()
fi

## HDF5
if [ "${HDF5}" == "true" ]; then
  echo
  echo "enabling HDF5"
  echo
  hdf=(--with-hdf5 HDF5_INC="/usr/include/hdf5/openmpi/" HDF5_LIBS="-L/usr/lib/x86_64-linux-gnu/hdf5/openmpi")
else
  hdf=()
fi

## HIP
if [ "${HIP}" == "true" ]; then
  echo
  echo "enabling HIP"
  echo
  hip=(--with-hip HIPCC=g++ HIP_FLAGS="-O2 -g -std=c++17" HIP_PLATFORM=cpu HIP_INC=./external_libs/ROCm-HIP-CPU/include HIP_LIBS="-ltbb -lpthread -lstdc++")
else
  hip=()
fi

## special testflags
if [ "${TESTFLAGS}" == "check-mcmodel-medium" ]; then
  # note: this is a work-around as using the 'env:' parameter in the workflow 'CI.yml' with TESTFLAGS: FLAGS_CHECK=".."
  #       won't work as the FLAGS_CHECK string will then get split up and ./configure .. complains about unknown parameters.
  #       here, we re-define TESTFLAGS with a single quote around FLAGS_CHECK=".." to avoid the splitting.
  # use FLAGS_CHECK
  flags=(FLAGS_CHECK="-O3 -mcmodel=medium -std=f2008 -Wall -Wno-do-subscript -Wno-conversion -Wno-maybe-uninitialized")
  TESTFLAGS=""  # reset
else
  flags=()
fi

# configuration
echo
echo "configuration:"
echo

# split TESTFLAGS into individual items
set -- ${TESTFLAGS}

###########################################################
# configuration & compilation
###########################################################
# configuration

if [ "${TESTCOV}" == "true" ]; then
  echo "configuration: for coverage"
  ./configure \
    "${adios[@]}" \
    "${hdf[@]}" \
    "${hip[@]}" \
    "${flags[@]}" \
    FLAGS_CHECK="-fprofile-arcs -ftest-coverage -O0" CFLAGS="-coverage -O0"
    FC=gfortran MPIFC=mpif90 CC=gcc "$@"
else
  if [ "$CUDA" == "true" ]; then
    echo "configuration: for cuda"
    ./configure \
      "${adios[@]}" \
      "${hdf[@]}" \
      "${hip[@]}" \
      "${flags[@]}" \
      CUDA_LIB="${CUDA_HOME}/lib64" CUDA_INC="${CUDA_HOME}/include" \
      CUDA_FLAGS="-Xcompiler -Wall,-Wno-unused-function,-Wno-unused-const-variable,-Wfatal-errors -g -G" \
      FC=${FC} MPIFC=${MPIFC} CC=${CC} "$@"
  else
    echo "configuration: default"
    ./configure \
      "${adios[@]}" \
      "${hdf[@]}" \
      "${hip[@]}" \
      "${flags[@]}" \
      FC=gfortran MPIFC=mpif90 CC=gcc "$@"
  fi
fi

# checks
if [[ $? -ne 0 ]]; then echo "configuration failed:"; cat config.log; echo ""; echo "exiting..."; exit 1; fi

# layered example w/ NGLL = 6
if [ "$TESTNGLL" == "6" ]; then
  sed -i "s:NGLLX =.*:NGLLX = 6:" setup/constants.h
fi

# we output to console
sed -i "s:IMAIN .*:IMAIN = ISTANDARD_OUTPUT:" setup/constants.h

# inversion example
if [ "$TESTDIR" == "EXAMPLES/applications/inversion_examples/fwi_test_acoustic/" ]; then
  sed -i "s:IMAIN .*:IMAIN = 42:" setup/constants.h
  sed -i "s:INVERSE_LOG_FILE .*:INVERSE_LOG_FILE = 6:" src/inverse_problem_for_model/inverse_problem_par.f90
fi

# compilation
echo
echo "clean"
echo
make clean

echo
echo "compilation"
echo
make -j4 all

# checks
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "done "
echo `date`
echo
