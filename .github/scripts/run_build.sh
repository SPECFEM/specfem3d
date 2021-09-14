#!/bin/bash
#
# builds all executables
#

# getting updated environment (CUDA_HOME, PATH, ..)
if [ -f $HOME/.tmprc ]; then source $HOME/.tmprc; fi

WORKDIR=`pwd`

# info
echo "work directory: $WORKDIR"
echo `date`
echo
echo "**********************************************************"
echo
echo "configuration test: TESTFLAGS=${TESTFLAGS} TESTNGLL=${TESTNGLL}"
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


# configuration
echo
echo "configuration: default"
echo

./configure FC=gfortran MPIFC=mpif90 CC=gcc ${TESTFLAGS}

# checks
if [[ $? -ne 0 ]]; then echo "configuration failed:"; cat config.log; echo ""; echo "exiting..."; exit 1; fi

# layered example w/ NGLL = 6
if [ "$TESTNGLL" == "6" ]; then
  sed -i "s:NGLLX =.*:NGLLX = 6:" setup/constants.h
fi

# we output to console
sed -i "s:IMAIN .*:IMAIN = ISTANDARD_OUTPUT:" setup/constants.h

# compilation
echo
echo "compilation:"
make clean; make -j2 all

# checks
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "done "
echo `date`
echo
