#!/bin/bash
#
# runs a test example case
#

# getting updated environment (CUDA_HOME, PATH, ..)
if [ -f $HOME/.tmprc ]; then source $HOME/.tmprc; fi

WORKDIR=`pwd`
dir=${TESTDIR}

# info
echo "work directory: $WORKDIR"
echo
echo "**********************************************************"
echo
echo "configuration test: TESTDIR=${TESTDIR} TESTFLAGS=${TESTFLAGS}"
echo
echo "    test directory: $dir"
echo
echo "**********************************************************"
echo

# bash function for checking seismogram output with reference solutions
my_test(){
  echo "testing seismograms:"
  ln -s $WORKDIR/utils/compare_seismogram_correlations.py
  ./compare_seismogram_correlations.py REF_SEIS/ OUTPUT_FILES/
  if [[ $? -ne 0 ]]; then exit 1; fi
  ./compare_seismogram_correlations.py REF_SEIS/ OUTPUT_FILES/ | grep min/max | cut -d \| -f 3 | awk '{print "correlation:",$1; if ($1 < 0.999 ){print $1,"failed"; exit 1;}else{ print $1,"good"; exit 0;}}'
  if [[ $? -ne 0 ]]; then exit 1; fi
  rm -rf OUTPUT_FILES/
}

# configuration
echo "configuration: default"
./configure FC=gfortran MPIFC=mpif90 CC=gcc ${TESTFLAGS}
if [[ $? -ne 0 ]]; then echo "configuration failed:"; cat config.log; echo ""; echo "exiting..."; exit 1; fi

# we output to console
sed -i "s:IMAIN .*:IMAIN = ISTANDARD_OUTPUT:" setup/constants.h

# layered example w/ NGLL = 6
if [ "$TESTDIR" == "EXAMPLES/layered_halfspace/" ]; then
  sed -i "s:NGLLX =.*:NGLLX = 6:" setup/constants.h
fi

# compilation
echo
echo "compilation:"
make clean; make -j2 all
if [[ $? -ne 0 ]]; then exit 1; fi

# test example
echo
echo "test directory: $dir"
echo
cd $dir

# default setup
# limit time steps for testing
sed -i "s:^NSTEP .*:NSTEP    = 200:" DATA/Par_file
# shortens output interval to avoid timeouts
sed -i "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO    = 50:" DATA/Par_file

# limit time steps for specific examples
# simple mesh example
if [ "$TESTDIR" == "EXAMPLES/meshfem3D_examples/simple_model/" ]; then
  sed -i "s:^NSTEP .*:NSTEP    = 800:" DATA/Par_file
fi
# tpv5 example
if [ "$TESTDIR" == "EXAMPLES/fault_examples/tpv5/" ]; then
  sed -i "s:^NSTEP .*:NSTEP    = 1500:" DATA/Par_file
fi
# layered halfspace example
if [ "$TESTDIR" == "EXAMPLES/layered_halfspace/" ]; then
  sed -i "s:^NSTEP .*:NSTEP    = 200:" DATA/Par_file
fi

# default script
./run_this_example.sh

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# simulation done
echo
echo "simulation done: `pwd`"
echo `date`
echo

# seismogram comparison
my_test

echo
echo "all good"
echo `date`
echo
