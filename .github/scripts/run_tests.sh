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
echo `date`
echo
echo "**********************************************************"
echo
echo "test directory: $dir"
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
}

# test example
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
# small adjoint example
if [ "$TESTDIR" == "EXAMPLES/small_adjoint_multiple_sources/" ]; then
  # full length
  sed -i "s:^NSTEP .*:NSTEP    = 1000:" DATA/Par_file
fi
# socal examples
if [ "$TESTDIR" == "EXAMPLES/meshfem3D_examples/socal1D/" ]; then
  # full length
  sed -i "s:^NSTEP .*:NSTEP    = 840:" DATA/Par_file
  # model setup
  if [ "$TESTID" == "1" ]; then
    # 1D_socal
    sed -i "s:^MODEL .*:MODEL    = 1d_socal:" DATA/Par_file
    rm -f REF_SEIS; ln -s REF_SEIS.1d_socal REF_SEIS
  elif [ "$TESTID" == "2" ]; then
    # 1D_prem
    sed -i "s:^MODEL .*:MODEL    = 1d_prem:" DATA/Par_file
    rm -f REF_SEIS; ln -s REF_SEIS.1d_prem REF_SEIS
  elif [ "$TESTID" == "3" ]; then
    # 1D_cascadia
    sed -i "s:^MODEL .*:MODEL    = 1d_cascadia:" DATA/Par_file
    rm -f REF_SEIS; ln -s REF_SEIS.1d_cascadia REF_SEIS
  else
    # default
    continue
  fi
fi
# coupling FK
if [ "$TESTDIR" == "EXAMPLES/small_example_coupling_FK_specfem/" ]; then
  sed -i "s:^NSTEP .*:NSTEP    = 1000:" DATA/Par_file
fi
# elastic halfspace, no absorbing
if [ "$TESTDIR" == "EXAMPLES/homogeneous_halfspace_HEX8_elastic_no_absorbing/" ]; then
  sed -i "s:^NSTEP .*:NSTEP    = 600:" DATA/Par_file
fi
# waterlayered_halfspace example
if [ "$TESTDIR" == "EXAMPLES/waterlayered_halfspace/" ]; then
  sed -i "s:^NSTEP .*:NSTEP    = 600:" DATA/Par_file
fi
# tomographic model
if [ "$TESTDIR" == "EXAMPLES/tomographic_model/" ]; then
  sed -i "s:^NSTEP .*:NSTEP    = 500:" DATA/Par_file
fi
# cavity example
if [ "$TESTDIR" == "EXAMPLES/meshfem3D_examples/cavity/" ]; then
  sed -i "s:^NSTEP .*:NSTEP    = 2000:" DATA/Par_file
fi
# SEP example
if [ "$TESTDIR" == "EXAMPLES/meshfem3D_examples/sep_bathymetry/" ]; then
  sed -i "s:^NSTEP .*:NSTEP    = 1000:" DATA/Par_file
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
if [ "${DEBUG}" == "true" ]; then
  # no comparisons
  continue
else
  my_test
fi

# cleanup
rm -rf OUTPUT_FILES/ DATABASES_MPI/

echo
echo "all good"
echo `date`
echo
