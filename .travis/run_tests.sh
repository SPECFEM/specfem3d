#!/bin/bash

# getting updated environment (CUDA_HOME, PATH, ..)
if [ -f $HOME/.tmprc ]; then source $HOME/.tmprc; fi

# checks if anything to do
echo "run checks: $RUN_CHECKS"
if [ "$RUN_CHECKS" == "0" ]; then
  echo "  no run checks required, exiting..."
  exit 0
else
  echo "  run checks required, start testing..."
fi
echo

###########################################################
# setup
###########################################################
# chooses example directory
case "$TESTDIR" in
  0) dir=./ ;;
  1) dir=EXAMPLES/homogeneous_halfspace_HEX8_elastic_absorbing_Stacey_5sides/ ;;
  2) dir=EXAMPLES/homogeneous_acoustic/ ;;
  3) dir=EXAMPLES/homogeneous_poroelastic/ ;;
  4) dir=EXAMPLES/meshfem3D_examples/simple_model/ ;;
  5) dir=EXAMPLES/homogeneous_acoustic/ ;;
  6) dir=EXAMPLES/homogeneous_acoustic/ ;;
  7) dir=EXAMPLES/homogeneous_halfspace/ ;;
  8) dir=EXAMPLES/CPML_examples/homogeneous_halfspace_HEX8_acoustic_absorbing_CPML_5sides/ ;;
  9) dir=EXAMPLES/CPML_examples/homogeneous_halfspace_HEX8_elastic_absorbing_CPML_5sides/ ;;
  10) dir=EXAMPLES/noise_tomography/ ;;
  11) dir=EXAMPLES/tomographic_model/ ;;
  12) dir=EXAMPLES/homogeneous_acoustic/ ;;
  13) dir=EXAMPLES/waterlayered_halfspace/ ;;
  14) dir=EXAMPLES/homogeneous_halfspace_HEX27_elastic_no_absorbing/ ;;
  15) dir=EXAMPLES/homogeneous_acoustic/ ;;
  16) dir=EXAMPLES/fault_examples/tpv5/ ;;
  17) dir=EXAMPLES/meshfem3D_examples/socal1D/ ;;
  18) dir=EXAMPLES/meshfem3D_examples/cavity/ ;;
  19) dir=EXAMPLES/meshfem3D_examples/sep_bathymetry/ ;;
  20) dir=EXAMPLES/small_example_coupling_FK_specfem/ ;;
  21) dir=EXAMPLES/layered_halfspace/ ;;
  22) dir=EXAMPLES/homogeneous_halfspace_HEX8_elastic_no_absorbing/ ;;
  23) dir=EXAMPLES/Gmsh_simple_lddrk/ ;;
  24) dir=EXAMPLES/decompose_mesh_MPI/ ;;
  25) dir=EXAMPLES/meshfem3D_examples/regular_element_mesh/ ;;
  26) dir=EXAMPLES/small_adjoint_multiple_sources/ ;;
  27) dir=EXAMPLES/Gmsh_simple_box_hex27/ ;;
  28) dir=EXAMPLES/waterlayered_poroelastic/ ;;
  29) dir=EXAMPLES/inversion_examples/fwi_test_acoustic/ ;;
  *) dir=EXAMPLES/homogeneous_halfspace/ ;;
esac

# info
#echo $TRAVIS_BUILD_DIR
echo $WORKDIR
echo `date`
echo
echo "**********************************************************"
echo
echo "run test: TESTID=${TESTID} TESTDIR=${TESTDIR} TESTCOV=${TESTCOV} "
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


###########################################################
# test examples
###########################################################

# testing internal mesher example (short & quick for all configuration)
echo 'Tests...' && echo -en 'travis_fold:start:tests\\r'

# runs test
echo "test directory: $dir"
echo

cd $dir

if [ "$TESTID" == "4" ]; then
  # runs default tests
  make tests
elif [ "$TESTID" == "14" ]; then
  # noise example
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  sed -i '10,$ d' NOISE_TOMOGRAPHY/S_squared  # truncates file, deletes all lines from line 10 till end
  ./run_this_example.sh
else
  # default setup
  # limit time steps for testing
  sed -i "s:^NSTEP .*:NSTEP    = 200:" DATA/Par_file
  # shortens output interval to avoid timeouts
  sed -i "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO    = 50:" DATA/Par_file

  # more specific directory setups
  # acoustic examples
  if [ "$TESTDIR" == "2" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 1000:" DATA/Par_file
  fi
  # serial/no mpi examples
  if [ "$TESTDIR" == "5" ]; then
    sed -i "s:^NPROC .*:NPROC    = 1:" DATA/Par_file
    sed -i "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO    = 20:" DATA/Par_file
  fi
  # waterlayered_halfspace example
  if [ "$TESTDIR" == "13" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 600:" DATA/Par_file
  fi
  # tpv5 example
  if [ "$TESTDIR" == "16" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 1200:" DATA/Par_file
  fi
  # socal examples
  if [ "$TESTDIR" == "17" ]; then
    sed -i "s:^NPROC .*:NPROC    = 1:" DATA/Par_file
    sed -i "s:^NSTEP .*:NSTEP    = 840:" DATA/Par_file
  fi
  # SEP example
  if [ "$TESTDIR" == "19" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 1000:" DATA/Par_file
  fi
  # small_adjoint example
  if [ "$TESTDIR" == "26" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 500:" DATA/Par_file
  fi

  # more specific test setups
  # simple mesh example
  if [ "$TESTID" == "8" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 800:" DATA/Par_file
  fi
  # debug+double precision example (acoustic)
  if [ "$TESTID" == "10" ]; then
    sed -i "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO    = 20:" DATA/Par_file
  fi
  # debug example (elastic)
  if [ "$TESTID" == "11" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 15:" DATA/Par_file
    sed -i "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO    = 5:" DATA/Par_file
  fi
  # tomographic model
  if [ "$TESTID" == "15" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 500:" DATA/Par_file
  fi
  # kernel example
  if [ "$TESTID" == "16" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 1000:" DATA/Par_file
  fi
  # elastic example w/ CUDA compilation
  if [ "$TESTID" == "18" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 500:" DATA/Par_file
  fi
  # acoustic example w/ CUDA compilation
  if [ "$TESTID" == "19" ]; then
    sed -i "s:^NPROC .*:NPROC    = 1:" DATA/Par_file
    sed -i "s:^NSTEP .*:NSTEP    = 200:" DATA/Par_file
  fi
  # socal example w/ 1d_socal
  if [ "$TESTID" == "22" ]; then
    sed -i "s:^MODEL .*:MODEL    = 1d_socal:" DATA/Par_file
    rm -f REF_SEIS; ln -s REF_SEIS.1d_socal REF_SEIS
  fi
  # socal example w/ 1d_prem
  if [ "$TESTID" == "23" ]; then
    sed -i "s:^MODEL .*:MODEL    = 1d_prem:" DATA/Par_file
    rm -f REF_SEIS; ln -s REF_SEIS.1d_prem REF_SEIS
  fi
  # socal example w/ 1d_cascadia
  if [ "$TESTID" == "24" ]; then
    sed -i "s:^MODEL .*:MODEL    = 1d_cascadia:" DATA/Par_file
    rm -f REF_SEIS; ln -s REF_SEIS.1d_cascadia REF_SEIS
  fi
  #
  # cavity example
  if [ "$TESTID" == "25" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 2000:" DATA/Par_file
  fi
  # coupled with FK
  if [ "$TESTID" == "27" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 1000:" DATA/Par_file
  fi
  # elastic, no absorbing
  if [ "$TESTID" == "29" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 600:" DATA/Par_file
  fi
  # Gmsh example w/ LDDRK
  if [ "$TESTID" == "30" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 2000:" DATA/Par_file
  fi
  # regular elements example
  if [ "$TESTID" == "32" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 1000:" DATA/Par_file
  fi
  # Gmsh example w/ hex27
  if [ "$TESTID" == "34" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 2000:" DATA/Par_file
  fi
  # waterlayered poroelastic
  if [ "$TESTID" == "35" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 1000:" DATA/Par_file
  fi
  # inversion example
  if [ "$TESTID" == "36" ]; then
    sed -i "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO    = 500:" DATA/Par_file
    sed -i "s/Niter .*/Niter       : 0/" DATA/inverse_problem/inversion_fwi.dat
  fi

  # coverage run
  if [ "$TESTCOV" == "1" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  fi

  if [ "$TESTID" == "16" ]; then
    # kernel script
    ./run_this_example_kernel.sh
    if [[ $? -ne 0 ]]; then exit 1; fi

    # reverse order of seismogram output for comparison
    mv OUTPUT_FILES/DB.X20.MXP.semp tmp
    tac tmp > OUTPUT_FILES/DB.X20.MXP.semp
  else
    # default script
    ./run_this_example.sh
  fi
  if [[ $? -ne 0 ]]; then exit 1; fi

  # simulation done
  echo
  echo "simulation done: `pwd`"
  echo `date`
  echo

  # seismogram comparison
  if [ "$TESTCOV" == "0" ] && [ ! "$TESTID" == "11" ] && [ ! "$TESTID" == "36" ]; then
    my_test
  fi
fi

# checks
if [[ $? -ne 0 ]]; then exit 1; fi

# simulation done
echo
echo "test done: `pwd`"
echo `date`
echo

echo -en 'travis_fold:end:tests\\r'
echo

# code coverage: https://app.codecov.io/gh/SPECFEM/specfem3d/
# additional runs for coverage
#
# note: log becomes too long, trying to fold each test output
cd $WORKDIR

##
## homogeneous halfspace examples
##
#echo 'Coverage...' && echo -en 'travis_fold:start:coverage.elastic-noabs\\r'
#if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "0" ]; then
#  ##
#  ## testing homogeneous halfspace
#  ##
#  echo "##################################################################"
#  echo "EXAMPLES/homogeneous_halfspace_HEX8_elastic_no_absorbing/"
#  echo
#  cd EXAMPLES/homogeneous_halfspace_HEX8_elastic_no_absorbing/
#  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
#  ./run_this_example.sh
#  if [[ $? -ne 0 ]]; then exit 1; fi
#  cd $WORKDIR
#fi
#echo -en 'travis_fold:end:coverage.elastic-noabs\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.elastic-hex27-noabs\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "0" ]; then
  ##
  ## testing hex27 example
  ##
  echo "##################################################################"
  echo "EXAMPLES/homogeneous_halfspace_HEX27_elastic_no_absorbing/"
  echo
  cd EXAMPLES/homogeneous_halfspace_HEX27_elastic_no_absorbing/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.elastic-hex27-noabs\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.poro\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "0" ]; then
  ##
  ## testing poroelastic
  ##
  echo "##################################################################"
  echo "EXAMPLES/homogeneous_poroelastic/"
  echo
  cd EXAMPLES/homogeneous_poroelastic/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.poro\\r'

## kernel example
echo 'Coverage...' && echo -en 'travis_fold:start:coverage.kernel\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "0" ]; then
  ##
  ## testing acoustic kernel simulation
  ##
  echo "##################################################################"
  echo "EXAMPLES/homogeneous_acoustic/"
  echo
  cd EXAMPLES/homogeneous_acoustic/
  cp -v DATA/Par_file DATA/Par_file.org
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  sed -i "s:300:5:" run_this_example_kernel.sh
  sed -i "s:^t_start.*:t_start=-6.0:" create_adjoint_sources.sh
  sed -i "s:^t_end.*:t_end=-5.55:" create_adjoint_sources.sh
  ./run_this_example_kernel.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cp -v DATA/Par_file.org DATA/Par_file
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.kernel\\r'

#echo 'Coverage...' && echo -en 'travis_fold:start:coverage.acoustic\\r'
#if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
#  ##
#  ## testing acoustic
#  ##
#  cd EXAMPLES/homogeneous_acoustic/
#  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
#  ./run_this_example.sh
#  cd $WORKDIR
#fi
#echo -en 'travis_fold:end:coverage.acoustic\\r'

##
## CPMLs
##
echo 'Coverage...' && echo -en 'travis_fold:start:coverage.pmlacoustic\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
  ##
  ## testing PML acoustic
  ##
  echo "##################################################################"
  echo "EXAMPLES/CPML_examples/homogeneous_halfspace_HEX8_acoustic_absorbing_CPML_5sides/"
  echo
  cd EXAMPLES/CPML_examples/homogeneous_halfspace_HEX8_acoustic_absorbing_CPML_5sides/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.pmlacoustic\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.pmlelastic\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
  ##
  ## testing PML elastic
  ##
  echo "##################################################################"
  echo "EXAMPLES/CPML_examples/homogeneous_halfspace_HEX8_elastic_absorbing_CPML_5sides/"
  echo
  cd EXAMPLES/CPML_examples/homogeneous_halfspace_HEX8_elastic_absorbing_CPML_5sides/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.pmlelastic\\r'

##
## noise
##
echo 'Coverage...' && echo -en 'travis_fold:start:coverage.noise\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
  ##
  ## testing noise
  ##
  echo "##################################################################"
  echo "EXAMPLES/noise_tomography/"
  echo
  cd EXAMPLES/noise_tomography/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  sed -i '10,$ d' NOISE_TOMOGRAPHY/S_squared  # truncates file, deletes all lines from line 10 till end
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.noise\\r'

##
## 3d model examples
##
echo 'Coverage...' && echo -en 'travis_fold:start:coverage.tomo\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
  ##
  ## testing tomographic model
  ##
  echo "##################################################################"
  echo "EXAMPLES/tomographic_model/"
  echo
  cd EXAMPLES/tomographic_model/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.tomo\\r'


echo 'Coverage...' && echo -en 'travis_fold:start:coverage.layered\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
  ##
  ## testing layered halfspace
  ##
  echo "##################################################################"
  echo "EXAMPLES/layered_halfspace/"
  echo
  cd EXAMPLES/layered_halfspace/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.layered\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.waterlayer\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
  ##
  ## testing waterlayered serial
  ##
  echo "##################################################################"
  echo "EXAMPLES/waterlayered_halfspace/"
  echo
  cd EXAMPLES/waterlayered_halfspace/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.waterlayer\\r'


##
## meshfem3D examples checks
##
echo 'Coverage...' && echo -en 'travis_fold:start:coverage.meshfem3D-simple\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "2" ]; then
  ##
  ## testing simple model
  ##
  echo "##################################################################"
  echo "EXAMPLES/meshfem3D_examples/simple_model/"
  echo
  cd EXAMPLES/meshfem3D_examples/simple_model/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.meshfem3D-simple\\r'


echo 'Coverage...' && echo -en 'travis_fold:start:coverage.meshfem3D-socal\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "2" ]; then
  ##
  ## testing socal example
  ##
  echo "##################################################################"
  echo "EXAMPLES/meshfem3D_examples/socal1D/"
  echo
  cd EXAMPLES/meshfem3D_examples/socal1D/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.meshfem3D-socal\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.meshfem3D-socal.1d_socal\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "2" ]; then
  ##
  ## testing socal example
  ##
  echo "##################################################################"
  echo "EXAMPLES/meshfem3D_examples/socal1D/"
  echo
  cd EXAMPLES/meshfem3D_examples/socal1D/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  sed -i "s:^MODEL .*:MODEL    = 1d_socal:" DATA/Par_file
  rm -f REF_SEIS; ln -s REF_SEIS.1d_socal REF_SEIS
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.meshfem3D-socal.1d_socal\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.meshfem3D-socal.1d_prem\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "2" ]; then
  ##
  ## testing socal example
  ##
  echo "##################################################################"
  echo "EXAMPLES/meshfem3D_examples/socal1D/"
  echo
  cd EXAMPLES/meshfem3D_examples/socal1D/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  sed -i "s:^MODEL .*:MODEL    = 1d_prem:" DATA/Par_file
  rm -f REF_SEIS; ln -s REF_SEIS.1d_prem REF_SEIS
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.meshfem3D-socal.1d_prem\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.meshfem3D-socal.1d_cascadia\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "2" ]; then
  ##
  ## testing socal example
  ##
  echo "##################################################################"
  echo "EXAMPLES/meshfem3D_examples/socal1D/"
  echo
  cd EXAMPLES/meshfem3D_examples/socal1D/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  sed -i "s:^MODEL .*:MODEL    = 1d_cascadia:" DATA/Par_file
  rm -f REF_SEIS; ln -s REF_SEIS.1d_cascadia REF_SEIS
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.meshfem3D-socal.1d_cascadia\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.meshfem3D-cavity\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "2" ]; then
  ##
  ## testing cavity example
  ##
  echo "##################################################################"
  echo "EXAMPLES/meshfem3D_examples/cavity/"
  echo
  cd EXAMPLES/meshfem3D_examples/cavity/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.meshfem3D-cavity\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.meshfem3D-sep\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "2" ]; then
  ##
  ## testing SEP model
  ##
  echo "##################################################################"
  echo "EXAMPLES/meshfem3D_examples/sep_bathymetry/"
  echo
  cd EXAMPLES/meshfem3D_examples/sep_bathymetry/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.meshfem3D-sep\\r'

## special examples
echo 'Coverage...' && echo -en 'travis_fold:start:coverage.fault\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "2" ]; then
  ##
  ## testing fault example
  ##
  echo "##################################################################"
  echo "EXAMPLES/fault_examples/tpv5/"
  echo
  cd EXAMPLES/fault_examples/tpv5/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.fault\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.coupleFK\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "2" ]; then
  ##
  ## testing coupling FK
  ##
  echo "##################################################################"
  echo "EXAMPLES/small_example_coupling_FK_specfem/"
  echo
  cd EXAMPLES/small_example_coupling_FK_specfem/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.coupleFK\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.Gmsh\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "2" ]; then
  ##
  ## testing Gmsh example
  ##
  echo "##################################################################"
  echo "EXAMPLES/Gmsh_simple_lddrk/"
  echo
  cd EXAMPLES/Gmsh_simple_lddrk/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.Gmsh\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.decompose_mpi\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "2" ]; then
  ##
  ## testing decompose_mesh_MPI example
  ##
  echo "##################################################################"
  echo "EXAMPLES/decompose_mesh_MPI/"
  echo
  cd EXAMPLES/decompose_mesh_MPI/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.decompose_mpi\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.regular_elements\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "2" ]; then
  ##
  ## testing regular elements example
  ##
  echo "##################################################################"
  echo "EXAMPLES/meshfem3D_examples/regular_element_mesh/"
  echo
  cd EXAMPLES/meshfem3D_examples/regular_element_mesh/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.regular_elements\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.small_adjoint\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "2" ]; then
  ##
  ## testing regular elements example
  ##
  echo "##################################################################"
  echo "EXAMPLES/small_adjoint_multiple_sources/"
  echo
  cd EXAMPLES/small_adjoint_multiple_sources/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.small_adjoint\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.waterlayered_poroelastic\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "2" ]; then
  ##
  ## testing regular elements example
  ##
  echo "##################################################################"
  echo "EXAMPLES/waterlayered_poroelastic/"
  echo
  cd EXAMPLES/waterlayered_poroelastic/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.waterlayered_poroelastic\\r'


##
## serial
##
echo 'Coverage...' && echo -en 'travis_fold:start:coverage.serial-elastic\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "3" ]; then
  ##
  ## testing elastic serial
  ##
  echo "##################################################################"
  echo "EXAMPLES/homogeneous_halfspace/"
  echo
  cd EXAMPLES/homogeneous_halfspace/
  cp -v DATA/Par_file DATA/Par_file.org
  sed -i "s:^NPROC .*:NPROC    = 1:" DATA/Par_file
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cp -v DATA/Par_file.org DATA/Par_file
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.serial-elastic\\r'

#echo 'Coverage...' && echo -en 'travis_fold:start:coverage.serial\\r'
#if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "3" ]; then
#  ##
#  ## testing acoustic serial
#  ##
#  cd EXAMPLES/homogeneous_acoustic/
#  cp -v DATA/Par_file DATA/Par_file.org
#  sed -i "s:^NPROC .*:NPROC    = 1:" DATA/Par_file
#  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
#  ./run_this_example.sh
#  cp -v DATA/Par_file.org DATA/Par_file
#  cd $WORKDIR
#fi
#echo -en 'travis_fold:end:coverage.serial\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.serial-meshfem\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "3" ]; then
  ##
  ## testing simple model serial
  ##
  echo "##################################################################"
  echo "EXAMPLES/meshfem3D_examples/simple_model"
  echo
  cd EXAMPLES/meshfem3D_examples/simple_model
  cp -v DATA/Par_file DATA/Par_file.org
  sed -i "s:^NPROC .*:NPROC    = 1:" DATA/Par_file
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  cp -v DATA/meshfem3D_files/Mesh_Par_file DATA/meshfem3D_files/Mesh_Par_file.org
  sed -i "s:^NPROC_XI .*:NPROC_XI    = 1:" DATA/meshfem3D_files/Mesh_Par_file
  sed -i "s:^NPROC_ETA .*:NPROC_ETA    = 1:" DATA/meshfem3D_files/Mesh_Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cp -v DATA/Par_file.org DATA/Par_file
  cp -v DATA/meshfem3D_files/Mesh_Par_file.org DATA/meshfem3D_files/Mesh_Par_file
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.serial-meshfem\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.Gmsh-hex27\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "3" ]; then
  ##
  ## testing Gmsh-hex27 example
  ##
  echo "##################################################################"
  echo "EXAMPLES/Gmsh_simple_box_hex27/"
  echo
  cd EXAMPLES/Gmsh_simple_box_hex27/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.Gmsh-hex27\\r'

# done
echo "all done"
echo `date`
echo

