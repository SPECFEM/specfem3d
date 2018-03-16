#!/bin/bash

# getting updated environment (CUDA_HOME, PATH, ..)
if [ -f $HOME/.tmprc ]; then source $HOME/.tmprc; fi

###########################################################
# setup
###########################################################
# chooses example directory
case "$TESTMAKE" in
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
  18) dir=EXAMPLES/meshfem3D_examples/socal1D/ ;;
  19) dir=EXAMPLES/meshfem3D_examples/socal1D/ ;;
  20) dir=EXAMPLES/meshfem3D_examples/socal1D/ ;;
  21) dir=EXAMPLES/meshfem3D_examples/cavity/ ;;
  22) dir=EXAMPLES/meshfem3D_examples/sep_bathymetry/ ;;
  23) dir=EXAMPLES/small_example_coupling_FK_specfem/ ;;
  24) dir=EXAMPLES/layered_halfspace/ ;;
  25) dir=EXAMPLES/homogeneous_halfspace_HEX8_elastic_no_absorbing/ ;;
  26) dir=EXAMPLES/Gmsh_simple_lddrk/ ;;
  *) dir=EXAMPLES/homogeneous_halfspace/ ;;
esac

# info
echo $TRAVIS_BUILD_DIR
echo $WORKDIR
echo
echo "**********************************************************"
echo
echo "configuration test: TESTMAKE=${TESTMAKE} TEST=${TEST} FLAGS=${TESTFLAGS}"
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
  ./compare_seismogram_correlations.py REF_SEIS/ OUTPUT_FILES/ | grep min/max | cut -d \| -f 3 | awk '{print "correlation:",$1; if ($1 < 0.9 ){print $1,"failed"; exit 1;}else{ print $1,"good"; exit 0;}}'
  if [[ $? -ne 0 ]]; then exit 1; fi
  rm -rf OUTPUT_FILES/
}

###########################################################
# configuration & compilation
###########################################################
# configuration
echo 'Configure...' && echo -en 'travis_fold:start:configure\\r'
echo "configuration:"
if [ "$TESTCOV" == "1" ]; then
  echo "configuration: for coverage"
  ./configure FC=${FC} MPIFC=${MPIFC} CC=${CC} ${TEST} FLAGS_CHECK="-fprofile-arcs -ftest-coverage -O0" CFLAGS="-coverage -O0"
else
  if [ "$CUDA" == "true" ]; then
    echo "configuration: for cuda"
    ./configure FC=${FC} MPIFC=${MPIFC} CC=${CC} ${TEST} CUDA_LIB="${CUDA_HOME}/lib64" CUDA_INC="${CUDA_HOME}/include" CUDA_FLAGS="-Xcompiler -Wall,-Wno-unused-function,-Wno-unused-const-variable,-Wfatal-errors -g -G"
  else
    echo "configuration: default"
    ./configure FC=${FC} MPIFC=${MPIFC} CC=${CC} ${TEST}
  fi
  # we output to console
  sed -i "s:IMAIN .*:IMAIN = ISTANDARD_OUTPUT:" setup/constants.h
fi
# layered example
if [ "$TESTMAKE" == "24" ]; then
  sed -i "s:NGLLX =.*:NGLLX = 6:" setup/constants.h
fi
echo -en 'travis_fold:end:configure\\r'

# compilation
echo 'Build...' && echo -en 'travis_fold:start:build\\r'
echo "compilation:"
make clean; make -j2 all
echo -en 'travis_fold:end:build\\r'

###########################################################
# test examples
###########################################################
# testing internal mesher example (short & quick for all configuration)
echo 'Tests...' && echo -en 'travis_fold:start:tests\\r'
# runs test
echo "test directory: $dir"
echo
cd $dir
if [ "$TESTMAKE" == "0" ]; then
  # runs default tests
  make tests
elif [ "$TESTMAKE" == "10" ]; then
  # noise example
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file_step1
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file_step2
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file_step3
  ./run_this_example.sh
else
  # limit time steps for testing
  sed -i "s:^NSTEP .*:NSTEP    = 200:" DATA/Par_file
  # shortens output interval to avoid timeouts
  sed -i "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO    = 100:" DATA/Par_file
  #
  # acoustic example
  if [ "$TESTMAKE" == "2" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 1000:" DATA/Par_file
  fi
  # simple mesh example
  if [ "$TESTMAKE" == "4" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 400:" DATA/Par_file
  fi
  # serial example
  if [ "$TESTMAKE" == "5" ]; then
    sed -i "s:^NPROC .*:NPROC    = 1:" DATA/Par_file
    sed -i "s:^NSTEP .*:NSTEP    = 200:" DATA/Par_file
    sed -i "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO    = 20:" DATA/Par_file
  fi
  # debug+double precision example (acoustic)
  if [ "$TESTMAKE" == "6" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 200:" DATA/Par_file
    sed -i "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO    = 20:" DATA/Par_file
  fi
  # debug example (elastic)
  if [ "$TESTMAKE" == "7" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 15:" DATA/Par_file
    sed -i "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO    = 5:" DATA/Par_file
  fi
  # PML example (acoustic)
  if [ "$TESTMAKE" == "8" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 200:" DATA/Par_file
    sed -i "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO    = 50:" DATA/Par_file
  fi
  # tomographic model
  if [ "$TESTMAKE" == "11" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 500:" DATA/Par_file
    sed -i "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO    = 50:" DATA/Par_file
  fi
  # kernel example
  if [ "$TESTMAKE" == "12" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 1000:" DATA/Par_file
  fi
  # elastic example w/ CUDA compilation
  if [ "$TESTMAKE" == "14" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 500:" DATA/Par_file
  fi
  # acoustic example w/ CUDA compilation
  if [ "$TESTMAKE" == "15" ]; then
    sed -i "s:^NPROC .*:NPROC    = 1:" DATA/Par_file
    sed -i "s:^NSTEP .*:NSTEP    = 500:" DATA/Par_file
  fi
  # socal example
  if [ "$TESTMAKE" == "17" ]; then
    sed -i "s:^NPROC .*:NPROC    = 1:" DATA/Par_file
    sed -i "s:^NSTEP .*:NSTEP    = 840:" DATA/Par_file
  fi
  # socal example w/ 1d_socal
  if [ "$TESTMAKE" == "18" ]; then
    sed -i "s:^NPROC .*:NPROC    = 1:" DATA/Par_file
    sed -i "s:^NSTEP .*:NSTEP    = 840:" DATA/Par_file
    sed -i "s:^MODEL .*:MODEL    = 1d_socal:" DATA/Par_file
    rm -f REF_SEIS; ln -s REF_SEIS.1d_socal REF_SEIS
  fi
  # socal example w/ 1d_prem
  if [ "$TESTMAKE" == "19" ]; then
    sed -i "s:^NPROC .*:NPROC    = 1:" DATA/Par_file
    sed -i "s:^NSTEP .*:NSTEP    = 840:" DATA/Par_file
    sed -i "s:^MODEL .*:MODEL    = 1d_prem:" DATA/Par_file
    rm -f REF_SEIS; ln -s REF_SEIS.1d_prem REF_SEIS
  fi
  # socal example w/ 1d_cascadia
  if [ "$TESTMAKE" == "20" ]; then
    sed -i "s:^NPROC .*:NPROC    = 1:" DATA/Par_file
    sed -i "s:^NSTEP .*:NSTEP    = 840:" DATA/Par_file
    sed -i "s:^MODEL .*:MODEL    = 1d_cascadia:" DATA/Par_file
    rm -f REF_SEIS; ln -s REF_SEIS.1d_cascadia REF_SEIS
  fi
  # cavity example
  if [ "$TESTMAKE" == "21" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 2000:" DATA/Par_file
  fi
  # SEP example
  if [ "$TESTMAKE" == "22" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 500:" DATA/Par_file
  fi
  # coupled with FK
  if [ "$TESTMAKE" == "23" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 500:" DATA/Par_file
  fi
  # elastic, no absorbing
  if [ "$TESTMAKE" == "25" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 600:" DATA/Par_file
  fi
  # Gmsh example w/ LDDRK
  if [ "$TESTMAKE" == "26" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 2000:" DATA/Par_file
  fi

  # coverage run
  if [ "$TESTCOV" == "1" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  fi

  if [ "$TESTMAKE" == "12" ]; then
    # kernel script
    ./run_this_example_kernel.sh

    # reverse order of seismogram output for comparison
    mv OUTPUT_FILES/DB.X20.MXP.semp tmp
    tac tmp > OUTPUT_FILES/DB.X20.MXP.semp
  else
    # default script
    ./run_this_example.sh
  fi

  # seismogram comparison
  if [ "$TESTCOV" == "0" ]; then
    my_test
  fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:tests\\r'


# code coverage: https://codecov.io/gh/geodynamics/specfem3d/
# additional runs for coverage
# note: log becomes too long, trying to fold each test output

##
## homogeneous halfspace examples
##
echo 'Coverage...' && echo -en 'travis_fold:start:coverage.elastic-noabs\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTMAKE" == "1" ]; then
  ##
  ## testing socal example
  ##
  cd EXAMPLES/homogeneous_halfspace_HEX8_elastic_no_absorbing/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.elastic-noabs\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.poro\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTMAKE" == "1" ]; then
  ##
  ## testing poroelastic
  ##
  cd EXAMPLES/homogeneous_poroelastic/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.poro\\r'

## kernel example
echo 'Coverage...' && echo -en 'travis_fold:start:coverage.kernel\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTMAKE" == "1" ]; then
  ##
  ## testing acoustic kernel simulation
  ##
  cd EXAMPLES/homogeneous_acoustic/
  cp -v DATA/Par_file DATA/Par_file.org
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  sed -i "s:^t_start.*:t_start=-6.0:" create_adjoint_sources.sh
  sed -i "s:^t_end.*:t_end=-5.55:" create_adjoint_sources.sh
  ./run_this_example_kernel.sh
  cp -v DATA/Par_file.org DATA/Par_file
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.kernel\\r'

#echo 'Coverage...' && echo -en 'travis_fold:start:coverage.acoustic\\r'
#if [ "$TESTCOV" == "1" ] && [ "$TESTMAKE" == "1" ]; then
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
if [ "$TESTCOV" == "1" ] && [ "$TESTMAKE" == "1" ]; then
  ##
  ## testing PML acoustic
  ##
  cd EXAMPLES/CPML_examples/homogeneous_halfspace_HEX8_acoustic_absorbing_CPML_5sides/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.pmlacoustic\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.pmlelastic\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTMAKE" == "1" ]; then
  ##
  ## testing PML elastic
  ##
  cd EXAMPLES/CPML_examples/homogeneous_halfspace_HEX8_elastic_absorbing_CPML_5sides/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.pmlelastic\\r'

##
## noise
##
echo 'Coverage...' && echo -en 'travis_fold:start:coverage.noise\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTMAKE" == "1" ]; then
  ##
  ## testing noise
  ##
  cd EXAMPLES/noise_tomography/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file_step1
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file_step2
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file_step3
  ./run_this_example.sh
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.noise\\r'

##
## 3d model examples
##
echo 'Coverage...' && echo -en 'travis_fold:start:coverage.tomo\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTMAKE" == "1" ]; then
  ##
  ## testing tomographic model
  ##
  cd EXAMPLES/tomographic_model/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.tomo\\r'


echo 'Coverage...' && echo -en 'travis_fold:start:coverage.layered\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTMAKE" == "1" ]; then
  ##
  ## testing socal example
  ##
  cd EXAMPLES/layered_halfspace/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.layered\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.waterlayer\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTMAKE" == "1" ]; then
  ##
  ## testing waterlayered serial
  ##
  cd EXAMPLES/waterlayered_halfspace/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.waterlayer\\r'


##
## meshfem3D examples checks
##
echo 'Coverage...' && echo -en 'travis_fold:start:coverage.meshfem3D-simple\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTMAKE" == "2" ]; then
  ##
  ## testing simple model
  ##
  cd EXAMPLES/meshfem3D_examples/simple_model/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.meshfem3D-simple\\r'


echo 'Coverage...' && echo -en 'travis_fold:start:coverage.meshfem3D-socal\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTMAKE" == "2" ]; then
  ##
  ## testing socal example
  ##
  cd EXAMPLES/meshfem3D_examples/socal1D/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.meshfem3D-socal\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.meshfem3D-socal.1d_socal\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTMAKE" == "2" ]; then
  ##
  ## testing socal example
  ##
  cd EXAMPLES/meshfem3D_examples/socal1D/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  sed -i "s:^MODEL .*:MODEL    = 1d_socal:" DATA/Par_file
  rm -f REF_SEIS; ln -s REF_SEIS.1d_socal REF_SEIS
  ./run_this_example.sh
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.meshfem3D-socal.1d_socal\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.meshfem3D-socal.1d_prem\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTMAKE" == "2" ]; then
  ##
  ## testing socal example
  ##
  cd EXAMPLES/meshfem3D_examples/socal1D/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  sed -i "s:^MODEL .*:MODEL    = 1d_prem:" DATA/Par_file
  rm -f REF_SEIS; ln -s REF_SEIS.1d_prem REF_SEIS
  ./run_this_example.sh
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.meshfem3D-socal.1d_prem\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.meshfem3D-socal.1d_cascadia\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTMAKE" == "2" ]; then
  ##
  ## testing socal example
  ##
  cd EXAMPLES/meshfem3D_examples/socal1D/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  sed -i "s:^MODEL .*:MODEL    = 1d_cascadia:" DATA/Par_file
  rm -f REF_SEIS; ln -s REF_SEIS.1d_cascadia REF_SEIS
  ./run_this_example.sh
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.meshfem3D-socal.1d_cascadia\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.meshfem3D-cavity\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTMAKE" == "2" ]; then
  ##
  ## testing socal example
  ##
  cd EXAMPLES/meshfem3D_examples/cavity/
  sed -i "s:^NSTEP .*:NSTEP    = 200:" DATA/Par_file
  ./run_this_example.sh
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.meshfem3D-cavity\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.meshfem3D-sep\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTMAKE" == "2" ]; then
  ##
  ## testing socal example
  ##
  cd EXAMPLES/meshfem3D_examples/sep_bathymetry/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.meshfem3D-sep\\r'

## special examples
echo 'Coverage...' && echo -en 'travis_fold:start:coverage.fault\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTMAKE" == "2" ]; then
  ##
  ## testing fault example
  ##
  cd EXAMPLES/fault_examples/tpv5/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.fault\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.coupleFK\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTMAKE" == "2" ]; then
  ##
  ## testing socal example
  ##
  cd EXAMPLES/small_example_coupling_FK_specfem/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.coupleFK\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.Gmsh\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTMAKE" == "2" ]; then
  ##
  ## testing Gmsh example
  ##
  cd EXAMPLES/Gmsh_simple_lddrk/
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.Gmsh\\r'



##
## serial
##
echo 'Coverage...' && echo -en 'travis_fold:start:coverage.serial-elastic\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTMAKE" == "5" ]; then
  ##
  ## testing elastic serial
  ##
  cd EXAMPLES/homogeneous_halfspace/
  cp -v DATA/Par_file DATA/Par_file.org
  sed -i "s:^NPROC .*:NPROC    = 1:" DATA/Par_file
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  ./run_this_example.sh
  cp -v DATA/Par_file.org DATA/Par_file
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.serial-elastic\\r'

#echo 'Coverage...' && echo -en 'travis_fold:start:coverage.serial\\r'
#if [ "$TESTCOV" == "1" ] && [ "$TESTMAKE" == "1" ]; then
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
if [ "$TESTCOV" == "1" ] && [ "$TESTMAKE" == "5" ]; then
  ##
  ## testing simple model serial
  ##
  cd EXAMPLES/meshfem3D_examples/simple_model
  cp -v DATA/Par_file DATA/Par_file.org
  sed -i "s:^NPROC .*:NPROC    = 1:" DATA/Par_file
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
  cp -v DATA/meshfem3D_files/Mesh_Par_file DATA/meshfem3D_files/Mesh_Par_file.org
  sed -i "s:^NPROC_XI .*:NPROC_XI    = 1:" DATA/meshfem3D_files/Mesh_Par_file
  sed -i "s:^NPROC_ETA .*:NPROC_ETA    = 1:" DATA/meshfem3D_files/Mesh_Par_file
  ./run_this_example.sh
  cp -v DATA/Par_file.org DATA/Par_file
  cp -v DATA/meshfem3D_files/Mesh_Par_file.org DATA/meshfem3D_files/Mesh_Par_file
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.serial-meshfem\\r'


# done
echo "done `pwd`"

