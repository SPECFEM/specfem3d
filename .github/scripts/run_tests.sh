#!/bin/bash
#
# runs a test example case
#

# getting updated environment (CUDA_HOME, PATH, ..)
if [ -f $HOME/.tmprc ]; then source $HOME/.tmprc; fi

WORKDIR=`pwd`
dir=${TESTDIR}
TESTID=${TESTID:-}
TESTCOV=${TESTCOV:-}

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
  echo "######################################################################################################################"
  echo "testing seismograms"
  ln -s $WORKDIR/utils/scripts/compare_seismogram_correlations.py
  ./compare_seismogram_correlations.py REF_SEIS/ OUTPUT_FILES/
  if [[ $? -ne 0 ]]; then exit 1; fi
  ./compare_seismogram_correlations.py REF_SEIS/ OUTPUT_FILES/ | grep min/max | cut -d \| -f 3 | awk '{print "correlation:",$1; if ($1 < 0.999 ){print $1,"failed"; exit 1;}else{ print $1,"good"; exit 0;}}'
  if [[ $? -ne 0 ]]; then exit 1; fi
  echo "######################################################################################################################"
}

my_kernel_test(){
  # kernel value test - checks rho/kappa/mu kernel value outputs
  echo "######################################################################################################################"
  echo "testing kernel values"
  file_ref=REF_KERNEL/output_solver.txt
  file_out=output.log        # captures the OUTPUT_FILES/output_solver.txt when running solver since IMAIN was set to standard out
  if [ ! -e $file_ref ]; then echo "Please check if file $file_ref exists..."; ls -alR ./; exit 1; fi
  if [ ! -e $file_out ]; then echo "Please check if file $file_out exists..."; ls -alR ./; exit 1; fi
  # gets reference expected kernel values from REF_KERNEL/ folder
  RHO=`grep -E 'maximum value of rho[[:space:]]+kernel' $file_ref | cut -d = -f 2 | tr -d ' '`
  KAPPA=`grep -E 'maximum value of kappa[[:space:]]+kernel' $file_ref | cut -d = -f 2 | tr -d ' '`
  MU=`grep -E 'maximum value of mu[[:space:]]+kernel' $file_ref | cut -d = -f 2 | tr -d ' '`
  # need at least rho & kappa (for acoustic kernels)
  if [ "$RHO" == "" ] || [ "$KAPPA" == "" ]; then
    echo "  missing reference kernel values: RHO=$RHO KAPPA=$KAPPA MU=$MU"
    echo
    exit 1
  else
    echo "  reference kernel values: RHO=$RHO KAPPA=$KAPPA MU=$MU"
  fi
  # compares with test output - using a relative tolerance of 0.001 (1 promille) with respect to expected value
  # final test result
  PASSED=0
  # checks rho kernel value
  if [ "$RHO" != "" ]; then
    VAL=`grep -E 'maximum value of rho[[:space:]]+kernel' $file_out | cut -d = -f 2 | tr -d ' '`
    echo "kernel rho   : $VAL"
    echo "" | awk '{diff=ex-val;diff_abs=(diff >= 0)? diff:-diff;diff_rel=diff_abs/ex;print "  value: expected = "ex" gotten = "val" - difference absolute = "diff_abs" relative = "diff_rel; if (diff_rel>0.001){print "  failed"; exit 1;}else{print "  good"; exit 0;} }' ex=$RHO val=$VAL
    if [[ $? -ne 0 ]]; then PASSED=1; fi
  fi
  # checks kappa kernel value
  if [ "$KAPPA" != "" ]; then
    VAL=`grep -E 'maximum value of kappa[[:space:]]+kernel' $file_out | cut -d = -f 2 | tr -d ' '`
    echo "kernel kappa : $VAL"
    echo "" | awk '{diff=ex-val;diff_abs=(diff >= 0)? diff:-diff;diff_rel=diff_abs/ex;print "  value: expected = "ex" gotten = "val" - difference absolute = "diff_abs" relative = "diff_rel; if (diff_rel>0.001){print "  failed"; exit 1;}else{print "  good"; exit 0;} }' ex=$KAPPA val=$VAL
    if [[ $? -ne 0 ]]; then PASSED=1; fi
  fi
  # checks mu kernel value (if available for elastic kernel)
  if [ "$MU" != "" ]; then
    VAL=`grep -E 'maximum value of mu[[:space:]]+kernel' $file_out | cut -d = -f 2 | tr -d ' '`
    echo "kernel mu    : $VAL"
    echo "" | awk '{diff=ex-val;diff_abs=(diff >= 0)? diff:-diff;diff_rel=diff_abs/ex;print "  value: expected = "ex" gotten = "val" - difference absolute = "diff_abs" relative = "diff_rel; if (diff_rel>0.001){print "  failed"; exit 1;}else{print "  good"; exit 0;} }' ex=$MU val=$VAL
    if [[ $? -ne 0 ]]; then PASSED=1; fi
  fi
  # overall pass
  if [[ $PASSED -ne 0 ]]; then
    echo "testing kernel values: failed"; exit 1;
  else
    echo "testing kernel values: all good"
  fi
  echo "######################################################################################################################"
}

# test example
cd $dir

# default setup
if [ -e DATA/Par_file ]; then
  # limit time steps for testing
  sed -i "s:^NSTEP .*:NSTEP    = 200:" DATA/Par_file
  # shortens output interval to avoid timeouts
  sed -i "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO    = 100:" DATA/Par_file
fi

# limit time steps for specific examples
# simple mesh example
if [ "$TESTDIR" == "EXAMPLES/applications/meshfem3D_examples/simple_model/" ]; then
  sed -i "s:^NSTEP .*:NSTEP    = 800:" DATA/Par_file
fi
# tpv5 example
if [ "$TESTDIR" == "EXAMPLES/applications/fault_examples/tpv5/" ]; then
  sed -i "s:^NSTEP .*:NSTEP    = 1500:" DATA/Par_file
fi
# layered halfspace example
if [ "$TESTDIR" == "EXAMPLES/applications/layered_halfspace/" ]; then
  sed -i "s:^NSTEP .*:NSTEP    = 200:" DATA/Par_file
fi
# small adjoint example
if [ "$TESTDIR" == "EXAMPLES/applications/small_adjoint_multiple_sources/" ]; then
  # full length
  sed -i "s:^NSTEP .*:NSTEP    = 1000:" DATA/Par_file
fi
# socal examples
if [ "$TESTDIR" == "EXAMPLES/applications/meshfem3D_examples/socal1D/" ]; then
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
    # just continue
    :
  fi
fi
# coupling FK
if [ "$TESTDIR" == "EXAMPLES/applications/small_example_coupling_FK_specfem/" ]; then
  sed -i "s:^NSTEP .*:NSTEP    = 1000:" DATA/Par_file
fi
# coupling SPECFEM
if [ "$TESTDIR" == "EXAMPLES/applications/small_example_coupling_SPECFEM_specfem/" ]; then
  # turning off mesh output to avoid using too much memory on the test nodes
  sed -i "s:^SAVE_MESH_FILES .*:SAVE_MESH_FILES    = .false.:" DATA.regional/Par_file
  sed -i "s:^SAVE_MESH_FILES .*:SAVE_MESH_FILES    = .false.:" DATA.local/Par_file
  sed -i "s:^CREATE_VTK_FILES .*:CREATE_VTK_FILES    = .false.:" DATA.regional/meshfem3D_files/Mesh_Par_file
  sed -i "s:^CREATE_VTK_FILES .*:CREATE_VTK_FILES    = .false.:" DATA.local/meshfem3D_files/Mesh_Par_file
fi
# elastic halfspace, no absorbing
if [ "$TESTDIR" == "EXAMPLES/applications/homogeneous_halfspace_HEX8_elastic_no_absorbing/" ]; then
  sed -i "s:^NSTEP .*:NSTEP    = 600:" DATA/Par_file
fi
# waterlayered_halfspace example
if [ "$TESTDIR" == "EXAMPLES/applications/waterlayered_halfspace/" ]; then
  sed -i "s:^NSTEP .*:NSTEP    = 600:" DATA/Par_file
fi
# tomographic model
if [ "$TESTDIR" == "EXAMPLES/applications/tomographic_model/" ]; then
  sed -i "s:^NSTEP .*:NSTEP    = 500:" DATA/Par_file
fi
# cavity example
if [ "$TESTDIR" == "EXAMPLES/applications/meshfem3D_examples/cavity/" ]; then
  sed -i "s:^NSTEP .*:NSTEP    = 2000:" DATA/Par_file
fi
# SEP example
if [ "$TESTDIR" == "EXAMPLES/applications/meshfem3D_examples/sep_bathymetry/" ]; then
  sed -i "s:^NSTEP .*:NSTEP    = 1000:" DATA/Par_file
fi
# serial/no-mpi coverage row
if [ "$TESTDIR" == "EXAMPLES/applications/homogeneous_acoustic/" ]; then
  sed -i "s:^NPROC .*:NPROC    = 1:" DATA/Par_file
  sed -i "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO    = 20:" DATA/Par_file
fi
# inversion example coverage setup (travis parity)
if [ "$TESTDIR" == "EXAMPLES/applications/inversion_examples/fwi_test_acoustic/" ]; then
  sed -i "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO    = 500:" DATA/Par_file
  # reduce mesh resolution for CI runtime
  if [ -e DATA/meshfem3D_files/Mesh_Par_file.INIT ]; then sed -i "s:20:16:g" DATA/meshfem3D_files/Mesh_Par_file.INIT; fi
  if [ -e DATA/meshfem3D_files/Mesh_Par_file.TRUE ]; then sed -i "s:20:16:g" DATA/meshfem3D_files/Mesh_Par_file.TRUE; fi
  if [ -e DATA/meshfem3D_files/interfaces.dat ]; then sed -i "s:20:16:g" DATA/meshfem3D_files/interfaces.dat; fi
  if [ -e DATA/inverse_problem/inversion_fwi.dat ]; then sed -i "s/Niter .*/Niter       : 1/" DATA/inverse_problem/inversion_fwi.dat; fi
  if [ -e DATA/inverse_problem/acquisition.dat ]; then sed -i "s/NSTEP .*/NSTEP         : 100/" DATA/inverse_problem/acquisition.dat; fi
fi

# coverage runs use short steps
if [ "$TESTCOV" == "true" ]; then
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
fi

## HDF5 - i/o example
if [ "${HDF5}" == "true" ]; then
  echo
  echo "test run w/ HDF5"
  echo
  # turns on HDF5
  echo "turning on HDF5"
  sed -i "s:^HDF5_ENABLED .*:HDF5_ENABLED    = .true.:" DATA/Par_file
  sed -i "s:^HDF5_FOR_MOVIES .*:HDF5_FOR_MOVIES    = .true.:" DATA/Par_file
  sed -i "s:^HDF5_IO_NODES .*:HDF5_IO_NODES    = 1:" DATA/Par_file
  # replaces run script
  cp -v run_this_example_HDF5_IO_server.sh run_this_example.sh
fi

## adios
if [ "${ADIOS2}" == "true" ]; then
  # turns on ADIOS
  echo "turning on ADIOS"
  sed -i "s:^ADIOS_ENABLED .*:ADIOS_ENABLED = .true.:" DATA/Par_file
fi

## GPU
if [ "${GPU}" == "true" ]; then
  # turns on GPU
  echo "turning on GPU"
  sed -i "s:^GPU_MODE .*:GPU_MODE = .true.:" DATA/Par_file
fi

# save Par_file state
if [ -e DATA/Par_file ]; then
  cp -v DATA/Par_file DATA/Par_file.bak
fi

# runs simulation
if [ "${RUN_KERNEL}" == "true" ]; then
  # use kernel script
  ./run_this_example_kernel.sh | tee output.log
else
  # default script
  ./run_this_example.sh
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# simulation done
echo
echo "simulation done: `pwd`"
echo `date`
echo

# seismogram comparison
RUN_COMPARE=true
# turn off for non-default runs
if [ "${TESTCOV}" == "true" ]; then RUN_COMPARE=false; fi
if [ "${DEBUG}" == "true" ]; then RUN_COMPARE=false; fi
if [ "${RUN_KERNEL}" == "true" ]; then RUN_COMPARE=false; fi
if [ "$TESTDIR" == "EXAMPLES/applications/inversion_examples/fwi_test_acoustic/" ]; then RUN_COMPARE=false; fi

if [ "${RUN_COMPARE}" == "true" ]; then
  my_test
else
  # no comparisons
  :     # do nothing
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# kernel test
if [ "${RUN_KERNEL}" == "true" ]; then
  # check kernel values
  my_kernel_test
  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi
  # clean up
  rm -rf OUTPUT_FILES/ SEM/ output.log

  # re-run kernel test w/ UNDO_ATT
  echo
  echo "*****************************************"
  echo "run kernel w/ UNDO_ATTENUATION_AND_OR_PML"
  echo "*****************************************"
  echo

  # turns on UNDO_ATTENUATION_AND_OR_PML
  echo "turning on UNDO_ATTENUATION_AND_OR_PML"
  sed -i "s:^UNDO_ATTENUATION_AND_OR_PML .*:UNDO_ATTENUATION_AND_OR_PML = .true.:" DATA/Par_file

  # use kernel script
  ./run_this_example_kernel.sh | tee output.log
  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi
  # kernel test
  my_kernel_test
  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi
fi

# restore original Par_file
if [ -e DATA/Par_file.bak ]; then
  cp -v DATA/Par_file.bak DATA/Par_file
fi

# cleanup
rm -rf OUTPUT_FILES/
if [ -e DATABASES_MPI ]; then rm -rf DATABASES_MPI/; fi
if [ -e SEM ]; then rm -rf SEM/; fi

echo
echo "all good"
echo `date`
echo
