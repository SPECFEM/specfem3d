#!/bin/bash
#
################################################################

# noise simulation step
step=$1

################################################################

if [ "$step" == "" ]; then echo "usage: ./run_single_step.sh step[1/2/3]"; exit 1; fi

echo "running step: $step"
echo `date`
echo
currentdir=`pwd`

# setup Par_file
#case $step in
#1) cp -v DATA/Par_file_step1 DATA/Par_file ;;
#2) cp -v DATA/Par_file_step2 DATA/Par_file ;;
#3) cp -v DATA/Par_file_step3 DATA/Par_file ;;
#*) echo "step not recognized: $step"; echo "please use as step number 1, 2 or 3"; exit 1 ;;
#esac

case $step in
1) sed -i "s:^SIMULATION_TYPE .*:SIMULATION_TYPE = 1:" DATA/Par_file
   sed -i "s:^NOISE_TOMOGRAPHY .*:NOISE_TOMOGRAPHY = 1:" DATA/Par_file
   sed -i "s:^SAVE_FORWARD .*:SAVE_FORWARD = .false.:" DATA/Par_file
   ;;
2) sed -i "s:^SIMULATION_TYPE .*:SIMULATION_TYPE = 1:" DATA/Par_file
   sed -i "s:^NOISE_TOMOGRAPHY .*:NOISE_TOMOGRAPHY = 2:" DATA/Par_file
   sed -i "s:^SAVE_FORWARD .*:SAVE_FORWARD = .true.:" DATA/Par_file
   ;;
3) sed -i "s:^SIMULATION_TYPE .*:SIMULATION_TYPE = 3:" DATA/Par_file
   sed -i "s:^NOISE_TOMOGRAPHY .*:NOISE_TOMOGRAPHY = 3:" DATA/Par_file
   sed -i "s:^SAVE_FORWARD .*:SAVE_FORWARD = .false.:" DATA/Par_file
   ;;
*) echo "step not recognized: $step"; echo "please use as step number 1, 2 or 3"; exit 1 ;;
esac
cp -v DATA/Par_file DATA/Par_file_step${step}
echo

# get the number of processors, ignoring comments in the Par_file
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

if [ "$step" == "1" ]; then
  # sets up directory structure in current example directory
  echo
  echo "   setting up example..."
  echo

  # checks if executables were compiled and available
  if [ ! -e ../../bin/xspecfem3D ]; then
    echo "Please compile first all binaries in the root directory, before running this example..."; echo
    exit 1
  fi

  # cleans output files
  mkdir -p OUTPUT_FILES
  rm -rf OUTPUT_FILES/*

  # links executables
  mkdir -p bin
  rm -f bin/*
  cd bin/
  ln -s ../../../bin/xdecompose_mesh
  ln -s ../../../bin/xgenerate_databases
  ln -s ../../../bin/xspecfem3D
  ln -s ../../../bin/xcombine_vol_data_vtk
  cd ../

  # DATABASES directory
  BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `
  mkdir -p $BASEMPIDIR

  # adjoint sources directory
  mkdir -p SEM
fi

# stores setup
cp DATA/Par_file OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/

# mesh setup
if [ "$step" == "1" ]; then
  # decomposes mesh using the pre-saved mesh files in MESH-default
  echo
  echo "  decomposing mesh..."
  echo
  ./bin/xdecompose_mesh $NPROC ./MESH-default $BASEMPIDIR
  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi

  # runs database generation
  if [ "$NPROC" -eq 1 ]; then
    # This is a serial simulation
    echo
    echo "  running database generation..."
    echo
    ./bin/xgenerate_databases
  else
    # This is a MPI simulation
    echo
    echo "  running database generation on $NPROC processors..."
    echo
    mpirun -np $NPROC ./bin/xgenerate_databases
  fi
  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi

fi

# runs simulation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo
  echo "  running solver..."
  echo
  ./bin/xspecfem3D
else
  # This is a MPI simulation
  echo
  echo "  running solver on $NPROC processors..."
  echo
  mpirun -np $NPROC ./bin/xspecfem3D
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# backup
mkdir -p OUTPUT_FILES/step_${step}
rm -rf OUTPUT_FILES/step_${step}/*
echo

mv -v OUTPUT_FILES/*.sem* OUTPUT_FILES/step_${step}/
mv -v OUTPUT_FILES/output_*.txt OUTPUT_FILES/step_${step}/
mv -v OUTPUT_FILES/STATIONS OUTPUT_FILES/step_${step}/
mv -v OUTPUT_FILES/Par_file OUTPUT_FILES/step_${step}/

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done: `date`"
echo


