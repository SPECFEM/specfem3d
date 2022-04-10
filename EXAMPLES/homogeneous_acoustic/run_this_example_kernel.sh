#!/bin/bash

echo "running example: `date`"
currentdir=`pwd`

# sets up directory structure in current example directory
echo
echo "setting up example..."
echo

rm -f change_simulation_type.pl
ln -s ../../utils/change_simulation_type.pl

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
cd bin/
rm -f *
ln -s ../../../bin/xdecompose_mesh
ln -s ../../../bin/xgenerate_databases
ln -s ../../../bin/xspecfem3D
ln -s ../../../bin/xcombine_vol_data_vtk
cd ../


## limits run time
cp DATA/Par_file DATA/Par_file.org
sed -i "s:^NSTEP .*:NSTEP = 300:" DATA/Par_file
sed -i "s:^USE_TRICK_FOR_BETTER_PRESSURE .*:USE_TRICK_FOR_BETTER_PRESSURE = .false.:" DATA/Par_file

# stores setup
cp DATA/Par_file OUTPUT_FILES/
cp DATA/CMTSOLUTION OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/

# get the number of processors, ignoring comments in the Par_file
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `
mkdir -p $BASEMPIDIR

# decomposes mesh using the pre-saved mesh files in MESH-default
echo
echo "decomposing mesh..."
echo
./bin/xdecompose_mesh $NPROC ./MESH-default $BASEMPIDIR
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# runs database generation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo
  echo "running database generation..."
  echo
  ./bin/xgenerate_databases
else
  # This is a MPI simulation
  echo
  echo "running database generation on $NPROC processors..."
  echo
  mpirun -np $NPROC ./bin/xgenerate_databases
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo "#########################################################"
echo "forward simulation"
echo "#########################################################"
echo "(running forward simulation with saving forward wavefield)"
echo
./change_simulation_type.pl -F

# backup
cp DATA/Par_file OUTPUT_FILES/Par_file.for

# runs simulation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo
  echo "running solver..."
  echo
  ./bin/xspecfem3D
else
  # This is a MPI simulation
  echo
  echo "running solver on $NPROC processors..."
  echo
  mpirun -np $NPROC ./bin/xspecfem3D
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo "#########################################################"
echo "adjoint sources"
echo "#########################################################"
echo "setting up adjoint sources"
echo
# setup adjoint sources directory
mkdir -p SEM
rm -rf SEM/*

./create_adjoint_sources.sh

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# backup forward outputs
mkdir -p OUTPUT_FILES/forward_run
rm -rf OUTPUT_FILES/forward_run/*

mv OUTPUT_FILES/timestamp* OUTPUT_FILES/forward_run/
mv OUTPUT_FILES/output_* OUTPUT_FILES/forward_run/
mv OUTPUT_FILES/DB.* OUTPUT_FILES/forward_run/
mv OUTPUT_FILES/plot_* OUTPUT_FILES/sr.vtk OUTPUT_FILES/forward_run/

echo "#########################################################"
echo "kernel simulation"
echo "#########################################################"
echo "(running kernel simulation: SIMULATION_TYPE == 3)"
echo
./change_simulation_type.pl -b
# stores output
cp DATA/Par_file DATA/Par_file.kernel
cp DATA/CMTSOLUTION DATA/STATIONS* OUTPUT_FILES/

# runs simulation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo
  echo "running solver (kernel run)..."
  echo
  ./bin/xspecfem3D
else
  # This is a MPI simulation
  echo
  echo "running solver (kernel run) on $NPROC processors..."
  echo
  mpirun -np $NPROC ./bin/xspecfem3D
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "see results in directory       : OUTPUT_FILES/"
echo "    kernel outputs in directory: $BASEMPIDIR"
echo
echo "done"
echo `date`


