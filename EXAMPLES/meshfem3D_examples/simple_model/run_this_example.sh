#!/bin/bash

echo "running example: `date`"
currentdir=`pwd`

# sets up directory structure in current example directory
echo
echo "   setting up example..."
echo

# checks if executables were compiled and available
if [ ! -e ../../../bin/xspecfem3D ]; then
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
ln -s ../../../../bin/xmeshfem3D
ln -s ../../../../bin/xgenerate_databases
ln -s ../../../../bin/xspecfem3D
cd ../

# stores setup
cp DATA/meshfem3D_files/Mesh_Par_file OUTPUT_FILES/
cp DATA/Par_file OUTPUT_FILES/
cp DATA/CMTSOLUTION OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/

# get the number of processors, ignoring comments in the Par_file
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

# HDF5 i/o server
# get the number of io server propcessors, ignoring comments in the Par_file
IO_NODES=`grep ^HDF5_IO_NODES DATA/Par_file | sed 's/#.*//' | cut -d = -f 2`
if [ "$IO_NODES" -gt 0 ]; then
  echo
  echo "Par_file setup:"
  echo "  NPROC         = $NPROC"
  echo "  HDF5_IO_NODES = $IO_NODES"
  echo
fi

BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `
mkdir -p $BASEMPIDIR

# runs in-house mesher
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo
  echo "  running mesher..."
  echo
  ./bin/xmeshfem3D
else
  # This is a MPI simulation
  echo
  echo "  running mesher on $NPROC processors..."
  echo
  mpirun -np $NPROC ./bin/xmeshfem3D
fi
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

# HDF5 i/o server
# setting new total number of MPI processes
# solver runs with (NPROC + IO_NODES) mpi processes
if [ "$IO_NODES" -gt 0 ]; then
  NPROC=$(($NPROC+$IO_NODES))
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

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
echo `date`


