#!/bin/bash
#
# single forward simulation
#

echo "running example: `date`"
currentdir=`pwd`

# sets up directory structure in current example directory
echo
echo "setting up example..."
echo

# backup
cp -vp DATA/Par_file DATA/Par_file.org

## single forward simulation
# turn off inverse flag
sed -i "s:^INVERSE_FWI_FULL_PROBLEM .*:INVERSE_FWI_FULL_PROBLEM = .false.:" DATA/Par_file

# takes true model
cp -v DATA/meshfem3D_files/Mesh_Par_file.TRUE DATA/meshfem3D_files/Mesh_Par_file

# takes 1. source
cp -v DATA/CMTSOLUTION.1 DATA/CMTSOLUTION

# checks if executables were compiled and available
if [ ! -e ../../../bin/xspecfem3D ]; then
  echo "Please compile first all binaries in the root directory, before running this example..."; echo
  exit 1
fi

# cleans output files
mkdir -p OUTPUT_FILES
rm -rf OUTPUT_FILES/*

# SPECFEM3D root directory
ROOT=../../../

# links executables
mkdir -p bin
cd bin/
rm -f *
exec_files=../$ROOT/bin/x*
for f in $exec_files; do
  ln -s $f
done
cd ../

# stores setup
cp DATA/Par_file OUTPUT_FILES/
cp DATA/CMTSOLUTION OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/

# get the number of processors, ignoring comments in the Par_file
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

mkdir -pv $BASEMPIDIR

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


# restore original Par_file
cp -vp DATA/Par_file.org DATA/Par_file

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
echo `date`


