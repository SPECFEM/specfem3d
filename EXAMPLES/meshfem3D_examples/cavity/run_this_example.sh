#!/bin/bash

echo "running example: `date`"
currentdir=`pwd`

# sets up directory structure in current example directory
echo
echo "   setting up example..."
echo

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
ln -s ../../../../bin/xcreate_movie_shakemap_AVS_DX_GMT
cd ../

# stores setup
cp DATA/meshfem3D_files/Mesh_Par_file OUTPUT_FILES/
cp DATA/Par_file OUTPUT_FILES/
cp DATA/FORCESOLUTION OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/

# get the number of processors, ignoring comments in the Par_file
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

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

# creates AVS movie files
start_timestep=1
end_timestep=6001
display_type=1     #norm velocity

rm -f OUTPUT_FILES/AVS_movie_*.inp

./bin/xcreate_movie_shakemap_AVS_DX_GMT <<EOF
2
${start_timestep}
${end_timestep}
1
${display_type}
EOF
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi


echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
echo `date`

