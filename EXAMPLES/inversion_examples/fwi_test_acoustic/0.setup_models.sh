#!/bin/bash
#
#
# script to create two models for the inversion setup:
# (1) as the true model
# (2) as the (homogeneous) initial starting model
#
echo "setup_models: `date`"
currentdir=`pwd`

# sets up directory structure in current example directory
echo
echo "setting up example..."
echo

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

cd $currentdir/

# backup original Par_file
cp -v DATA/Par_file DATA/Par_file.org

# get the number of processors, ignoring comments in the Par_file
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

mkdir -p $BASEMPIDIR
rm -rf $BASEMPIDIR/*

# cleans output files
mkdir -p OUTPUT_FILES
rm -rf OUTPUT_FILES/*

# setup STATIONS
if [[ 1 == 0 ]]; then
  echo "# top" > tmp.stations
  echo "" | awk '{i=0;for(k=1;k<10;k++){for(j=1;j<10;j++){i++;x=k*100.0;y=j*100.0; print "A"i,"DB",x,y,0.0,10.0}}}' >> tmp.stations
  echo "# bottom" >> tmp.stations
  echo "" | awk '{i=0;for(k=1;k<10;k++){for(j=1;j<10;j++){i++;x=k*100.0;y=j*100.0; print "B"i,"DB",x,y,0.0,990.0}}}' >> tmp.stations
  echo "# left" >> tmp.stations
  echo "" | awk '{i=0;for(k=1;k<10;k++){for(j=1;j<10;j++){i++;x=k*100.0;y=j*100.0; print "C"i,"DB",10.0,x,0.0,y}}}' >> tmp.stations
  echo "# right" >> tmp.stations
  echo "" | awk '{i=0;for(k=1;k<10;k++){for(j=1;j<10;j++){i++;x=k*100.0;y=j*100.0; print "D"i,"DB",990.0,x,0.0,y}}}' >> tmp.stations
  echo "# front" >> tmp.stations
  echo "" | awk '{i=0;for(k=1;k<10;k++){for(j=1;j<10;j++){i++;x=k*100.0;y=j*100.0; print "E"i,"DB",x,10.0,0.0,y}}}' >> tmp.stations
  echo "# back" >> tmp.stations
  echo "" | awk '{i=0;for(k=1;k<10;k++){for(j=1;j<10;j++){i++;x=k*100.0;y=j*100.0; print "F"i,"DB",x,990.0,0.0,y}}}' >> tmp.stations
  mv -v tmp.stations DATA/STATIONS
fi

# setup dummy CMTSOLUTION (and STATIONS) file needed to run mesher
cp -v DATA/CMTSOLUTION.1 DATA/CMTSOLUTION

# cleans previous runs in case
if [ -e run0001/DATA/Par_file ]; then rm -f run0001/DATA/Par_file; fi

##
## meshing
##

# model
models=( TRUE INIT )

for model in ${models[@]}; do

MODEL_DIR=MODEL_$model
echo ""
echo "model: $MODEL_DIR"
echo ""

# cleans output files
mkdir -pv OUTPUT_FILES
rm -rf OUTPUT_FILES/*

## 1. true model
# Par_file for meshing
cp -v DATA/Par_file.org DATA/Par_file
cp -v DATA/meshfem3D_files/Mesh_Par_file.$model DATA/meshfem3D_files/Mesh_Par_file

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

# backup files for solver
cp OUTPUT_FILES/values_from_mesher.h $BASEMPIDIR/

# moves output to model directory
mkdir -p $MODEL_DIR
rm -rf $MODEL_DIR/*

mv -v $BASEMPIDIR/* $MODEL_DIR/
mv -v OUTPUT_FILES $MODEL_DIR/

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# cleanup
#rm -f DATA/Par_file
rm -f DATA/meshfem3D_files/Mesh_Par_file

echo
echo "  created model in directory: $MODEL_DIR/"
echo

done

echo
echo "setup done"
echo




