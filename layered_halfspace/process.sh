#!/bin/bash
#
# script runs decomposition,database generation and solver
# using this example setup
#
# prior to running this script, you must create the mesh files
# in directory MESH/
#

###################################################

# number of processes
NPROC=4

##################################################

echo "running example: `date`"
currentdir=`pwd`

echo
echo "(will take about 27 minutes)"
echo


# sets up directory structure in current example directoy
echo
echo "   setting up example..."
echo

mkdir -p bin
mkdir -p in_out_files/OUTPUT_FILES
mkdir -p in_out_files/DATABASES_MPI

rm -rf in_out_files/OUTPUT_FILES/*
rm -rf in_out_files/DATABASES_MPI/*

# compiles executables in root directory
cd ../../

# compiles with flag for a point force (with a ricker source time function)
cd src/shared/
cp constants.h constants.h.org
sed -e "s:USE_FORCE_POINT_SOURCE.*:USE_FORCE_POINT_SOURCE = .true. :" constants.h.org > constants.h.tmp
sed -e "s:FACTOR_FORCE_SOURCE.*:FACTOR_FORCE_SOURCE = -1.d15 :" constants.h.tmp > constants.h
rm -f constants.h.tmp
cd ../../
make > tmp.log

# backup & restores original file again
cp src/shared/constants.h $currentdir/in_out_files/OUTPUT_FILES/constants.h.bak
mv src/shared/constants.h.org src/shared/constants.h

cd $currentdir

# links executables
cd bin/
ln -s ../../../bin/xdecompose_mesh_SCOTCH
ln -s ../../../bin/xgenerate_databases
ln -s ../../../bin/xspecfem3D
cd ../

# decomposes mesh
echo
echo "  decomposing mesh..."
echo
./bin/xdecompose_mesh_SCOTCH $NPROC MESH/ in_out_files/DATABASES_MPI/

# stores setup
cp in_data_files/Par_file in_out_files/OUTPUT_FILES/
cp in_data_files/CMTSOLUTION in_out_files/OUTPUT_FILES/
cp in_data_files/STATIONS in_out_files/OUTPUT_FILES/

# runs database generation
echo
echo "  running database generation..."
echo
cd bin/
mpirun -np $NPROC ./xgenerate_databases
cd ../

# runs simulation
echo
echo "  running solver..."
echo
cd bin/
mpirun -np $NPROC ./xspecfem3D
cd ../

echo
echo "see results in directory: in_out_files/OUTPUT_FILES/"
echo
echo "done"
echo `date`


