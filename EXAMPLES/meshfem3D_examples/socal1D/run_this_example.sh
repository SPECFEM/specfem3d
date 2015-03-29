#!/bin/bash
#
# script runs mesher,database generation and solver
# using this example setup
#

###################################################

# number of processes
NPROC=1

##################################################

echo "running example: `date`"
currentdir=`pwd`

echo
echo "(will take about 5 minutes)"
echo

# sets up directory structure in current example directory
echo
echo "   setting up example..."
echo

mkdir -p bin
mkdir -p OUTPUT_FILES/DATABASES_MPI

rm -f OUTPUT_FILES/*
rm -rf OUTPUT_FILES/DATABASES_MPI/*

cd $currentdir

# links executables
rm -f bin/*
cp ../../../bin/* bin/
if [ ! -e bin/xspecfem3D ]; then echo "compilation failed, please check..."; exit 1; fi

# stores setup
cp DATA/meshfem3D_files/Mesh_Par_file OUTPUT_FILES/
cp DATA/Par_file OUTPUT_FILES/
cp DATA/CMTSOLUTION OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/

# decomposes mesh
echo
echo "running mesher..."
echo
mpirun -np $NPROC ./bin/xmeshfem3D
mv OUTPUT_FILES/output_mesher.txt OUTPUT_FILES/output_meshfem3D.txt


# runs database generation
echo
echo "running database generation..."
echo
mpirun -np $NPROC ./bin/xgenerate_databases

# exit here if you want to do mesher only
#exit

# runs simulation
echo
echo "running solver..."
echo
mpirun -np $NPROC ./bin/xspecfem3D

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
date

# To create a full mesh using combine_vol_data:
# cd bin
# xcombine_vol_data 0 0 vs ../OUTPUT_FILES/DATABASES_MPI/ ../OUTPUT_FILES/DATABASES_MPI/ 1
# cd ../OUTPUT_FILES/DATABASES_MPI/
# (check that mesh2vtu.pl is working)
# ../../../../../utils/Visualization/Paraview/mesh2vtu.pl -i vs.mesh -o vs.vtu
#


