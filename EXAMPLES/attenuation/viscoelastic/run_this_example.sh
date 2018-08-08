#!/bin/bash
#
# script runs mesher and solver (in serial)
# using this example setup
#

echo "running example: `date`"
currentdir=`pwd`

# sets up directory structure in current example directoy
echo
echo "   setting up example..."
echo

rm -r -f OUTPUT_FILES
mkdir -p OUTPUT_FILES
mkdir -p OUTPUT_FILES/DATABASES_MPI

mkdir -p DATA DATA/meshfem3D_files bin

# sets up local DATA/ directory
cd DATA/
cp ../Par_file_attenuation_3D Par_file
cp ../FORCESOLUTION_attenuation_3D FORCESOLUTION
cp ../STATIONS_attenuation_3D STATIONS
cp ../interface1_attenuation_3D.dat ./meshfem3D_files/interface1.dat
cp ../interfaces_attenuation_3D.dat ./meshfem3D_files/interfaces.dat
cp ../Mesh_Par_file_attenuation_3D ./meshfem3D_files/Mesh_Par_file
cd ../

cd $currentdir

# links executables
cd bin
rm -f xmeshfem3D xgenerate_databases xspecfem3D
ln -s ../../../../bin/xmeshfem3D
ln -s ../../../../bin/xgenerate_databases
ln -s ../../../../bin/xspecfem3D
cd ..

# Get the number of processors
NPROC=`grep ^NPROC DATA/Par_file | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`

# runs database generation
echo
echo "  running mesher..."
echo
mpirun -n $NPROC ./bin/xmeshfem3D
mpirun -n $NPROC ./bin/xgenerate_databases

# runs simulation
echo
echo "  running solver..."
echo
mpirun -n $NPROC ./bin/xspecfem3D

echo
echo "see results in directory: OUTPUT_FILES/"
echo "seismograms can be plotted using the gnuplot script present in this directory"
echo
echo "done"
date
