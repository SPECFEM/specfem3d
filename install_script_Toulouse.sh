#!/bin/bash

###########################
# PARAMETERS AND SETTINGS #
###########################

# Get absolute path of installation repository
rootdir=`pwd`

# Define local parameters (for Toulouse)
FC='ifort'
CC='icc'
MPIFC='mpiifort'
MPICC='mpiicc'
CCFLAGS='-O3'
FFLAGS='-O3 -xHost -assume byterecl'
MPI_INC='/opt/intel/compilers_and_libraries_2017.2.174/linux/mpi/include64'

###################
# COMPILE SPECFEM #
###################

# Go to specfem directory
cd specfem3d_Git_devel

# Launch configure
./configure FC=$FC MPIFC=$MPIFC CC=$CC MPI_INC=$MPI_INC

# Modify specfem Makefile
cp Makefile Makefile.old
line="FLAGS_CHECK = ${FFLAGS}"
sed -i "/^FLAGS_CHECK/s/.*/${line}/" Makefile
idline=`grep -n "^FC =" Makefile | cut -d : -f 1`
line='DEBUG_COUPLED_FLAG=-DDEBUG_COUPLED -DUSE_VTK_INSTEAD_OF_MESH'
sed -i "${idline} i DEBUG_COUPLED_FLAG=-DDEBUG_COUPLED -DUSE_VTK_INSTEAD_OF_MESH" Makefile

# Mofifications of files constants.h and save_databases.f90
cp setup/constants.h setup/constants.h.old
line="  logical, parameter :: USE_SOURCES_RECEIVERS_Z = .true."
sed -i "/:: USE_SOURCES_RECEIVERS_Z/s/.*/${line}/" setup/constants.h

cp src/meshfem3D/save_databases.F90 src/meshfem3D/save_databases.F90.old
line="  logical, parameter :: SAVE_MESH_AS_CUBIT = .true."
sed -i "/:: SAVE_MESH_AS_CUBIT/s/.*/${line}/" src/meshfem3D/save_databases.F90

# Launch compilation
make realclean
make all -j 8

# Go back to rootdir
cd $rootdir

##########################
# COMPILE AXISEM SPECFEM #
##########################

# Go to repository of utils
cd EXTERNAL_CODES_coupled_with_SPECFEM3D/AxiSEM_for_SPECFEM3D/UTILS_COUPLING_SpecFEM

# Modify config.h
cp config.h config.h.old
sed -i "s/^CC=.*/CC = $CC/" config.h
sed -i "s/^FC=.*/FC = $MPIFC/" config.h
sed -i "s/^CCFLAGS.*/CCFLAGS = $CCFLAGS/" config.h
sed -i "s/^FFLAGS.*/FFLAGS = $FFLAGS/" config.h

# Compile
make clean
make all

# Go to AxiSEM to change macro file
cd ../AxiSEM_modif_for_coupling_with_specfem
cp make_axisem.macros make_axisem.macros.old
echo -e "#Intel compiler for Toulouse station
CC\t= $CC
FC\t= $MPIFC
FFLAGS\t= $FFLAGS
CFLAGS\t= $CCFLAGS
LDFLAGS\t= $CCFLAGS\n" >> make_axisem.macros

# Go back to root directory
cd $rootdir


#######################
# COMPILE SPECFEM_INV #
#######################

# Go to specfem repository
cd inverse_problem

# Change makefile
cp Makefile Makefile.old
sed -i "s/^FC =.*/FC = $FC/" Makefile
sed -i "s/^MPIFC =.*/MPIFC = $MPIFC/" Makefile
sed -i "s/^FLAGS_CHECK =.*/FLAGS_CHECK = $FFLAGS/" Makefile

# Compile
make clean
make

# Go back to root directory
cd $rootdir


##############################
# CHANGE ROOTDIR OF EXAMPLES #
##############################
./setup_hardcoded_path_in_examples.sh

