#!/bin/bash

date
echo "running directory: `pwd`"
echo


# get the number of processors
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`


# remove trash
rm -rf bin
rm -rf DATABASES_MPI
rm -rf OUTPUT_FILES
rm -rf generating_field
rm -rf correlation_field
rm -f DATABASES_MPI/*
rm -f OUTPUT_FILES/*
rm -f NOISE_TOMOGRAPHY/S_squared
rm -f convert_utm2geo
rm -f NOISE_TOMOGRAPHY.py
rm -f PetersonNoiseModel.py
rm -f convert_utm2geo.py

# make directories
mkdir bin
mkdir DATABASES_MPI
mkdir OUTPUT_FILES

# link executables
cd bin/
ln -s ../../../bin/xmeshfem3D
ln -s ../../../bin/xgenerate_databases
ln -s ../../../bin/xspecfem3D
ln -s ../../../bin/xcreate_movie_shakemap_AVS_DX_GMT
cd ../


# link utilities
ln -s ../../utils/convert_utm2geo.py .
ln -s ../noise_tomography/NOISE_TOMOGRAPHY.py .
ln -s ../noise_tomography/PetersonNoiseModel.py


# generate noise source time function
echo
echo " generating noise source time function..."
echo
./run_generate_S_squared.sh
if [[ $? -ne 0 ]]; then exit 1; fi

# run mesher
if [ "$NPROC" -eq 1 ]; then
  echo
  echo " running mesh generation..."
  echo
  ./bin/xmeshfem3D
else
  echo
  echo " running mesh generation on $NPROC processors..."
  echo
  mpirun -np $NPROC ./bin/xmeshfem3D
fi
if [[ $? -ne 0 ]]; then exit 1; fi


# run database generation
if [ "$NPROC" -eq 1 ]; then
  echo
  echo " running database generation..."
  echo
  ./bin/xgenerate_databases
else
  echo
  echo " running database generation on $NPROC processors..."
  echo
  mpirun -np $NPROC ./bin/xgenerate_databases
fi
if [[ $? -ne 0 ]]; then exit 1; fi


# generate noise distribution and direction
echo
echo " generating noise distribution and direction..."
echo
./generate_noise_distribution_direction.py
if [[ $? -ne 0 ]]; then exit 1; fi


# run noise simulation step 1
cp DATA/Par_file_step1 DATA/Par_file
if [ "$NPROC" -eq 1 ]; then
  echo
  echo " running noise simulation step 1..."
  echo
  ./bin/xspecfem3D
else
  echo
  echo " running noise simulation step 1 on $NPROC processors..."
  echo
  mpirun -np $NPROC ./bin/xspecfem3D
fi
if [[ $? -ne 0 ]]; then exit 1; fi

# save generating field movie
echo
echo " saving generating field movie..."
echo
./bin/xcreate_movie_shakemap_AVS_DX_GMT < generate_movie_input.txt
if [[ $? -ne 0 ]]; then exit 1; fi
mkdir generating_field
mv OUTPUT_FILES/AVS_movie* generating_field


# run noise simulation step 2
cp DATA/Par_file_step2 DATA/Par_file
if [ "$NPROC" -eq 1 ]; then
  echo
  echo " running noise simulation step 2..."
  echo
  ./bin/xspecfem3D
else
  echo
  echo " running noise simulation step 2 on $NPROC processors..."
  echo
  mpirun -np $NPROC ./bin/xspecfem3D
fi
if [[ $? -ne 0 ]]; then exit 1; fi

# save correlation field  movie
echo
echo " saving generating field movie..."
echo
./bin/xcreate_movie_shakemap_AVS_DX_GMT < generate_movie_input.txt
if [[ $? -ne 0 ]]; then exit 1; fi
mkdir correlation_field
mv OUTPUT_FILES/AVS_movie* correlation_field

echo
date
echo "done"
