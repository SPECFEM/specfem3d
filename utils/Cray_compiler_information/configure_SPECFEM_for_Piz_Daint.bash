#!/bin/bash

mpif90=ftn
mpicc=cc
f90=ftn
cc=cc

# run flags
#flags="-eF -em -rm"
#cflags="-h list=m"

## memory > 2GB
flags="-eF -em -rm -O3,fp3 -hpic -dynamic"
cflags="-h list=m -hpic -dynamic"

# debug flags
#flags="-g -Rb -eF -rm -eC -eD" # -hfp2"
#cflags="-g -h list=m"

###
### CUDA
###
CUDA_INC="${CRAY_CUDATOOLKIT_DIR}/include"
CUDA_LIB="${CRAY_CUDATOOLKIT_DIR}/lib64"

###
### mpi.h / mpif.h
###
MPI_INC="${CRAY_MPICH2_DIR}/include"

# DK DK June 2018: on Pascal P100 (for instance on Piz Daint in Switzerland) do *NOT* change --with-cuda=cuda8 below to --with-cuda=cuda9;
# DK DK June 2018: only change that on machines that have Volta V100, not Pascal P100

./configure \
--with-cuda=cuda8 CUDA_INC="$CUDA_INC" CUDA_LIB="$CUDA_LIB" MPI_INC="$MPI_INC" \
MPIFC=$mpif90 MPICC=$mpicc FC=$f90 CC=$cc CXX=$cc FLAGS_CHECK="$flags" CFLAGS="$cflags"

