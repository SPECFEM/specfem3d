#!/bin/bash

mpif90=ftn
mpicc=cc
f90=ftn
cc=cc

# run flags
#flags="-eF -em -rm"
#cflags="-h list=m"

## memory > 2GB
#flags="-eF -em -rm -O3,fp3 -hpic -dynamic"
#cflags="-h list=m -hpic -dynamic"

# debug flags
flags="-g -Rb -eF -rm -eC -eD" # -hfp2"
cflags="-g -h list=m"

###
### CUDA
###
CUDA_INC="${CRAY_CUDATOOLKIT_DIR}/include"
CUDA_LIB="${CRAY_CUDATOOLKIT_DIR}/lib64"

###
### mpi.h / mpif.h
###
MPI_INC="${CRAY_MPICH2_DIR}/include"

./configure \
--with-cuda=cuda5 CUDA_INC="$CUDA_INC" CUDA_LIB="$CUDA_LIB" MPI_INC="$MPI_INC" \
MPIFC=$mpif90 MPICC=$mpicc FC=$f90 CC=$cc CXX=$cc FLAGS_CHECK="$flags" CFLAGS="$cflags"
