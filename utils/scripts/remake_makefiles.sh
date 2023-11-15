#!/bin/bash
rm -r autom4te.cache
aclocal -I ./m4/
echo "autoconf!"
autoconf configure.ac > test_conf
if [ $? -eq 0 ]
then
    echo "configure!"
    chmod +x ./test_conf
    # example on guinan, which has "normal" defaults for a workstation or server
    ./test_conf MPIFC=mpif90 FC=mpif90 CUDA_LIB="-L/usr/local/cuda/lib64/" MPI_INC="-I/usr/include/mpich2/" --with-cuda CUDA_INC="-I/usr/local/cuda/include"
fi
