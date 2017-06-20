#!/bin/bash

export CC=gcc
export FC=gfortran
export INSTALL_DIR=$HOME/local
export LDFLAGS=-L$INSTALL_DIR/lib
export CPPFLAGS="-I$INSTALL_DIR/include -DgFortran"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$INSTALL_DIR/lib
set -e


#zlib 1.2.8
rm -rf zlib-1.2.8 zlib-1.2.8.tar.gz
wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4/zlib-1.2.8.tar.gz
tar -xvf zlib-1.2.8.tar.gz
cd zlib-1.2.8
./configure --prefix=$INSTALL_DIR
make check
make install
cd ..


#HDF 1.8.12
rm -rf hdf5-1.8.12 hdf5-1.8.12.tar.bz2
wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4/hdf5-1.8.12.tar.gz ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4/hdf5-1.8.9.tar.gz
tar -xvf hdf5-1.8.12.tar.gz
cd hdf5-1.8.12

./configure --prefix=$INSTALL_DIR --with-zlib=$INSTALL_DIR

# -j parallelizes make;  -s reduces output
make -sj
make check
make install
cd ..


# netcdf-4.3.1.1
rm -rf netcdf-4.3.1.1.tar.gz netcdf-4.3.1.1
wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.3.1.1.tar.gz
tar -xvf netcdf-4.3.1.1.tar.gz
cd netcdf-4.3.1.1
./configure --enable-netcdf-4 --enable-dap --enable-shared --prefix=$INSTALL_DIR
make -sj
make check # all tests should succeed.
make install
cd ..


#netcdf-fortran-4.4$
rm -rf netcdf-fortran-4.4-beta1.tar.gz netcdf-fortran-4.4-beta1
export LDFLAGS="-L$INSTALL_DIR/lib -lnetcdf "
wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-fortran-4.4-beta1.tar.gz
tar -xvf netcdf-fortran-4.4-beta1.tar.gz
cd netcdf-fortran-4.4-beta1
./configure --prefix=$INSTALL_DIR --enable-netcdf-4
make -sj
make check
make install
cd ..


echo 'You might want to put these to lines in your Makefile to make sure, the'
echo 'libraries just compiled are used:'
echo ''
echo 'LIBS = -lm  -L $(INSTALL_DIR)/lib -lnetcdff -Wl,-rpath,$(INSTALL_DIR)/lib'
echo 'INCLUDE = -I $(INSTALL_DIR)/include'

