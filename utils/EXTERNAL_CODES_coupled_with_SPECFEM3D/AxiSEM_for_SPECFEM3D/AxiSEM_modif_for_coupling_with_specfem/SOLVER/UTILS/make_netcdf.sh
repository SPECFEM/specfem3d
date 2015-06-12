#!/bin/bash

export CC=gcc
export FC=gfortran
export INSTALL_DIR=$HOME/local
export LDFLAGS=-L$INSTALL_DIR/lib
export CPPFLAGS="-I$INSTALL_DIR/include -DgFortran"
export LD_LIBRARY_PATH=$INSTALL_DIR/lib 
set -e


#zlib 1.2.7
rm -rf zlib-1.2.7 zlib-1.2.7.tar.gz
wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4/zlib-1.2.7.tar.gz
tar -xvf zlib-1.2.7.tar.gz
cd zlib-1.2.7
./configure --prefix=$INSTALL_DIR
make check
make install
cd ..


#HDF 1.8.9
rm -rf hdf5-1.8.9 hdf5-1.8.9.tar.bz2
wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4/hdf5-1.8.9.tar.gz
tar -xvf hdf5-1.8.9.tar.gz
cd hdf5-1.8.9

./configure --prefix=$INSTALL_DIR --with-zlib=$INSTALL_DIR

# -j parallelizes make;  -s reduces output
make -sj
make check 
make install
cd ..


# netcdf-4.3.0
rm -rf netcdf-4.3.0.tar.gz netcdf-4.3.0
wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.3.0.tar.gz 
tar -xvf netcdf-4.3.0.tar.gz
cd netcdf-4.3.0
./configure --enable-netcdf-4 --enable-dap --enable-shared --prefix=$INSTALL_DIR
make check # all tests should succeed.
make install
cd ..


#netcdf-fortran-4.4$
rm -rf v4.4.0-rc1.tar.gz netcdf-fortran-4.4.0-rc1
export LDFLAGS="-L$INSTALL_DIR/lib -lnetcdf "
wget https://github.com/Unidata/netcdf-fortran/archive/v4.4.0-rc1.tar.gz
tar -xvf v4.4.0-rc1.tar.gz
cd netcdf-fortran-4.4.0-rc1
cmake -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR -DCMAKE_BUILD_TYPE=Release .
make install
cd ..


echo 'You might want to put these to lines in your Makefile to make sure, the'
echo 'libraries just compiled are used:'
echo ''
echo 'LIBS = -lm  -L $(INSTALL_DIR)/lib -lnetcdff -Wl,-rpath,$(INSTALL_DIR)/lib'
echo 'INCLUDE = -I $(INSTALL_DIR)/include'

