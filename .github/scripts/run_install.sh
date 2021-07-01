#!/bin/bash
#
# script to install needed packages
#

# updates repository
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 6B05F25D762E3157
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 78BD65473CB3BD13
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 762E3157
sudo apt-get update

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# fortran/openMPI compiler
sudo apt-get install -yq --no-install-recommends gfortran g++ openmpi-bin libopenmpi-dev

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi
echo 

# python script needs numpy
sudo apt-get install -qq python-numpy

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi
echo

# version info
echo "Python on path: $(which python)"
python --version
echo
echo "pip on path   : $(which pip)"
pip --version
echo
echo "numpy version : "
python -c "import numpy; print(numpy.__version__)"
echo

# compiler infos
echo "compiler versions:"
echo "gcc --version"
gcc --version
echo "gfortran --version"
gfortran --version
echo "mpif90 --version"
mpif90 --version
echo 

# exports
export TERM=xterm
echo "exports:"
export
echo 

