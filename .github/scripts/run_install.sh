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
#sudo apt-get install -qq python-numpy # not working, likely installed on older python version
pip install --user --upgrade pip setuptools wheel
pip install --user --only-binary=numpy numpy

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

# MPI
# github actions uses for Linux virtual machines a 2-core CPU environment
# see: https://docs.github.com/en/actions/using-github-hosted-runners/about-github-hosted-runners#supported-runners-and-hardware-resources
#
# job issue with mpirun -np 4 .. command:
#
#    There are not enough slots available in the system to satisfy the 4
#    slots that were requested by the application:
#
echo "MPI environment:"
echo "mpif90 on path: $(which mpif90)"
echo "ompi_info on path: $(which ompi_info)"
echo
echo "ompi_info:"
ompi_info
echo
echo "ompi_info all:"
ompi_info --all
echo
echo "ompi_info param all:"
ompi_info --param all
echo
echo "ompi_info param_files:"
ompi_info --param param_files
echo

#export OMPI_MCA_orte_default_hostfile
# mca param_files points to: /home/runner/.openmpi/mca-params.conf
echo "home: $HOME"
mkdir -p $HOME/.openmpi
echo "localhost slots=8" >> $HOME/.openmpi/mca-params.conf
echo
#echo "OMPI_MCA_param_files=$HOME/.openmpi/mca-params.conf" >> $GITHUB_ENV
#echo "OMPI_MCA_orte_set_default_slots=8" >> $GITHUB_ENV
#
#echo "export OMPI_MCA_btl_vader_single_copy_mechanism=none" >> $HOME/.tmprc
#echo "export OMPI_MCA_btl=^openib" >> $HOME/.tmprc


# exports
export TERM=xterm
echo "exports:"
export
echo 

