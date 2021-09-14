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
echo

## gets MPI setting
#echo "ompi_info on path: $(which ompi_info)"
#echo "ompi_info:"
#ompi_info
#echo
#echo "ompi_info all:"
#ompi_info --all
#echo
#echo "ompi_info param all:"
#ompi_info --param all all --level 9
#echo
## allow for more MPI processes than cores (2-core VM nodes)
# 1. option: oversubscribe
#echo "home: $HOME"
# mca param_files points to: /home/runner/.openmpi/mca-params.conf
#mkdir -p $HOME/.openmpi
#echo "rmaps_base_oversubscribe = 1" >> $HOME/.openmpi/mca-params.conf
#echo "rmaps_base_inherit = 1" >> $HOME/.openmpi/mca-params.conf
# 2. option: increase number of slots
#echo "orte_set_default_slots = 8">> $HOME/.openmpi/mca-params.conf
#echo "orte_default_hostfile = $HOME/.openmpi/openmpi-default-hostfile" >> $HOME/.openmpi/mca-params.conf
#echo "localhost slots=8" >> $HOME/.openmpi/openmpi-default-hostfile
#or: export OMPI_MCA_orte_default_hostfile=$HOME/.openmpi/openmpi-default-hostfile

# storing updated environment parameters for following bash-script
echo "export PATH=${PATH}" > $HOME/.tmprc
echo "export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}" >> $HOME/.tmprc
echo "export CUDA_HOME=${CUDA_HOME}" >> $HOME/.tmprc

## avoids MPI issue with number of slots
#echo "export OMPI_MCA_rmaps_base_oversubscribe=1" >> $HOME/.tmprc
#echo "export OMPI_MCA_rmaps_base_inherit=1" >> $HOME/.tmprc
# uses github environment to store (for subsequent steps)
echo "OMPI_MCA_rmaps_base_oversubscribe=1" >> $GITHUB_ENV
echo "OMPI_MCA_rmaps_base_inherit=1" >> $GITHUB_ENV

## avoid MPI warnings when running in container
#echo "export OMPI_MCA_btl_vader_single_copy_mechanism=none" >> $HOME/.tmprc
#echo "export OMPI_MCA_btl=^openib" >> $HOME/.tmprc

# exports for xterm output (for make tests)
echo "TERM=xterm" >> $GITHUB_ENV

echo
echo "exports:"
export
echo

