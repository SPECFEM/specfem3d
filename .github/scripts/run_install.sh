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

## parallel HDF5
if [ "${HDF5}" == "true" ]; then
  echo
  echo "HDF5 additional installation:"
  echo
  sudo apt-get install -yq --no-install-recommends libhdf5-mpi-dev
  ## checks installation paths
  #echo
  #dpkg -L libhdf5-mpi-dev
  #echo
  #dpkg -L libhdf5-openmpi-dev
  #echo
  #echo "hdf5 module paths:"
  #find /usr/ -iname 'hdf5.mod'
  #echo "hdf5 library paths:"
  #find /usr/ -iname 'libhdf5hl_fortran*'
  #echo
fi

## HIP
if [ "${HIP}" == "true" ]; then
  echo
  echo "HIP additionals installation:"
  echo
  sudo apt-get install -yq --no-install-recommends libtbb-dev
fi

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi
echo

# python3 pip upgrade might complain: "ERROR: launchpadlib 1.10.13 requires testresources"
sudo apt-get install -yq --no-install-recommends python3-testresources
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

## ADIOS2
if [ "${ADIOS2}" == "true" ]; then
  echo
  echo "ADIOS2 installation:"
  echo
  # installs cmake wget
  sudo apt-get install -yq --no-install-recommends cmake wget
  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi
  # uses /opt as installation directory
  mkdir -p /opt; cd /opt
  # download source
  wget https://github.com/ornladios/ADIOS2/archive/refs/tags/v2.10.1.tar.gz
  tar zxf v2.10.1.tar.gz
  cd ADIOS2-2.10.1/
  # build source
  mkdir -p build; cd build/
  CC=gcc CXX=g++ FC=gfortran cmake -DADIOS2_USE_Fortran=ON \
    -DADIOS2_USE_HDF5=OFF -DADIOS2_BUILD_EXAMPLES=OFF -DBUILD_TESTING=OFF \
    -DCMAKE_INSTALL_PREFIX=/opt/ADIOS2 ../
  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi
  make -j4
  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi
  make install
  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi
  # environment for directory
  echo "ADIOS2_DIR=/opt/ADIOS2" >> $GITHUB_ENV
  echo; echo "done ADIOS2"; echo
fi

# installs the CUDA toolkit
if [ "${CUDA}" == "true" ]; then
  # Linux environment
  ## distribution from ubuntu 22.04
  UBUNTU_VERSION=ubuntu2204

  # CUDA_VERSION - specifies CUDA toolkit version
  # http://developer.download.nvidia.com/compute/cuda/repos/
  CUDA_VERSION=13.1.1-1

  # default architecture amd64
  CUDA_OS=x86_64
  CUDA_ARCH=amd64
  if [ "${RUNNER_ARCH}" == "arm64" ]; then
    CUDA_OS=sbsa
    CUDA_ARCH=arm64 # ARM
  fi

  echo "Installing CUDA library"
  echo "CUDA version  : ${CUDA_VERSION}"
  echo "UBUNTU version: ${UBUNTU_VERSION}"
  echo "CUDA OS       : ${CUDA_OS}"
  echo "CUDA arch     : ${CUDA_ARCH}"

  # package needs key
  # see: https://developer.nvidia.com/blog/updating-the-cuda-linux-gpg-repository-key/
  # old:
  #sudo apt-key adv --fetch-keys http://developer.download.nvidia.com/compute/cuda/repos/${UBUNTU_VERSION}/${CUDA_OS}/7fa2af80.pub
  # new:
  # manually add new key (not recommended):
  #sudo apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/${UBUNTU_VERSION}/${CUDA_OS}/3bf863cc.pub
  #echo
  # gets packages
  #INSTALLER=cuda-repo-${UBUNTU_VERSION}_${CUDA_VERSION}_${CUDA_ARCH}.deb
  #wget http://developer.download.nvidia.com/compute/cuda/repos/${UBUNTU_VERSION}/${CUDA_OS}/${INSTALLER}
  #sudo dpkg -i ${INSTALLER}
  #echo
  # (preferred) w/ new keyring package:
  # see https://forums.developer.nvidia.com/t/notice-cuda-linux-repository-key-rotation/212772
  # if it doesn't work yet with error:
  #   E:Conflicting values set for option Signed-By regarding source
  # remove outdated key:
  sudo apt-key del 7fa2af80
  sudo sed -i '/developer\.download\.nvidia\.com\/compute\/cuda\/repos/d' /etc/apt/sources.list
  sudo rm -f /etc/apt/sources.d/cuda*.list
  sudo rm -f /etc/apt/sources.list.d/cuda.list
  sudo rm -f /etc/apt/sources.list.d/nvidia-ml.list
  # for ubuntu1804/ppc64el ../$distro/$arch/.. becomes ../${UBUNTU_VERSION}/${CUDA_OS}/..
  wget https://developer.download.nvidia.com/compute/cuda/repos/${UBUNTU_VERSION}/${CUDA_OS}/cuda-keyring_1.0-1_all.deb
  sudo dpkg -i cuda-keyring_1.0-1_all.deb
  echo

  # update
  echo "Updating libraries"
  sudo apt-get update -qq
  dpkg -l | grep cuda
  export CUDA_APT=${CUDA_VERSION:0:4}  # version 13.1
  export CUDA_APT=${CUDA_APT/./-}
  echo "CUDA: ${CUDA_APT}"  # apt version 13-1 -> package name: cuda-compiler-13-1

  # installs packages
  CUDA_PACKAGES="cuda-drivers cuda-compiler-${CUDA_APT} cuda-cudart-dev-${CUDA_APT}"
  echo "Installing ${CUDA_PACKAGES}"
  sudo apt-get install -y --no-install-recommends ${CUDA_PACKAGES}
  sudo apt-get clean
  export CUDA_HOME=/usr/local/cuda-${CUDA_VERSION:0:4}    # version 13.1
  export LD_LIBRARY_PATH=${CUDA_HOME}/lib64:${LD_LIBRARY_PATH}
  export PATH=${CUDA_HOME}/bin:${PATH}
  echo ""
  nvcc --version
  echo ""

  ## OpenCL additionals
  if [ "${OPENCL}" == "true" ]; then
    echo "OpenCL installation"
    #echo "dpkg toolkit:"
    #dpkg -l | grep toolkit
    #echo ""
    #echo "dpkg opencl:"
    #dpkg -l | grep opencl
    #echo ""
    #echo "apt-cache opencl:"
    #apt-cache search opencl
    # possible packages for OpenCL:
    #sudo apt-get install -y --no-install-recommends cuda-toolkit-${CUDA_APT}
    #sudo apt-get install opencl-headers
    # for ppc64 architecture: to be able to compile/link OpenCL version
    sudo apt-get install nvidia-opencl-dev
    echo ""
  fi
else
  export CUDA_HOME=""
fi

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

