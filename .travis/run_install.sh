#!/bin/bash

# fortran/openMPI compiler
sudo apt-get install gfortran libgomp1 openmpi-bin libopenmpi-dev

# python script needs numpy
echo "Python on path: `which python`"
echo "pip on path: $(which pip)"
# for distribution precise
#sudo apt-get install -qq python-numpy #python-scipy
# trusty:
# (Sep2017 update: adding flag --user : see https://github.com/travis-ci/travis-ci/issues/8382)
pip install --user --upgrade pip setuptools wheel
pip install --user --only-binary=numpy numpy
# version info
python --version

# version infos
echo "compiler versions:" ${FC} ${MPIFC} ${CC}
${FC} --version
${MPIFC} --version
${CC} --version

# installs the CUDA toolkit
if [ "$CUDA" == "true" ]; then
  echo "Installing CUDA library"
  # note: travis could stall and time out here
  #       one could try to add: travis_retry sudo dpgk -i ..
  #       https://docs.travis-ci.com/user/common-build-problems/#travis_retry
  #
  # remove old nvidia-cuda packages
  #sudo apt-get remove nvidia-cuda-* ;
  # gets packages
  ## distribution trusty: from ubuntu 14.04
  wget http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1404/x86_64/cuda-repo-ubuntu1404_${CUDA_VERSION}_amd64.deb
  sudo dpkg -i cuda-repo-ubuntu1404_${CUDA_VERSION}_amd64.deb
  ## distribution precise: from ubuntu 12.04
  #wget http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1204/x86_64/cuda-repo-ubuntu1204_${CUDA_VERSION}_amd64.deb;
  #sudo dpkg -i cuda-repo-ubuntu1204_${CUDA_VERSION}_amd64.deb;
  # update
  echo "Updating libraries"
  sudo apt-get update -qq
  dpkg -l | grep cuda
  export CUDA_APT=${CUDA_VERSION:0:3}
  export CUDA_APT=${CUDA_APT/./-}
  # installs packages
  # CUDA_PACKAGES="cuda-drivers cuda-core-${CUDA_APT} cuda-cudart-dev-${CUDA_APT} cuda-cufft-dev-${CUDA_APT}";
  CUDA_PACKAGES="cuda-drivers cuda-core-${CUDA_APT} cuda-cudart-dev-${CUDA_APT}"
  echo "Installing ${CUDA_PACKAGES}"
  sudo apt-get install -y ${CUDA_PACKAGES}
  sudo apt-get clean
  export CUDA_HOME=/usr/local/cuda-${CUDA_VERSION:0:3}
  export LD_LIBRARY_PATH=${CUDA_HOME}/lib64:${LD_LIBRARY_PATH}
  export PATH=${CUDA_HOME}/bin:${PATH}
  echo ""
  nvcc --version
fi

# storing updated environment parameters for following bash-script
echo "export PATH=${PATH}" > $HOME/.tmprc
echo "export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}" >> $HOME/.tmprc
echo "export CUDA_HOME=${CUDA_HOME}" >> $HOME/.tmprc

echo "exports:"
export
echo ""
