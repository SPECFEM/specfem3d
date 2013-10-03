#!/bin/bash

# this script must be executed from the GEOCUBIT base directory with the command
#   source setpaths.sh

PWD=`pwd`
export GEOCUBITHOME=$PWD
echo base directory for GEOCUBIT is $GEOCUBITHOME

# check paths to CUBIT
# note: CUBIT 12.2 will not work with hex27
#       CUBIT 13.0 has a bug associated with merging meshes
#       CUBIT 14.0 is under a new license and name (trelis)
#                  and has not been fully tested
# we recommend using CUBIT 12.2
echo checking paths for CUBIT
echo CUBITHOME = $CUBITHOME
echo CUBITDIR = $CUBITDIR
echo PYTHONPATH = $PYTHONPATH
echo LD_LIBRARY_PATH = $LD_LIBRARY_PATH
echo PATH = $PATH

# set paths
# note: if you need to load a python module, then it should be done beforehand
echo setting paths for GEOCUBIT
export PYTHONPATH=$PYTHONPATH:$GEOCUBITHOME
export PATH=$PATH:$GEOCUBITHOME

