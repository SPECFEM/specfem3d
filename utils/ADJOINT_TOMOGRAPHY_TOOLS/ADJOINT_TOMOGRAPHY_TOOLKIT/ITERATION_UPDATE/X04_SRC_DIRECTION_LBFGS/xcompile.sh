#!/bin/sh
# Author: Hejun Zhu, hejunzhu@princeton.edu
# Princeton University, New Jersey, USA
# Last modified: Tue Jan 25 17:19:32 EST 2011

if [ ! -f ../../SHARE_FILES/HEADER_FILES/constants.h ]; then
  echo WRONG! NO constants.h
  exit
fi
if [ ! -f ../../SHARE_FILES/HEADER_FILES/values_from_mesher.h ]; then
  echo WRONG! NO values_from_mesher.h
  exit
fi
if [ ! -f ../../SHARE_FILES/HEADER_FILES/precision.h ]; then
  echo WRONG! NO precision.h
  exit
fi

mpif90 -O3 -o xcompute_direction_lbfgs compute_direction_lbfgs.f90 exit_mpi.f90
