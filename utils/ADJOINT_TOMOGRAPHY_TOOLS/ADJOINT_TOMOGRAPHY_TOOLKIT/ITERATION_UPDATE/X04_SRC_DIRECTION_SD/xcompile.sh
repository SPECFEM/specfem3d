#!/bin/sh
# Author: Hejun Zhu, hejunzhu@princeton.edu
# Princeton University, New Jersey, USA
# Last modified: Tue Jan 25 17:19:32 EST 2011

if [ ! -f ../../SHARE_FILES/HEADER_FILES/constants.h ]; then
  echo WRONG! NO constants.h in SHARE_FILES
  exit
fi
if [ ! -f ../../SHARE_FILES/HEADER_FILES/values_from_mesher.h ]; then
  echo WRONG! NO values_from_mesher.h in SHARE_FILES
  exit
fi
if [ ! -f ../../SHARE_FILES/HEADER_FILES/precision.h ]; then
  echo WRONG! NO precision.h in SHARE_FILES
  exit
fi

mpif90 -O3 -o xcompute_direction_sd compute_direction_sd.f90 exit_mpi.f90
