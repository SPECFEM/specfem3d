#!/bin/sh
# Author: Hejun Zhu, hejunzhu@princeton.edu
# Princeton University, New Jersey, USA
# Last modified: Mon Jan 24 20:34:14 EST 2011

if [ -f xcopy_local_forward ]; then
  echo RM xcopy_local_forward
  rm xcopy_local_forward
fi

mpif90 -o xcopy_local_forward -O3  copy_local.f90

