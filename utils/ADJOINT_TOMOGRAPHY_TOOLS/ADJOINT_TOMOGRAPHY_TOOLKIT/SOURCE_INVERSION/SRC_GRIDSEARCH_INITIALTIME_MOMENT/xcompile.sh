#!/bin/sh
# Author: Hejun Zhu, hejunzhu@princeton.edu
# Princeton University, New Jersey, USA
# Last modified: Wed Sep 12 11:17:06 EDT 2012


if [ -f xgridsearch_time_moment ]; then
  echo RM xgridsearch_time_moment
  rm xgridsearch_time_moment
fi

gfortran -o xgridsearch_time_moment -O3 gridsearch_time_moment.f90 -L/home/hejunzhu/BIN/sac-101.4/lib -lsacio
#gfortran -o xgridsearch_time_moment -O3 gridsearch_time_moment.f90 -L${SACHOME}/lib -lsacio
