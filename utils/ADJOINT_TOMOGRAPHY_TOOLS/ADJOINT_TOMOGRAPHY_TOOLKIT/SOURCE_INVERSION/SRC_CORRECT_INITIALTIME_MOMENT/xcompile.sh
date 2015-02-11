#!/bin/sh
# Author: Hejun Zhu, hejunzhu@princeton.edu
# Princeton University, New Jersey, USA
# Last modified: Mon Jan 24 20:34:14 EST 2011

if [ -f xcorrect_syn_time_moment ]; then
  echo RM xcorrect_syn_time_moment
  rm xcorrect_syn_time_moment
fi


gfortran -o xcorrect_syn_time_moment -O3 correct_syn_time_moment.f90 -L/home/hejunzhu/BIN/sac-101.4/lib -lsacio

