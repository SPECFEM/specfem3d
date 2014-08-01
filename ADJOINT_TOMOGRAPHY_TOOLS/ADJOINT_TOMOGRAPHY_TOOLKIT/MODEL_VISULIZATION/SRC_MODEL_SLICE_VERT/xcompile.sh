#!/bin/sh
# Author: Hejun Zhu, hejunzhu@princeton.edu
# Princeton University, New Jersey, USA
# Last modified: Tue Jan 25 17:19:32 EST 2011


mpif90 -O3 -o xsem_model_slice sem_model_slice.f90 rthetaphi_xyz.f90 exit_mpi.f90 get_value_parameters.f90
