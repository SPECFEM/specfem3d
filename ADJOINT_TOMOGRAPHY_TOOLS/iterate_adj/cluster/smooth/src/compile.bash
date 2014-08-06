#!/bin/bash

flags="-O3  -implicitnone -warn argument_checking -warn unused -warn declarations  -check nobounds"

mpif90 $flags -o smooth_sem_fun smooth_sem_fun.f90 smooth_sub.f90 gll_library.f90 exit_mpi.f90

#ifort $flags -o smooth_sem_fun smooth_specfem_function.f90 smooth_sub.f90 gll_library.f90 exit_mpi.f90

rm *.o

