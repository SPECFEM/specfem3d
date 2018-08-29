#!/bin/sh

module unload intel
module load intel/17.2
module load intelmpi/2017.1.132
module load hdf5/1.8.17

h5pfc combine_signal_files_into_a_single_one.f90 -o combine
