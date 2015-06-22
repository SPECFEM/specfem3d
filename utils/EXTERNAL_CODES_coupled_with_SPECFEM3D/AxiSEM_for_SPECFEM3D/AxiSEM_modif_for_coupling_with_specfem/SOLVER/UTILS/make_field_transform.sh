#!/bin/bash

gfortran -Wall -pedantic -fbacktrace -c -I $HOME/local/include/ -I /usr/include field_transform.f90
gfortran field_transform.o -o field_transform -L $HOME/local/lib -lnetcdff -Wl,-rpath,$HOME/local/lib -L /usr/lib -lm -lfftw3_threads -lfftw3 -lfftw3f
