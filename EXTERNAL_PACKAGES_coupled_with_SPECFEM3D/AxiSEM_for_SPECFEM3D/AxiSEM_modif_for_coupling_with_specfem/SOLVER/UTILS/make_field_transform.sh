#!/bin/bash

gfortran -Dunc -Wall -pedantic -fbacktrace -c  -I /usr/include field_transform.F90
gfortran field_transform.o -o field_transform -lnetcdff -lm -lfftw3_threads -lfftw3 -lfftw3f
