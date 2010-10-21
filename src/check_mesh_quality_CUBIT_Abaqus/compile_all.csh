#!/bin/sh

rm -f *.o  ./xcheck_mesh_quality_CUBIT_Abaqus ./xconvert_skewness_to_angle ./xmultiply_CUBIT_Abaqus_mesh_by_1000

gfortran -Wall -o xcheck_mesh_quality_CUBIT_Abaqus check_mesh_quality_CUBIT_Abaqus.f90
gfortran -Wall -o xconvert_skewness_to_angle convert_skewness_to_angle.f90
gfortran -Wall -o xmultiply_CUBIT_Abaqus_mesh_by_1000 multiply_CUBIT_Abaqus_mesh_by_1000.f90

