#!/bin/bash

mkdir -p bin

###cd Part1_create_SPECFEM3D_Cartesian_mesh_for_DSM;make;cd ..
cd Part2_compute_DSM_coefficients_with_FEMs_SH; make;cd ../
cd Part2_compute_DSM_coefficients_with_FEMs_PSV; make;cd ../
cd Part3_read_DSM_coefficients_back_SH/TraPSV_MPI_read_vertical_faces; make;cd ..
cd TraPSV_MPI_read_Zmin; make;cd ../..
cd Part3_read_DSM_coefficients_back_PSV/TraPSV_MPI_read_vertical_faces; make;cd ..
cd TraPSV_MPI_read_Zmin; make;cd ../..
cd Part4_modify_DSM_results_for_SPECFEM/FFT_MPI_FACES_ZMIN_SH; make;cd ..
cd FFT_MPI_FACES_VERT_SH; make;cd ..
cd FFT_MPI_FACES_ZMIN_PSV; make;cd ..
cd FFT_MPI_FACES_VERT_PSV; make;cd ..
cd FFT_MPI_FACES_ZMIN_FULL; make;cd ..
cd FFT_MPI_FACES_VERT_FULL; make;cd ..
cd ChangeFormat_zmin_disp; make;cd ..
cd ChangeFormat_zmin; make;cd ..
cd ChangeFormat_disp; make;cd ..
cd ChangeFormat; make;cd ..
cd Interf_SPECFEM3D_DSM;make;cd ../..
cd UTILS;make;cd ..

