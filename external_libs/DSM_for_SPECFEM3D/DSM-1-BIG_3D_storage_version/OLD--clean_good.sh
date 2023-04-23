#!/bin/bash

cd Part1_create_SPECFEM3D_Cartesian_mesh_for_DSM;make clean;cd ..
cd Part2_compute_DSM_coefficients_with_FEMs_SH; make clean;cd ../
cd Part2_compute_DSM_coefficients_with_FEMs_PSV; make clean;cd ../
cd Part3_read_DSM_coefficients_back_SH/TraPSV_MPI_read_vertical_faces; make clean;cd ..
cd TraPSV_MPI_read_Zmin; make clean;cd ../..
cd Part3_read_DSM_coefficients_back_PSV/TraPSV_MPI_read_vertical_faces; make clean;cd ..
cd TraPSV_MPI_read_Zmin; make clean;cd ../..
cd Part4_modify_DSM_results_for_SPECFEM/FFT_MPI_FACES_ZMIN_SH; make clean;cd ..
cd FFT_MPI_FACES_VERT_SH; make clean;cd ..
cd FFT_MPI_FACES_ZMIN_PSV; make clean;cd ..
cd FFT_MPI_FACES_VERT_PSV; make clean;cd ..
cd FFT_MPI_FACES_ZMIN_FULL; make clean;cd ..
cd FFT_MPI_FACES_VERT_FULL; make clean;cd ..
cd ChangeFormat_zmin_disp; make clean;cd ..
cd ChangeFormat_zmin; make clean;cd ..
cd ChangeFormat_disp; make clean;cd ..
cd ChangeFormat; make clean;cd ..
cd Interf_SPECFEM3D_DSM;make clean;cd ../..
cd UTILS;make clean;cd ..

