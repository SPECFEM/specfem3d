#!/bin/bash

rm *~

echo " "
echo "Cleaning from the HOME_DSM_MAIN_DIR directory"
echo " "

cd Part-1a-compute_DSM_coefficients_with_FEMs_SH/ ; rm *~ ; make clean

cd Part-1b-compute_DSM_coefficients_with_FEMs_PSV/ ; rm *~ ; make clean

cd Part-2a-read_DSM_coefficients_back_SH/TraPSV_MPI_read_vertical_faces/ ; rm *~ ; make clean
cd Part-2a-read_DSM_coefficients_back_SH/TraPSV_MPI_read_Zmin/ ; rm *~ ; make clean

cd Part-2b-read_DSM_coefficients_back_PSV/TraPSV_MPI_read_vertical_faces/ ; rm *~ ; make clean
cd Part-2b-read_DSM_coefficients_back_PSV/TraPSV_MPI_read_Zmin/ ; rm *~ ; make clean

cd Part-3-modify_DSM_results_for_SPECFEM/FFT_MPI_FACES_ZMIN_SH ; rm *~ ; make clean
cd Part-3-modify_DSM_results_for_SPECFEM/FFT_MPI_FACES_VERT_SH ; rm *~ ; make clean
cd Part-3-modify_DSM_results_for_SPECFEM/FFT_MPI_FACES_ZMIN_PSV ; rm *~ ; make clean
cd Part-3-modify_DSM_results_for_SPECFEM/FFT_MPI_FACES_VERT_PSV ; rm *~ ; make clean
cd Part-3-modify_DSM_results_for_SPECFEM/FFT_MPI_FACES_ZMIN_FULL ; rm *~ ; make clean
cd Part-3-modify_DSM_results_for_SPECFEM/FFT_MPI_FACES_VERT_FULL ; rm *~ ; make clean
cd Part-3-modify_DSM_results_for_SPECFEM/ChangeFormat_zmin_disp ; rm *~ ; make clean
cd Part-3-modify_DSM_results_for_SPECFEM/ChangeFormat_zmin ; rm *~ ; make clean
cd Part-3-modify_DSM_results_for_SPECFEM/ChangeFormat_disp ; rm *~ ; make clean
cd Part-3-modify_DSM_results_for_SPECFEM/ChangeFormat ; rm *~ ; make clean
cd Part-3-modify_DSM_results_for_SPECFEM/Interf_SPECFEM3D_DSM ; rm *~ ; make clean

cd UTILS/ ; rm *~ ; make clean 

cd ..
rm bin/*

