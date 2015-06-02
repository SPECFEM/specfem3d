#!/bin/bash

vert='\e[0;32m'
blanc='\e[0;37m'
jaune='\e[1;33m'
neutre='\e[0;m'

rm *~

HOME_DSM_MAIN_DIR=/ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/EXTERNAL_CODES_coupled_with_SPECFEM3D/DSM_for_SPECFEM3D/
export HOME_DSM_MAIN_DIR

echo " "
echo "Cleaning this directory, and remove the directory bin/ in HOME_DSM_MAIN_DIR"
echo -e "Check the path of HOME_DSM_MAIN_DIR in clean.sh, currently defined as: ${jaune}\033[1m${HOME_DSM_MAIN_DIR}\033[0m${neutre}"
echo " "

cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-1a-compute_DSM_coefficients_with_FEMs_SH/ ; rm *~ ; make clean

cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-1b-compute_DSM_coefficients_with_FEMs_PSV/ ; rm *~ ; make clean

cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-2a-read_DSM_coefficients_back_SH/TraPSV_MPI_read_vertical_faces/ ; rm *~ ; make clean
cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-2a-read_DSM_coefficients_back_SH/TraPSV_MPI_read_Zmin/ ; rm *~ ; make clean

cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-2b-read_DSM_coefficients_back_PSV/TraPSV_MPI_read_vertical_faces/ ; rm *~ ; make clean
cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-2b-read_DSM_coefficients_back_PSV/TraPSV_MPI_read_Zmin/ ; rm *~ ; make clean

cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-3-modify_DSM_results_for_SPECFEM/FFT_MPI_FACES_ZMIN_SH ; rm *~ ; make clean
cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-3-modify_DSM_results_for_SPECFEM/FFT_MPI_FACES_VERT_SH ; rm *~ ; make clean
cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-3-modify_DSM_results_for_SPECFEM/FFT_MPI_FACES_ZMIN_PSV ; rm *~ ; make clean
cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-3-modify_DSM_results_for_SPECFEM/FFT_MPI_FACES_VERT_PSV ; rm *~ ; make clean
cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-3-modify_DSM_results_for_SPECFEM/FFT_MPI_FACES_ZMIN_FULL ; rm *~ ; make clean
cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-3-modify_DSM_results_for_SPECFEM/FFT_MPI_FACES_VERT_FULL ; rm *~ ; make clean
cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-3-modify_DSM_results_for_SPECFEM/ChangeFormat_zmin_disp ; rm *~ ; make clean
cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-3-modify_DSM_results_for_SPECFEM/ChangeFormat_zmin ; rm *~ ; make clean
cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-3-modify_DSM_results_for_SPECFEM/ChangeFormat_disp ; rm *~ ; make clean
cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-3-modify_DSM_results_for_SPECFEM/ChangeFormat ; rm *~ ; make clean
cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-3-modify_DSM_results_for_SPECFEM/Interf_SPECFEM3D_DSM ; rm *~ ; make clean

cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/UTILS/ ; rm *~ ; make clean

rm -rf ${HOME_DSM_MAIN_DIR}/bin/

cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version

