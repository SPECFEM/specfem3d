#!/bin/bash

vert='\e[0;32m'
blanc='\e[0;37m'
jaune='\e[1;33m'
neutre='\e[0;m'

echo " "

HOME_DSM_MAIN_DIR=/ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d/utils/EXTERNAL_CODES_coupled_with_SPECFEM3D/DSM_for_SPECFEM3D/
export HOME_DSM_MAIN_DIR

echo -e "!! CAUTION !! Verify the definition of your ${jaune}\033[1mHOME_DSM_MAIN_DIR\033[0m${neutre} in install.sh, currently defined as: "
echo " "
echo -e "${jaune}\033[1m${HOME_DSM_MAIN_DIR}\033[0m${neutre}"

echo " "
echo "The file install.sh is written for the DSM 3D big storage version"

mkdir -p bin

echo " "
echo "============================================================================================================================"
echo "============================================================================================================================"
echo " "
echo "The mesher part is now generated in meshfem3d (as a particular case)"

echo " "
cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-1a-compute_DSM_coefficients_with_FEMs_SH/ ; make

cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-1b-compute_DSM_coefficients_with_FEMs_PSV/ ; make

cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-2a-read_DSM_coefficients_back_SH/TraPSV_MPI_read_vertical_faces/ ; make
cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-2a-read_DSM_coefficients_back_SH/TraPSV_MPI_read_Zmin/ ; make

cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-2b-read_DSM_coefficients_back_PSV/TraPSV_MPI_read_vertical_faces/ ; make
cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-2b-read_DSM_coefficients_back_PSV/TraPSV_MPI_read_Zmin/ ; make

cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-3-modify_DSM_results_for_SPECFEM/FFT_MPI_FACES_ZMIN_SH ; make
cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-3-modify_DSM_results_for_SPECFEM/FFT_MPI_FACES_VERT_SH ; make
cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-3-modify_DSM_results_for_SPECFEM/FFT_MPI_FACES_ZMIN_PSV ; make
cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-3-modify_DSM_results_for_SPECFEM/FFT_MPI_FACES_VERT_PSV ; make
cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-3-modify_DSM_results_for_SPECFEM/FFT_MPI_FACES_ZMIN_FULL ; make
cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-3-modify_DSM_results_for_SPECFEM/FFT_MPI_FACES_VERT_FULL ; make
cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-3-modify_DSM_results_for_SPECFEM/ChangeFormat_zmin_disp ; make
cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-3-modify_DSM_results_for_SPECFEM/ChangeFormat_zmin ; make
cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-3-modify_DSM_results_for_SPECFEM/ChangeFormat_disp ; make
cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-3-modify_DSM_results_for_SPECFEM/ChangeFormat ; make
cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/Part-3-modify_DSM_results_for_SPECFEM/Interf_SPECFEM3D_DSM ; make

cd ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/UTILS/ ; make

mv ${HOME_DSM_MAIN_DIR}/DSM-1-BIG_3D_storage_version/bin ${HOME_DSM_MAIN_DIR}

cd ${HOME_DSM_MAIN_DIR}

