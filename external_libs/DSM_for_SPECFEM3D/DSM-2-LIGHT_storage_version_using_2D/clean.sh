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
echo "To adapt for the LIGHT storage version of DSM using 2D"
