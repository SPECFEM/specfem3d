#!/bin/bash

###
### Absolute paths :
###

HOME_SPECFEM3D=/ccc/work/cont003/gen6351/vmont/codes/git_my_specfem3d_21_09_2015/specfem3d-devel
export HOME_SPECFEM3D

SPECFEM3D_BINARY_PATH=${HOME_SPECFEM3D}/bin
export SPECFEM3D_BINARY_PATH
HOME_AxiSEM=${HOME_SPECFEM3D}/utils/EXTERNAL_CODES_coupled_with_SPECFEM3D/AxiSEM_for_SPECFEM3D
export HOME_AxiSEM

###
### Relative path :
###

AxiSEM_FILE_1D_MODEL=${HOME_SPECFEM3D}/EXAMPLES/coupling_with_EXTERNAL_CODES/with_AxiSEM/example_1st_for_validation/INPUT_AxiSEM/prem_dsm
export AxiSEM_FILE_1D_MODEL



