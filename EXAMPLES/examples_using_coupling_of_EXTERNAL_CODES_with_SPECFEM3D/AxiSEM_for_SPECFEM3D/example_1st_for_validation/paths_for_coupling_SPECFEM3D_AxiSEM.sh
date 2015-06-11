#!/bin/bash

###
### Absolute paths : 
###

HOME_SPECFEM3D=/ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d
export HOME_SPECFEM3D

SPECFEM3D_BINARY_PATH=${HOME_SPECFEM3D}/bin
export SPECFEM3D_BINARY_PATH
HOME_AxiSEM=${HOME_SPECFEM3D}/utils/EXTERNAL_CODES_coupled_with_SPECFEM3D/AxiSEM_for_SPECFEM3D
export HOME_AxiSEM
AxiSEM_BINARY_PATH=${HOME_SPECFEM3D}/utils/EXTERNAL_CODES_coupled_with_SPECFEM3D/AxiSEM_for_SPECFEM3D/bin
export AxiSEM_BINARY_PATH

###
### Relative path :
###

AxiSEM_FILE_1D_MODEL=${HOME_SPECFEM3D}/EXAMPLES/examples_using_coupling_of_EXTERNAL_CODES_with_SPECFEM3D/AxiSEM_for_SPECFEM3D/example_1st_for_validation/INPUT_AxiSEM/iasp91_dsm
export AxiSEM_FILE_1D_MODEL



