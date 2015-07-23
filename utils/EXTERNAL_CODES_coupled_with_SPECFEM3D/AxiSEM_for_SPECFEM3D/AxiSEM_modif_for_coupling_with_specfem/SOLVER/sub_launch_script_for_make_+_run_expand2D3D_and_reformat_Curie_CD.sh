#!/bin/bash

current_dir=$(pwd)

cd ../../../../../EXAMPLES/coupling_with_EXTERNAL_CODES/with_AxiSEM/example_1st_for_validation/DATA/
mkdir AxiSEM_tractions/
mkdir AxiSEM_tractions/1/
cd ${current_dir}


dir_for_run=$(pwd)/$1

####
#### We compile expand_2D_3D and reformat for interface
####

cd ../../UTILS_COUPLING_SpecFEM/

make clean

make all

####
#### We submit expand_2D_3D from here, which will submit reformat at the end
####

cd ${dir_for_run}

ccc_msub -q standard ../sub_called_batch_for_expand2D3D_Curie_CD.sh

