#!/bin/bash


#             ------------ BACTH AND SPECIFIC CLUSTER DIRECTIVES  ------

#MSUB -r Interface_expand_2D_3D        # Nom du job
#MSUB -n 160
#MSUB -N 10
#MSUB -T 3000
#MSUB -q standard
#MSUB -e interface_expand_2D_3D_run.e
#MSUB -o interface_expand_2D_3D_run.o
#MSUB -A gen7165

set -x
cd ${BRIDGE_MSUB_PWD}

#
## ------------------ INPUTS -----------------------------

# NBPROC is declared as integer (important do not change)
declare -i NPROC NPROC_MINUS_ONE CPUS CHOICE MIDDLE

# NUMBER OF MPI PROCESSES
NPROC=160
CPUS=160

# MPIRUN COMMAND
MPIRUN=ccc_mprun

# ENTER OPTION FOR MPIRUN
OPTION=

### Define relative path
UTILS_COUPLING=../../../UTILS_COUPLING_SpecFEM

### Copy meshes and params files in result dir
cp ../input_box* .
cp ../inparam_advanced ../inparam_basic ../inparam_source ../inparam_hetero .

# do not change
NPROC_MINUS_ONE="$NPROC-1"

$MPIRUN ${UTILS_COUPLING}/xexpand_2D_3D > OUTPUT_expand_2D_3D

### We submit reformat from here when expand is finished

ccc_msub -q standard ../sub_called_batch_for_reformat_Curie_CD.sh

