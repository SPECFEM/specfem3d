#!/bin/bash


#             ------------ BACTH AND SPECIFIC CLUSTER DIRECTIVES  ------

#MSUB -r Interface_expand_2D_3D_160p        # Nom du job
#MSUB -n 160
#MSUB -N 10
#MSUB -T 7200
#MSUB -q standard
#MSUB -e interface_expand_2D_3D_160_run.e
#MSUB -o interface_expand_2D_3D_160_run.o
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

### Copy meshes file in result dir
cp ../input_box* .

# do not change
NPROC_MINUS_ONE="$NPROC-1"

$MPIRUN ${UTILS_COUPLING}/xexpand_2D_3D > OUTPUT_expand_2D_3D

### We submit reformat from here when expand is finished

####ccc_msub -q standard ../sub_called_batch_for_reformat_Curie_CD.sh

