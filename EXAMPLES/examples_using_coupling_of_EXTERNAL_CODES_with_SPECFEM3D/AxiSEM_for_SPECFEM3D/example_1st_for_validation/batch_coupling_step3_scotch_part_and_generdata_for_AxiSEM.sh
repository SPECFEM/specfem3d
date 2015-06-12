#!/bin/bash


#             ------------ BACTH AND SPECIFIC CLUSTER DIRECTIVES  ------

#MSUB -r Coupling_step3_part_and_generdata_for_AxiSEM       # Nom du job
#MSUB -n 32
#MSUB -N 2
#MSUB -T 5400
#MSUB -q standard
#MSUB -e Test_generdata_32p.e
#MSUB -o Test_generdata_32p.o
#MSUB -A ra2410

set -x
cd ${BRIDGE_MSUB_PWD}

######################################################################################################################

# NBPROC is declared as integer (important do not change)
declare -i NPROC NPROC_MINUS_ONE CPUS CHOICE MIDDLE

# NUMBER OF MPI PROCESSES
NPROC=32
CPUS=32

# MPIRUN COMMAND
MPIRUN=ccc_mprun

# ENTER OPTION FOR MPIRUN
OPTION=

# do not change
NPROC_MINUS_ONE="$NPROC-1"

# Define different paths and folders
source paths_for_coupling_SPECFEM3D_AxiSEM.sh

# ----------------------------------------------------
# 2 / ------- AxiSEM mesher and solver
# ----------------------------------------------------

### Cf AxiSEM directory

# ----------------------------------------------------
# 3 / ------- create specfem3D data base
# ----------------------------------------------------

source ${HOME_AxiSEM}/scripts_SPECFEM3D_for_AxiSEM.sh
run_create_partitioning_and_specfem_databases

