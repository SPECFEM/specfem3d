#!/bin/bash


#             ------------ BACTH AND SPECIFIC CLUSTER DIRECTIVES  ------

#MSUB -r Coupling_step1_meshfem3d_for_AxiSEM        # Nom du job
#MSUB -n 1
#MSUB -N 1
#MSUB -T 1000
#MSUB -q standard
#MSUB -e Coupling_step1_meshfem3d_for_AxiSEM_run.e
#MSUB -o Coupling_step1_meshfem3d_for_AxiSEM_run.o
#MSUB -A ra2410

set -x
cd ${BRIDGE_MSUB_PWD}

######################################################################################################################
#
#
## ------------------ INPUTS -----------------------------

# NBPROC is declared as integer (important do not change)
declare -i NPROC NPROC_MINUS_ONE CPUS CHOICE MIDDLE

# NUMBER OF MPI PROCESSES
NPROC=1
CPUS=1

# MPIRUN COMMAND
MPIRUN=ccc_mprun

# ENTER OPTION FOR MPIRUN
OPTION=

###OPTION=" -np "${NPROC}
###OPTION_SIMU=" -np "${CPUS}
###OPTION=" -np "${CPUS}"  -machinefile "${OAR_NODEFILE}" -bysocket -bind-to-core"
###OPTION_SIMU=" -np "${CPUS}"  -machinefile "${OAR_NODEFILE}" -bysocket -bind-to-core"

# do not change
NPROC_MINUS_ONE="$NPROC-1"

# Define different paths and folders
mkdir -p MESH/
source paths_for_coupling_SPECFEM3D_AxiSEM.sh

# ----------------------------------------------------
# 1 / ------- create mesh
# ----------------------------------------------------
source ${HOME_SPECFEM3D}/scripts_SPECFEM3D_for_AxiSEM.sh
run_create_mesh

