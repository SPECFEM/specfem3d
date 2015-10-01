#!/bin/bash


#             ------------ BACTH AND SPECIFIC CLUSTER DIRECTIVES  ------

#MSUB -r Coupling_step8_xspecfem3d       # Nom du job
#MSUB -n 32
#MSUB -N 2
#MSUB -T 16000
#MSUB -q standard
#MSUB -e Xspecfem3D_32p_run.e
#MSUB -o Xspecfem3D_32p_run.o
#MSUB -A gen7165

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

### Check AxiSEM directory for previous steps

# ----------------------------------------------------
# 8 / ------- Run simulation with specfem3D
# ----------------------------------------------------

source ${HOME_AxiSEM}/scripts_SPECFEM3D_for_AxiSEM.sh

run_specfem3d_simulation

