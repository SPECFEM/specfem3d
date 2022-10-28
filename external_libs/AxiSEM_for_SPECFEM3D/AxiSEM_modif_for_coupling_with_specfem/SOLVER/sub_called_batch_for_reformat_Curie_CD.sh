#!/bin/bash


#             ------------ BACTH AND SPECIFIC CLUSTER DIRECTIVES  ------

#MSUB -r Reformat_for_AxiSEM_224p         # Nom du job
#MSUB -n 224
#MSUB -N 14
#MSUB -T 3600
#MSUB -q standard
#MSUB -e reformat_for_AxiSEM_224p_run.e
#MSUB -o reformat_for_AxiSEM_224p_run.o
#MSUB -A gen7165

set -x
cd ${BRIDGE_MSUB_PWD}

#
## ------------------ INPUTS -----------------------------

# NBPROC is declared as integer (important do not change)
declare -i NPROC NPROC_MINUS_ONE CPUS CHOICE MIDDLE

# NUMBER OF MPI PROCESSES
NPROC=224
CPUS=224

# MPIRUN COMMAND
MPIRUN=ccc_mprun

# ENTER OPTION FOR MPIRUN
OPTION=

### Define relative path
UTILS_COUPLING=../../../UTILS_COUPLING_SpecFEM

# do not change
NPROC_MINUS_ONE="$NPROC-1"

$MPIRUN ${UTILS_COUPLING}/xreformat > OUTPUT_reformat

