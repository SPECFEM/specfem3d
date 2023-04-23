#!/bin/bash


#             ------------ BACTH AND SPECIFIC CLUSTER DIRECTIVES  ------

#MSUB -r Reformat_for_AxiSEM_32p         # Nom du job
#MSUB -n 32
#MSUB -N 2
#MSUB -T 8000
#MSUB -q standard
#MSUB -e reformat_for_AxiSEM_32_run.e
#MSUB -o reformat_for_AxiSEM_32_run.o
#MSUB -A gen7165

set -x
cd ${BRIDGE_MSUB_PWD}

#
## ------------------ INPUTS -----------------------------

# NBPROC is declared as integer (important do not change)
declare -i NPROC NPROC_MINUS_ONE CPUS CHOICE MIDDLE

# NUMBER OF MPI PROCESSES
NPROC=32
CPUS=32

# MPIRUN COMMAND
MPIRUN=ccc_mprun

# ENTER OPTION FOR MPIRUN
OPTION=

### Define relative path
UTILS_COUPLING=../../../UTILS_COUPLING_SpecFEM

# do not change
NPROC_MINUS_ONE="$NPROC-1"

$MPIRUN ${UTILS_COUPLING}/xreformat > OUTPUT_reformat

