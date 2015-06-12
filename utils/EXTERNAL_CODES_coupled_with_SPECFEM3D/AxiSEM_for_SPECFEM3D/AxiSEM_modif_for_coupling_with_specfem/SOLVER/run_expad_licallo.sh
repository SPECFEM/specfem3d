#!/bin/bash 
#OAR -l /nodes=4,walltime=2:00:00
#OAR -n expand_24_cores
##OAR -p gpu='NO' AND ibpool='FDR'

source /softs/env-intel15-impi50.sh

SEMUTILS=/home/monteill/codes/GIT/specfem3d/branches/devel/utils/Coupling_with_AxiSEM/UTILS_COUPLING_SpecFEM
#20
mpiexec.hydra -f $OAR_FILE_NODES -bootstrap-exec oarsh -perhost 12 $SEMUTILS/xexpand_2D_3D

