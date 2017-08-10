#!/bin/bash


#             ------------ BACTH AND SPECIFIC CLUSTER DIRECTIVES  ------

#MSUB -r Convert_stations        # Nom du job
#MSUB -n 1
#MSUB -N 1
#MSUB -T 3600
#MSUB -q standard
#MSUB -e Convet_stations_run.e
#MSUB -o Convert_stations_run.o
#MSUB -A gen7165


# convetion station

source  ./paths_for_coupling_SPECFEM3D_AxiSEM.sh


#
z_bottom=6131000.00000000 # must be take in MESH/model_1D.in

${HOME_AxiSEM}/UTILS_COUPLING_SpecFEM/create_stations.x<<EOF
$z_bottom
EOF
