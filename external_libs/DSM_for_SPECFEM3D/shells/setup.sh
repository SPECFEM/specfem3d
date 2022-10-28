#!/bin/bash


$BIN_DSM/xcreate_inputs_files<<EOF
parfile_for_benchmark
EOF

source params.in
source $SCRIPTS/scripts_specfem3D.sh
source $SCRIPTS/scripts_dsm.sh


