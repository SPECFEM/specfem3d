#!/bin/bash


$BIN_DSM/xcreate_input<<EOF
parfile_for_benchmark
EOF

source params.in
source $SCRIPTS/scrpits_specfem3D.sh
source $SCRIPTS/scripts_dsm.sh


