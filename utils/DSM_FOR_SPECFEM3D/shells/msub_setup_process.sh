#!/bin/bash


#             ------------ BACTH AND SPECIFIC CLUSTER DIRECTIVES  ------

#MSUB -r run_mesh
#MSUB -n 1
#MSUB -T 10800
#MSUB -q standard
#MSUB -o run_mesh.o
#MSUB -e run_mesh.e

set -x
cd $BRIDGE_MSUB_PWD

#

######################################################################################################################
#
#
#      BENCHMARK FOR HYBRID DSM/SPECFEM3D METHOD
#
# INPUTS :
#
#   1/ input directoy : ./input_dsm
#      containts
#             -- Double_para.txt
#             -- FrqsMpi.txt
#             -- iasp91
#             -- iasp91_dsm
#             -- st
#
#   2/ input file : parfile_for_benchmark
#
#
#   3/ SPECFEM3D input directory : ./DATA
#      containts
#             -- Par_file
#             -- STATIONS
#             -- CMTSOLUTION
#
#
# the script runs :
#
#   1/ MESHER
#   2/ DSM to compute tractions on the chunk boundary
#   3/ SCHOTCH + CREATE DATABASE FOR SPECFEM3D
#   4/ ADD DSM TRACTION TO THE SPECFEM3D DATABASE
#   5/ RUN SPECFEM3D
#   6/ MAKE MOVIE
#
#
#
#  Vadim Monteiller April 2013.
#
# reference :
# "A hybrid method to compute short-period synthetic seismograms of teleseismic body waves in a 3-D regional model"
# Monteiller, V; Chevrot, S; Komatitsch, D; Fuji, N
# GEOPHYSICAL JOURNAL INTERNATIONAL Volume:192 Issue:1 Pages:230-247 DOI:10.1093/gji/ggs006 Published: JAN 2013
#
#####################################################################################################################


## ------------------ INPUTS -----------------------------

# DSM BINARY : (to do supprimer peut-etre ca de params.in??)
BIN_DSM=$HOME_SPECFEM3D/utils/DSM_FOR_SPECFEM3D/bin
# NBPROC is declared as integer (important do not change)
declare -i NPROC NPROC_MINUS_ONE

# NUMBER OF MPI PROCESSES
NPROC=32

# ENTER OPTION FOR MPIRUN
OPTION=

# do not change
NPROC_MINUS_ONE="$NPROC-1"

# log file for output
flog_file=$(pwd)/log.benchmark

# choose the movie
PREFIX_MOVIE=velocity_Z_it

# directory where SPECFEM3D writes outputs
IN_MOVIE=$(pwd)/OUTPUT_FILES/DATABASES_MPI/

# output movie directory
OUT_MOVIE=$(pwd)/movie

#------- input files creation
# you must write the absolute path for : xcreate_input
# you must edit and complete : parfile_for_benchmark
$BIN_DSM/xcreate_inputs_files<<EOF
parfile_for_benchmark
EOF


#
# ------------------------ FROM HERE DO NOT CHANGE ANYTHING --------------------

# ----- load script and path ---
source params.in
source $SCRIPTS/scripts_specfem3D.sh
source $SCRIPTS/scripts_dsm.sh

# clean and make directories SPECFEM3D
clean_and_make_dir

# clean and make directories DSM
clean_and_make_dir_dsm

# mv some input files in rigth place
mv input_dsm_for_write_coef $IN_DSM/inputIASP.infTra_for_coef
mv input_dsm_for_read_xmin  $IN_DSM/inputIASP.infTra_stxmin
mv input_dsm_for_read_xmax  $IN_DSM/inputIASP.infTra_stxmax
mv input_dsm_for_read_ymin  $IN_DSM/inputIASP.infTra_stymin
mv input_dsm_for_read_ymax  $IN_DSM/inputIASP.infTra_stymax
mv input_dsm_for_read_zmin  $IN_DSM/inputIASP.infTra_stzmin
# copy model file
cp $IN_DSM/ak135 $MESH/.


## open the log file
echo > $flog_file
echo " set up  " >> $flog_file
echo >> $flog_file
echo $(date) >> $flog_file


# 1 / ------- create mesh
run_create_mesh

echo $(date) >> $flog_file

