#!/bin/bash


#             ------------ BACTH AND SPECIFIC CLUSTER DIRECTIVES  ------
#MSUB -r run_write
#MSUB -n 2000
#MSUB -x
#MSUB -T 10800
#MSUB -q standard
#MSUB -o run_write.o
#MSUB -e run_write.e

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
# NBPROC is declared as integer (important do not change)
declare -i NPROC NPROC_MINUS_ONE

# NUMBER OF MPI PROCESSES
NPROC=32

# ENTER OPTION FOR MPIRUN
OPTION=

# do not change
NPROC_MINUS_ONE="$NPROC-1"

# log file for output
flog_file=$(pwd)/log.compute_coef

# choose the movie
PREFIX_MOVIE=velocity_Z_it

# directory where SPECFEM3D writes outputs
IN_MOVIE=$(pwd)/OUTPUT_FILES/DATABASES_MPI/

# output movie directory
OUT_MOVIE=$(pwd)/movie
#
# ------------------------ FROM HERE DO NOT CHANGE ANYTHING --------------------

# ----- load script and path ---
source params.in
source $SCRIPTS/scripts_specfem3D.sh
source $SCRIPTS/scripts_dsm.sh

## open the log file
echo > $flog_file
echo " WRITE COEF " >> $flog_file
echo >> $flog_file
echo $(date) >> $flog_file

# 2 / ----- compute DSM tractions


mkdir $DSM_tractions
cd $DSM_tractions

make_dir_exp
copy_input_files_exp
compute_exp_coeff

echo $(date) >> $flog_file

#


