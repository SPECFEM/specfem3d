#!/bin/bash


#             ------------ BACTH AND SPECIFIC CLUSTER DIRECTIVES  ------ 

#MSUB -r run_benchmark
#MSUB -n 12
#MSUB -N 1
#MSUB -x 
#MSUB -T 10800
#MSUB -q standard
#MSUB -o run_benchmark.o
#MSUB -e run_benchmark.e 

# set -x
# cd $BRIDGE_MSUB_PWD


#OAR -n bench_hy
#OAR -l nodes=1,walltime=2:00:00
#OAR -p cluster
##OAR -q development
#OAR -O bench_hy.%jobid%.out 
#OAR -E bench_hy.%jobid%.err


## Chargement des modules module load ... 
#module load intel/13.0
#module load openmpi/intel/1.6.3
#module list
#echo ${OAR_NODEFILE} 
#CPUS=$(wc -l ${OAR_NODEFILE} | awk '{print $1}')
#echo $CPUS


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


# NBPROC is declared as integer (important do not change)
declare -i NPROC NPROC_MINUS_ONE CPUS

# NUMBER OF MPI PROCESSES
NPROC=24
CPUS=24

# MPIRUN COMMAND 
MPIRUN="mpirun"

# ENTER OPTION FOR MPIRUN 
OPTION=" -n "${CPUS}
OPTION_SIMU=" -np "${CPUS}
#OPTION=" -np "${CPUS}"  -machinefile "${OAR_NODEFILE}" -bysocket -bind-to-core"
#OPTION_SIMU=" -np "${CPUS}"  -machinefile "${OAR_NODEFILE}" -bysocket -bind-to-core"

# do not change
NPROC_MINUS_ONE="$NPROC-1"

# log file for output 
flog_file=$(pwd)/log.benchmark

# choose the movie
PREFIX_MOVIE=velocity_Y_it

# directory where SPECFEM3D writes outputs  
IN_MOVIE=$(pwd)/OUTPUT_FILES/DATABASES_MPI/

# output movie directory 
OUT_MOVIE=$(pwd)/movie


# ----- load script and path --- 
source params.in
source $SCRIPTS/scrpits_specfem3D.sh
source $SCRIPTS/scripts_dsm.sh



# 6 / ----------------- make movie
echo "" >> $flog_file
echo " MAKE movie" >> $flog_file
echo $(date) >> $flog_file

create_movie $PREFIX_MOVIE $IN_MOVIE $OUT_MOVIE 25 8500 

# to do chane 25 and 8500
echo $(date) >> $flog_file

                               
