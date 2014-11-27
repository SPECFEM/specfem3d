#!/bin/bash


#             ------------ BACTH AND SPECIFIC CLUSTER DIRECTIVES  ------ 

#MSUB -r Benchmark_SPECFEM3D_DSM_HEX27_32p        # Nom du job  
#MSUB -n 32
#MSUB -N 2
#MSUB -T 5400
#MSUB -q standard
#MSUB -e Benchmark_SPECFEM3D_DSM_HEX27_32p_run.e
#MSUB -o Benchmark_SPECFEM3D_DSM_HEX27_32p_run.o
#MSUB -A ra2410 
      
set -x
cd ${BRIDGE_MSUB_PWD}


######################################################################################################################
#
#
#      BENCHMARK FOR HYBRID DSM/SPECFEM3D METHOD
#      
# INPUTS :
#
#   1/ input directory : ./input_dsm
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
#
### ----- First thing to do : the home of SPECFEM3D is in parmans.in ---------
##source params.in

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
 
###OPTION=" -np "${NPROC}
###OPTION_SIMU=" -np "${CPUS}
###OPTION=" -np "${CPUS}"  -machinefile "${OAR_NODEFILE}" -bysocket -bind-to-core"
###OPTION_SIMU=" -np "${CPUS}"  -machinefile "${OAR_NODEFILE}" -bysocket -bind-to-core"

# do not change
NPROC_MINUS_ONE="$NPROC-1"

# log file for output 
flog_file=$(pwd)/log.benchmark

# Define the home of specfem3d
HOME_SPECFEM3D=/ccc/scratch/cont003/gen7165/durochtc/Codes/SPECFEM3Ds/specfem3d

# choose the movie
PREFIX_MOVIE=velocity_Z_it

# directory where SPECFEM3D writes outputs  
IN_MOVIE=$(pwd)/OUTPUT_FILES/DATABASES_MPI/

# output movie directory 
OUT_MOVIE=$(pwd)/movie/

#------- input files creation 
# you must write the absolute path for : xcreate_input
# you must edit and complete : parfile_for_benchmark  
${HOME_SPECFEM3D}/utils/DSM_FOR_SPECFEM3D/bin/xcreate_inputs_files<<EOF
parfile_for_benchmark
EOF

###echo '!!!!!!!!!!!!!!!!!! SHELLS STEP1 : fin de lecture parfile_for_benchmark !!!!!!!!!!!!!!!!'

# CHOOSE the computation type. CHOICE =1/2/3 means SH/PSV/FULL wavefield computation
CHOICE=3
if [ $CHOICE -gt 3  ]
  then
  MIDDLE=3
  echo 'Wrong CHOICE definition, Reset the CHOICE to 3'
elif [ $CHOICE -lt 1  ]
  then
  MIDDLE=1
  echo 'Wrong CHOICE definition, Reset the CHOICE to 1'
else
  MIDDLE=$CHOICE
fi
CHOICE=$MIDDLE
echo 'The value of CHOICE variable is' $CHOICE
echo 'The value of CHOICE variable is' $CHOICE >  $flog_file


#
# ------------------------ FROM HERE DO NOT CHANGE ANYTHING --------------------

# ----- load script and path --- 
source params.in
source $SCRIPTS/scripts_specfem3D.sh
###echo '!!!!!!!!!!!!!!!!!! SHELLS STEP2 : fin de lecture scripts_specfem3D.sh !!!!!!!!!!!!!!!!'
if [ $CHOICE -eq 1 ]
 then
 source $SCRIPTS/scripts_dsm_SH.sh
elif [ $CHOICE -eq 2 ]
 then
 source $SCRIPTS/scripts_dsm_PSV.sh
else
 source $SCRIPTS/scripts_dsm_full.sh
fi

# clean and make directories SPECFEM3D
clean_and_make_dir

# clean and make directories DSM
clean_and_make_dir_dsm

# mv some input files in right place
mv input_dsm_for_write_coef $IN_DSM/inputIASP.infTra_for_coef
mv input_dsm_for_read_xmin  $IN_DSM/inputIASP.infTra_stxmin
mv input_dsm_for_read_xmax  $IN_DSM/inputIASP.infTra_stxmax
mv input_dsm_for_read_ymin  $IN_DSM/inputIASP.infTra_stymin
mv input_dsm_for_read_ymax  $IN_DSM/inputIASP.infTra_stymax
mv input_dsm_for_read_zmin  $IN_DSM/inputIASP.infTra_stzmin
# copy model file 
cp $IN_DSM/iasp91 $MESH/.


## open the log file 
echo >> $flog_file
echo " BENCHMARK RUN  " >> $flog_file
echo >> $flog_file
echo $(date) >> $flog_file

# ----------------------------------------------------
# 1 / ------- create mesh
# ----------------------------------------------------
run_create_mesh


# ----------------------------------------------------
# 2 / ----- compute DSM tractions
# ----------------------------------------------------
cd $DSM_tractions
make_dir_exp
copy_input_files_exp
compute_exp_coeff
run_dsm_traction


# ----------------------------------------------------
# 3 / ------- create specfem3D data base
# ----------------------------------------------------
echo "" >> $flog_file
echo $(date) >> $flog_file
run_create_specfem_databases
echo " create specfem3D data base" >> $flog_file
echo $(date) >> $flog_file


# ----------------------------------------------------
# 4 / -------- create tractions for specfem3D from DSM
# ----------------------------------------------------
echo "" >> $flog_file
echo " create traction database" >> $flog_file
echo $(date) >> $flog_file

run_create_tractions_for_specfem

echo $(date) >> $flog_file


# ----------------------------------------------------
# 5 / --------------- run simulation
# ----------------------------------------------------
echo "" >> $flog_file
echo " simulation" >> $flog_file
echo $(date) >> $flog_file

run_simu

echo $(date) >> $flog_file




#### 6 / ----------------- make movie
###echo "" >> $flog_file
###echo " MAKE movie" >> $flog_file
###echo $(date) >> $flog_file
###create_movie $PREFIX_MOVIE $IN_MOVIE $OUT_MOVIE 25 8500 

#### to do chane 25 and 8500
###echo $(date) >> $flog_file

                               
