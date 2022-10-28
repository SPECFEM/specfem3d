#             ------------ BACTH AND SPECIFIC CLUSTER DIRECTIVES  ------
#=========== Global directives ===========
# @ job_name = benchmark_hybrid
# @ output = $(job_name).$(step_name).$(jobid)
# @ error = $(output)

#=========== Step 1 directives ===========
#======= Sequential preprocessing ========
# @ step_name = sequential_preprocessing
# @ job_type = serial
# @ wall_clock_limit = 3600
# @ queue

#=========== Step 2 directives ===========
#============= Parallel step =============
# @ step_name = parallel_step_compute_coef
# @ dependency = (sequential_preprocessing == 0)
#   (executed only if previous step completed without error)
# @ job_type = parallel
# @ total_tasks = 12
# @ wall_clock_limit = 3600
# @ queue

#=========== Step 3 directives ===========
#============= Parallel step_2 =============
# @ step_name = parallel_step_read_coef
# @ dependency = (sequential_preprocessing == 0)
#   (executed only if previous step completed without error)
# @ job_type = parallel
# @ total_tasks = 12
# @ wall_clock_limit = 36000
# @ queue

#=========== Step 4 directives ===========
#======= Sequential postprocessing =======
# @ step_name = sequential_postprocessing
# @ dependency = (parallel_step == 0)
#   (executed only if previous step completed without error)
# @ job_type = serial
# @ wall_clock_limit = 36000
# @ queue



#
# Repertoire temporaire de travail
cd $TMPDIR

# La variable LOADL_STEP_INITDIR est automatiquement positionnee par
# LoadLeveler au repertoire dans lequel on tape la commande llsubmit
cp  $LOADL_STEP_INITDIR/input_data.tar.bz2 .
tar -jxvf input_data.tar.bz2

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
$HOME_SPECFEM3D/utils/DSM_FOR_SPECFEM3D/bin/xcreate_inputs_files<<EOF
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
cp $IN_DSM/iasp91 $MESH/.


## open the log file
echo > $flog_file
echo " BENCHMARK RUN  " >> $flog_file
echo >> $flog_file
echo $(date) >> $flog_file


# 1 / ------- create mesh
run_create_mesh




# 2 / ----- compute DSM tractions
run_dsm_traction




# 3 / ------- create specfem3D data base
run_create_specfem_databases




# 4 / -------- create tractions for specfem3D from DSM
echo "" >> $flog_file
echo " create traction database" >> $flog_file
echo $(date) >> $flog_file

run_create_tractions_for_specfem

echo $(date) >> $flog_file





# 5 / --------------- run simulation
echo "" >> $flog_file
echo " simulation" >> $flog_file
echo $(date) >> $flog_file

run_simu

echo $(date) >> $flog_file




# 6 / ----------------- make movie
echo "" >> $flog_file
echo " MAKE movie" >> $flog_file
echo $(date) >> $flog_file

create_movie $PREFIX_MOVIE $IN_MOVIE $OUT_MOVIE 25 8500

# to do chane 25 and 8500
echo $(date) >> $flog_file


