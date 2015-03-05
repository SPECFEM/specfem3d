#=========== Global directives ===========
# @ job_name = ADA-test2
# @ output = $(job_name).$(step_name).$(jobid)
# @ error = $(output)

#=========== Step 1 directives ===========
#======= Sequential preprocessing ========
# @ step_name = sequential_setup
# @ job_type = serial
# @ wall_clock_limit = 3600
# @ queue

#=========== Step 2 directives ===========
#============= Parallel step =============
# @ step_name = parallel_write_coef
# @ dependency = (sequential_setup == 0)
#   (executed only if previous step completed without error)
#@ job_type   = parallel
#@ total_tasks = 2000
# @ wall_clock_limit = 72000
# @ queue

#=========== Step 3 directives ===========
#============= Parallel step =============
# @ step_name = parallel_read_coef_xmin
# @ dependency = (parallel_write_coef == 0)
#   (executed only if previous step completed without error)
#@ job_type   = parallel
#@ total_tasks = 480
# @ wall_clock_limit = 72000
# @ queue


#=========== Step 3 directives ===========
#============= Parallel step =============
# @ step_name = parallel_read_coef_xmax
# @ dependency = (parallel_write_coef == 0)
#   (executed only if previous step completed without error)
#@ job_type   = parallel
#@ total_tasks    = 480
# @ wall_clock_limit = 72000
# @ queue


#=========== Step 3 directives ===========
#============= Parallel step =============
# @ step_name = parallel_read_coef_ymin
# @ dependency = (parallel_write_coef == 0)
#   (executed only if previous step completed without error)
#@ job_type    = parallel
#@ total_tasks = 480
# @ wall_clock_limit = 72000
# @ queue


#=========== Step 3 directives ===========
#============= Parallel step =============
# @ step_name = parallel_read_coef_ymax
# @ dependency = (parallel_write_coef == 0)
#   (executed only if previous step completed without error)
#@ job_type   = parallel
#@ total_tasks = 480
# @ wall_clock_limit = 72000
# @ queue


#=========== Step 3 directives ===========
#============= Parallel step =============
# @ step_name = parallel_read_coef_zmin
# @ dependency = (parallel_write_coef == 0)
#   (executed only if previous step completed without error)
#@ job_type    = parallel
#@ total_tasks = 480
# @ wall_clock_limit = 72000
# @ queue



#=========== Step 4 directives ===========
#======= Sequential postprocessing =======
# @ step_name = sequential_postprocessing
# @ dependency = (parallel_read_coef_zmin == 0 ) && (parallel_read_coef_ymin == 0 ) && (parallel_read_coef_ymax == 0 ) && (parallel_read_coef_xmin == 0 ) && (parallel_read_coef_xmax == 0 )
#   (executed only if previous step completed without error)
# @ job_type = serial
# @ class = archive
# @ queue

# loading setup fonctions

# loading setup fonctions
SHELL_DSM=/smphome/rech/ubv/rubv002/progs/shells_specfem3D_hybrid

declare -i NPROC NPROC_MINUS_ONE

# NUMBER OF MPI PROCESSES
NPROC=32

# MPI COMMAND
MPIRUN=poe

# ENTER OPTION FOR MPIRUN
OPTION=

# do not change
NPROC_MINUS_ONE="$NPROC-1"

# log file for output
flog_file=$(pwd)/log.benchmark

BIN_DSM=/smphome/rech/ubv/rubv002/progs/DSM_FOR_SPECFEM3D/bin/

export LC_ALL=C
export MP_DEBUG_TIMEOUT_SECONDS=84000







case ${LOADL_STEP_NAME} in

  #============ Step 1 commands ============
  #======= Sequential preprocessing ========
  sequential_setup )
    set -x
    cd $TMPDIR
    printenv | sort
    # on recupere les datas
    cp ${LOADL_STEP_INITDIR}/input.data.tar.bz2 .
    tar -jxvf input.data.tar.bz2
    source $SHELL_DSM/setup.sh

    # clean and make directories SPECFEM3D
    clean_and_make_dir
    # clean and make directories DSM
    clean_and_make_dir_dsm
    # move inputs
    mv input_dsm_for_write_coef $IN_DSM/inputIASP.infTra_for_coef
    mv input_dsm_for_read_xmin  $IN_DSM/inputIASP.infTra_stxmin
    mv input_dsm_for_read_xmax  $IN_DSM/inputIASP.infTra_stxmax
    mv input_dsm_for_read_ymin  $IN_DSM/inputIASP.infTra_stymin
    mv input_dsm_for_read_ymax  $IN_DSM/inputIASP.infTra_stymax
    mv input_dsm_for_read_zmin  $IN_DSM/inputIASP.infTra_stzmin
    # create mesh
    run_create_mesh
    # create dsm dir
    mkdir $DSM_tractions
    cd $DSM_tractions
    make_dir_exp
    copy_input_files_exp
    make_dir_faces
    copy_input_files_faces
  ;;

  #============ Step 2 commands ============
  #============= Parallel step =============
  parallel_write_coef )
    set -x
    cd $TMPDIR
    source params.in
    source $SCRIPTS/scripts_specfem3D.sh
    source $SCRIPTS/scripts_dsm.sh
    cd $DSM_tractions
    compute_exp_coeff
  ;;

  #============ Step 3 commands ============
  #============= Parallel step =============
  parallel_read_coef_xmin )
    set -x
    cd $TMPDIR
    source params.in
    source $SCRIPTS/scripts_specfem3D.sh
    source $SCRIPTS/scripts_dsm.sh
    cd $DSM_tractions
    run_dsm_traction_xmin
  ;;

  #============ Step 3 commands ============
  #============= Parallel step =============
  parallel_read_coef_xmax )
    set -x
    cd $TMPDIR
    source params.in
    source $SCRIPTS/scripts_specfem3D.sh
    source $SCRIPTS/scripts_dsm.sh
    cd $DSM_tractions
    run_dsm_traction_xmax
  ;;

  #============ Step 3 commands ============
  #============= Parallel step =============
  parallel_read_coef_ymin )
    set -x
    cd $TMPDIR
    source params.in
    source $SCRIPTS/scripts_specfem3D.sh
    source $SCRIPTS/scripts_dsm.sh
    cd $DSM_tractions
    run_dsm_traction_ymin
  ;;

  #============ Step 3 commands ============
  #============= Parallel step =============
  parallel_read_coef_ymax )
    set -x
    cd $TMPDIR
    source params.in
    source $SCRIPTS/scripts_specfem3D.sh
    source $SCRIPTS/scripts_dsm.sh
    cd $DSM_tractions
    run_dsm_traction_ymax
  ;;


  #============ Step 3 commands ============
  #============= Parallel step =============
  parallel_read_coef_zmin )
    set -x
    cd $TMPDIR
    source params.in
    source $SCRIPTS/scripts_specfem3D.sh
    source $SCRIPTS/scripts_dsm.sh
    cd $DSM_tractions
    run_dsm_traction_zmin
  ;;




  #============ Step 4 commands ============
  #======= Sequential postprocessing =======
  sequential_postprocessing )
    set -x
    cd $TMPDIR
    source params.in
    source $SCRIPTS/scripts_specfem3D.sh
    source $SCRIPTS/scripts_dsm.sh
    cd $DSM_tractions


    move_output_wihtout_change

    # on ecrit les resumtats
    mfput XMIN.tar.bz2 EV5/.

    mfput XMAX.tar.bz2 EV5/.

    mfput YMIN.tar.bz2 EV5/.

    mfput YMAX.tar.bz2 EV5/.

    mfput ZMIN.tar.bz2 EV5/.



  ;;

esac
