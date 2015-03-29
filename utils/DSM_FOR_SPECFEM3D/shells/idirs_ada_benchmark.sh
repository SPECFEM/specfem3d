#=========== Global directives ===========
# @ job_name = multi-steps-benchmark
# @ output = $(job_name).$(step_name).$(jobid)
# @ error = $(output)

#=========== Step 1 directives ===========
#======= Sequential preprocessing ========
# @ step_name = sequential_stetup
# @ job_type = serial
# @ wall_clock_limit = 3600
# @ queue

#=========== Step 2 directives ===========
#============= Parallel step =============
# @ step_name = parallel_write_coef
# @ dependency = (sequential_setup == 0)
#   (executed only if previous step completed without error)
# @ job_type = parallel
# @ total_tasks = 12
# @ wall_clock_limit = 7200
# @ queue

#=========== Step 3 directives ===========
#============= Parallel step =============
# @ step_name = parallel_read_coef
# @ dependency = (parallel_write_coef == 0)
#   (executed only if previous step completed without error)
# @ job_type = parallel
# @ total_tasks = 12
# @ wall_clock_limit = 7200
# @ queue

#=========== Step 4 directives ===========
#======= Sequential postprocessing =======
# @ step_name = sequential_postprocessing
# @ dependency = (parallel_write_coef == 0)
#   (executed only if previous step completed without error)
# @ job_type = serial
# @ wall_clock_limit = 3600
# @ queue

# loading setup fonctions
SHELL_DSM=/smphome/rech/ubv/rubv002/progs/shells_specfem3D_hybrid
source $SHELL_DSM/setup.sh

case ${LOADL_STEP_NAME} in

  #============ Step 1 commands ============
  #======= Sequential preprocessing ========
  sequential_preprocessing )
    set -x
    cd $TMPDIR
    cp ${LOADL_STEP_INITDIR}/input.data.tar.bz2 .
    tar -jxvf input.data.tar.bz2
    setup_process
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
    compute_exp_coeff
  ;;

  #============ Step 3 commands ============
  #============= Parallel step =============
  parallel_read_coef )
    set -x
    cd $TMPDIR
    cd $DSM_tractions
    run_dsm_traction
  ;;

  #============ Step 4 commands ============
  #======= Sequential postprocessing =======
  sequential_postprocessing )
    set -x
    cd $TMPDIR

    mfput $REP/tractxmin.bin .
    mfput $REP/velxmin.bin .

    mfput $REP/tractxmax.bin .
    mfput $REP/velxmax.bin .

    mfput $REP/tractymin.bin .
    mfput $REP/velymin.bin .

    mfput $REP/tractymax.bin .
    mfput $REP/velymax.bin .

    mfput $REP/tractzmin.bin .
    mfput $REP/velzmin.bin .

  ;;

esac
