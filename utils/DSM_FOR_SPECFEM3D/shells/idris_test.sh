#=========== Global directives ===========
# @ job_name = multi-steps
# @ output = $(job_name).$(step_name).$(jobid)
# @ error = $(output)

#=========== Step 1 directives ===========
#======= Sequential preprocessing ========
# @ step_name = sequential_preprocessing
# @ job_type = serial
# @ wall_clock_limit = 600
# @ queue

#=========== Step 2 directives ===========
#============= Parallel step =============
# @ step_name = parallel_step
# @ dependency = (sequential_preprocessing == 0)
#   (executed only if previous step completed without error)
# @ job_type = parallel
# @ total_tasks = 2
# @ wall_clock_limit = 3600
# @ queue

#=========== Step 3 directives ===========
#======= Sequential postprocessing =======
# @ step_name = sequential_postprocessing
# @ dependency = (parallel_step == 0)
#   (executed only if previous step completed without error)
# @ job_type = serial
# @ wall_clock_limit = 1200
# @ queue

source /smphome/rech/ubv/rubv002/progs/shells_specfem3D_hybrid/shells_tests.sh

case ${LOADL_STEP_NAME} in

  #============ Step 1 commands ============
  #======= Sequential preprocessing ========
  sequential_preprocessing )
    set -x
    cd $TMPDIR
    echo "Step 1" > step_1.out
    echo $TMPDIR >> step_1.out
    mkdir rep_test
    ls >> step_1.out
    pwd >> step_1.out
    test >> step_1.out
  ;;

  #============ Step 2 commands ============
  #============= Parallel step =============
  parallel_step )
    set -x
    cd $TMPDIR
    echo "parallel step" > step_2.out
    echo " pwd " >> step_2.out
    pwd >> step_2.out
    echo "ls " >> step_2.out
    ls -altrh >> step_2.out
    cd rep_test
    echo "pwd " >> ../step_2.out
    pwd >> ../step_2.out
    test >> ../step_2.out
  ;;

  #============ Step 3 commands ============
  #======= Sequential postprocessing =======
  sequential_postprocessing )
    set -x
    cd $TMPDIR
    cp step_1.out ${LOADL_STEP_INITDIR}/.
    cp step_2.out ${LOADL_STEP_INITDIR}/.
  ;;

esac
