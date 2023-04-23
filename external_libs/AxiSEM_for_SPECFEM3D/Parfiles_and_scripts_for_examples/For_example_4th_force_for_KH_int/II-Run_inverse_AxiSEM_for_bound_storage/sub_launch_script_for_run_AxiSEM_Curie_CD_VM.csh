#!/bin/csh -f

set pwd_noblank = `echo $PWD | sed 's/ //g'`
test "$pwd_noblank" != "$PWD" && echo "ERROR: your path contains a blank, please rename" && exit
set homedir = $PWD

#
##### For coupling with specfem3d #####

./add_line_number_in_points_lists.sh
echo " "
echo "Add line number in lists of points DONE"
echo " "

#######################################
#

if ( ${#argv} < 1 || "$1" == "-h" ) then
    echo ""
    echo "=================================================="
    echo " Argument 1:  directory name for the simulation"
    echo ""
    echo " Optional arguments after directory name: "
    echo " -q <queue>, where <queue> can be:"
    echo "      'lsf': submit to lsf queue using bsub"
    echo "       Default: submit locally"
    echo "=================================================="
    echo""
    exit
else if ( -d $1) then
    echo " ERROR: Run directory" $1 "exists....... its content:"
    ls $1
    exit
endif

#@TODO What happens if not defined? Checkint these in parameters.F90 is quite late then...
set datapath = `grep "^DATA_DIR" inparam_advanced  |awk '{print $2}'| sed 's/\"//g'`
set infopath = `grep "^INFO_DIR" inparam_advanced |awk '{print $2}'| sed 's/\"//g'`
set meshdir = "MESHES/"`grep "^MESHNAME" inparam_basic | awk '{print $2}'`
set mpiruncmd = `grep "^MPIRUN" ../make_axisem.macros | awk '{print $3}'`
set serial = `grep "^SERIAL" ../make_axisem.macros | awk '{print $3}'`

#Check whether NetCDF is requested and whether the code is compiled with it
set netcdf_compiled = `grep "^USE_NETCDF" ../make_axisem.macros | awk '{print $3}'`
set netcdf_requested = `grep "^USE_NETCDF" inparam_advanced |awk '{print $2}'| sed 's/\"//g'`
if ( $netcdf_requested == 'true' && $netcdf_compiled != 'true') then
  echo "NetCDF compiled  (../make_axisem.macros): " $netcdf_compiled
  echo "NetCDF requested (inparam_advanced):      " $netcdf_requested
  echo "ERROR: NetCDF is requested in inparam_advanced, but disabled in ../make_axisem.macros"
  exit
endif

set gitversion = `git describe --dirty --abbrev=4 --always --tags`
echo $gitversion "GIT_VERSION"  > runinfo
set username = `whoami`
echo $username "USER_NAME" >> runinfo
set hostname = `hostname`
echo $hostname "HOST_NAME" >> runinfo
set FFLAGS = `grep "^FFLAGS" ../make_axisem.macros`
echo $FFLAGS  >> runinfo
set CFLAGS = `grep "^CFLAGS" ../make_axisem.macros`
echo $CFLAGS >> runinfo
set LDFLAGS = `grep "^LDFLAGS" ../make_axisem.macros`
echo $LDFLAGS >> runinfo

if ( -d $meshdir) then
  echo "Using mesh " $meshdir
else
  echo "ERROR: Mesh " $meshdir " not found."
  echo "Available meshes:"
  ls MESHES
  exit
endif

if ( -d $datapath) then
    if (`ls $datapath -1 | wc -l` == 0) then
        echo " Saving data into $datapath"
    else
        echo "ERROR: $datapath is not empty. Its content:"
        ls $datapath
        exit
    endif
endif

set bgmodel = `grep ^BACKGROUND_MODEL $meshdir/inparam_mesh | awk '{print $2}'`

# Since the compiling does not depend on mesh_params.h anymore, we just copy it here anyway.
# actually mesh_params.h is not needed anymore by the solver, just keeping it for
# informational purposes
echo 'copying mesh_params.h from ' $meshdir
cp $meshdir/mesh_params.h .

# Check arguments: source types and submission queues
set newqueue = 'false'
if ( "$2" == '-q') then
    set queue = $3
    set newqueue = 'true'
  echo "Submitting to queue type" $queue
endif

set multisrc = 'false'
set simtype = `grep "^SIMULATION_TYPE" inparam_basic |awk '{print $2}'`
set srcfile = 'inparam_source'

if ( $simtype == 'single') then
    set multisrc = 'false'
else if ( $simtype == 'force') then
    set multisrc = 'true'
else if ( $simtype == 'moment') then
    set multisrc = 'true'
else
    echo "ERROR: unknown simulation type: " $simtype
    exit
endif

# Run make to see whether the code has to be rebuilt and if so, do it.
# If 'make' returns an Error (something >0), then exit.

if ! { make -j }  then
  echo "ERROR: Compilation failed, please check the errors."
  exit
endif


if ( ! -f $homedir/$srcfile ) then
    echo "ERROR: Source file $srcfile does not exist"
    exit
endif

# identify receiver input file
set rec_file_type = `grep "^RECFILE_TYPE" inparam_basic |awk '{print $2}'`
echo "Receiver file type:" $rec_file_type

if ( $rec_file_type == 'colatlon' ) then
    set recfile = 'receivers.dat'
else if ( $rec_file_type == 'stations' ) then
    set recfile = 'STATIONS'
else if ( $rec_file_type == 'database' ) then
    set recfile = 'database'
    echo "this is a dummy database receiver file" >! $homedir/$recfile
endif

if ( ! -f $homedir/$recfile ) then
    echo "ERROR: Receiver file $recfile does not exist"
    exit
endif
echo "Source file:" $srcfile, "Receiver file:" $recfile

if ( $multisrc == 'true' ) then
    # multiple simulations
    echo "setting up multiple simulations for full" $simtype "source type"
    if ( $simtype == 'moment' ) then
        set srcapp = ( MZZ MXX_P_MYY MXZ_MYZ MXY_MXX_M_MYY )
        set srctype  = ( "mrr" "mtt_p_mpp" "mtr" "mtp" )
        set srcdepth = `grep "depth: " $homedir/CMTSOLUTION  |awk '{print $2}'`
        set srclat   = `grep "latitude: " $homedir/CMTSOLUTION  |awk '{print $2}'`
        set srclon   = `grep "longitude: " $homedir/CMTSOLUTION  |awk '{print $2}'`

    else if ( $simtype == 'force' ) then
        set srcapp   = ( PZ PX )

        set forcetype = `grep "^SOURCE_TYPE" $srcfile |awk '{print $2}'`
        if ( $forcetype == 'thetaforce' ) then
            set srctype  = ( "vertforce" "thetaforce" )
        else if ( $forcetype == 'phiforce' ) then
            set srctype  = ( "vertforce" "phiforce" )
        endif

        # TODO hardcoded for testing. need to define an input file for force sources!
##        set srcdepth = '0.0'
##        set srclat   = '90.0'
##        set srclon   = '0.0'
        set srcdepth = `grep "^SOURCE_DEPTH" $srcfile  |awk '{print $2}'`
        set srclat   = `grep "^SOURCE_LAT" $srcfile  |awk '{print $2}'`
        set srclon   = `grep "^SOURCE_LON" $srcfile  |awk '{print $2}'`
        set srcampl  = `grep "^SOURCE_AMPLITUDE" $srcfile  |awk '{print $2}'`

    else
        echo " ERROR: Unrecognized source type" $srctype
        echo " Choose either 'moment', or leave blank for one simulation as in inparam_source"
        exit
    endif

else if ( $multisrc == 'false' ) then
    # one simulation
    set srctype = `grep "^SOURCE_TYPE" $srcfile  |awk '{print $2}'`
    set srcapp = ( "./"  )
endif

echo 'source names:' $srcapp
echo 'source components:' $srctype

echo 'Create the run directory ' $1
mkdir $1


# Copy the make_axisem.macros file, in which the exact compiler settings are stored
echo 'copying make_axisem.macros from ../'
cp ../make_axisem.macros $1

# Copy inparam_mesh, just for archival purposes
cp $meshdir/inparam_mesh $1

cd $1
set mainrundir = $PWD

# make sure moment tensor is copied correctly
cp -p $homedir/$srcfile $mainrundir/


# Prepare and copy relevant files for each simulation
set i = 0
foreach isim  (${srcapp})

    @ i ++

    echo ""
    echo "Setting up simulation" $isim
    # construct different source file for each simulation
    if  ( $multisrc == 'true' ) then
        echo "constructing separate source files for" $isim

        echo 'SOURCE_TYPE'  $srctype[$i]  >  $srcfile.$isim
        echo 'SOURCE_DEPTH' $srcdepth     >> $srcfile.$isim
        echo 'SOURCE_LAT'   $srclat       >> $srcfile.$isim
        echo 'SOURCE_LON'   $srclon       >> $srcfile.$isim
        if ( $simtype == 'moment' ) then
            echo 'SOURCE_AMPLITUDE  1.E20'    >> $srcfile.$isim
        else if ( $simtype == 'force' ) then
            echo 'SOURCE_AMPLITUDE' $srcampl    >> $srcfile.$isim
        endif

        mkdir $isim
        cd $isim
    endif


    if ( $datapath == './Data' ) then
        mkdir $datapath
    else
        if ( $multisrc == 'true' ) then
            set datapath_isim = $datapath/$isim
            echo "creating $datapath_isim"
        else
            set datapath_isim = $datapath
        endif
        mkdir -p $datapath_isim
        ln -s $datapath_isim ./Data
    endif

    if ( -d $infopath) then
        echo " saving info into $infopath"
    else
        echo "creating $infopath"
        mkdir $infopath
    endif

    mkdir Code
    cp -Lp $homedir/*.c   Code
    cp -Lp $homedir/*.f90 Code
    cp -Lp $homedir/*.F90 Code
    cp -Lp $homedir/Makefile Code

    echo "copying crucial files for the simulation..."

    if ( $multisrc == 'true' ) then
        mv ../$srcfile.$isim $srcfile
    else
        cp $homedir/$srcfile $srcfile
    endif

    cp $homedir/axisem .
    cp $homedir/mesh_params.h .
    cp $homedir/runinfo .
    cp $homedir/$recfile .
    cp $homedir/inparam_basic .
    cp $homedir/inparam_advanced .
    cp $homedir/inparam_hetero .

    if ( $multisrc == 'false' ) then
        ln -s ../$meshdir/ Mesh
    else
        ln -s ../../$meshdir/ Mesh
    endif

    if ( $bgmodel == 'external' ) then
        cp Mesh/external_model.bm .
    endif
    cd $mainrundir

end

cp $homedir/mesh_params.h .
cp $homedir/inparam_basic .
cp $homedir/inparam_advanced .
cp $homedir/inparam_hetero .
cp $homedir/input_box.txt .       ### VM VM ici on copie en dur input_box
cp $homedir/input_box.txt ../.


if ( $simtype == 'moment' ) then
    # TODO could also be forces
    cp $homedir/CMTSOLUTION .
endif

# write a script that runs fieldtransform in all rundirs
if ( $netcdf_requested == 'true') then
    if ( $simtype == 'moment' ) then
        cp ../UTILS/fieldtransform_moment.sh fieldtransform.sh
        chmod +x fieldtransform.sh
    else if ( $simtype == 'force' ) then
        cp ../UTILS/fieldtransform_force.sh fieldtransform.sh
        chmod +x fieldtransform.sh
    endif
endif

########################################################
######### submit the jobs ##############################
########################################################

set nodnum = `grep nproc_mesh $homedir/mesh_params.h |awk '{print $6}'`
echo "preparing job on $nodnum nodes..."

foreach isim (${srcapp})
    cd $isim

    if ( $multisrc == 'true' ) then
        set outputname = "OUTPUT_"`echo $isim |sed 's/\//_/g'`
    else
        set outputname = "OUTPUT_"`echo $1 |sed 's/\//_/g'`
    endif

    if ( $newqueue == 'true' ) then

        set jobname = `echo $1 |sed 's/\//_/g'`"_"`echo $isim |sed 's/\//_/g'`


        ########## LSF SCHEDULER ######################
        if ( $queue == 'lsf' ) then
            # for Brutus: http://brutuswiki.ethz.ch/brutus/OpenMPI#Issues_when_Using_Many_Cores
            #unset OMPI_MCA_btl_openib_receive_queues
            #bsub -R "rusage[mem=2048]" -n $nodnum -W 167:59 $mpiruncmd -n $nodnum ./axisem > $outputname &
            bsub -R "rusage[mem=2048]" -I -n $nodnum $mpiruncmd -n $nodnum ./axisem 2>&1 > $outputname &

        ######## slurm  #######
        else if ( $queue == 'slurmlocal' ) then
          aprun -n $nodnum ./axisem >& $outputname &

        else if ( $queue == 'slurm' ) then

            set ntaskspernode = 32
            echo "ntaskspernode = $ntaskspernode"

            echo '#\!/bin/bash -l'                          >  sbatch.sh
            echo "#SBATCH --ntasks=$nodnum"                 >> sbatch.sh
            echo "#SBATCH --ntasks-per-node=$ntaskspernode" >> sbatch.sh
            echo "#SBATCH --time=00:59:00"                  >> sbatch.sh

            echo "module load slurm"                        >> sbatch.sh

            echo 'echo "The current job ID is $SLURM_JOB_ID"'           >> sbatch.sh
            echo 'echo "Running on $SLURM_JOB_NUM_NODES nodes"'         >> sbatch.sh
            echo 'echo "Using $SLURM_NTASKS_PER_NODE tasks per node"'   >> sbatch.sh
            echo 'echo "A total of $SLURM_NTASKS tasks is used"'        >> sbatch.sh

            echo  'aprun -n $SLURM_NTASKS ./axisem >& '$outputname      >> sbatch.sh

            sbatch sbatch.sh

        ######## TORQUE/MAUI SCHEDULER #######
        else if ( $queue == 'torque' ) then
          # this is a crazy line, but with pure integer division its hard to handle.
            #set nodes = `echo ${nodnum} | awk '{printf "%.0f\n", $1/16+0.49}'`

            echo "# Sample PBS for parallel jobs" > run_solver.pbs
            echo "#PBS -l nodes=$nodnum,walltime=7:59:00" >> run_solver.pbs
            #echo "#PBS -l nodes=${nodes}:ppn=16" >> run_solver.pbs
            echo "ulimit -s unlimited " >> run_solver.pbs
            echo "cd $PWD " >> run_solver.pbs
            echo "$mpiruncmd -n ${nodnum} $PWD/axisem  > $outputname " >> run_solver.pbs
            qsub run_solver.pbs

        ############### SuperMUC ###################
        else if ( $queue == 'SuperMUC') then
            set current_dir=$PWD
            @ nnodes = ($nodnum / 16)
            echo "# Job file for AxiSEM run, followed by field_transform" > job.cmd
            echo "#@ job_name = $jobname"                                >> job.cmd
            echo " "                                                     >> job.cmd
            echo "# JOB STEP SOLVER"                                     >> job.cmd
            echo "#@ step_name = SOLVER "                                >> job.cmd
            echo '#@ output = job_$(jobid).out '                         >> job.cmd
            echo '#@ error = job_$(jobid).err  '                         >> job.cmd
            echo "#@ job_type = parallel "                               >> job.cmd
            echo "#@ class = general "                                   >> job.cmd
            echo "#@ total_tasks=$nodnum "                               >> job.cmd
            echo "#@ node = $nnodes "                                    >> job.cmd
            echo "#@ island_count = 1"                                   >> job.cmd
            echo "#@ network.MPI = sn_all,not_shared,us "                >> job.cmd
            echo "#@ wall_clock_limit =00:10:00"                         >> job.cmd
            echo "#@ initialdir = $current_dir"                          >> job.cmd
            echo "#@ executable = $current_dir/exe_solver.sh"            >> job.cmd
            echo "#@ notification=always"                                >> job.cmd
            echo "#@ notify_user = MAILADRESS"                           >> job.cmd
            echo "#@ energy_policy_tag = Axisem_Solver  "                >> job.cmd
            echo "#@ minimize_time_to_solution = yes    "                >> job.cmd
            echo "#@ queue "                                             >> job.cmd
            echo " "                                                     >> job.cmd
            echo "# JOB STEP FIELD_TRANSFORM"                            >> job.cmd
            echo "#@ step_name = FIELD_TRANSFORM"                        >> job.cmd
            echo "#@ dependency = (SOLVER == 0)"                         >> job.cmd
            echo '#@ output = job_$(jobid).out '                         >> job.cmd
            echo '#@ error = job_$(jobid).err  '                         >> job.cmd
            echo "#@ class = micro   "                                   >> job.cmd
            echo "#@ total_tasks=1 "                                     >> job.cmd
            echo "#@ node = 1 "                                          >> job.cmd
            echo "#@ wall_clock_limit =48:00:00"                         >> job.cmd
            echo "#@ initialdir = $current_dir"                          >> job.cmd
            echo "#@ executable = $current_dir/exe_FT.sh"                >> job.cmd
            echo "#@ notification=always"                                >> job.cmd
            echo "#@ notify_user = MAILADRESS"                           >> job.cmd
            echo "#@ energy_policy_tag = Axisem_FT  "                    >> job.cmd
            echo "#@ minimize_time_to_solution = yes    "                >> job.cmd
            echo "#@ queue "                                             >> job.cmd

            # Create Solver executable script
            echo ". /etc/profile"                                         > exe_solver.sh
            echo ". /etc/profile.d/modules.sh"                           >> exe_solver.sh
            echo "module load mpi.ibm"                                   >> exe_solver.sh
            echo "module load netcdf/mpi/4.3"                            >> exe_solver.sh
            echo "module load fortran/intel"                             >> exe_solver.sh
            echo "poe ./axisem > $outputname "                           >> exe_solver.sh

            # Create Field transform executable script
            echo ". /etc/profile"                                         > exe_FT.sh
            echo ". /etc/profile.d/modules.sh"                           >> exe_FT.sh
            echo "module load mpi.ibm"                                   >> exe_FT.sh
            echo "module load netcdf/mpi/4.3"                            >> exe_FT.sh
            echo "module load fortran/intel"                             >> exe_FT.sh
            echo "../xfield_transform > OUTPUT_FT "                      >> exe_FT.sh
            llsubmit job.cmd
        endif

    ######## SUBMIT LOCALLY #######
    else
        #ulimit -s unlimited
        #setenv OMP_NUM_THREADS 4
        unlimit stacksize
        if ( $serial == 'true' ) then
           ./axisem >& $outputname &
        else if ( $serial == 'false' ) then ## VM VM add batch for Curie
            #$mpiruncmd -n $nodnum ./axisem >& $outputname &
             cp $homedir/sub_called_batch_for_AxiSEM_Curie_CD_VM.sh .
             echo " launch command"
             echo  \$"MPIRUN ./axisem > $outputname "
             echo  \$"MPIRUN ./axisem > $outputname " >> sub_called_batch_for_AxiSEM_Curie_CD_VM.sh
             ccc_msub -q standard ./sub_called_batch_for_AxiSEM_Curie_CD_VM.sh
        else
            echo 'ERROR: value for SERIAL in make_axisem.macros should be either "true" or "false"'
            echo "SERIAL = $serial"
            exit
        endif
    endif

    echo "Job running in directory $isim"
    cd $mainrundir
end


######## post processing ##################################################

cd $homedir
cd $1

cp -p $homedir/UTILS/nc_postroutines.F90 .
cp -p $homedir/UTILS/post_processing.csh .
cp -p $homedir/UTILS/post_processing.F90 .
cp -p $homedir/UTILS/xpost_processing .

cp -p $homedir/UTILS/xfield_transform .
cp -p $homedir/UTILS/field_transform.F90 .

cp -p $homedir/UTILS/plot_recfile_seis.csh .
cp -p $homedir/UTILS/plot_recs.plot .
cp -p $homedir/UTILS/taup_allrec.csh .
cp -p $homedir/UTILS/plot_record_section.m .

echo "To convolve and sum seismograms, run ./post_processing.csh after the simulations in:"
echo $mainrundir
echo ".... the post-processing input file param_post_processing is generated in the solver"
echo ".... based on guesses. Edit please."
echo " ~ ~ ~ ~ ~ ~ ~ h a n g   o n   &   l o o s e ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~"


