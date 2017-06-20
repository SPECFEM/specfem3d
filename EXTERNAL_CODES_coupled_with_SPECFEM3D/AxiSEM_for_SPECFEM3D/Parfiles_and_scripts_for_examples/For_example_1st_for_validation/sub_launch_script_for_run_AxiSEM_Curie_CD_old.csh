#!/bin/csh -f

set pwd_noblank = `echo $PWD | sed 's/ //g'`
test "$pwd_noblank" != "$PWD" && echo "ERROR: your path contains a blank, please rename" && exit
set homedir = $PWD

### Go to l. 369 to see command for Curie ###
#
##### For coupling with specfem3d #####

./add_line_number_in_points_lists.sh
echo " "
echo "Add line number in lists of points DONE"
echo " "

#######################################
#

if ( ${#argv} < 1 || "$1" == "-h" ) then

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
    echo " Run or directory" $1 "exists....... its content:"
    ls $1
    exit
endif

#@TODO What happens if not defined? Checkint these in parameters.F90 is quite late then...
set datapath = `grep "^DATA_DIR" inparam_advanced  |awk '{print $2}'| sed 's/\"//g'`
set infopath = `grep "^INFO_DIR" inparam_advanced |awk '{print $2}'| sed 's/\"//g'`
set meshdir = "MESHES/"`grep "^MESHNAME" inparam_basic | awk '{print $2}'`
set mpiruncmd = `grep "^MPIRUN" ../make_axisem.macros | awk '{print $3}'`

#Check whether NetCDF is requested and whether the code is compiled with it
set netcdf_compiled = `grep "^USE_NETCDF" ../make_axisem.macros | awk '{print $3}'`
set netcdf_requested = `grep "^USE_NETCDF" inparam_advanced |awk '{print $2}'| sed 's/\"//g'`
if ( $netcdf_requested == 'true' && $netcdf_compiled != 'true') then
  echo "NetCDF compiled  (../make_axisem.macros): " $netcdf_compiled
  echo "NetCDF requested (inparam_advanced):      " $netcdf_requested
  echo "ERROR: NetCDF is requested in inparam_advanced, but disabled in ../make_axisem.macros"
  exit
endif

set svnrevision = `svnversion`
echo $svnrevision "SVN_VERSION      " > runinfo
set username = `whoami`
echo $username "USER_NAME        " >> runinfo
set hostname = `hostname`
echo $hostname "HOST_NAME        " >> runinfo
set FFLAGS = `grep "^FFLAGS" ../make_axisem.macros`
echo $FFLAGS  >> runinfo
set CFLAGS = `grep "^CFLAGS" ../make_axisem.macros`
echo $CFLAGS >> runinfo
set LDFLAGS = `grep "^LDFLAGS" ../make_axisem.macros`
echo $LDFLAGS >> runinfo

if ( -d $meshdir) then
  echo "Using mesh " $meshdir
else
  echo "Mesh " $meshdir " not found."
  echo "Available meshes:"
  ls MESHES
  exit
endif

set bgmodel = `grep ^BACKGROUND_MODEL $meshdir/inparam_mesh | awk '{print $2}'`

if ( ! -f inparam_hetero) then
  cp inparam_hetero.TEMPLATE inparam_hetero
endif

# if the mesh has different mesh_params.h, copy here
if ( ! -f mesh_params.h || `diff mesh_params.h $meshdir/mesh_params.h | wc -l` != "0" ) then
  echo 'copying mesh_params.h from ' $meshdir
  cp $meshdir/mesh_params.h .
endif

# if the mesh has different background_models.f90, copy over
if ( `diff background_models.f90 $meshdir/background_models.f90 | wc -l` != "0" ) then
  echo 'copying background_models.f90 from ' $meshdir
  cp $meshdir/background_models.f90 .
endif

# Check arguments: source types and submission queues
set newqueue = 'false'
if ( "$2" == '-q') then
    set queue = $3
    set newqueue = 'true'
endif

set multisrc = 'false'

# @TODO grep is not stable if SIMULATION_TYPE is there twice, e.g. in a comment line!!

set srctype = `grep "^SIMULATION_TYPE" inparam_basic |awk '{print $2}'`
set src_file_type = 'sourceparams'
set srcfile = 'inparam_source'

if ( $srctype == 'single') then
    set multisrc = 'false'
else if ( $srctype == 'force') then
    set multisrc = 'true'
else if ( $srctype == 'moment') then
    set multisrc = 'true'
endif

if ( $newqueue == 'true' ) then
  echo "Submitting to queue type" $queue
endif

# Run make to see whether the code has to be rebuilt and if so, do it.
# If 'make' returns an Error (something >0), then exit.

# MvD: I do not get this: it checks == 0 where 0 is the status when exited without
# problems???
if ( { make -j } == 0 ) then
  echo "Compilation failed, please check the errors."
  exit
endif


if ( ! -f $homedir/$srcfile ) then
    echo "file $srcfile does not exist"
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
    echo "file $recfile does not exist"
    exit
endif
echo "Source file:" $srcfile, "Receiver file:" $recfile

set num_src = 1
set num_src_arr = ( 1 )
if ( $multisrc == 'true' ) then
    # multiple simulations
    echo "setting up multiple simulations for full" $srctype "source type"
    if ( $srctype == 'moment' ) then
        set mij_sourceparams = ( 0. 0. 0. 0. 0. 0. )
        set map_mij = ( 1 2 4 6 )
        set numsim = 4
        set srcapp = ( MZZ MXX_P_MYY MXZ_MYZ MXY_MXX_M_MYY )
        set srctype  = ( "mrr" "mtt_p_mpp" "mtr" "mtp" )
        set srcdepth = `grep "depth: " $homedir/CMTSOLUTION  |awk '{print $2}'`
        set srclat   = `grep "latitude: " $homedir/CMTSOLUTION  |awk '{print $2}'`
        set srclon   = `grep "longitude: " $homedir/CMTSOLUTION  |awk '{print $2}'`

    else if ( $srctype == 'force' ) then
        set numsim   = 2
        set srcapp   = ( PZ PX )
        set srctype  = ( "vertforce" "xforce" )

    else
        echo " Unrecognized source type" $srctype
        echo " Choose either 'moment', 'force', or leave blank for one simulation as in inparam_source"
        exit
    endif

else if ( $multisrc == 'false' ) then
    # one simulation
    set numsim = 1;
    set srctype = `grep "^SOURCE_TYPE" $srcfile  |awk '{print $2}'`
    set srcapp = ( "./"  )
endif

echo 'source names:' $srcapp
echo 'source components:' $srctype

mkdir $1
cd $1
set mainrundir = $PWD

# make sure moment tensor is copied correctly
cp -p $homedir/$srcfile $mainrundir/



# Prepare and copy relevant files for each simulation
foreach isrc (${num_src_arr})
    set i = 0
    foreach isim  (${srcapp})

        @ i ++

        set num = 6
        echo ""
        echo "Setting up simulation" $isim
        # construct different source file for each simulation
        if  ( $multisrc == 'true' ) then
            echo "constructing separate source files for" $isim

            echo 'SOURCE_TYPE'  $srctype[$i]  >  $srcfile.$isrc.$isim
            echo 'SOURCE_DEPTH' $srcdepth     >> $srcfile.$isrc.$isim
            echo 'SOURCE_LAT'   $srclat       >> $srcfile.$isrc.$isim
            echo 'SOURCE_LON'   $srclon       >> $srcfile.$isrc.$isim
            echo 'SOURCE_AMPLITUDE  1.E20'    >> $srcfile.$isrc.$isim
        endif

        if ( $multisrc == 'false' ) then
            set simdir = './'
        else
            if ( $num_src == 1 ) then
                set simdir = $isim
                mkdir $simdir
                cd $simdir
            else
                set simdir = $isrc"_"$isim
                mkdir $simdir
                cd $simdir
            endif
        endif

        if ( -d $datapath) then
            echo " Saving data into $datapath"
        else
            echo "creating $datapath"
            mkdir $datapath
        endif

        if ( -d $infopath) then
            echo " saving info into $infopath"
        else
            echo "creating $infopath"
            mkdir $infopath
        endif

        mkdir Code
        cp -p $homedir/*.f90 Code
        cp -p $homedir/*.F90 Code
        cp -p $homedir/Makefile Code

        echo "copying crucial files for the simulation..."

        if ( $multisrc == 'true' ) then
            mv ../$srcfile.$isrc.$isim $srcfile
        else
            cp $homedir/$srcfile $srcfile
        endif

        cp $homedir/axisem .
        cp $homedir/mesh_params.h .
        #cp $homedir/mesh_params.dat .
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

        cp $homedir/mesh_params.h .
        cp $homedir/inparam_basic .
        cp $homedir/inparam_advanced .
        cp $homedir/inparam_hetero .

        if ( $multisrc == 'true' ) then
            cp $homedir/CMTSOLUTION .
        endif
    end
end


########################################################
######### submit the jobs ##############################
########################################################

set nodnum = `grep nproc_mesh $homedir/mesh_params.h |awk '{print $6}'`
echo "preparing job on $nodnum nodes..."

foreach isrc (${num_src_arr})
    foreach isim  (${srcapp})
        if ( $num_src == 1) then
            cd $isim
        else
            cd $isrc"_"$isim
        endif

        if ( $multisrc == 'true' ) then
            set outputname = "OUTPUT_"`echo $isim |sed 's/\//_/g'`
        else
            set outputname = "OUTPUT_"`echo $1 |sed 's/\//_/g'`
        endif

        if ( $newqueue == 'true' ) then

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

            endif

        else
            #ulimit -s unlimited
            #setenv OMP_NUM_THREADS 4
            #$mpiruncmd -n $nodnum ./axisem >& $outputname &

            ccc_msub -q standard ../sub_called_batch_for_AxiSEM_Curie_CD.sh
        endif

        echo "Job running in directory $isim"
        cd $mainrundir
    end
end

######## post processing ##################################################

cd $homedir
cd $1

#cp -p $homedir/UTILS/xpost_processing .
cp -p $homedir/UTILS/post_processing.F90 .
cp -p $homedir/UTILS/field_transform.F90 .
cp -p $homedir/UTILS/nc_postroutines.F90 .

echo "Compiling postprocessing routines"
make -f $homedir/UTILS/Makefile -sj
rm *.o *.mod

cp -p $homedir/UTILS/post_processing.csh .
cp -p $homedir/UTILS/plot_recfile_seis.csh .
cp -p $homedir/UTILS/plot_recs.plot .
cp -p $homedir/UTILS/taup_allrec.csh .
cp -p $homedir/UTILS/plot_record_section.m .

echo "To convolve and sum seismograms, run ./post_processing.csh after the simulations in:"
echo $mainrundir
echo ".... the post-processing input file param_post_processing is generated in the solver"
echo ".... based on guesses. Edit please."
echo " ~ ~ ~ ~ ~ ~ ~ h a n g   o n   &   l o o s e ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~"


