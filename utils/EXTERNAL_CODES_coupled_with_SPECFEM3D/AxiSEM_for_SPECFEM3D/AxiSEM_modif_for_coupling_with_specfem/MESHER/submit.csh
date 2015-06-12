#!/bin/csh -f

if ( $1 == '-h' ) then 
  echo "Argument options:"
  echo "   default (no argument): submit xmesh on current machine"
  echo "   lsf: submit to a lsf queue using bsub"
  echo "   torque: submit to a Torque/Mauo queue using qsub"
  exit
endif

# tidy up:
rm -rf OUTPUT
rm -rf meshdb.dat*
rm -rf mesh_params.h*
# create emtpy output file
touch OUTPUT

if ( ! -d Diags ) then
  mkdir Diags
else
  rm -rf Diags/*
endif

if ( { make -j 5 all } == 0 ) then
  echo "Compilation failed, please check the errors."
  exit
endif

set bgmodel = `grep "^BACKGROUND_MODEL" inparam_mesh | awk '{print $2}'`
if ( $bgmodel == 'external') then
  set fnam_extmodel = `grep "^EXT_MODEL" inparam_mesh | awk '{print $2}'`
  echo "Using external mesh file " $fnam_extmodel
  if ( ! -f $fnam_extmodel ) then
    echo "External mesh " $fnam_extmodel " does not exist!"
    exit
  endif
  cp $fnam_extmodel external_model.bm
endif

if ( $1 == 'lsf' ) then 
  ########## LSF SCHEDULER ######################
  bsub -R "rusage[mem=2048]" -I -n 1 ./xmesh > OUTPUT &

else if ( $1 == 'torque' ) then 
    ######## TORQUE/MAUI SCHEDULER #######
    echo "# Sample PBS for serial jobs" > run_mesh.pbs
    echo "#PBS -l nodes=1,walltime=2:00:00" >> run_mesh.pbs
    echo "ulimit -s unlimited " >> run_mesh.pbs
    echo "cd $PWD " >> run_mesh.pbs
    echo "./xmesh > OUTPUT " >> run_mesh.pbs
    qsub run_mesh.pbs

else if ( $1 == 'slurm' ) then 
    echo '#\!/bin/bash -l' > sbatch.sh
    echo "#SBATCH --ntasks=1" >> sbatch.sh
    echo "#SBATCH --nodes=1" >> sbatch.sh
    echo "#SBATCH --time=00:59:00" >> sbatch.sh
    echo "export OMP_NUM_THREADS=4" >> sbatch.sh
    echo "aprun -n 1 -d 4 ./xmesh > OUTPUT" >> sbatch.sh

    sbatch sbatch.sh 

else if ( $1 == 'slurmlocal' ) then 
    ######## slurm #######
    aprun -n 1 ./xmesh > OUTPUT &
else
    ######## SUBMIT LOCALLY #######
    #setenv OMP_NUM_THREADS 4
    nohup ./xmesh > OUTPUT &
    # uncomment the following three lines to monitor memory usage of the mesher
    #cd UTILS
    #python monitor_memory.py > ../memory_output &
    #cd ..
endif

echo 'xmesh submitted, output in "OUTPUT"'
echo "After the run, move the mesh to a new directory <meshdir> via:"
echo "./movemesh <meshdir>"
echo "This will be located in ../SOLVER/MESHES/<meshdir>"

