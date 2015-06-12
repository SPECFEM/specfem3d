#!/bin/csh -f

source /etc/profile.d/modules.csh
        module purge
#        module load CINES/20110112
        module load intel/12.1.3
        module load intelmpi/4.0.3
#	module load openmpi/1.2.7-11

which mpif90

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

	echo "#PBS -S /bin/bash" > submit_axisemmesh.sh
	echo "#PBS -N AxSEMesh" >> submit_axisemmesh.sh
	echo "#PBS -e AxSEMesh.err" >> submit_axisemmesh.sh
	echo "#PBS -o OUTPUT">> submit_axisemmesh.sh
	echo "#PBS -l select=1:ncpus=8:mpiprocs=1">> submit_axisemmesh.sh
	echo "#PBS -l walltime=1:59:00">> submit_axisemmesh.sh
	echo 'cat $PBS_NODEFILE'>> submit_axisemmesh.sh
	echo 'cd $PBS_O_WORKDIR'>> submit_axisemmesh.sh
	echo ". /etc/profile.d/modules.sh">> submit_axisemmesh.sh
	echo "module purge">> submit_axisemmesh.sh
	echo "module load intel/12.1.3">> submit_axisemmesh.sh
	echo "module load intelmpi/4.0.3">> submit_axisemmesh.sh
#	echo "module load CINES/20110112">> submit_axisemmesh.sh
#	echo "module load openmpi/1.4.4" >> submit_axisemmesh.sh
	echo 'cat $PBS_NODEFILE | uniq > mpd.hosts'>> submit_axisemmesh.sh
	echo "mpdboot --rsh=ssh -v -n `cat mpd.hosts|wc -l` -f mpd.hosts">> submit_axisemmesh.sh
	echo 'nprocs=$(cat $PBS_NODEFILE|wc -l)'>> submit_axisemmesh.sh
	echo "ulimit -s unlimited">> submit_axisemmesh.sh
#	echo 'TMPDIR="./"' >> submit_axisemmesh.sh
	echo 'nohup $PBS_O_WORKDIR/xmesh '>> submit_axisemmesh.sh
	echo "mpdallexit">> submit_axisemmesh.sh

	chmod +x submit_axisemmesh.sh
	qsub ./submit_axisemmesh.sh
#!    nohup ./xmesh > OUTPUT &
    # uncomment the following three lines to monitor memory usage of the mesher
    #cd UTILS
    #python monitor_memory.py > ../memory_output &
    #cd ..

endif

echo 'xmesh submitted, output in "OUTPUT"'
echo "After the run, move the mesh to a new directory <meshdir> via:"
echo "./movemesh <meshdir>"
echo "This will be located in ../SOLVER/MESHES/<meshdir>"

