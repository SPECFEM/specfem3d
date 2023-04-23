#!/bin/bash
#original file name: script_to_copy_all_the_initial_input_files_correctly_for_the_NUMBER_OF_SIMULTANEOUS_RUNS_option.bash
#BATCH -N 1
#SBATCH -t 5:09:00
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4

#cd $SBATCH_O_WORKDIR
cd /scratch/gpfs/etienneb/specfem3d-devel/
rm OUTPUT_FILES/*

echo `pwd`
module purge
module load cudatoolkit/8.0
module load openmpi/intel-16.0/1.10.2/64
module load intel/16.0/64/16.0.4.258

module list
NPROC=4
rm output*
rm OUTPUT_FILES/*
rm OUTPUT_FILES/DATABASES_MPI/*
rm run*/OUTPUT_FILES/*
rm run*/OUTPUT_FILES/DATABASES_MPI/*
rm -rf run0001
cp p DATA/Par_file

# creates and decomposes mesh
echo
echo "running mesher..."
echo
mpirun -n 1 ./bin/xmeshfem3D >> outputmesher
mpirun -n 1 ./bin/xgenerate_databases >> outputdatabases

cp -R m_run0001 run0001
rm DATA/Par_file
cp OUTPUT_FILES/values_from_mesher.h run0001/OUTPUT_FILES
cp OUTPUT_FILES/values_from_mesher.h run0002/OUTPUT_FILES
cp OUTPUT_FILES/values_from_mesher.h run0003/OUTPUT_FILES
cp OUTPUT_FILES/values_from_mesher.h run0004/OUTPUT_FILES

cp OUTPUT_FILES/DATABASES_MPI/* run0001/OUTPUT_FILES/DATABASES_MPI
cp OUTPUT_FILES/DATABASES_MPI/* run0002/OUTPUT_FILES/DATABASES_MPI
cp OUTPUT_FILES/DATABASES_MPI/* run0003/OUTPUT_FILES/DATABASES_MPI
cp OUTPUT_FILES/DATABASES_MPI/* run0004/OUTPUT_FILES/DATABASES_MPI




# runs simulation
echo
echo "  running inverse problem..."
echo
mpirun -n $NPROC ./bin/xinverse_problem_for_model l-bfgs >> outputinverseproblem
#mpirun -n $NPROC ./bin/xinverse_problem_for_model forward >> outputinverseproblem

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
date

