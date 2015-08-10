#!/bin/bash
#SBATCH -J monnom
#SBATCH -D /scratch/cnt0023/git6091/beller/EGU_2015/Real_Data/AxiSEM_simulations/axisem-master/SOLVER/mondossier
#SBATCH --get-user-env
#SBATCH --nodes=20
#SBATCH --ntasks=480
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1 
#SBATCH --time=00:30:00
#SBATCH --output INTERP
module purge
module load intel/15.0.0.090
module load bullxmpi/1.2.8.3
ulimit -s unlimited
time srun --mpi=pmi2 -K1 --resv-ports -n $SLURM_NTASKS ./interpolate_3D_wavefield.x > OUTPUT_interp

