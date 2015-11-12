#!/bin/bash
#SBATCH -J MesherAxisem
#SBATCH -D /scratch/cnt0023/git6091/beller/EGU_2015/Real_Data/AxiSEM_simulations/axisem-master/MESHER
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --time=00:30:00
#SBATCH --output AxiMesh
module purge
module load intel/15.0.0.090
module load bullxmpi/1.2.8.3
export OMP_NUM_THREADS=24
ulimit -s unlimited
time srun --mpi=pmi2 -K1 --resv-ports -n $SLURM_NTASKS ./xmesh > OUTPUT

