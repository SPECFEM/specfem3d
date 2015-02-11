#!/bin/sh

#PBS -q tromp
#PBS -N XPRECOND_KERNELS
#PBS -l nodes=13:ppn=8
#PBS -l walltime=15:00:00
#PBS -j oe
#PBS -k oe
#PBS -o precond_kernels.log


echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

iter=M00

input_dir=../SUMMED_KERNELS_$iter

if [ ! -d $input_dir ]; then
  echo WRONG! NO $input_dir
  exit
fi


echo submit precondition kernels
mpiexec -np 100 ./xprecond_kernels $input_dir
echo done successfully


