#!/bin/sh

#PBS -q tromp
#PBS -N XCOMPUTE_CG_DIRECTION
#PBS -l nodes=13:ppn=8
#PBS -l walltime=15:00:00
#PBS -j oe
#PBS -k oe
#PBS -o xcompute_direction_cg.log


echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

iter_new=M01
iter_old=M00

xoutput_tag=XTAG_$iter_new


direction_0=../DIRECTION_CG_$iter_old
direction_1=../DIRECTION_CG_$iter_new

gradient_0=../SUMMED_KERNEL_$iter_old
gradient_1=../SUMMED_KERNEL_$iter_new

if [ ! -d $direction_0 ]; then
  echo WRONG! NO $direction_0
  exit
fi
if [ ! -d $gradient_0 ]; then
  echo WRONG! NO $gradient_0
  exit
fi
if [ ! -d $gradient_1 ]; then
  echo WRONG! NO $gradient_1
  exit
fi

if [ ! -d $direction_1 ]; then
  echo MKDIR $direction_1
  mkdir $direction_1
fi

echo submit compute cg direction
mpiexec -np 100 ./xcompute_direction_cg $direction_0 $direction_1 $gradient_0 $gradient_1  > $xoutput_tag
echo done successfully


