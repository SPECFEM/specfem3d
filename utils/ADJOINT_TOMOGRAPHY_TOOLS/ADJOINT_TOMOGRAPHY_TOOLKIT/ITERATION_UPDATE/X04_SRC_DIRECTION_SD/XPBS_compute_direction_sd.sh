#!/bin/sh

#PBS -q tromp
#PBS -N XCOMPUTE_CG_DIRECTION
#PBS -l nodes=13:ppn=8
#PBS -l walltime=15:00:00
#PBS -j oe
#PBS -k oe
#PBS -o xcompute_direction_sd.log


echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

iter=M00

xoutput_tag=XTAG_$iter

direction_dir=../DIRECTION_SD_$iter
gradient_dir=../SUMMED_KERNEL_$iter


if [ ! -d $gradient_dir ]; then
  echo WRONG! NO $gradient_dir
  exit
fi

if [ ! -d $direction_dir ]; then
  echo MKDIR $direction_dir
  mkdir $direction_dir
fi

echo submit compute sd direction
mpiexec -np 100 ./xcompute_direction_sd $direction $gradient  > $xoutput_tag
echo done successfully


