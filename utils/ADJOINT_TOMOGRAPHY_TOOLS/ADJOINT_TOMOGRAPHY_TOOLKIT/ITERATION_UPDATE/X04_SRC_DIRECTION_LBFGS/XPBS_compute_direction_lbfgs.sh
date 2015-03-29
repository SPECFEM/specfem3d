#!/bin/sh

#PBS -q tromp
#PBS -N XCOMPUTE_LBFGS_DIRECTION
#PBS -l nodes=13:ppn=8
#PBS -l walltime=15:00:00
#PBS -j oe
#PBS -k oe
#PBS -o job_src2.log


echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

iter_start=35
iter_current=36

xoutput_tag=XTAG_$iter_current


# check directories
model_0=../MODEL_M$iter_start
model_1=../MODEL_M$iter_current

gradient_0=../SUMMED_KERNEL_M$iter_start
gradient_1=../SUMMED_KERNEL_M$iter_current

direction=../DIRECTION_LBFGS_M$iter_current

if [ ! -d $model_0 ]; then
  echo WRONG! NO $model_0
  exit
fi
if [ ! -d $model_1 ]; then
  echo WRONG! NO $model_1
  exit
fi
if [ ! -d $gradient_1 ]; then
  echo WRONG! NO $gradient_1
  exit
fi
if [ ! -d $gradient_0 ]; then
  echo WRONG! NO $gradient_0
  exit
fi

if [ ! -d $direction ]; then
  echo MKDIR $direction
  mkdir $direction
fi


# submit job
echo submit compute direction lbfgs
mpiexec -np 100 ./xcompute_direction_lbfgs $iter_start $iter_current  > $xoutput_tag
echo done successfully


