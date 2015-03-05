#!/bin/sh

#PBS -q tromp
#PBS -N XUPDATE_MODEL
#PBS -l nodes=13:ppn=8
#PBS -l walltime=15:00:00
#PBS -j oe
#PBS -k oe
#PBS -o job_src2.log


echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

iter_old=M17
iter_new=M18
step_len=0.01

xoutput_tag=XTAG_$iter_new


input_model=../MODEL_$iter_old
input_kernel=../DIRECTION_CG_$iter_old
output_model=../MODEL_$iter_new

if [ ! -d $input_model ]; then
  echo WRONG! NO $input_model
  exit
fi
if [ ! -d $input_kernel ]; then
  echo WRONG! NO $input_kernel
  exit
fi
if [ ! -d $output_model ]; then
  echo MKDIR $output_model
  mkdir $output_model
fi

echo submit updata model
mpiexec -np 100 ./xadd_model_globe $step_len $input_model $input_kernel $output_model > $xoutput_tag
echo done successfully


