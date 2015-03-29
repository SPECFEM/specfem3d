#!/bin/sh

#PBS -q tromp
#PBS -N XSUM_KERNELS
#PBS -l nodes=13:ppn=8
#PBS -l walltime=15:00:00
#PBS -j oe
#PBS -k oe
#PBS -o sum_kernels.log


echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

iter=M00

input_dir=../../FORWARD_ADJOINT/ADJOINT_$iter
output_dir=../SUMMED_KERNELS_$iter
eventid=../../SHARE_FILES/EVENTID_CENTER/XEVENTID


# checking directories
if [ ! -d $input_dir ]; then
  echo WRONG! NO $inputdir
  exit
fi

if [ ! -f $eventid ]; then
  echo WRONG! NO $eventid
  exit
fi

# making directories
if [ ! -d $output_dir ]; then
  echo MKDIR $output_dir
  mkdir $output_dir
fi


echo submit summing kernels
mpiexec -np 100 ./xsum_kernels $input_dir $output_dir $eventid
echo done successfully


