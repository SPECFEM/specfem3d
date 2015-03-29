#!/bin/sh
#PBS -q tromp
#PBS -N XXCORRECT_201105192015A
#PBS -l nodes=1:ppn=1
#PBS -l walltime=15:00:00
#PBS -j oe
#PBS -k oe

echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR


input=XCORRECT_INPUT_M42/CORRECT_201105192015A

./xcorrect_syn_time_moment < $input
