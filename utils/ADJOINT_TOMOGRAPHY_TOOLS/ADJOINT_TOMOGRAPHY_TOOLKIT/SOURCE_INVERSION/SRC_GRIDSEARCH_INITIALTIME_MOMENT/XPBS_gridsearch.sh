#!/bin/sh
#PBS -q tromp
#PBS -N XOUTPUT_GRID_201105192015A
#PBS -l nodes=1:ppn=1
#PBS -l walltime=30:00:00
#PBS -j oe
#PBS -k oe

echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR


input=XGRIDSEARCH_INPUT_M42/PAR_201105192015A

./xgridsearch_time_moment < $input
