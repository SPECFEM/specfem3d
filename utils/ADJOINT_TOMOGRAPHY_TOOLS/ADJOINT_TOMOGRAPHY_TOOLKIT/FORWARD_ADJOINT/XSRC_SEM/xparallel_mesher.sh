#!/bin/sh

#PBS -q tromp
#PBS -N XEUMESHER
#PBS -l nodes=13:ppn=8
#PBS -l walltime=2:00:00
#PBS -j oe
#PBS -k oe
#PBS -o job_src2.log

#module load openmpi

echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

local_dir=`grep LOCAL_PATH DATA/Par_file | cut -d '=' -f 2`

if [ ! -d $local_dir ]; then
  mkdir $local_dir
  echo create $local_dir
fi

mpiexec -np 100 ./xmeshfem3D
echo mesher done successfully


