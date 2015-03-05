#!/bin/sh

#PBS -q tromp
#PBS -N XEUSOLVER
#PBS -l nodes=13:ppn=8
#PBS -l walltime=15:00:00
#PBS -j oe
#PBS -k oe
#PBS -o job_src2.log


echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

local_dir=`grep LOCAL_PATH DATA/Par_file | cut -d '=' -f 2`

# obtain lsf job information
cat $PBS_NODEFILE > OUTPUT_FILES/compute_nodes
echo "$PBS_JOBID" > OUTPUT_FILES/jobid


if [ -d /scratch/hejunzhu ]; then
  echo rm /scratch/hejunzhu
  pbsdsh -u rm -rf /scratch/hejunzhu
fi
echo mkdir local dir...
pbsdsh -u mkdir -p /scratch/hejunzhu

echo copying mesh file...
mpiexec -np 100 ./xcopy_local_forward

echo submit job...
mpiexec -np 100 ./xspecfem3D
echo solver done successfully

#echo clean the folder
#mv OUTPUT_FILES/*.sac SYN/tmp/


echo cleaning nodes...
pbsdsh -u rm -rf /scratch/hejunzhu

echo done sucessfully
