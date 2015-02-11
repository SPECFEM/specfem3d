#!/bin/sh

#PBS -q tromp
#PBS -N XSLICE
#PBS -l nodes=12:ppn=8+1:ppn=4
#PBS -l walltime=15:00:00
#PBS -j oe
#PBS -k oe
#PBS -o job_src2.log


depthslice=XYZ_FILE/VERT_SLICE_England.xyz
topo_path=/scratch/lustre/hejunzhu/2011EUROPE_SHOW_MODEL/EUROPE_TOPOLOGY_FILE_NO_TOPO_ELL
model_path=/scratch/lustre/hejunzhu/2011EUROPE_SHOW_MODEL/MODEL_M30_PERT_STW
tag=dvsv
gmtout=/scratch/lustre/hejunzhu/2011EUROPE_SHOW_MODEL/XSRC_MODEL_SLICE_VERT/MODEL_M30_PERT_STW/VERT_SLICE_England_dvsv.xyz


echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

echo submit slicing
mpiexec -np 100 ./xsem_model_slice $depthslice $topo_path $model_path $tag $gmtout
echo done successfully


