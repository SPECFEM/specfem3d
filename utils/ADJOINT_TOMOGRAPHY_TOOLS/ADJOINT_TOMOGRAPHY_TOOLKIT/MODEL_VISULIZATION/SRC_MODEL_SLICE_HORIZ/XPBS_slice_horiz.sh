#!/bin/sh


#PBS -q tromp
#PBS -N XSLICE
#PBS -l nodes=12:ppn=8+1:ppn=4
#PBS -l walltime=15:00:00
#PBS -j oe
#PBS -k oe
#PBS -o job_src2.log


depthslice=XYZ_FILE_HORIZ/DEPTH_SLICE_075.xyz
topo_path=/scratch/lustre/hejunzhu/2012SHEAR_ATTENUATION_ITERATION_UPDATE/EUROPE_TOPOLOGY_FILE
model_path=/scratch/lustre/hejunzhu/2012SHEAR_ATTENUATION_ITERATION_UPDATE/MODEL_M42_PERT_STW
tag=dQmu
gmtout=/scratch/lustre/hejunzhu/2012SHEAR_ATTENUATION_ITERATION_UPDATE/XSRC_MODEL_SLICE_HORIZ/MODEL_M42_PERT_STW/DEPTH_SLICE_075_dQmu.xyz


echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

echo submit slicing
mpiexec -np 100 ./xsem_model_slice $depthslice $topo_path $model_path $tag $gmtout
echo done successfully


