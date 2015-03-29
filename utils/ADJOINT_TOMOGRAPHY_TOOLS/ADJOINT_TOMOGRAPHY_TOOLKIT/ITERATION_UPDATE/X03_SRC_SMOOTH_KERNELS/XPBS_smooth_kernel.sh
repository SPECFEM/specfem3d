#!/bin/sh

#PBS -q tromp
#PBS -N XSMOOTH_eta_kernel_precond
#PBS -l nodes=13:ppn=8
#PBS -l walltime=15:00:00
#PBS -j oe
#PBS -k oe
#PBS -o job_src2.log


echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

iter=M18
sigma_h=50
sigma_v=5
topo_path=../EUROPE_TOPOLOGY_FILE
kernel_path=../SUMMED_KERNEL_$iter
tag=eta_kernel_precond

output_tag="XTAG_SMOOTH_"$iter"_"$tag

echo submit smoothing $tag kernel
mpiexec -np 100 ./xsmooth_sem_globe $sigma_h $sigma_v $tag $kernel_path $topo_path  > $output_tag
echo smoothing $tag kernel done successfully
