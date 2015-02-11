#!/bin/sh
#PBS -q tromp
#PBS -N XCOMBINE_KERNEL_200704090832A
#PBS -l nodes=1:ppn=1
#PBS -l walltime=5:00:00
#PBS -j oe
#PBS -k oe
#PBS -o job_src2.log

echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

slicefile=XSLICE_FILE
topo_path=EUROPE_TOPOLOGY_FILE
local_path=../MODEL_INVERSION/ADJOINT_M18/CMTSOLUTION_200704090832A/KERNEL/
out_path=VTU_SINGLE_KERNEL_M18/CMTSOLUTION_200704090832A

#for tag in bulk_c_kernel bulk_betav_kernel bulk_betah_kernel eta_kernel rho_kernel hess_kernel
for tag in bulk_betav_kernel
do
  ./xcombine_vol_data $slicefile $tag $topo_path $local_path $out_path 0 1
done

echo combine kernels successfully

