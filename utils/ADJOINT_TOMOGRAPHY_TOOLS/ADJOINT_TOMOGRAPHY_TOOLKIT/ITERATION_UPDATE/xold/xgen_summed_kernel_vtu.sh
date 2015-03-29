#!/bin/sh
# Author: Hejun Zhu, hejunzhu@princeton.edu
# Princeton University, New Jersey, USA
# Last modified: Tue Feb 22 09:34:05 EST 2011



iter=M18

slicefile=XSLICE_FILE
topo_path=/tigress-hsm/hejunzhu/2011EUROPE_ITERATION_UPDATE/EUROPE_TOPOLOGY_FILE
local_path=/tigress-hsm/hejunzhu/2011EUROPE_ITERATION_UPDATE/SUMMED_KERNEL_$iter
out_path=/tigress-hsm/hejunzhu/2011EUROPE_ITERATION_UPDATE/VTU_SUMMED_KERNEL_$iter

if [ ! -f $slicefile ] ;then
  echo WRONG! NO $slicefile
  exit
fi
if [ ! -d $topo_path ]; then
  echo WRONG! NO $topo_path
  exit
fi
if [ ! -d $local_path ]; then
  echo WROGN! NO $local_path
  exit
fi
if [ ! -d $out_path ]; then
  echo MKDIR $out_path
  mkdir $out_path
fi


#for tag in bulk_c_kernel bulk_betav_kernel bulk_betah_kernel eta_kernel rho_kernel hess_kernel
#for tag in bulk_betav_kernel bulk_betah_kernel hess_kernel
#for tag in bulk_betah_kernel
#for tag in bulk_c_kernel_precond bulk_betav_kernel_precond bulk_betah_kernel_precond eta_kernel_precond rho_kernel_precond #hess_kernel_precond
#for tag in bulk_betav_kernel_precond bulk_betah_kernel_precond hess_kernel_precond
#for tag in bulk_betah_kernel_precond
#for tag in bulk_betav_kernel_precond_smooth bulk_betah_kernel_precond_smooth bulk_c_kernel_precond_smooth eta_kernel_precond_smooth
for tag in bulk_betav_kernel_precond_smooth bulk_betah_kernel_precond_smooth
#for tag in bulk_betav_kernel_precond_smooth
do
  ./xcombine_vol_data $slicefile $tag $topo_path $local_path $out_path 0 1
done

echo combine kernels successfully

