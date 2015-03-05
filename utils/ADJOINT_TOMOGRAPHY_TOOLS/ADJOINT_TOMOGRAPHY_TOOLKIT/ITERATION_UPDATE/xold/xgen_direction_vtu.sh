#!/bin/sh
# Author: Hejun Zhu, hejunzhu@princeton.edu
# Princeton University, New Jersey, USA
# Last modified: Tue Feb 22 09:34:05 EST 2011



iter=M00
method=CG

slicefile=XSLICE_FILE
topo_path=EUROPE_TOPOLOGY_FILE

if [ $method == SD ]; then
  local_path=DIRECTION_SD_$iter
  out_path=VTU_DIRECTION_SD_$iter
elif [ $method == CG ]; then
  local_path=DIRECTION_CG_$iter
  out_path=VTU_DIRECTION_CG_$iter
elif [ $method == LBFGS ]; then
  local_path=DIRECTION_LBFGS_$iter
  out_path=VTU_DIRECTION_LBFGS_$iter
else
  echo WRONG! NO $method
  exit
fi



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


for tag in bulk_c_kernel_precond bulk_betav_kernel_precond bulk_betah_kernel_precond eta_kernel_precond rho_kernel_precond hess_kernel_precond
do
  ./xcombine_vol_data $slicefile $tag $topo_path $local_path $out_path 0 1
done

echo combine kernels successfully

