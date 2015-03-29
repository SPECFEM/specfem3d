#!/bin/sh
# Author: Hejun Zhu, hejunzhu@princeton.edu
# Princeton University, New Jersey, USA
# Last modified: Fri Sep 14 15:08:23 EDT 2012


# This script is used to generate vtu file from bin file

iter=M00
vtu=KERNEL
if [ $vtu == DIRECTION ]; then
  method=CG
fi


topo_path=EUROPE_TOPOLOGY_FILE
slicefile=XSLIEC_FILE



if [ $vtu == KERNEL ]; then
  local_path=SUMMED_KERNEL_$iter
  out_path=VTU_SUMMED_KERNEL_$iter
elif [ $vtu == DIRECTION ]; then
  if [ $method == SD ]; then
    local_path=DIRECTION_SD_$iter
    output_path=VTU_DIRECTION_SD_$iter
  elif [ $method == CG ]; then
    local_path=DIRECTION_CG_$iter
    output_path=VTU_DIRECTION_CG_$iter
  elif [ $method == LBFGS ]; then
    local_path=DIRECTION_LBFGS_$iter
    output_path=VTU_DIRECTION_LBFGS_$iter
  else
    echo WRONG! NO $method
  fi
elif [ $vtu == MODEL ]; then
  local_path=MODEL_$iter
  output_path=VTU_MODEL_$iter
elif [ $vtu == MODEL_PERT ]; then
  local_path=MODEL_$iter"_PERT_STW"
  output_path=VTU_MODEL_$iter"_PERT_STW"
else
  echo WRONG! NO $vtu
  exit
fi

if [ $vtu == KERNEL ]; then
  for tag in bulk_betav_kernel_precond_smooth bulk_betah_kernel_precond_smooth bulk_c_kernel_precond_smooth eta_kernel_precond_smooth
  do
    ./xcombine_vol_data $slicefile $tag $topo_path $local_path $out_path 0 1
  done
elif [ $vtu == DIRECTION ]; then
  for tag in bulk_c_kernel_precond bulk_betav_kernel_precond bulk_betah_kernel_precond eta_kernel_precond rho_kernel_precond hess_kernel_precond
  do
    ./xcombine_vol_data $slicefile $tag $topo_path $local_path $out_path 0 1
  done
elif [ $vtu == MODEL ]; then
  for tag in vpv vph vsv vsh eta rho
  do

    ./xcombine_vol_data $slicefile $tag $topo_path $local_path $out_path 0 1
  done
elif [ $vtu == MODEL_PERT ]; then
  for tag in dvsh dvpv dvph deta drho
  do
    ./xcombine_vol_data $slicefile $tag $topo_path $local_path $out_path 0 1
  done

else
  echo WRONG! NO $vtu
  exit
fi




