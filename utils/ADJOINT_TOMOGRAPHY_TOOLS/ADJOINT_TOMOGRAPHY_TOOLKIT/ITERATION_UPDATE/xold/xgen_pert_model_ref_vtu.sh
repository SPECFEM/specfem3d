#!/bin/sh
# Author: Hejun Zhu, hejunzhu@princeton.edu
# Princeton University, New Jersey, USA
# Last modified: Thu Feb 17 21:56:34 EST 2011


iter=M18_PERT_STW
slicefile=XSLICE_FILE

topo_path=/tigress-hsm/hejunzhu/2011EUROPE_ITERATION_UPDATE/EUROPE_TOPOLOGY_FILE
local_path=/tigress-hsm/hejunzhu/2011EUROPE_ITERATION_UPDATE/MODEL_$iter
out_path=/tigress-hsm/hejunzhu/2011EUROPE_ITERATION_UPDATE/VTU_MODEL_$iter


if [ ! -f $slicefile ]; then
  echo WRONG! NO $slicefile
  exit
fi
if [ ! -d $topo_path ]; then
  echo WRONG! NO $topo_path
  exit
fi
if [ ! -d $local_path ]; then
  echo WRONG! NO $local_path
  exit
fi
if [ ! -d $out_path ]; then
  echo MKDIR $out_path
  mkdir $out_path
fi



#for tag in dvsv dvsh dvpv dvph deta drho
for tag in dvsh dvpv dvph deta drho
do
  ./xcombine_vol_data $slicefile $tag $topo_path $local_path $out_path 0 1
done

echo combine kernels successfully

