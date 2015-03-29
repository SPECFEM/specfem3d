#!/bin/sh
# Author: Hejun Zhu, hejunzhu@princeton.edu
# Princeton University, New Jersey, USA
# Last modified: Tue Jan 25 20:58:19 EST 2011

iter=M42_PERT_STW

topo_path="\/scratch\/lustre\/hejunzhu\/2012SHEAR_ATTENUATION_ITERATION_UPDATE\/EUROPE_TOPOLOGY_FILE"
topo_path1=/scratch/lustre/hejunzhu/2012SHEAR_ATTENUATION_ITERATION_UPDATE/EUROPE_TOPOLOGY_FILE

model_path="\/scratch\/lustre\/hejunzhu\/2012SHEAR_ATTENUATION_ITERATION_UPDATE\/MODEL_$iter"
model_path1=/scratch/lustre/hejunzhu/2012SHEAR_ATTENUATION_ITERATION_UPDATE/MODEL_$iter

out_path="\/scratch\/lustre\/hejunzhu\/2012SHEAR_ATTENUATION_ITERATION_UPDATE\/XSRC_MODEL_SLICE_HORIZ\/MODEL_$iter"
out_path1=/scratch/lustre/hejunzhu/2012SHEAR_ATTENUATION_ITERATION_UPDATE/XSRC_MODEL_SLICE_HORIZ/MODEL_$iter


if [ ! -f xsem_model_slice ]; then
  echo WRONG! NO xsem_model_slice
  exit
fi

if [ ! -d $topo_path1 ]; then
  echo WRONG !NO $topo_path1
  exit
fi

if [ ! -d $model_path1 ]; then
  echo WRONG! NO $model_path1
  exit
fi


if [ ! -d $out_path1 ] ;then
  echo MKDIR $out_path1
  mkdir $out_path1
fi

#for depth in 005 025 050 075 125 175 225 275 325 375 425 475 525 575 625 675
for depth in 075
do
  depthslice="XYZ_FILE_HORIZ\/DEPTH_SLICE_$depth.xyz"
  depthslice1="XYZ_FILE_HORIZ/DEPTH_SLICE_$depth.xyz"
  if [ ! -f $depthslice1 ]; then
    echo WRONG! NO $depthslice1
    exit
  fi

  #for tag in dvsh_vsv dbulk dpossion
  #for tag in dvsh dvsh_vsv dbulk dpossion
  #for tag in dvp_vs
  for tag in dvsv dQmu
  do
    gmtout="\/scratch\/lustre\/hejunzhu\/2012SHEAR_ATTENUATION_ITERATION_UPDATE\/XSRC_MODEL_SLICE_HORIZ\/MODEL_$iter\/DEPTH_SLICE_"$depth"_"$tag".xyz"

    echo $depth $tag

    sed -e "s/^depthslice=.*$/depthslice=$depthslice/g"\
          -e "s/^topo_path=.*$/topo_path=$topo_path/g"\
          -e "s/^model_path=.*$/model_path=$model_path/g"\
          -e "s/^tag=.*$/tag=$tag/g"\
          -e "s/^gmtout=.*$/gmtout=$gmtout/g"\
    xgen_slice_pbs.sh > xgen_slice_pbs.sh.out
    mv xgen_slice_pbs.sh.out xgen_slice_pbs.sh
    qsub xgen_slice_pbs.sh
    sleep 5

  done
done


