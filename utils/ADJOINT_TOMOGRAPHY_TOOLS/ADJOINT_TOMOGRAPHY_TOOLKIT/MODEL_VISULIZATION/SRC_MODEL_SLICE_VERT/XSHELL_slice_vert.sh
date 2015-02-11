#!/bin/sh
# Author: Hejun Zhu, hejunzhu@princeton.edu
# Princeton University, New Jersey, USA
# Last modified: Tue Jan 25 20:58:19 EST 2011

iter=M30_PERT_STW

topo_path="\/scratch\/lustre\/hejunzhu\/2011EUROPE_SHOW_MODEL\/EUROPE_TOPOLOGY_FILE_NO_TOPO_ELL"
topo_path1=/scratch/lustre/hejunzhu/2011EUROPE_SHOW_MODEL/EUROPE_TOPOLOGY_FILE_NO_TOPO_ELL

model_path="\/scratch\/lustre\/hejunzhu\/2011EUROPE_SHOW_MODEL\/MODEL_$iter"
model_path1=/scratch/lustre/hejunzhu/2011EUROPE_SHOW_MODEL/MODEL_$iter

#model_path="\/scratch\/lustre\/hejunzhu\/2011EUROPE_ITERATION_UPDATE\/MODEL_$iter"
#model_path1=/scratch/lustre/hejunzhu/2011EUROPE_ITERATION_UPDATE/MODEL_$iter

out_path="\/scratch\/lustre\/hejunzhu\/2011EUROPE_SHOW_MODEL\/XSRC_MODEL_SLICE_VERT\/MODEL_$iter"
out_path1=/scratch/lustre/hejunzhu/2011EUROPE_SHOW_MODEL/XSRC_MODEL_SLICE_VERT/MODEL_$iter



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


#for sliceid in A B C D E F G H
#for sliceid in SINGLE
for sliceid in England
do
  echo submitting slice $sliceid

  depthslice="XYZ_FILE\/VERT_SLICE_$sliceid.xyz"
  depthslice1="XYZ_FILE/VERT_SLICE_$sliceid.xyz"


  if [ ! -f $depthslice1 ]; then
    echo WRONG! NO $depthslice1
    exit
  fi

  #for tag in dvsh_vsv
  #for tag in dvsv dvsh dvsh_vsv
  #for tag in dvsh_vsv dvs dvp  #dvsh dbulk # dvp
  for tag in dvsv
  do
    gmtout="\/scratch\/lustre\/hejunzhu\/2011EUROPE_SHOW_MODEL\/XSRC_MODEL_SLICE_VERT\/MODEL_$iter\/VERT_SLICE_"$sliceid"_"$tag".xyz"


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

