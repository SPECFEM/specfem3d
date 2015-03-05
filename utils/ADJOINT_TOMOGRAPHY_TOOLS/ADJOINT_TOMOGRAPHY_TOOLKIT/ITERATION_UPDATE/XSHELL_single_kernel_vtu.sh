#!/bin/sh
# Author: Hejun Zhu, hejunzhu@princeton.edu
# Princeton University, New Jersey, USA
# Last modified: Sat Jan 22 18:47:49 EST 2011


eventfile=../SHARE_FILES/EVENTID_CENTER/XEVENTID
iter=M18
slicefile=XSLICE_FILE

savedir=VTU_SINGLE_KERNEL_$iter
topo_path1=EUROPE_TOPOLOGY_FILE
topo_path2="EUROPE_TOPOLOGY_FILE"


if [ ! -f $eventfile ]; then
  echo WRONG! NO $eventfile
  exit
fi
if [ ! -f $slicefile ]; then
  echo WRONG! NO $slicefile
  exit
fi
if [ ! -f xcombine_vol_data ]; then
  echo WRONG! NO xcombine_vol_data
  exit
fi
if [ ! -d $savedir ]; then
  echo MKDIR $savedir
  mkdir $savedir
fi
if [ ! -d $topo_path1 ]; then
  echo WRONG !NO $topo_path 1
  exit
fi


while read line
do
  cmtid=`echo $line | awk -F"_" '{print $2}'`
  local_path1=../FORWARD_ADJOINT/ADJOINT_$iter/$line/KERNEL/
  local_path2="..\/FORWARD_ADJOINT\/ADJOINT_$iter\/$line\/KERNEL\/"

  out_path1=VTU_SINGLE_KERNEL_$iter/$line
  out_path2="VTU_SINGLE_KERNEL_$iter\/$line"

  proctag="#PBS -N XCOMBINE_KERNEL_$cmtid"

  if [ ! -d $out_path1 ] ;then
    echo MKDIR $out_path1
    mkdir $out_path1
  fi

  if [ ! -d $local_path1 ]; then
    echo WRONG! NO $local_path1
    exit
  fi


  sed -e "s/^#PBS -N.*$/$proctag/g" \
      -e "s/^slicefile=.*$/slicefile=$slicefile/g" \
      -e "s/^topo_path=.*$/topo_path=$topo_path2/g" \
      -e "s/^local_path=.*$/local_path=$local_path2/g" \
      -e "s/^out_path=.*$/out_path=$out_path2/g" \
      XPBS_single_kernel_vtu.sh > XPBS_single_kernel_vtu.sh.out
      mv XPBS_single_kernel_vtu.sh.out XPBS_single_kernel_vtu.sh

  echo qsub $line
  qsub XPBS_single_kernel_vtu.sh
  sleep 3

done < $eventfile


