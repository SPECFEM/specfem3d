#!/bin/sh
# generate and clean code dir for several simulations
# Input: cmtid: CMTSOLUTION id
#        iter: iteration number
#    opt: 1-generate code
#             2-clean code
#             3-clean /scratch/lustre/hejunzhu/*** directory
# Hejun Zhu, Mar 23, 2010

eventfile="../SHARE_FILES/EVENTID_CENTER/XEVENTID"
iter=M00

parent="XSRC_SEM"
cmtcenter="../SHARE_FILES/CMTSOLUTION_CENTER"
stacenter="../SHARE_FILES/STATION_CENTER"
xcopy="COPY_LOCAL/xcopy_local_forward"

pnm=FORWARD_$iter

if [ ! -f $eventfile ]; then
  echo WRONG! NO $eventfile
  exit
fi
if [ ! -d $pnm ]; then
  echo mkdir $pnm ...
  mkdir $pnm
fi
if [ ! -d $cmtcenter ]; then
  echo WRONG! NO $cmtcenter
  exit
fi
if [ ! -d $stacenter ]; then
  echo WRONG! NO $stacenter
  exit
fi
if [ ! -f $xcopy ]; then
  echo WRONG! NO $xcopy
  exit
fi



while read line
do
  echo "*** processing $line ***"

  dir=$pnm"/"$line
  eid=`echo $line | awk -F"_" '{print $NF}'`
  stafile=STATION_$eid

  if [ ! -d $dir ]; then
    echo "   make dir...   "
    mkdir $dir
  fi

  if [ ! -d $dir/DATA ]; then
    echo " make DATA"
    mkdir $dir/DATA
  fi

  # only copy useful part of the code to sub directory
  cp -r $parent/DATA/Par_file  $dir/DATA/
  cp -r $parent/DATA/STATIONS  $dir/DATA/
  cp -r $parent/DATA/CMTSOLUTION  $dir/DATA/
  cp -r $parent/OUTPUT_FILES $dir/
  cp $parent/xparallel_forward_solver.sh $dir/
  cp $parent/xspecfem3D $dir/
  cp $xcopy $dir/

        cat $cmtcenter/$line > $dir/DATA/CMTSOLUTION
  cat $stacenter/$stafile > $dir/DATA/STATIONS

  cd $dir/DATA

  if [ ! -f Par_file ]; then
    echo WRONG! NO Par_file
    exit
  fi
  echo "   change Par_file...   "
  path="\/scratch\/hejunzhu"
  sed -e "s/^LOCAL_PATH.*$/LOCAL_PATH               = $path/g" Par_file >Par_file_out
  mv Par_file_out Par_file

  cd ..
  if [ ! -f xparallel_forward_solver.sh ]; then
    echo WRONG! NO xparallel_forward_solver.sh
    exit
  fi
  echo "    changing pbs file   "
  tag="#PBS -N  $line"
  sed -e "s/^#PBS -N.*$/$tag/g" xparallel_forward_solver.sh > xparallel_forward_solver_out.sh
  mv xparallel_forward_solver_out.sh xparallel_forward_solver.sh

  if [ ! -f xspecfem3D ] ;then
    echo WRONG! NO xspecfem3D
    exit
  fi
  if [ ! -f xcopy_local_forward ]; then
    echo WRONG! NO xcopy_local_forward
    exit
  fi

  mimic=`grep USE_ATTENUATION_MIMIC OUTPUT_FILES/values_from_mesher.h | awk -F"=" '{print $NF}'`
  if [ $mimic != ".false." ]; then
    echo WRONG! NEED reompile solver with MIMIC attenuation
    exit
  fi


  echo submitting job
  qsub xparallel_forward_solver.sh

  sleep 3

  cd ../../

done < $eventfile
