#!/bin/sh

# This script is used to submit adjoint simulation
# require adjoint source and station file first
# Hejun Zhu, Aug 25, 2010


eventfile=../SHARE_FILES/EVENTID_CENTER/XEVENTID
iter=M00

#fmin=15
#fmax=50
#ext1=`printf "%03i\n" $fmin`
#ext2=`printf "%03i\n" $fmax`
#ext="T"$ext1"_"$ext2
ext="comb"


parent="XSRC_SEM"
cmtcenter="../SHARE_FILES/CMTSOLUTION_CENTER"
stacenter="../SHARE_FILES/STATION_CENTER"
xcopy="COPY_LOCAL/xcopy_local_forward"
xchange="PERL_SRC/change_simulation_type.pl"

pnm=ADJOINT_$iter

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
if [ ! -f $xcopy ]; then
  echo WRONG! NO $xcopy
  exit
fi
if [ ! -f $xchange ]; then
  echo WRONG! NO $xchange
  exit
fi

while read line
do
  dir=$pnm"/"$line
  eid=`echo $line | awk -F"_" '{print $NF}'`
  stafile=STATION_$eid

  if [ ! -d $dir ]; then
    echo "   make dir...   "
    mkdir $dir
  fi
  if [ ! -d $dir/DATA ]; then
    echo " make DATA ... "
    mkdir $dir/DATA
  fi


  adjdir="SYN_"$iter"/"$line"_ADJ_"$ext
  if [ ! -d $adjdir ]; then
    echo WRONG! NO $adjdir
    exit
  fi


  # only copy useful part of the code to sub directory
  cp -r $parent/DATA/Par_file  $dir/DATA
  cp -r $parent/DATA/STATIONS  $dir/DATA
  cp -r $parent/DATA/CMTSOLUTION  $dir/DATA
  cp -r $parent/OUTPUT_FILES $dir/
  cp -r $parent/SEM $dir/
  cp -r $parent/KERNEL $dir/
  cp $parent/xparallel_adjoint_solver.sh $dir/
  cp $parent/xspecfem3D $dir/
  cp $xcopy $dir/
  cp $xchange $dir/

  # distribute adjoint source and stations
  echo distribute adjoint sources and stations
  if [ ! -f $adjdir/STATIONS_ADJOINT ]; then
    echo WRONG! NO $adjdir/STATIONS_ADJOINT
    exit
  fi
  cp $adjdir/*.LH[ZEN].adj $dir/SEM/
  cp $adjdir/STATIONS_ADJOINT $dir/DATA

  cat $cmtcenter/$line > $dir/DATA/CMTSOLUTION
  cat $stacenter/$stafile > $dir/DATA/STATIONS

  cd $dir/DATA

  if [ ! -f Par_file ]; then
    echo WRONG! NO Par_file
    exit
  fi
  echo "   change Par_file...   "
  path="\/scratch\/hejunzhu"
  sed -e "s/^LOCAL_PATH.*$/LOCAL_PATH               = $path/g" \
      Par_file > Par_file_out
  mv Par_file_out Par_file


  cd ..
  if [ ! -f xparallel_adjoint_solver.sh ]; then
    echo WRONG! NO xparallel_adjoint_solver.sh
    exit
  fi
  echo "    changing pbs file   "
  tag="#PBS -N  $line"
  sed -e "s/^#PBS -N.*$/$tag/g" xparallel_adjoint_solver.sh > xparallel_adjoint_solver_out.sh
  mv xparallel_adjoint_solver_out.sh xparallel_adjoint_solver.sh


  # check code exist
  if [ ! -f xspecfem3D ] ;then
    echo WRONG! NO xspecfem3D
    exit
  fi
  if [ ! -f xcopy_local_forward ]; then
    echo WRONG! NO xcopy_local_forward
    exit
  fi
  if [ ! -f change_simulation_type.pl ]; then
    echo WRONG! NO change_simulation_type.pl
    exit
  fi
  mimic=`grep USE_ATTENUATION_MIMIC OUTPUT_FILES/values_from_mesher.h | awk -F"=" '{print $NF}'`
  if [ $mimic != ".true." ]; then
    echo WRONG! NEED reompile solver with MIMIC attenuation
    exit
  fi


  echo submitting job
  qsub xparallel_adjoint_solver.sh

  sleep 3
  cd ../../

done < $eventfile
