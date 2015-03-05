#!/bin/sh
# Author: Hejun Zhu, hejunzhu@princeton.edu
# Princeton University, New Jersey, USA
# Last modified: Thu Sep  6 07:09:24 EDT 2012


iter=M42

ext1=T015_040
ext2=T030_100

eventfile=../EVENTID_CENTER/XEVENTID

meadir=XMEASUREMENT_INPUT_$iter
grddir=XGRIDSEARCH_INPUT_$iter

if [ ! -d $meadir ]; then
  echo WRONG! NO $meadir
  exit
fi
if [ ! -d $grddir ]; then
  echo MKDIR $grddir
  mkdir $grddir
fi
if [ ! -f $eventfile ]; then
  echo WRONG! NO $eventfile
  exit
fi

while read line
do
  echo $line
  cmtid=`echo $line | awk -F"_" '{print $NF}'`
  m1=$meadir/MEASUREMENT_$cmtid'_'$ext1
  m2=$meadir/MEASUREMENT_$cmtid'_'$ext2
  ofile=$grddir/GRIDSEARCH_$cmtid
  if [ ! -f $m1 ]; then
    echo WRONG! NO $m1
    exit
  fi
  if [ ! -f $m2 ]; then
    echo WRONG! NO $m2
    exit
  fi
  if [ -f $ofile ]; then
    echo DELETE $ofile
    rm $ofile
  fi

  n1=`sed -n '5p' $m1`
  n2=`sed -n '5p' $m2`

        bd_z=`sed -n '1p' $m1 | awk '{print $2}'`
  bd_r=`sed -n '1p' $m1 | awk '{print $3}'`
  bd_t=`sed -n '1p' $m1 | awk '{print $4}'`
  sw_z=`sed -n '1p' $m2 | awk '{print $2}'`
  sw_r=`sed -n '1p' $m2 | awk '{print $3}'`
  sw_t=`sed -n '1p' $m2 | awk '{print $4}'`

  ntotal=`echo $n1 + $n2 | bc -l`

  echo $ntotal $bd_z $bd_r $bd_t $sw_z $sw_r $sw_t > $ofile
  awk 'FNR>5' $m1 >> $ofile
  awk 'FNR>5' $m2 >> $ofile


done < $eventfile


