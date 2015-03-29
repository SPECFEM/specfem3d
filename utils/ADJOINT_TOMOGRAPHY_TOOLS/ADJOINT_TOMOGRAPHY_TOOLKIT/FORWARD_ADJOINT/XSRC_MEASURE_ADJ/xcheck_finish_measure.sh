#!/bin/sh
# Author: Hejun Zhu, hejunzhu@princeton.edu
# Princeton University, New Jersey, USA
# Last modified: Fri Feb 11 14:03:25 EST 2011

eventfile=../EVENTID_CENTER/XEVENTID
iter=M18
fmin=25
fmax=150

ext1=`printf "%03i\n" $fmin`
ext2=`printf "%03i\n" $fmax`
ext="T"$ext1"_"$ext2




if [ ! -f $eventfile ]; then
  echo WRONG! NO $eventfile
  exit
fi
while read line
do

  cmtid=`echo $line | awk -F"_" '{print $NF}'`
  mea_input="XMEASUREMENT_INPUT_"$iter"/MEASUREMENT_"$cmtid"_"$ext
  mea_output="../SYN_"$iter"/"$line"_MT_"$ext

  nwin=`sed -n '5p' $mea_input`

  nadj=`ls $mea_output/*.adj | wc -l`

  echo $line $nwin $nadj
  if [ $nwin -ne $nadj ]; then
    echo WRONG!!!!
  fi

done < $eventfile
