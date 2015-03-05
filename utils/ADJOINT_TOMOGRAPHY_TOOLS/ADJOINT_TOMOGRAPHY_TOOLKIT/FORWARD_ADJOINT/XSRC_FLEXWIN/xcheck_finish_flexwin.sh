#!/bin/sh
# Author: Hejun Zhu, hejunzhu@princeton.edu
# Princeton University, New Jersey, USA
# Last modified: Mon Oct  4 15:48:11 EDT 2010


# This script is used to check whether FLEXWIN finish completely


eventfile=../EVENTID_CENTER/XEVENTID_TMP
iter=M18
fmin=30
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
  flexwin_output="XFLEXWIN_OUTPUT_"$iter"/OUTPUT_FLEXWIN_"$cmtid"_"$ext
  mefile="../MEASUREMENT_CENTER/MEASUREMENT_"$cmtid

  if [ ! -f $flexwin_output ]; then
    echo WRONG! NO $flexwin_output
    exit
  fi
  if [ ! -f $mefile ]; then
    echo WRONG! NO $mefile
    exit
  fi


  lastline=`tail -1 $mefile`
  tag=`echo $lastline| awk -F" " '{print $1}'`

  exist=`grep $tag $flexwin_output | awk '{print $1}'`
  echo $line $exist


done < $eventfile
