#!/bin/sh
# Author: Hejun Zhu, hejunzhu@princeton.edu
# Princeton University, New Jersey, USA
# Last modified: Thu Nov 18 15:31:44 EST 2010


eventfile=../SHARE_FILES/EVENTID_CENTER/XEVENTID
iter=M18

# Check files
if [ ! -f $eventfile ]; then
  echo WRONG! NO $eventfile
  exit
fi

if [ ! -d SYN_$iter ]; then
  echo mkdir SYN_$iter
  mkdir SYN_$iter
fi


while read line
do
  echo $line
  cmtid=`echo $line | awk -F"_" '{print $2}'`
  mefile=MEASUREMENT_CENTER/MEASUREMENT_$cmtid
  dir="SYN_"$iter"/"$line
  if [ ! -d $dir ]; then
    echo mkdir $dir
    mkdir $dir
  fi
  if [ ! -f $mefile ]; then
    echo WRONG! NO $mefile
    exit
  fi

  while read name
  do
    if [ ! -f FORWARD_$iter/$line/OUTPUT_FILES/$name.LHZ.sem.sac ]; then
      echo WRONG! NO FORWARD_$iter/$line/OUTPUT_FILES/$name.LHZ.sem.sac
      exit
    fi
    if [ ! -f FORWARD_$iter/$line/OUTPUT_FILES/$name.LHE.sem.sac ]; then
      echo WRONG! NO FORWARD_$iter/$line/OUTPUT_FILES/$name.LHE.sem.sac
      exit
    fi
    if [ ! -f FORWARD_$iter/$line/OUTPUT_FILES/$name.LHN.sem.sac ]; then
      echo WRONG! NO FORWARD_$iter/$line/OUTPUT_FILES/$name.LHN.sem.sac
      exit
    fi


    cp FORWARD_$iter/$line/OUTPUT_FILES/$name.*.sac SYN_$iter/$line
  done < $mefile

done < $eventfile
