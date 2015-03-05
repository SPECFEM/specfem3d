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
  mea_output="../SYN_"$iter"/"$line"_ADJ_"$ext
  adjstafile="../SYN_"$iter"/"$line"_ADJ_"$ext"/STATIONS_ADJOINT"

  if [ ! -f $adjstafile ] ;then
    echo WRONG!!!
  fi

  n=`wc -l $adjstafile | awk '{print $1}'`
  n1=`echo "($n*3)" | bc -l`

  nadj=`ls $mea_output/*LH[ENZ].adj | wc -l`

  echo $line $n1 $nadj

  if [ $n1 -ne $nadj ]; then
    echo WRONG!!!!
  fi

done < $eventfile
