#!/bin/sh
# Author: Hejun Zhu, hejunzhu@princeton.edu
# Princeton University, New Jersey, USA
# Last modified: Fri Jan 21 14:54:09 EST 2011

# This script is used to generate adjoint source for ENZ component

eventfile=../EVENTID_CENTER/XEVENTID
iter=M18
fmin=25
fmax=150

ext1=`printf "%03i\n" $fmin`
ext2=`printf "%03i\n" $fmax`
ext="T"$ext1"_"$ext2

if [ ! -f $eventfile ] ;then
  echo WRONG! NO $eventfile
  exit
fi


while read line
do

  echo $line
  cmtid=`echo $line | awk -F"_" '{print $NF}'`
  cmtfile="../CMTSOLUTION_CENTER/CMTSOLUTION_"$cmtid
  stafile="../STATIONS_CENTER/STATIONS.ORFEUS"
  meaoutput="../SYN_"$iter"/CMTSOLUTION_"$cmtid"_MT_"$ext
  adjoutput="../SYN_"$iter"/CMTSOLUTION_"$cmtid"_ADJ_"$ext
# datdir="../DATASET_EUROPE/CMTSOLUTION_"$cmtid
# syndir="../SYN_"$iter"/CMTSOLUTION_"$cmtid
# mefile="XMEASUREMENT_INPUT_"$iter"/MEASUREMENT_"$cmtid"_"$ext

  if [ ! -f $cmtfile ]; then
    echo WRONG! NO $cmtfile
    exit
  fi
  if [ ! -f $stafile ]; then
    echo WRONG! NO $stafile
    exit
  fi
  if [ ! -d $meaoutput ]; then
    echo WRONG! NO $meaoutput
    exit
  fi
  if [ ! -d $adjoutput ]; then
    echo MKDIR $adjoutput
    mkdir $adjoutput
  fi
# if [ ! -d $datdir ] ;then
#   echo WRONG! NO $datdir
#   exit
# fi
# if [ ! -d $syndir ]; then
#   echo WRONG! NO $syndir
#   exit
# fi
# if [ ! -f $mefile ]; then
#   echo WRONG! NO $mefile
#   exit
# fi


  # rotate adjoint source
  ./prepare_adj_src.pl -m $cmtfile -s $stafile -z LH -o $adjoutput $meaoutput/*.adj

  # plot adjoint source
  while read name
  do

    sta=`echo $name | awk '{print $1}'`
    net=`echo $name | awk '{print $2}'`


    for comp in Z E N R T
    do
      #if [ $ext == "T050_150" ]; then
      #if [ $ext == "T040_150" ]; then
        fileold=$adjoutput"/"$sta"."$net".LH"$comp".iker07.adj"
      #elif [ $ext == "T015_050" ]; then
      #elif [ $ext == "T015_040" ]; then
      # fileold=$adjoutput"/"$sta"."$net".LH"$comp".iker05.adj"
      #else
      # echo WRONG! $ext
      # exit
      #fi
      #if [ ! -f $fileold ]; then
      # echo WRONG! NO $fileold
      # exit
      #fi
      filenew=$adjoutput"/"$sta"."$net".LH"$comp".adj"
      if [ -f $fileold ]; then
        mv $fileold $filenew
      fi
    done


  done < STATIONS_ADJOINT

  # remove station file
  mv STATIONS_ADJOINT $adjoutput/

done < $eventfile


