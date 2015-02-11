#!/bin/sh
# Author: Hejun Zhu, hejunzhu@princeton.edu
# Princeton University, New Jersey, USA
# Last modified: Mon Jun  4 16:25:38 EDT 2012


# input parameters
iter1=M42
iter2=M43

mstart=0.2
mend=1.8
dm=0.01

tstart=-5
tend=5
dt=0.1

fact_am=0.857
fact_tt=0.143
linesearch_mbest=0.02

body_bandpass=T015_040
surf_bandpass=T040_100

criteria=phase_amplitude

eventfile=../../SHARE_FILES/EVENTID_CENTER/XEVENTID
dir1=../../SHARE_FILES/CMTSOLUTION_CENTER_$iter1
dir2=../../SHARE_FILES/CMTSOLUTION_CENTER_$iter2

# check directories and files
if [ ! -f $eventfile ]; then
  echo WRONG! NO $eventfile
  exit
fi

if [ ! -d $dir1 ]; then
  echo wrong! no $dir1
  exit
fi

if [ ! -d $dir2 ]; then
  echo mkdir $dir2
  mkdir $dir2
fi

if [ ! -d XGRIDSEARCH_OUTPUT_$iter1 ];then
  echo MKDIR XGRIDSEARCH_OUTPUT_$iter1
  mkdir XGRIDSEARCH_OUTPUT_$iter1
fi

while read line
do

  cmtid=`echo $line | awk -F"_" '{print $NF}'`
  grd_tag="#PBS -N XOUTPUT_GRID_$cmtid"
  parinput=XGRIDSEARCH_INPUT_$iter1/PAR_$cmtid
  parinput_tag="XGRIDSEARCH_INPUT_"$iter1"\/PAR_"$cmtid


  cmtorig="../../SHARE_FILES/CMTSOLUTION_CENTER_$iter1/CMTSOLUTION_$cmtid"
  cmtnew="../../SHARE_FILES/CMTSOLUTION_CENTER_$iter2/CMTSOLUTION_$cmtid"
  gridinput="XGRIDSEARCH_INPUT_$iter1/GRIDSEARCH_"$cmtid
  gridoutput1="XGRIDSEARCH_OUTPUT_$iter1/XMINMAX_$cmtid"
  gridoutput2="XGRIDSEARCH_OUTPUT_$iter1/XMISFITFUNCTION_$cmtid"


  if [ -f $parinput ] ;then
    echo $parinput exist, delete it
    rm $parinput
  fi

  if [ ! -f $cmtorig ]; then
    echo WRONG! NO $cmtorig
    exit
  fi

  if [ ! -f $gridinput ]; then
    echo WRONG! NO $gridinput
    exit
  fi

  echo GERENATING... $parinput
  echo $cmtorig > $parinput
  echo $cmtnew >> $parinput
  echo $gridinput >> $parinput
  echo $gridoutput1 >> $parinput
  echo $gridoutput2 >> $parinput
  echo "$mstart    $mend    $dm" >> $parinput
  echo "$tstart    $tend    $dt"  >> $parinput
  echo "$fact_am" >> $parinput
  echo "$fact_tt" >> $parinput
  echo "$linesearch_mbest" >> $parinput
  echo "$body_bandpass" >> $parinput
  echo "$surf_bandpass" >> $parinput
        echo "$criteria" >> $parinput


  sed -e "s/^#PBS -N.*$/$grd_tag/g" \
  -e "s/^input=.*$/input=$parinput_tag/g" \
  XPBS_gridsearch.sh > XPBS_gridsearch.sh.out
  mv XPBS_gridsearch.sh.out XPBS_gridsearch.sh


  echo RUNNING GRIDSEARCH FOR EVENT $cmtid
  qsub XPBS_gridsearch.sh
  sleep 3

done < $eventfile
