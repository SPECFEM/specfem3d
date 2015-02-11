#!/bin/sh
# Author: Hejun Zhu, hejunzhu@princeton.edu
# Princeton University, New Jersey, USA
# Last modified: Wed Aug 18 10:00:14 EDT 2010

# This script is used to run flexwin
# using flexwin.input in XFLEXWIN_INPUT_OUTPUT
# and generate flexwin.output in XFLEXWIN_INPUT_OUTPUT


eventfile=../EVENTID_CENTER/XEVENTID_00
iter=M13
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
  cmtid=`echo $line | awk -F"_" '{print $NF}'`
  flexwin_input="XFLEXWIN_INPUT_"$iter"\/INPUT_FLEXWIN_"$cmtid"_"$ext
  flexwin_output="XFLEXWIN_OUTPUT_"$iter"\/OUTPUT_FLEXWIN_"$cmtid"_"$ext
  flexwin_tag="#PBS -N XFLEXWIN_$cmtid"
  syndir="..\/SYN_$iter"
  dir="CMTSOLUTION_"$cmtid"_WIN_"$ext

  dir2=../SYN_$iter/$dir
  if [ ! -d $dir2 ]; then
    echo MKDIR $dir2
    mkdir $dir2
  fi


  sed -e "s/^#PBS -N.*$/$flexwin_tag/g" \
      -e "s/^input=.*$/input=$flexwin_input/g" \
      -e "s/^output=.*$/output=$flexwin_output/g" \
      -e "s/^syndir=.*$/syndir=$syndir/g" \
      -e "s/^dir=.*$/dir=$dir/g" \
      xrun_flexwin_pbs.sh > xrun_flexwin_pbs.sh.out
      mv xrun_flexwin_pbs.sh.out xrun_flexwin_pbs.sh

  # change PAR_FILE
  #echo replace PAR_FILE with fmax = $fmax and fmin = $fmin
  #tag_min="WIN_MIN_PERIOD      = $fmin.0"
  #tag_max="WIN_MAX_PERIOD      = $fmax.0"
  #sed -e "24s/^.*$/$tag_min/g" PAR_FILE > PAR_FILE.sed
  #mv PAR_FILE.sed PAR_FILE
  #sed -e "25s/^.*$/$tag_max/g" PAR_FILE > PAR_FILE.sed
  #mv PAR_FILE.sed PAR_FILE

  echo running flexwin for fmax = $fmax and fmin = $fmin for $cmtid
  #./flexwin < $flexwin_input > $flexwin_output
  qsub xrun_flexwin_pbs.sh
  sleep 10
done < $eventfile

