#!/bin/sh
# Author: Hejun Zhu, hejunzhu@princeton.edu
# Princeton University, New Jersey, USA
# Last modified: Fri Jan 21 13:18:09 EST 2011


# this script is used to run measure_adj parallelly
# INPUT :  measurement from FLEXWIN as, XMEASUREMENT_INPUT_M00/MEASUREMENT_012202A_T050_120
# OUTPUT:  adjoint source saved in as, ../SYN_M00/CMTSOLUTION_012204A_MT_T050_120
#          window_chi and window_index in XWINDOW_OUTPUT_M00/WINDOW_CHI_012202A_T050_120 and XWINDOW_OUTPUT_M00/WINDOW_INDEX_012202A_T050_120
#          outputtag in XOUTPUT_TAG_M00/OUTPUT_012202A_T050_120


eventfile=../EVENTID_CENTER/XEVENTID
iter=M18
fmin=15
fmax=40

ext1=`printf "%03i\n" $fmin`
ext2=`printf "%03i\n" $fmax`
ext="T"$ext1"_"$ext2

if [ ! -f $eventfile ] ;then
  echo WRONG! NO $eventfile
  exit
fi

if [ ! -d XOUTPUT_TAG_$iter ]; then
  echo mkdir XOUTPUT_TAG_$iter
  mkdir XOUTPUT_TAG_$iter
fi
if [ ! -d XWINDOW_OUTPUT_$iter ]; then
  echo mkdir XWINDOW_OUTPUT_$iter
  mkdir XWINDOW_OUTPUT_$iter
fi


while read line
do
  echo $line
  cmtid=`echo $line | awk -F"_" '{print $NF}'`
  # input for measure_adj
  inputfile="XMEASUREMENT_INPUT_"$iter"/MEASUREMENT_"$cmtid"_"$ext
  inputfile_tag="XMEASUREMENT_INPUT_"$iter"\/MEASUREMENT_"$cmtid"_"$ext
  if [ ! -f $inputfile ]; then
    echo WRONG! NO $inputfile
    exit
  fi
  outdir=`sed -n '2p' $inputfile`
  if [ ! -d $outdir ]; then
    echo mkdir $outdir
    mkdir $outdir
  fi
  output="XOUTPUT_TAG_"$iter"\/OUTPUT_"$cmtid"_"$ext
  meatag="#PBS -N XMEASURE_$cmtid"
  syndir="..\/SYN_"$iter
  mtdir="CMTSOLUTION_"$cmtid"_MT_"$ext

  sed -e "s/^#PBS -N.*$/$meatag/g" \
      -e "s/^inputfile=.*$/inputfile=$inputfile_tag/g" \
      -e "s/^output=.*$/output=$output/g" \
      -e "s/^syndir=.*$/syndir=$syndir/g" \
      -e "s/^mtdir=.*$/mtdir=$mtdir/g" \
      xrun_measure_adj_pbs.sh > xrun_measure_adj_pbs.sh.out
      mv xrun_measure_adj_pbs.sh.out xrun_measure_adj_pbs.sh

  qsub xrun_measure_adj_pbs.sh
  sleep 10
done < $eventfile
