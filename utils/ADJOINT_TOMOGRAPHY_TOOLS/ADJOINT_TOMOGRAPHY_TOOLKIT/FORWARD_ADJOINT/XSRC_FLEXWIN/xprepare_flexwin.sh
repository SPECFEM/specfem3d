#!/bin/sh
# Author: Hejun Zhu, hejunzhu@princeton.edu
# Princeton University, New Jersey, USA
# Last modified: Wed Aug 18 10:02:23 EDT 2010

# This script is used to prepare flexwin input file
# flexwin.input
# INPUT:
#   cmtid:    cmtid from Harvard Solution
#   iter:     iteration number
#   fmin:     minimum period
#   fmax:     maximum period (Tmin and Tmax are used to generate file extension)
# Name:
#   datdir:   data directory
#   syndir:   syn directory
#   windir:   win directory (save all output from FLEXWIN)



eventfile=../EVENTID_CENTER/XEVENTID_00
iter=M01
fmin=25
fmax=150

datpnm=../DATASET_EUROPE
synpnm="../SYN_"$iter
ext1=`printf "%03i\n" $fmin`
ext2=`printf "%03i\n" $fmax`
ext="T"$ext1"_"$ext2

if [ ! -f $eventfile ]; then
  echo WRONG! NO $eventfile
  exit
fi

if [ -f t1 ] ; then
  echo deleting t1
  rm t1
fi

flexwin_input_dir="XFLEXWIN_INPUT_"$iter
flexwin_output_dir="XFLEXWIN_OUTPUT_"$iter

if [ ! -d $flexwin_input_dir ]; then
  echo MKDIR $flexwin_input_dir
  mkdir $flexwin_input_dir
fi
if [ ! -d $flexwin_output_dir ]; then
  echo MKDIR $flexwin_output_dir
  mkdir $flexwin_output_dir
fi


while read name
do

  cmtid=`echo $name | awk -F"_" '{print $NF}'`

  mefile="../MEASUREMENT_CENTER/MEASUREMENT_"$cmtid
  datdir=$datpnm"/CMTSOLUTION_"$cmtid
  syndir=$synpnm"/CMTSOLUTION_"$cmtid
  windir=$synpnm"/CMTSOLUTION_"$cmtid"_WIN_"$ext
  flexwin_input=$flexwin_input_dir"/INPUT_FLEXWIN_"$cmtid"_"$ext

  echo checking INPUT
  echo DATDIR:      $datdir
  echo SYNDIR:      $syndir
  echo WINDIR:      $windir
  echo MEASURE:           $mefile


  if [ ! -f $mefile ]; then
    echo WRONG! NO $mefile
    exit
  fi

  nsta=`wc -l $mefile | awk '{print $1}'`
  echo found $nsta stations


  n=0
  while read line
  do
    sta=`echo $line | awk '{print $1}'| awk -F"." '{print $1}'`
    net=`echo $line | awk '{print $1}'| awk -F"." '{print $2}'`

    datz=$datdir"/"$sta"."$net".BHZ.sac."$ext
    datr=$datdir"/"$sta"."$net".BHR.sac."$ext
    datt=$datdir"/"$sta"."$net".BHT.sac."$ext

    synz=$syndir"/"$sta"."$net".LHZ.sem.sac."$ext
    synr=$syndir"/"$sta"."$net".LHR.sem.sac."$ext
    synt=$syndir"/"$sta"."$net".LHT.sem.sac."$ext


    if [ ! -f $datz ]; then
      echo WRONG! NO $datz
      exit
    fi
    if [ ! -f $datr ]; then
      echo WRONG! NO $datr
      exit
    fi
    if [ ! -f $datt ]; then
      echo WRONG! NO $datt
      exit
    fi
    if [ ! -f $synz ]; then
      echo WRONG! NO $synz
      exit
    fi
    if [ ! -f $synr ]; then
      echo WRONG! NO $synr
      exit
    fi
    if [ ! -f $synt ]; then
      echo WRONG! NO $synt
      exit
    fi


    winz=$windir"/"$sta"."$net".LHZ"
    winr=$windir"/"$sta"."$net".LHR"
    wint=$windir"/"$sta"."$net".LHT"

    n=`echo $n | awk '{print $n+1}'`
    echo $datz >> t1
    echo $synz >> t1
    echo $winz >> t1

    n=`echo $n | awk '{print $n+1}'`
    echo $datr >> t1
    echo $synr >> t1
    echo $winr >> t1

    n=`echo $n | awk '{print $n+1}'`
    echo $datt >> t1
    echo $synt >> t1
    echo $wint >> t1

  done < $mefile

  echo $n > $flexwin_input
  cat t1 >> $flexwin_input
  rm t1
done < $eventfile

