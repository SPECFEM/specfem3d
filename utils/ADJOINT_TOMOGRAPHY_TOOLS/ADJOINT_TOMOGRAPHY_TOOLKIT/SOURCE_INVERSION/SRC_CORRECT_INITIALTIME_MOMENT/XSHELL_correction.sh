#!/bin/sh
# Author: Hejun Zhu, hejunzhu@princeton.edu
# Princeton University, New Jersey, USA
# Last modified: #!/bin/sh
# Author: Hejun Zhu, hejunzhu@princeton.edu
# Princeton University, New Jersey, USA
# Last modified: Sun Jun  3 17:19:13 EDT 2012

eventfile=../../SHARE_FILES/EVENTID_CENTER/XEVENTID
iter1=M42
iter2=M43

ext1=T015_040
ext2=T030_100


dirold=../SYN_$iter1
dirnew=../SYN_$iter2
dirinput=XCORRECT_INPUT_$iter1

# check directories
if [ ! -f $eventfile ]; then
  echo WRONG! NO $eventfile
  exit
fi
if [ ! -d $dirold ]; then
  echo WRONG! NO $dirold
  exit
fi
if [ ! -d $dirinput ]; then
  echo MKDIR $dirinput
  mkdir $dirinput
fi
if [ ! -d $dirnew ]; then
  echo MKDIR $dirnew
  mkdir $dirnew
fi

while read line
do
  echo running correction for $line...
  cmtid=`echo $line | awk -F"_" '{print $NF}'`
  minmax=../SRC_GRIDSEARCH_INITIALTIME_MOMENT/XGRIDSEARCH_OUTPUT_$iter1/XMINMAX_$cmtid
  input=$dirinput/CORRECT_$cmtid
  synold=$dirold/CMTSOLUTION_$cmtid
  synnew=$dirnew/CMTSOLUTION_$cmtid

  tag="#PBS -N XXCORRECT_$cmtid"
  input_tag="$dirinput\/CORRECT_$cmtid"

  # check directories
  if [ ! -f $minmax ]; then
    echo WRONG! NO $minmax
    exit
  fi
  if [ -f $input ]; then
    echo $input exist, rm it
    rm $input
  fi
  if [ ! -d $synold ]; then
    echo WRONG! NO $synold
    exit
  fi
  if [ ! -d $synnew ]; then
    echo NO $synnew, make it
    mkdir $synnew
  fi
  if [ -f xtmp ]; then
    echo xtmp exist, rm it
    rm xtmp
  fi

  dt0=`cat $minmax | awk '{print $1}'`
  dm0=`cat $minmax | awk '{print $2}'`

  ls $synold/*.$ext1 > xtmp
  ls $synold/*.$ext2 >> xtmp

  nfile=`wc -l xtmp | awk '{print $1}'`

  echo "$dt0    $dm0 " > $input
  echo $nfile >> $input
  while read line
  do
    name=`echo $line | awk -F"/" '{print $NF}'`
    name_new=$synnew/$name
    echo $line >> $input
    echo $name_new >> $input
  done < xtmp

  rm xtmp


  sed -e "s/^#PBS -N.*$/$tag/g" \
  -e "s/^input=.*$/input=$input_tag/g" \
  XPBS_correction.sh > XPBS_correction.sh.out
  mv XPBS_correction.sh.out XPBS_correction.sh

  echo RUNNING CORRECTION FOR EVENTID $cmtid
  qsub XPBS_correction.sh
  sleep 3


done < $eventfile

