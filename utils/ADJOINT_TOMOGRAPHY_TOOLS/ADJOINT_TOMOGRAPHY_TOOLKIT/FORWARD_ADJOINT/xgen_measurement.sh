#!/bin/sh
# This script is used to generate a file which contain the same
# data and synthetics
# format for MEASUREMENT file is sta.net  stlo   stla   dist  az
# Hejun Zhu, Jun 24, 2010

eventfile=../SHARE_FILES/EVENTID_CENTER/XEVENTID
iter=M00

while read line
do
  echo $line
  eid=`echo $line | awk -F"_" '{print $NF}'`
  datdir=DATASET_EUROPE/$line
  syndir="SYN_"$iter"/"$line
  mefile="MEASUREMENT_CENTER/MEASUREMENT_"$eid

  if [ -f $mefile ]; then
    echo deleting $mefile
    rm $mefile
  fi
  if [ ! -d $datdir ]; then
    echo WRONG! NO $datdir
    exit
  fi
  if [ ! -d $syndir ]; then
    echo WRONG! NO $syndir
    exit
  fi


  for file1 in $datdir/*BHE.sac; do
    sta=`echo $file1 | awk -F"/" '{print $NF}' | awk -F"." '{print $1}'`
    net=`echo $file1 | awk -F"/" '{print $NF}' | awk -F"." '{print $2}'`

    file2=$datdir"/"$sta"."$net".BHE.sac"


    file4=$syndir"/"$sta"."$net".LHE.sem.sac"
    if [ -f $file4 ] ; then
      #echo found $sta $net for $line
      #echo $sta"."$net $stlo $stla   $dist $az >> $mefile
      echo $sta"."$net  >> $mefile
    fi
  done
done < $eventfile



