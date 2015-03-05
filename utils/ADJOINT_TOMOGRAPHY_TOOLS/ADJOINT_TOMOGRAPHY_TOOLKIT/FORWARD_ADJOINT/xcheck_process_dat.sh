#!/bin/sh

ntot=0
for line in XPROC_DAT_*
do

  n=`wc -l $line | awk '{print $1}'`
  ntot=`echo $ntot +1 | bc -l`
  tag="successfully"
  tag1="DONE"
  text=`grep $tag $line`
  text1=`grep $tag1 $line`
  echo $line   $text  $text1 $n
done
echo found $ntot finshed
