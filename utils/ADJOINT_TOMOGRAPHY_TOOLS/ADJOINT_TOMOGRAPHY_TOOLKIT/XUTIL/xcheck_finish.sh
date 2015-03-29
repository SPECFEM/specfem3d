#!/bin/sh
# Author: Hejun Zhu, hejunzhu@princeton.edu
# Princeton University, New Jersey, USA
# Last modified: Wed Feb 22 20:58:13 EST 2012


ntot=0
for dir in CMTSOLUTION_*
do

  file=$dir/job_src2.log
  if [ -f $file ]; then
    ntot=`echo $ntot + 1| bc -l`
    tag="successfully"
    n=`wc -l $file | awk '{print $1}'`
    text=`grep $tag $file`
    echo $dir   $text    $n
  fi

done

echo found total $ntot events
