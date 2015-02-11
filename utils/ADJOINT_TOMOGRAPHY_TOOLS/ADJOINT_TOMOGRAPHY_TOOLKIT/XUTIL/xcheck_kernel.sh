#!/bin/sh
# Author: Hejun Zhu, hejunzhu@princeton.edu
# Princeton University, New Jersey, USA
# Last modified: Fri Feb 17 16:39:24 EST 2012


for dir in CMTSOLUTION_*
do
  n=`ls $dir/KERNEL/*.bin | wc -l`
  echo $dir $n

done
