#!/bin/sh
# Author: Hejun Zhu, hejunzhu@princeton.edu
# Princeton University, New Jersey, USA
# Last modified: Thu Feb 23 17:27:34 EST 2012

for dir in CMTSOLUTION_*
do
  n=`ls $dir/OUTPUT_FILES/*.sem.sac | wc -l`
  echo $dir $n

done
