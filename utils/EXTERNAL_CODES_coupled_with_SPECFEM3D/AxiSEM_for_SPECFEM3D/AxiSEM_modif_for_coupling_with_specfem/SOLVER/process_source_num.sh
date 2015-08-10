#!/bin/bash

numsrc=$1

lala=`squeue -u beller | grep beller | awk '{printf("%d:",$1)}'`
job2wait=`echo ${lala%:}#job1=$2`
#job2=$3
#job3=$4
#job4=$5

cp SAVE_CMTSOLUTIONS/CMTSOLUTION_$1 CMTSOLUTION

./submit_sb2.csh Source_$1


