#!/bin/bash

run_pwd=$PWD

#sleep 120m

for eid in `ls -1 -d [0-9]*`;do
#for eid in `ls -1 -d 13935988 9703873 9818433 9735129 9105672 10223765`; do

 cd $run_pwd
 cd $eid
 echo "==== $eid ====="

 #run.lsf
 run.lsf_hero
 #sleep 15m
 sleep 5s

done
