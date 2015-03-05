#!/bin/sh
# Author: Hejun Zhu, hejunzhu@princeton.edu
# Princeton University, New Jersey, USA
# Last modified: Thu Aug 26 15:08:23 EDT 2010

# This script uses combine_adj_src.pl to combine two
# sets of adjoint sources



eventfile=../EVENTID_CENTER/XEVENTID
iter=M18

while read line
do
  echo $line

  cmtid=`echo $line | awk -F"_" '{print $NF}'`
  adjdir1="../SYN_"$iter"/CMTSOLUTION_"$cmtid"_ADJ_T015_040"
  adjdir2="../SYN_"$iter"/CMTSOLUTION_"$cmtid"_ADJ_T025_150"
  adjdir_new="../SYN_"$iter"/CMTSOLUTION_"$cmtid"_ADJ_comb"

  stafile1="../SYN_"$iter"/CMTSOLUTION_"$cmtid"_ADJ_T015_040/STATIONS_ADJOINT"
  stafile2="../SYN_"$iter"/CMTSOLUTION_"$cmtid"_ADJ_T025_150/STATIONS_ADJOINT"
  stafile_new="../SYN_"$iter"/CMTSOLUTION_"$cmtid"_ADJ_comb/STATIONS_ADJOINT"

  if [ ! -d $adjdir1 ]; then
    echo WRONG! NO $adjdir1
    exit
  fi
  if [ ! -d $adjdir2 ]; then
    echo WRONG! NO $adjdir2
    exit
  fi
  if [ ! -f $stafile1 ]; then
    echo WRONG! NO $stafile1
    exit
  fi
  if [ ! -f $stafile2 ]; then
    echo WRONG! NO $stafile2
    exit
  fi
  if [ ! -d $adjdir_new ]; then
    echo mkdir $adjdir_new
    mkdir $adjdir_new
  fi

  ./combine_adj_src.pl $adjdir1 $adjdir2 $adjdir_new $stafile1 $stafile2 $stafile_new

done < $eventfile
