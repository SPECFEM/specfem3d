#!/bin/sh

# This script is used to process data and syn
# It used process_data.pl process_syn.pl
# Hejun Zhu Jun 07 2010

# INPUT:
# fmin :   minimum period
# fmax :   maximum period
# tstart:  cut start time
# tend:    cut end time
# cmtid:   event id
# opt:     1=process dat; 2=process syn; 3=cut data and synthetics
# Bandpass 15-25sec, 20-34sec, 26-56sec, 50-120sec (15-35sec, 50-120sec)


# Input
opt=2
iter=M13
fmin=25
fmax=150
eventfile=../SHARE_FILES/EVENTID_CENTER/XEVENTID
tstart=0
tend=1800


# Get names
datpnm=DATASET_EUROPE
synpnm="SYN_"$iter
ext1=`printf "%03i\n" $fmin`
ext2=`printf "%03i\n" $fmax`
ext="T"$ext1"_"$ext2


# Check files
if [ ! -f $eventfile ]; then
  echo WRONG! NO $eventfile
  exit
fi


while read line
do
  cmtfile="../SHARE_FILES/CMTSOLUTION_CENTER/$line"
  cmtfiletag="..\/SHARE_FILES\/CMTSOLUTION_CENTER\/$line"
  cmtid=`echo $line | awk -F"_" '{print $2}'`
  proctag="#PBS -N XPROC_DAT_$cmtid"
  datdir="$datpnm\/CMTSOLUTION_$cmtid"
  syndir="$synpnm\/$line"

  if [ ! -f $cmtfile ]; then
    echo WRONG! NO $cmtfile
    exit
  fi

  sed -e "s/^#PBS -N.*$/$proctag/g" \
      -e "s/^cmtfile=.*$/cmtfile=$cmtfiletag/g" \
      -e "s/^tstart=.*$/tstart=$tstart/g" \
      -e "s/^tend=.*$/tend=$tend/g" \
      -e "s/^fmin=.*$/fmin=$fmin/g" \
      -e "s/^fmax=.*$/fmax=$fmax/g" \
      -e "s/^ext=.*$/ext=$ext/g" \
      -e "s/^datdir=.*$/datdir=$datdir/g" \
      -e "s/^syndir=.*$/syndir=$syndir/g" \
      -e "s/^opt=.*$/opt=$opt/g" \
      XPBS_process_dat.sh > XPBS_process_dat.sh.out
      mv XPBS_process_dat.sh.out XPBS_process_dat.sh

  echo qsub $line
  qsub XPBS_process_dat.sh
  sleep 5
done < $eventfile


