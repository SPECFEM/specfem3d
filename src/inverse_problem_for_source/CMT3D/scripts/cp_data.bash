#!/bin/bash
# this is an old program used to copy tar'ed derivative synthetics from caltech machines to local machines.

if [ $# != 1 ]; then
  echo "Usage: cp_data.bash evid-file"; exit
fi

dir=/net/sierra/raid1/carltape/socal/socal_3D/DATA/FINAL

for evid in `cat $1`; do
  mkdir -p $evid/data $evid/backup; cd $evid/data
  echo " == $evid ==="
  ssh -l lqy wagholi.gps.caltech.edu "cd $dir/$evid; tar cjvf ~/data_${evid}.tar.bz *.sac >/dev/null "
  scp lqy@wagholi.gps.caltech.edu:data_${evid}.tar.bz .
#  ssh wagholi.gps.caltech.edu "rm -f data_${evid}.tar.bz"
  mv ../../data/data_${evid}.tar.bz .
  tar xjvf data_${evid}.tar.bz > /dev/null; mv data_${evid}.tar.bz ../backup
  cd ../..
done
