#!/bin/bash
# this script checks if all the files listed in MEASURE*ALL file exist or not
if [ $# != 1 ]; then
  echo "check_ds_exist.bash event-list"
  exit
fi

for evid in `cat $1 `; do
cd $evid
for file in ` cat MEASURE*ALL | grep '^[ds]' `; do
  if [ ! -f $file ]; then
    echo $file
  fi
done
cd ..
done
