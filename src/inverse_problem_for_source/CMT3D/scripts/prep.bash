#!/bin/bash

# this is the program that uses process_data_and_syn.pl to pre-process data and synthetics for cmt3d/grid3d inversions

# add data and syn headers
../scripts/process_data_and_syn.pl -d data,sac -s syn, -0 -m CMTSOLUTION -a STATIONS -x
# pad zeros to data and syn and cut accordingly
../scripts/process_data_and_syn.pl -d data,sac -s syn, -1 -b -26 -x

# filter at different frequency bands
for freq in 2 3 6; do
  ../scripts/process-data_and_syn.pl -d data,sac -s syn, -2 -t ${freq}/30 -x
done
