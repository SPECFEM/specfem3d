# add data and syn headers
../scripts/process_data_and_syn.pl -d data,sac  -s syn, -0 -m CMTSOLUTION -a STATIONS -x
# pad zeros to data and syn and cut accordingly
../scripts/process_data_and_syn.pl -d data,sac -s syn, -1 -x -b -26

# filter at different frequency bands
for freq in 2 3 6; do
  ../scripts/process_data_and_syn.pl -d data,sac -s syn, -2 -t ${freq}/30 -x
done
