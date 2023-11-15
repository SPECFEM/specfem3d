#!/bin/bash
#
# creates an external source time function file
#

# based on a copy from OUTPUT_FILES/plot_source_time_function.txt
if [ ! -e OUTPUT_FILES/plot_source_time_function.txt ]; then
  echo "no file OUTPUT_FILES/plot_source_time_function.txt, please run forward simulation first with the flag USE_EXTERNAL_SOURCE_FILE = .false."
  exit 1
fi

stf_file=stf.dat

echo "creating external source time function file..."
echo
echo "based on OUTPUT_FILES/plot_source_time_function.txt"
cp -v OUTPUT_FILES/plot_source_time_function.txt tmp.txt

nstep=`wc -l tmp.txt | awk '{print $1}'`
echo
echo "number of STF time steps: $nstep"
echo

# format for external:

# val_1
# val_2
# ..


# gets DT from Par_file
DT=`grep ^DT DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

echo "valid for simulation DT = $DT"

# comment lines
echo "# Source time function" > ${stf_file}
echo "#   simulation DT       : $DT" >> ${stf_file}
echo "#   number of time steps: $nstep" >> ${stf_file}
# data values
echo "# STF values" >> ${stf_file}
awk '{print $2}' tmp.txt >> ${stf_file}

mv -v ${stf_file} DATA/${stf_file}
rm -f tmp.txt

echo
echo "done"
echo "written to: DATA/${stf_file}"
echo


