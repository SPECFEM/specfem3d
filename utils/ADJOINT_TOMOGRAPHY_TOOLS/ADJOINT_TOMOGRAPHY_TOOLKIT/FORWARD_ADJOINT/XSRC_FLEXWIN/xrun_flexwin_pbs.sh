#!/bin/sh
#PBS -q tromp
#PBS -N XFLEXWIN_200604121652A
#PBS -l nodes=1:ppn=1
#PBS -l walltime=5:00:00
#PBS -j oe
#PBS -k oe
#PBS -o job_src2.log

echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR


input=XFLEXWIN_INPUT_M13/INPUT_FLEXWIN_200604121652A_T025_150
output=XFLEXWIN_OUTPUT_M13/OUTPUT_FLEXWIN_200604121652A_T025_150
syndir=../SYN_M13
dir=CMTSOLUTION_200604121652A_WIN_T025_150

./flexwin < $input > $output
cd $syndir
if [ -f $dir.tar.gz ]; then
  echo rm $dir.tar.gz
  rm $dir.tar.gz
fi
tar -czvf $dir.tar.gz $dir
rm -rf $dir
echo flexwin done successfully
