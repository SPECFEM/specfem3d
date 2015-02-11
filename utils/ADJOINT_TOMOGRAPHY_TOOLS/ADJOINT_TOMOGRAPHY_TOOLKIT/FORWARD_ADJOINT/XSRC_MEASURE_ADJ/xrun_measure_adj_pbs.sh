#!/bin/sh
#PBS -q tromp
#PBS -N XMEASURE_201006030432A
#PBS -l nodes=1:ppn=1
#PBS -l walltime=10:00:00
#PBS -j oe
#PBS -k oe
#PBS -o job_src2.log

echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

# input and output for mt_measure_adj
inputfile=XMEASUREMENT_INPUT_M18/MEASUREMENT_201006030432A_T015_040
output=XOUTPUT_TAG_M18/OUTPUT_201006030432A_T015_040
syndir=../SYN_M18
mtdir=CMTSOLUTION_201006030432A_MT_T015_040

echo run mt_measure_adj...
./mt_measure_adj < $inputfile > $output

echo delete useless files
cd $syndir"/"$mtdir
rm *.err_*

#cd $syndir
#if [ -f $mtdir.tar.gz ]; then
# echo rm $mtdir.tar.gz
# rm $mtdir.tar.gz
#fi
#tar -czvf $mtdir.tar.gz $mtdir
#rm $mtdir/*
echo measurement done successfully

