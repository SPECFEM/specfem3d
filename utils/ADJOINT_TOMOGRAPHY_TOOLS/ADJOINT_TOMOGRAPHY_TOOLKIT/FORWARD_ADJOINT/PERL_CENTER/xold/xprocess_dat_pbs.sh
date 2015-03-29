#!/bin/sh
#PBS -q tromp
#PBS -N XPROC_DAT_201006030432A
#PBS -l nodes=1:ppn=1
#PBS -l walltime=5:00:00
#PBS -j oe
#PBS -k oe
#PBS -o job_src2.log


echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

opt=3
tstart=0
tend=4800
fmin=80
fmax=150
ext=T080_150
datdir=DATASET_GLOBE/CMTSOLUTION_201006030432A
syndir=SYN_S00/CMTSOLUTION_201006030432A
cmtfile=CMTSOLUTION_CENTER/CMTSOLUTION_201006030432A

if [ $opt == 1 ]; then
  ./PERL_SRC/process_data.pl -m $cmtfile -s 1.0 -l $tstart/$tend -t $fmin/$fmax -f -i -p -x $ext $datdir/*.sac
  ./PERL_SRC/rotate.pl -l $tstart -L $tend $datdir/*.BHE.sac.$ext
elif [ $opt == 2 ]; then
  ./PERL_SRC/process_syn.pl -S -m $cmtfile -s 1.0 -l $tstart/$tend -t $fmin/$fmax -f -x $ext $syndir/*.sac
  ./PERL_SRC/rotate.pl -l $tstart -L $tend $syndir/*.LHE.sem.sac.$ext
elif [ $opt == 3 ]; then
  ./PERL_SRC/cut_data_syn.pl -d $syndir $datdir/*BH[ZRT]*.$ext
fi

echo PROCESSING DATA DONE
