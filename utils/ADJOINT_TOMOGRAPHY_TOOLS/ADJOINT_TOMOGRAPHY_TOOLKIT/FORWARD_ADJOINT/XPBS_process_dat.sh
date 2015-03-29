#!/bin/sh
#PBS -q tromp
#PBS -N XPROC_DAT_200604121652A
#PBS -l nodes=1:ppn=1
#PBS -l walltime=10:00:00
#PBS -j oe
#PBS -k oe
#PBS -o job_src2.log


echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

opt=2
tstart=0
tend=1800
fmin=25
fmax=150
ext=T025_150
datdir=DATASET_EUROPE/CMTSOLUTION_200604121652A
syndir=SYN_M13/CMTSOLUTION_200604121652A
cmtfile=CMTSOLUTION_CENTER/CMTSOLUTION_200604121652A

if [ $opt == 1 ]; then
  ./PERL_CENTER/process_data.pl -m $cmtfile -s 10.0 -l $tstart/$tend -t $fmin/$fmax -f -i -p -x $ext $datdir/*.sac
  ./PERL_CENTER/rotate.pl -l $tstart -L $tend $datdir/*.BHE.sac.$ext
  echo delete N E bandpass component
  cd $datdir
  rm *.BHN.sac.$ext
  rm *.BHE.sac.$ext

elif [ $opt == 2 ]; then
  ./PERL_CENTER/process_syn.pl -S -m $cmtfile -s 10.0 -l $tstart/$tend -t $fmin/$fmax -f -x $ext $syndir/*.sac
  ./PERL_CENTER/rotate.pl -l $tstart -L $tend $syndir/*.LHE.sem.sac.$ext
  ./PERL_CENTER/cut_data_syn.pl -d $syndir $datdir/*BH[ZRT]*.$ext
  echo delete N E bandpass component
  cd $syndir
  rm *.LHN.sem.sac.$ext
  rm *.LHE.sem.sac.$ext
fi

echo PROCESSING DATA successfully
