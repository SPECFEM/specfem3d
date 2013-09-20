#!/bin/sh -f
#PBS -N GEOCUBIT
#PBS -o run.log
#PBS -l select=1:ncpus=1
#PBS -l place=free
#PBS -l walltime=46:0:0
#PBS -m abe -M emanuele.casarotti@ingv.it
#PBS -j oe


cd $PBS_O_WORKDIR
python2.5 ./GEOCUBIT.py --collect --merge --meshfiles='/lscratch/users/casarotti/cal/SCAL_LR/mesh*.e' --cpux=15 --cpuy=15
