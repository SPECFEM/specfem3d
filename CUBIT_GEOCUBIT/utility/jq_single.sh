#!/bin/sh -f
#PBS -N GEOCUBIT
#PBS -o joba.log
#PBS -l walltime=2:0:0
#PBS -m abe -M emanuele.casarotti@ingv.it
#PBS -j oe
cd /home/casarotti/geocubit4
python GEOCUBIT.py --mesh --build_volume --cfg='/home/casarotti/geocubit4/example/abruzzo_single.cfg' --id_proc=${PBS_ARRAY_INDEX} > ${PBS_ARRAY_INDEX}.log
