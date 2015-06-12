#!/bin/bash 
#OAR -l /nodes=2,walltime=3:00:00
#OAR -n simu_2_nodes
#OAR -p gpu='NO' AND ibpool='FDR'

source /softs/env-intel15-impi50.sh

echo ${OAR_NODEFILE}

nb_cpu=$(grep "physical id" /proc/cpuinfo | sort -u | wc -l)


ulimit -s unlimited

mpiexec.hydra -f $OAR_FILE_NODES  -bootstrap-exec oarsh -perhost 20 ./axisem


