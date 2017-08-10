#!/usr/bin/env bash
#OAR -n axisemmesh
#OAR -l /nodes=1/cpu=2/core=10,walltime=00:30:00
#OAR -p ibpool='FDR'
#OAR -p gpu='NO'

source /softs/env_default.sh


ulimit -s unlimited

export OMP_NUM_THREAD=20


#echo $OAR_FILE_NODES
#cat $OAR_FILE_NODES > lalala
#uniq lalala > hostsfi

#make clean
#make
#time mpiexec.hydra  \
#-machinefile $OAR_FILE_NODES \
#-hostfile hostsfi \
#-bootstrap ssh \
#-bootstrap-exec /usr/bin/oarsh \
#-envall  \
#-perhost 1 ./lithos_fwi_time.x

./xmesh > OUTPUT


#time mpiexec.hydra  \
#-machinefile $OAR_FILE_NODES \
#-hostfile hostsfi \
#-bootstrap ssh \
#-bootstrap-exec /usr/bin/oarsh \
#-envall  \
#-perhost 1 -np 4 ./lithos_fwi_time.x

