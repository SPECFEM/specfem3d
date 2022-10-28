#!/usr/bin/env bash
#OAR -n axisem
#OAR -l /nodes=6/cpu=2/core=10,walltime=05:00:00
#OAR -p ibpool='FDR'
#OAR -p gpu='NO'

source /softs/env_default.sh
#export OMP_NUM_THREAD=20

#echo $OAR_FILE_NODES
#cat $OAR_FILE_NODES > lalala
#uniq lalala > hostsfi

#time mpiexec.hydra  \
#-machinefile $OAR_FILE_NODES \
#-hostfile hostsfi \
#-bootstrap ssh \
#-bootstrap-exec /usr/bin/oarsh \
#-envall  \
#-perhost 1 ./lithos_fwi_time.x
#time mpiexec.hydra  \
#-machinefile $OAR_FILE_NODES \
#-hostfile hostsfi \
#-bootstrap ssh \
#"-bootstrap-exec /usr/bin/oarsh \
#-envall  \
#-perhost 1 -np 4 ./lithos_fwi_time.x
ulimit -s unlimited

time mpiexec.hydra  \
-machinefile $OAR_FILE_NODES \
-bootstrap ssh \
-bootstrap-exec /usr/bin/oarsh \
-envall   ./axisem > OUTPUT_AxiSEMi


