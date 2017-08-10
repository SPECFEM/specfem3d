#!/bin/bash
#MSUB -r run_inversion
#MSUB -n 128
#MSUB -x
#MSUB -T 84000
#MSUB -q standard
#MSUB -o run_inversion.o
#MSUB -e run_inversion.e
#
#
set -x
cd $BRIDGE_MSUB_PWD

#
## Chargement des modules module load ...

# load parameters
source ./global_parameters.in
# load scripts
source $SHELL_SCRIPTS/functions_set_up.sh
source $SHELL_SCRIPTS/functions_general.sh
source $SHELL_SCRIPTS/functions_optimisation.sh
source $SHELL_SCRIPTS/functions_simu_mpi.sh
source $SHELL_SCRIPTS/functions_inversion.sh

#####-------------------------------------------

declare -i iter niter isrc nsrc
nsrc=$nb_eqks

iter=0;
niter=50;
long_liss=5.
long_lissv=0.1
sigma_beta=1.
sigma_gamma=1.

flog_file=log.inversion

echo $long_liss > longueur_lissage.txt
echo $long_lissv > longueur_lissage_vertical.txt
echo $sigma_beta > sigma_beta.txt
echo $sigma_gamma > sigma_gamma.txt

initialise_process
echo $(date) > $flog_file

#
while [ "$iter" -le "$niter" ]; do

# iteration complete
iteration_fwi

iter="$iter+1"
done



