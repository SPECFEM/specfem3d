source ./global_parameters.in
WORKING_DIRECTORY=$(pwd)
declare -i isrc nsrc
nsrc=$nb_eqks
source $SHELL_SCRIPTS/functions_general.sh
clean_directories
