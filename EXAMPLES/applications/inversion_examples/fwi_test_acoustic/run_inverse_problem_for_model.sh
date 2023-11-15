#!/bin/bash
#
#
# script to run inversion framework
#
#############################################
# USER PARAMETERS

# simulation mode
mode=$1

#############################################

if [ "$mode" == "" ]; then echo "usage: ./run_inverse_problem_for_model.sh mode[==forward or l-bfgs] (type==exploration or teleseismic, optional)"; exit 1; fi
echo
echo "run_inverse_problem_for_model: `date`"
echo
currentdir=`pwd`

# get the number of processors, ignoring comments in the Par_file
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
NSIM=`grep ^NUMBER_OF_SIMULTANEOUS_RUNS DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

# total number of processes to run with simultaneous solver
NPROC_TOT=$(($NPROC * $NSIM))

echo ""
echo "simulation parameters: using mode $mode"
echo "                       NPROC         = $NPROC"
echo "                       NSIMULTANEOUS = $NSIM"
echo
echo "total number of MPI processes for run: $NPROC_TOT"
echo

# runs inversion
if [ "$NPROC_TOT" -eq 1 ]; then
  # This is a serial simulation
  ./bin/xinverse_problem_for_model $mode
else
  # This is a MPI simulation
  mpirun -np $NPROC_TOT ./bin/xinverse_problem_for_model $mode
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "run done: `date`"
echo




