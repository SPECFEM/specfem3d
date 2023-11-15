#!/bin/bash

# get the number of time steps, ignoring comments in the Par_file
NSTEP=`grep ^NSTEP DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
DT=`grep ^DT DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
TMIN=12.0
TMAX=60.0


noise_nstep=$((2*NSTEP - 1))
noise_model=NLNM

echo "simulation setup:"
echo "  NSTEP = $NSTEP"
echo "  DT    = $DT"
echo
echo "noise:"
echo "  number of steps = $noise_nstep"
echo "  model           = $noise_model"
echo

# generates S_squared file
./NOISE_TOMOGRAPHY.py $noise_nstep $DT $TMIN $TMAX $noise_model

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# move to noise setup folder
mv -v S_squared NOISE_TOMOGRAPHY/
