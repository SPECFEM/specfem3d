#!/bin/bash

# get the number of time steps, ignoring comments in the Par_file
NSTEP=`grep ^NSTEP DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
DT=`grep ^DT DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`


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

mkdir -p NOISE_TOMOGRAPHY/

# generates S_squared file
if [ 0 == 1 ]; then
  ## matlab
  matlab -nosplash -nodisplay -r "NOISE_TOMOGRAPHY($noise_nstep,$DT,1.,15.,\'$noise_model\');exit"
else
  ## python & numpy
  ./NOISE_TOMOGRAPHY.py $noise_nstep $DT 1.0 15.0 $noise_model
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# move to noise setup folder
mv -v S_squared NOISE_TOMOGRAPHY/

