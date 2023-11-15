#!/bin/bash
#
#
#  CONVERSION OF ALL SNAPSHOTS COMPUTED BY SOLVER IN VTK FORMAT
#
#

NTSTEP_BETWEEN_FRAMES=`grep ^NTSTEP_BETWEEN_FRAMES DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
NSTEP=`grep ^NSTEP DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

for IT in `seq $NTSTEP_BETWEEN_FRAMES $NTSTEP_BETWEEN_FRAMES $NSTEP`;
do
./create_one_snapshot.sh $IT
done
