#!/bin/csh
#
# script to output travel times 
# based on iaspei-tau package
# intended to work for P and S phases
#
# usage: ./phtimes.csh depth(km) distance(degrees) phase
#
#  e.g.  ./phtimes.csh 611 52.474 P
#        ./phtimes.csh 611 52.474 S
# 
###########################################
# MODEL PARAMETER 

# locates your IASPEI-TAU installation
set dir="/opt/seismo-util/source/iaspei-tau"

# uses IASP91 tables
set model="$dir/tables/iasp91"

###########################################

# input arguments
set depth="$1"
set delta="$2"
set branch="$3"

# checks arguments
if( $depth =~ "" ) then
  echo "Usage: phtimes.csh depth(km) distance(degrees) phase"
  exit 1
endif
if( $delta =~ "" ) then
  echo "Usage: phtimes.csh depth(km) distance(degrees) phase"
  exit 1
endif
if( $branch =~ "" ) then
  echo "Usage: phtimes.csh depth(km) distance(degrees) phase"
  exit 1
endif


# travel times
$dir/bin/ttimes <<EOF  > tmp.log
Y
$model
$branch

N
$depth
$delta
-1
-1
EOF

# outputs single line with phase name and traveltime 
fgrep -A 1 "#" tmp.log | tail -n 1 | awk '{print $3,$4}'


# in case of problems, see the created tmp.log for further informations