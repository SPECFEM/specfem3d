#!/bin/bash
#################################################

# reference adjoint station
network="XX"
station="STA0008"
comp="BXZ"
en="semd"

# window start/end time
t_start=$1   # e.g. 10.0
t_end=$2     #      34.0

#################################################

# usage: ./create_adjoint_sources.sh t_start t_end

if [ "$t_start" == "" ] || [ "$t_end" == "" ]; then echo "usage: ./create_adjoint_sources.sh t_start t_end"; exit 1; fi

echo
echo "adjoint sources:"
echo "  window start/end = $t_start / $t_end"
echo

# adjoint sources will be in folder SEM/
currentdir=`pwd`
mkdir -p SEM

# needs traces
sta=$network.$station
if [ ! -e OUTPUT_FILES/$sta.$comp.$en ]; then echo "please make sure trace OUTPUT_FILES/$sta.$comp.$en is available"; exit 1; fi

rm -f SEM/$sta.*
cp -v OUTPUT_FILES/$sta.* SEM/

# compile adjoint_source tool
if [ ! -e xcreate_adjsrc_traveltime ]; then
  # creates adjoint sources
  cd ../../utils/adjoint_sources/traveltime

  # fortran compiler (as specified in Makefile)
  FC=`grep '^FC .*' ../../../Makefile | cut -d = -f 2 | sed "s/^[ \t]*//"`
  if [ "$FC" == "" ]; then echo "fortran compiler not found, exiting..."; exit 1; fi
  CC=`grep '^CC .*' ../../../Makefile | cut -d = -f 2 | sed "s/^[ \t]*//"`
  if [ "$CC" == "" ]; then echo "C compiler not found, exiting..."; exit 1; fi

  echo "compiling xcreate_adjsrc_traveltime:"
  echo "  using fortran compiler = $FC"
  echo "  using C compiler       = $CC"
  echo

  cp Makefile Makefile.host
  sed -i "s:F90 .*:F90 = $FC:" Makefile.host
  sed -i "s:CC .*:CC = $CC:" Makefile.host

  rm -rf xcreate_adjsrc_traveltime
  make -f Makefile.host
  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi

  cp -v xcreate_adjsrc_traveltime $currentdir/SEM/
  cd $currentdir
fi
if [ ! -e SEM/xcreate_adjsrc_traveltime ]; then echo "please make xcreate_adjsrc_traveltime and copy to SEM/"; exit 1; fi


echo
echo "running adjoint source creation"
echo
# creates adjoint sources
cd SEM/

./xcreate_adjsrc_traveltime $t_start $t_end 3 $sta.*
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

if [ ! -e $sta.$comp.adj ]; then echo "error creating adjoint sources, please check..."; exit 1; fi
echo

# renames
#rename .semd.adj .adj *semd.adj

# create STATIONS_ADJOINT file with adjoint source location
fgrep $station ../DATA/STATIONS > ./STATIONS_ADJOINT
cp -v ./STATIONS_ADJOINT ../DATA/

cd ../

echo
echo

