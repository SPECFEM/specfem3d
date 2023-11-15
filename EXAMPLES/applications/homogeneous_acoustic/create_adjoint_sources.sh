#!/bin/bash
#################################################

# reference adjoint station
network="DB"
station="X20"
comp="MXP"
en="semp"

# window start/end time
t_start=$1  # 9.0  # 32.0
t_end=$2    # 26.0   # 30.0   # 38.0

#################################################

# set defaults
if [ "$t_start" == "" ]; then t_start=9.0; fi
if [ "$t_end" == "" ]; then t_end=26.0; fi

# adjoint sources will be in folder SEM/
currentdir=`pwd`
mkdir -p SEM

# needs traces
sta=$network.$station
if [ ! -e OUTPUT_FILES/$sta.$comp.$en ]; then echo "please make sure trace OUTPUT_FILES/$sta.$comp.$en is available"; exit 1; fi

# replaces pressure component with X/Y/Z
stap=$sta.$comp.$en
compx=${comp/P/X}; stax=$network.$station.$compx.$en
compy=${comp/P/Y}; stay=$network.$station.$compy.$en
compz=${comp/P/Z}; staz=$network.$station.$compz.$en
rm -f SEM/$network.$station*
cp -v OUTPUT_FILES/$stap SEM/$stax
cp -v OUTPUT_FILES/$stap SEM/$stay
cp -v OUTPUT_FILES/$stap SEM/$staz
echo

# compile adjoint_source tool
if [ ! -e xcreate_adjsrc_traveltime ]; then
  # creates adjoint sources
  cd ../../../utils/adjoint_sources/traveltime

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

./xcreate_adjsrc_traveltime $t_start $t_end 1 $sta.*
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

if [ ! -e $sta.$compx.adj ]; then echo "error creating adjoint sources, please check..."; exit 1; fi
echo

# renames
# (acoustic adjoint sources are read from component 1 -> MXX.adj)
#mv -v $sta.$compx.$en.adj $sta.$compx.adj
#mv -v $sta.$compy.$en.adj $sta.$compy.adj
#mv -v $sta.$compz.$en.adj $sta.$compz.adj

# create STATIONS_ADJOINT file with adjoint source location
fgrep $station ../DATA/STATIONS > ./STATIONS_ADJOINT
cp -v ./STATIONS_ADJOINT ../DATA/

cd ../

echo
echo

