#!/bin/bash
testdir=`pwd`
me=`basename "$0"`

# sets source directory
cd $ROOT/
srcdir=`pwd`

cd $testdir/

# title
echo >> $testdir/results.log
echo "$me in: $testdir" >> $testdir/results.log
echo >> $testdir/results.log

echo "directory: `pwd`" >> $testdir/results.log

# default configuration
#$srcdir/configure >> $testdir/results.log 2>&1

# checks if all supplement executables exist
echo "checking if supplement binaries exist" >> $testdir/results.log
exec=( xcombine_surf_data \
       xcombine_vol_data \
       xcombine_vol_data_vtk \
       xconvolve_source_timefunction \
       xcreate_movie_shakemap_AVS_DX_GMT \
       xcheck_mesh_quality \
       xconvert_skewness_to_angle \
       xmodel_update \
       xsum_kernels \
       xsum_preconditioned_kernels \
      )

for var in ${exec[@]};
do
  # single compilation
  echo "compilation: $var" >> $testdir/results.log
  make clean >> $testdir/results.log 2>&1
  make -j 4 $var >> $testdir/results.log 2>&1

  # checks exit code
  if [[ $? -ne 0 ]]; then
    echo >> $testdir/results.log
    echo "compilation failed, please check..." >> $testdir/results.log
    exit 1
  fi

  echo "" >> $testdir/results.log
  # check
  if [ ! -e bin/$var ]; then
    echo "compilation of $var failed, please check..." >> $testdir/results.log
    exit 1
  else
    echo "binary exists: $var" >> $testdir/results.log
  fi
done
echo "" >> $testdir/results.log

echo "successful compilation" >> $testdir/results.log

