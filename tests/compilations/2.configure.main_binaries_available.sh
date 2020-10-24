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

# checks if all main executables exist
echo "checking if main binaries exist" >> $testdir/results.log
exec=( xdecompose_mesh \
       xmeshfem3D \
       xgenerate_databases \
       xspecfem3D \
      )

for var in ${exec[@]};
do
  # single compilation
  #echo "compilation: $var" >> $testdir/results.log
  #make clean >> $testdir/results.log 2>&1
  #make -j 4 $var >> $testdir/results.log 2>&1

  # check
  if [ ! -e bin/$var ]; then
    echo "compilation of $var failed, please check..." >> $testdir/results.log
    exit 1
  else
    echo "binary exists: $var" >> $testdir/results.log
  fi
done
echo "" >> $testdir/results.log

cd $testdir/
echo "successful compilation" >> $testdir/results.log

