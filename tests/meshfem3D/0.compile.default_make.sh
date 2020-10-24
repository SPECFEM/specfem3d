#!/bin/bash
testdir=`pwd`
me=`basename "$0"`

#checks if ROOT valid
if [ -z "${ROOT}" ]; then export ROOT=../../ ; fi

# sets source directory
cd $ROOT/
srcdir=`pwd`

cd $testdir/

# title
echo >> $testdir/results.log
echo "$me in: $testdir" >> $testdir/results.log
echo >> $testdir/results.log

echo "directory: `pwd`" >> $testdir/results.log

#cleanup
rm -rf config.log config.status
rm -rf ./bin ./obj ./setup ./OUTPUT_FILES

# default configuration
$srcdir/configure >> $testdir/results.log 2>&1

# executable
var=xmeshfem3D

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

#cleanup
rm -rf ./bin/*

echo "successful compilation" >> $testdir/results.log

