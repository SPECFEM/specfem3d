#!/bin/bash
testdir=`pwd`

#checks if ROOT valid
if [ -z "${ROOT}" ]; then export ROOT=../../ ; fi

# sets source directory
cd $ROOT/
srcdir=`pwd`

cd $testdir/

# title
echo >> $testdir/results.log
echo "compile default make" >> $testdir/results.log
echo >> $testdir/results.log

echo "directory: `pwd`" >> $testdir/results.log

#cleanup
rm -rf config.log config.status
rm -rf ./bin ./obj ./setup ./OUTPUT_FILES

# default configuration
$srcdir/configure >> $testdir/results.log 2>&1

# single compilation
echo "compilation: auxiliaries" >> $testdir/results.log
make clean >> $testdir/results.log 2>&1

make -j 4 aux >> $testdir/results.log 2>&1

# check
if [ ! -e bin/xsum_kernels ]; then
  echo "compilation of auxiliaries failed, please check..." >> $testdir/results.log
  exit 1
else
  echo "binary exists: xsum_kernels" >> $testdir/results.log
fi

#cleanup
rm -rf ./bin/* 

echo "successful compilation" >> $testdir/results.log

