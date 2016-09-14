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
echo "compilation: $testdir" >> $testdir/results.log
make clean >> $testdir/results.log 2>&1
make -j 4 tomography >> $testdir/results.log 2>&1
# checks exit code
if [[ $? -ne 0 ]]; then
  echo >> $testdir/results.log
  echo "compilation failed, please check..." >> $testdir/results.log
  exit 1
fi

echo "" >> $testdir/results.log

# checks executable availability
exec=( xadd_model_iso \
       xmodel_update \
       xsum_preconditioned_kernels \
      )

for var in ${exec[@]};
do
  # check
  if [ ! -e bin/$var ]; then
    echo "binary does not exist! $var " >> $testdir/results.log
    echo "compilation of $testdir failed, please check..." >> $testdir/results.log
    exit 1
  else
    echo "binary exists: $var " >> $testdir/results.log
  fi
done
echo "" >> $testdir/results.log

#cleanup
rm -rf ./bin/*

echo "successful compilation" >> $testdir/results.log

