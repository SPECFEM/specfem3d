#!/bin/sh
# Author: Hejun Zhu, hejunzhu@princeton.edu
# Princeton University, New Jersey, USA
# Last modified: Mon Mar  7 15:56:47 EST 2011

iter=M00

model=MODEL_$iter
dest=../FORWARD_ADJOINT/XSRC_SEM/DATA/GLL/

if [ ! -d $model ]; then
  echo WRONG! NO $model
  exit
fi
if [ ! -d $dest ]; then
  echo WRONG! NO $dest
  exit
fi

echo copy vpv from $model ...
cp $model/proc*_reg1_vpv.bin  $dest
echo copy vph from $model ...
cp $model/proc*_reg1_vph.bin  $dest
echo copy vsv from $model ...
cp $model/proc*_reg1_vsv.bin  $dest
echo copy vsh from $model ...
cp $model/proc*_reg1_vsh.bin  $dest
echo copy rho from $model ...
cp $model/proc*_reg1_rho.bin  $dest
echo copy eta from $model ...
cp $model/proc*_reg1_eta.bin  $dest
