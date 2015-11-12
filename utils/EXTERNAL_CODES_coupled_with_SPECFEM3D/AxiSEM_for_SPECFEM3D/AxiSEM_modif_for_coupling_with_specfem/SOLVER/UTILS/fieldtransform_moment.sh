#!/bin/bash
set -e
cd MZZ
echo "Transforming MZZ"
../xfield_transform
cd ..
cd MXX_P_MYY
echo "Transforming MXX_P_MYY"
../xfield_transform
cd ..
cd MXZ_MYZ
echo "Transforming MXZ_MYZ"
../xfield_transform
cd ..
cd MXY_MXX_M_MYY
echo "Transforming MXY_MXX_M_MYY"
../xfield_transform
cd ..
