#!/bin/bash
set -e
cd PX
echo "Transforming PX"
../xfield_transform
cd ..
cd PZ
echo "Transforming PZ"
../xfield_transform
cd ..
