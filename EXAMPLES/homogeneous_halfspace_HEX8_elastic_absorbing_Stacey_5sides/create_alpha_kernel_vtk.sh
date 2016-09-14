#!/bin/bash

mkdir -p bin.xcombine

# checks for executable
if [ ! -e bin.xcombine/xcombine_vol_data_vtk  ]; then echo "please make xcombine_vol_data_vtk and copy executable to bin.xcombine"; exit 1; fi

# checks for kernel files
if [ ! -e OUTPUT_FILES/DATABASES_MPI/proc000000_alpha_kernel.bin ]; then echo "proc***alpha_kernel.bin missing in OUTPUT_FILES/DATABASES_MPI/"; exit 1; fi


# creates kernel as vtk-file
cd bin.xcombine/
./xcombine_vol_data_vtk 0 3 alpha_kernel ../OUTPUT_FILES/DATABASES_MPI/ ../OUTPUT_FILES/ 1
cd ../

echo
echo

