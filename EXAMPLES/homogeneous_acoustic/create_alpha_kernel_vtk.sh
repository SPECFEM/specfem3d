#!/bin/bash

# checks for executable
if [ ! -e bin/xcombine_vol_data_vtk  ]; then echo "please make xcombine_vol_data_vtk and copy executable to bin/"; exit 1; fi

# checks for kernel files
if [ ! -e OUTPUT_FILES/DATABASES_MPI/proc000000_alpha_acoustic_kernel.bin ]; then echo "proc***alpha_acoustic_kernel.bin missing in OUTPUT_FILES/DATABASES_MPI/"; exit 1; fi


# creates kernel as vtk-file
./bin/xcombine_vol_data_vtk 0 3 alpha_acoustic_kernel OUTPUT_FILES/DATABASES_MPI/ OUTPUT_FILES/ 1

echo
echo

