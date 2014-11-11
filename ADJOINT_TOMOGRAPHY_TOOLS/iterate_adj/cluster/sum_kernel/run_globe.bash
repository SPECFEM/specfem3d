#!/bin/bash

###########################################################
# USER PARAMETERS

# set this to your directory which contains all event kernel directories
outputs_kernel_directory=my_kernels

###########################################################

# clean up
mkdir -p OUTPUT_FILES
rm -f OUTPUT_FILES/*
rm -f OUTPUT_SUM/*

# update kernel links
cd INPUT_KERNELS/
rm -rf ./*.*
ln -s ${outputs_kernel_directory}/*.* ./
cd ..

# update input kernel list
ls -1 INPUT_KERNELS/ > kernels_run.globe

# compiles
#cd src/
#make -f Makefile_globe
#cd ..

# runs job
date
qsub go_globe_pbs.bash

