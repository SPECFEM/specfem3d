#!/bin/bash

# cleans outputs
mkdir -p OUTPUT_FILES
rm -f OUTPUT_FILES/*
rm -f OUTPUT_MODEL/*.*

# compiles executable
#cd src/
#make -f Makefile_globe
#cd ..

# submits job
date
qsub go_globe_pbs.bash

