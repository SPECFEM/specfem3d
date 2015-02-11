#!/bin/bash

# compiles executable in src/ directory
#current_pwd=$PWD
#cd src/
#make -f Makefile_globe
#cd ../

# cleans up outputs
rm -f OUTPUT_FILES/*

# runs job
date

# creates job scripts from template:
template=go_globe_pbs.bash

for tag in "bulk_c" "bulk_beta" "rho"
do

  echo "kernel: $tag"

  # prepares job script from template
  sed -e "s:bulk_c:$tag:" $template > go_globe_pbs.${tag}.bash

  # smooths: bulk_c, bulk_beta and rho kernels
  qsub go_globe_pbs.${tag}.bash
  sleep 2

done


