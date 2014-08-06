#!/bin/bash

for eid in `ls -1 -d [0-9]*`; do

 echo "============ $eid ============"
 ls -1 $eid/inout_smooth/*mu_kernel_smooth.bin | wc
 ls -1 $eid/inout_smooth/*kappa_kernel_smooth.bin | wc

done
