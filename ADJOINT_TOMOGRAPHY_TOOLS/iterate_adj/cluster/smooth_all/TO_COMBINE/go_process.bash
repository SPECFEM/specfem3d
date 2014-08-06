#!/bin/bash -v

dirbase=$PWD

# modified in the setup script
eid=9114763

current_pwd=$dirbase/$eid

cd $current_pwd

$current_pwd/xcombine_vol_data slice_file mu_kernel_smooth topo inout_smooth . 0
sleep 5s
mv mu_kernel_smooth.mesh mu_kernel_smooth_h006km_v001km.mesh
sleep 5s

$current_pwd/xcombine_vol_data slice_file kappa_kernel_smooth topo inout_smooth . 0
sleep 5s
mv kappa_kernel_smooth.mesh kappa_kernel_smooth_h006km_v001km.mesh
sleep 5s
