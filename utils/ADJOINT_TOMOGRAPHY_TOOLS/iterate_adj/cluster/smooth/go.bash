#!/bin/bash -v
#BSUB -o OUTPUT_FILES/%J.o
#BSUB -a openmpi
#BSUB -J smooth_model

echo "$LSB_MCPU_HOSTS" > OUTPUT_FILES/lsf_machines
echo "$LSB_JOBID" > OUTPUT_FILES/jobid
remap_lsf_machines.pl OUTPUT_FILES/lsf_machines >OUTPUT_FILES/machines

# number of slices 14x12, 0 stands for basin code (1 would be 1-chunk global,
# and 6 would be 6-chunk global)
# element size on the surface would be 2 km,
# the smoothing 'radius' sigma = 5 km

# the neigboring points within 3*sigma+element_size are used for the smoothing
# the 3D Gaussian function is defined as
#   G(x,y,z) = 1/(sqrt(2*pi)*sigma)**3 * exp[-r**2/(2*sigma**2)]
# which is parallel to the 1D Gaussian function

# the last two entries are the directories for (1) the unsmoothed and smoothed files and (2) topology

#--------------------------------------------------

#mpirun.lsf --gm-no-shmem --gm-copy-env smooth_sem_fun 14 12 0 2 6 1 mu_kernel inout_smooth topo
#sleep 5s

#mpirun.lsf --gm-no-shmem --gm-copy-env smooth_sem_fun 14 12 0 2 6 1 kappa_kernel inout_smooth topo
#sleep 5s

#xcombine_vol_data slice_file mu_kernel_smooth topo inout_smooth . 0
#sleep 5s
#mv mu_kernel_smooth.mesh mu_kernel_smooth_h006km_v001km.mesh
#sleep 5s

#xcombine_vol_data slice_file kappa_kernel_smooth topo inout_smooth . 0
#sleep 5s
#mv kappa_kernel_smooth.mesh kappa_kernel_smooth_h006km_v001km.mesh
#sleep 5s

#--------------------------------------------

# new command for Harvard cluster
#mpirun smooth_sem_fun 14 12 0 2 12 2 vs_m16 inout_smooth topo
#sleep 5s
#xcombine_vol_data slice_file vs_m16_smooth topo inout_smooth . 0
#sleep 5s
#mv vs_m16_smooth.mesh vs_m16_smooth_h012km_v002km.mesh
#sleep 5s

# new command for Harvard cluster
#mpirun smooth_sem_fun 14 12 0 2 12 2 vp_m16 inout_smooth topo
#sleep 5s
#xcombine_vol_data slice_file vp_m16_smooth topo inout_smooth . 0
#sleep 5s
#mv vp_m16_smooth.mesh vp_m16_smooth_h012km_v002km.mesh
#sleep 5s

#--------------------------------------------

# new command for Harvard cluster
mpirun smooth_sem_fun 14 12 0 2 4 1 rho_cbr_kernel inout_smooth topo
sleep 5s
xcombine_vol_data slice_file rho_cbr_kernel_smooth topo inout_smooth . 0
sleep 5s
mv rho_cbr_kernel_smooth.mesh rho_cbr_kernel_smooth_h004km_v001km.mesh
sleep 5s

