#!/bin/bash -v
#BSUB -o OUTPUT_FILES/%J.o
#BSUB -a mpich_gm
#BSUB -J smooth_model_update

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

mpirun.lsf --gm-no-shmem --gm-copy-env smooth_sem_fun 14 12 0 2 10 dm_mu_kernel_smooth_p068 inout_smooth topo
sleep 5s
xcombine_vol_data slice_file dm_mu_kernel_smooth_p068_smooth topo inout_smooth . 0
sleep 5s
mv dm_mu_kernel_smooth_p068_smooth.mesh dm_mu_kernel_smooth_p068_smooth_10km.mesh
sleep 5s

mpirun.lsf --gm-no-shmem --gm-copy-env smooth_sem_fun 14 12 0 2 10 dm_kappa_kernel_smooth_p068 inout_smooth topo
sleep 5s
xcombine_vol_data slice_file dm_kappa_kernel_smooth_p068_smooth topo inout_smooth . 0
sleep 5s
mv dm_kappa_kernel_smooth_p068_smooth.mesh dm_kappa_kernel_smooth_p068_smooth_10km.mesh
sleep 5s
