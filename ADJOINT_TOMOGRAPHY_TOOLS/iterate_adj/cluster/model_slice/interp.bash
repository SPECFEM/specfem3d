#!/bin/bash -v
#BSUB -o OUTPUT_FILES/%J.o
#BSUB -a mpich_gm
#BSUB -J Model_slice_1

echo "$LSB_MCPU_HOSTS" > OUTPUT_FILES/lsf_machines
echo "$LSB_JOBID" > OUTPUT_FILES/jobid
remap_lsf_machines.pl OUTPUT_FILES/lsf_machines >OUTPUT_FILES/machines

for ftag in `cat ftags`; do

#--------------------
# vertical cross sections

mpirun.lsf --gm-no-shmem --gm-copy-env sem_model_slice vert_xc_${ftag}_input /ibrixfs1/home/carltape/ADJOINT_TOMO_EXTRA/topo /ibrixfs1/home/carltape/ADJOINT_TOMO_EXTRA/KERNELS_MODELS/models/m16 vs_m16 vert_xc_vs_m16_${ftag}.gmt
sleep 5s

mpirun.lsf --gm-no-shmem --gm-copy-env sem_model_slice vert_xc_${ftag}_input /ibrixfs1/home/carltape/ADJOINT_TOMO_EXTRA/topo /ibrixfs1/home/carltape/ADJOINT_TOMO_EXTRA/KERNELS_MODELS/models/m00 vs_m00 vert_xc_vs_m00_${ftag}.gmt
sleep 5s

mpirun.lsf --gm-no-shmem --gm-copy-env sem_model_slice vert_xc_${ftag}_input /ibrixfs1/home/carltape/ADJOINT_TOMO_EXTRA/topo /ibrixfs1/home/carltape/ADJOINT_TOMO_EXTRA/KERNELS_MODELS/models/m16 vb_m16 vert_xc_vb_m16_${ftag}.gmt
sleep 5s

mpirun.lsf --gm-no-shmem --gm-copy-env sem_model_slice vert_xc_${ftag}_input /ibrixfs1/home/carltape/ADJOINT_TOMO_EXTRA/topo /ibrixfs1/home/carltape/ADJOINT_TOMO_EXTRA/KERNELS_MODELS/models/m00 vb_m00 vert_xc_vb_m00_${ftag}.gmt
sleep 5s

done
