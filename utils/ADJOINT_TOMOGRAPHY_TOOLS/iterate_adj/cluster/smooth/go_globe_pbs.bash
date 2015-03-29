#!/bin/bash
#PBS -S /bin/bash

## job name and output file
#PBS -N smooth_kernel
#PBS -j oe
#PBS -o OUTPUT_FILES/job_bulk_c.o

###########################################################
# USER PARAMETERS

## 144 CPUs ( 18*8  ), walltime 10 hour
#PBS -l nodes=18:ppn=8,walltime=10:00:00

numnodes=144

## horizontal: sigmah (in km)
## vertical:     sigmav (in km)
## e.g. period 40 s: wavelength lambda = 40 s * 4 km/s ~ 160 km
sigmah=160
sigmav=40

# kernel to smooth
kernel=bulk_c_kernel

###########################################################

# (takes about: 1h 48 min for sigmah=160/sigmav=40 ...)
#                            25 min for sigmah=40/sigmav=40 ... )

date
echo "kernel: $kernel"
cd $PBS_O_WORKDIR

# obtain lsf job information
cat $PBS_NODEFILE > OUTPUT_FILES/compute_nodes
echo "$PBS_JOBID" > OUTPUT_FILES/jobid


# the neigboring points within 3*sigma+element_size are used for the smoothing
# the 3D Gaussian function is defined as
#   G(x,y,z) = 1/(sqrt(2*pi)*sigma)**3 * exp[-r**2/(2*sigma**2)]
# which is parallel to the 1D Gaussian function

# the horizontal smoothing radius would be sigma_h:
#  160 km (~about wavelength 4km * 40 s )
# the vertical smoothing radius would be sigma_v:
#  40 km ( ~about quater wavelength 4km * 40 s / 4 )

# usage: smooth_sem sigma_h sigma_v kernel_name kernel_input_dir/ topo_dir/

echo "kernel smoothing: $kernel"
mpiexec -np $numnodes $PWD/xsmooth_sem $sigmah $sigmav $kernel $PWD/OUTPUT_SUM $PWD/topo
echo


echo "done successfully"
date
