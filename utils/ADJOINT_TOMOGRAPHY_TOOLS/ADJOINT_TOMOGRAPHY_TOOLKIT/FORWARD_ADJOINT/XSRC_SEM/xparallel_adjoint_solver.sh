#!/bin/sh

#PBS -q tromp
#PBS -N XEUSOLVER
#PBS -l nodes=13:ppn=8
#PBS -l walltime=15:00:00
#PBS -j oe
#PBS -k oe
#PBS -o job_src2.log


echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

local_dir=`grep LOCAL_PATH DATA/Par_file | cut -d '=' -f 2`

# obtain lsf job information
cat $PBS_NODEFILE > OUTPUT_FILES/compute_nodes
echo "$PBS_JOBID" > OUTPUT_FILES/jobid

if [ -d /scratch/hejunzhu ]; then
  echo rm /scratch/hejunzhu
  pbsdsh -u rm -rf /scratch/hejunzhu
fi
echo mkdir local dir...
pbsdsh -u mkdir -p /scratch/hejunzhu

echo change simulation type...
./change_simulation_type.pl -F

echo copying mesh file...
mpiexec -np 100 ./xcopy_local_forward

echo submit job...
mpiexec -np 100 ./xspecfem3D
echo solver done

echo tar the synthetics
cd OUTPUT_FILES
#tar -czvf SEM.forward.tar.gz *.sac
rm *.sac
cd ../
echo end tar

echo change simulation type again...
./change_simulation_type.pl -b

echo sumit adjoint simulation
mpiexec -np 100 ./xspecfem3D


echo tar the synthetics
cd OUTPUT_FILES
#tar -czvf SEM.backward.tar.gz *.sac
rm *.sac
cd ../
echo end tar

echo collect kernels...
pbsdsh -u bash -c 'cp -p /scratch/hejunzhu/proc*_reg1_bulk_betav_kernel.bin $PBS_O_WORKDIR/KERNEL/'
pbsdsh -u bash -c 'cp -p /scratch/hejunzhu/proc*_reg1_bulk_betah_kernel.bin $PBS_O_WORKDIR/KERNEL/'
pbsdsh -u bash -c 'cp -p /scratch/hejunzhu/proc*_reg1_bulk_c_kernel.bin $PBS_O_WORKDIR/KERNEL/'
pbsdsh -u bash -c 'cp -p /scratch/hejunzhu/proc*_reg1_eta_kernel.bin $PBS_O_WORKDIR/KERNEL/'
pbsdsh -u bash -c 'cp -p /scratch/hejunzhu/proc*_reg1_rho_kernel.bin $PBS_O_WORKDIR/KERNEL/'
pbsdsh -u bash -c 'cp -p /scratch/hejunzhu/proc*_reg1_hess_kernel.bin $PBS_O_WORKDIR/KERNEL/'



echo cleaning nodes...
pbsdsh -u rm -rf /scratch/hejunzhu

echo done sucessfully
