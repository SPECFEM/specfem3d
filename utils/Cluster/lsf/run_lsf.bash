#!/bin/bash

# this script launches the specfem simulation
# Qinya Liu, May 2007, Caltech

# use the normal queue unless otherwise directed

queue="-q normal"
if [ $# -eq 1 ]; then
	echo "Setting the queue to $1"
	queue="-q $1"
fi

d=`date`
echo "Starting compilation $d"
# regular forward/adjoint simulation
change_simulation_type.pl -f
# kernel simulation
#change_simulation_type.pl -b
make clean
make meshfem3D
make generate_databases
make specfem3D
d=`date`
echo "Finished compilation $d"

# get total number of nodes needed for solver
NPROC=`grep NPROC DATA/Par_file | cut -d = -f 2 `

# compute total number of nodes needed for mesher
NPROC_XI=`grep NPROC_XI DATA/meshfem3D_files/Mesh_Par_file | cut -d = -f 2 `
NPROC_ETA=`grep NPROC_ETA DATA/meshfem3D_files/Mesh_Par_file | cut -d = -f 2 `
# total number of nodes is the product of the values read
numnodes=$(( $NPROC_XI * $NPROC_ETA ))

# checks total number of nodes
if [ $numnodes -neq $NPROC ]; then
  echo "error number of procs mismatch: ",$NPROC_XI,$NPROC_ETA," with ",$NPROC
  exit
fi

echo "Submitting job"
bsub $queue -n $numnodes -W 60 -K < go_mesher_solver_lsf_basin.forward
#kernel simulation
#bsub $queue -n $numnodes -W 60 -K < go_mesher_solver_lsf_basin.kernel
