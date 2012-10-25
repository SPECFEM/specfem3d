#!/bin/bash -eu

#################### pre-simulation ###########################################################################
# in this part, we make preparations for the simulations

# now we are in this directory:   (prefix)SPECFEM3D/examples/noise_tomography
# save this path as "script_dir", later we will copy default files from this folder
script_dir=`pwd`

echo `date`
echo "running directory: $script_dir"
echo

echo
echo "(will take about 3 h 30 min)"
echo

# compile the package in root directory
cd ../../
make clean
make > $script_dir/tmp.log
make combine_vol_data >> $script_dir/tmp.log
cd $script_dir

# specify directories for executables, input files and output files
# those values are default for SPECFEM3D
bin="$script_dir/bin"
OUTPUT_FILES="$script_dir/OUTPUT_FILES"
DATA="$script_dir/DATA"

# specify which kernel we want to visualize
# since the Rayleigh wave is dominantly dependent on shear wave speed, we choose shear wave speed kernels
# you may visualize other kernels if you would like to
kernel="beta_kernel"

# create directories for noise simulations and adjoint simulations
# they are also default in SPECFEM3D
mkdir -p $OUTPUT_FILES
mkdir -p $OUTPUT_FILES/SEM
mkdir -p $OUTPUT_FILES/NOISE_TOMOGRAPHY
mkdir -p $OUTPUT_FILES/DATABASES_MPI
mkdir -p $OUTPUT_FILES

# create directories for storing kernels (first contribution and second contribution)
mkdir -p $OUTPUT_FILES/NOISE_TOMOGRAPHY/1st
mkdir -p $OUTPUT_FILES/NOISE_TOMOGRAPHY/2nd

# copy noise input files
cp $script_dir/NOISE_TOMOGRAPHY/S_squared                $OUTPUT_FILES/NOISE_TOMOGRAPHY/
cp $script_dir/NOISE_TOMOGRAPHY/irec_master_noise*       $OUTPUT_FILES/NOISE_TOMOGRAPHY/
cp $script_dir/NOISE_TOMOGRAPHY/nu_master                $OUTPUT_FILES/NOISE_TOMOGRAPHY/

# copy model information
cp $script_dir/DATABASES_MPI/proc*                       $OUTPUT_FILES/DATABASES_MPI/

# copy simulation parameter files
#cp $script_dir/DATA/Par_file*                   $DATA/
#cp $script_dir/DATA/CMTSOLUTION                 $DATA/
#cp $script_dir/DATA/STATIONS*                   $DATA/

# copy and compile subroutine for adjoint source calculation
#cp $script_dir/bin/adj_traveltime_filter.f90             $bin/
cd $bin
ifort adj_traveltime_filter.f90 > tmp.log
ln -s ../../../bin/xgenerate_databases
ln -s ../../../bin/xspecfem3D
ln -s ../../../bin/xcombine_vol_data

#****************************************************************************************************************************************************
#////////////////////////////// SIMULATION IS STARTING //////////////////////////////////////////////////////////////////////////////////////////////
# as theory states, one noise sensitivity kernel contains two contributions
# both the 1st and the 2nd contributions may be obtained through THREE steps
# each contribution requires a distinct 'master' reicever, as shown in Tromp et al., 2010, GJI
# each step requires slightly different Par_file, as documented in the Manual

# if you don't understand above sentences, you will probably get confused later
# please STOP now and go back to the paper & Manual
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#****************************************************************************************************************************************************

#################### first contribution ###########################################################################
# in this part, we start noise simulations for 1st contribution of the noise sensitivity kernels
echo `date`
echo "1. contribution..."
echo
# the master receiver is receiver 1
cp $OUTPUT_FILES/NOISE_TOMOGRAPHY/irec_master_noise_contribution1  $OUTPUT_FILES/NOISE_TOMOGRAPHY/irec_master_noise

# step 1 of noise simulation
cp $DATA/Par_file_step1                         $DATA/Par_file
mpirun -np 4 ./xgenerate_databases
mpirun -np 4 ./xspecfem3D

# step 2 of noise simulation
cp $DATA/Par_file_step2                         $DATA/Par_file
mpirun -np 4 ./xspecfem3D
mv $OUTPUT_FILES/X2.DB.BXZ.semd             $OUTPUT_FILES/SEM/

# calculating adjoint source
# note that "a.out" is compiled from "ifort adj_traveltime_filter.f90"
# this program produces two traces --- adj_sources_contribution1 & adj_sources_contribution2
./a.out
# since it's 1st contribution, we inject adjoint source 1 at receiver 2
# pay attention to "adj_sources_contribution1" & "X2.DB.BXZ.adj"
# we will be using "adj_sources_contribution2" & "X1.DB.BXZ.adj" for the 2nd contribution in next part
rm -f $OUTPUT_FILES/SEM/*.adj
cp $OUTPUT_FILES/SEM/adj_sources_contribution1           $OUTPUT_FILES/SEM/X2.DB.BXZ.adj
awk '{print $1,0.0}' $OUTPUT_FILES/SEM/X2.DB.BXZ.adj > $OUTPUT_FILES/SEM/X2.DB.BXX.adj
awk '{print $1,0.0}' $OUTPUT_FILES/SEM/X2.DB.BXZ.adj > $OUTPUT_FILES/SEM/X2.DB.BXY.adj

awk '{print $1,0.0}' $OUTPUT_FILES/SEM/X2.DB.BXZ.adj > $OUTPUT_FILES/SEM/X1.DB.BXX.adj
awk '{print $1,0.0}' $OUTPUT_FILES/SEM/X2.DB.BXZ.adj > $OUTPUT_FILES/SEM/X1.DB.BXY.adj
awk '{print $1,0.0}' $OUTPUT_FILES/SEM/X2.DB.BXZ.adj > $OUTPUT_FILES/SEM/X1.DB.BXZ.adj


# step 3 of noise simulation
cp $DATA/Par_file_step3                         $DATA/Par_file
mpirun -np 4 ./xspecfem3D

# store kernels
cp $OUTPUT_FILES/DATABASES_MPI/*kernel*                  $OUTPUT_FILES/NOISE_TOMOGRAPHY/1st/

# visualization (refer to other examples, if you don't know the visualization process very well)
# note that "xcombine_vol_data" is compiled by "make combine_vol_data"
# this program generates a file "$OUTPUT_FILES/$kernel.mesh"
./xcombine_vol_data 0 3 $kernel $OUTPUT_FILES/DATABASES_MPI/  $OUTPUT_FILES/ 1
# you may need to install "mesh2vtu" package first, before you can use "mesh2vtu.pl"
# convert "$OUTPUT_FILES/$kernel.mesh" to "$OUTPUT_FILES/NOISE_TOMOGRAPHY/1st_$kernel.vtu"
# which can be loaded and visualized in Paraview
mesh2vtu.pl -i $OUTPUT_FILES/$kernel.mesh -o $OUTPUT_FILES/NOISE_TOMOGRAPHY/1st_$kernel.vtu

# at the end of this part, we obtain the 1st contribution of the noise sensitivity kernel, stored as:
# $OUTPUT_FILES/NOISE_TOMOGRAPHY/1st_$kernel.vtu

echo

#################### second contribution ###########################################################################
# in this part, we start noise simulations for 2nd contribution of the noise sensitivity kernels

echo `date`
echo "2. contribution..."
echo

# the master receiver is receiver 2
cp $OUTPUT_FILES/NOISE_TOMOGRAPHY/irec_master_noise_contribution2  $OUTPUT_FILES/NOISE_TOMOGRAPHY/irec_master_noise

# step 1 of noise simulation
cp $DATA/Par_file_step1                         $DATA/Par_file
mpirun -np 4 ./xspecfem3D

# step 2 of noise simulation
cp $DATA/Par_file_step2                         $DATA/Par_file
mpirun -np 4 ./xspecfem3D

# calculating adjoint source
# since it's 2nd contribution, we inject adjoint source 2 at receiver 1
# pay attention to "adj_sources_contribution2" & "X1.DB.BXZ.adj"
# we have been using "adj_sources_contribution1" & "X2.DB.BXZ.adj" for the 1st contribution in previous part
rm -f $OUTPUT_FILES/SEM/*.adj
cp $OUTPUT_FILES/SEM/adj_sources_contribution2           $OUTPUT_FILES/SEM/X1.DB.BXZ.adj
awk '{print $1,0.0}' $OUTPUT_FILES/SEM/X1.DB.BXZ.adj > $OUTPUT_FILES/SEM/X1.DB.BXX.adj
awk '{print $1,0.0}' $OUTPUT_FILES/SEM/X1.DB.BXZ.adj > $OUTPUT_FILES/SEM/X1.DB.BXY.adj

awk '{print $1,0.0}' $OUTPUT_FILES/SEM/X1.DB.BXZ.adj > $OUTPUT_FILES/SEM/X2.DB.BXX.adj
awk '{print $1,0.0}' $OUTPUT_FILES/SEM/X1.DB.BXZ.adj > $OUTPUT_FILES/SEM/X2.DB.BXY.adj
awk '{print $1,0.0}' $OUTPUT_FILES/SEM/X1.DB.BXZ.adj > $OUTPUT_FILES/SEM/X2.DB.BXZ.adj


# step 3 of noise simulation
cp $DATA/Par_file_step3                         $DATA/Par_file
mpirun -np 4 ./xspecfem3D

# store kernels
cp $OUTPUT_FILES/DATABASES_MPI/*kernel*                  $OUTPUT_FILES/NOISE_TOMOGRAPHY/2nd/

# visualization (refer to other examples, if you don't know the visualization process very well)
# this program generates a file "$OUTPUT_FILES/$kernel.mesh"
./xcombine_vol_data 0 3 $kernel $OUTPUT_FILES/DATABASES_MPI/  $OUTPUT_FILES/ 1
# you may need to install "mesh2vtu" package first, before you can use "mesh2vtu.pl"
# convert "$OUTPUT_FILES/$kernel.mesh" to "$OUTPUT_FILES/NOISE_TOMOGRAPHY/1st_$kernel.vtu"
# which can be loaded and visualized in Paraview
mesh2vtu.pl -i $OUTPUT_FILES/$kernel.mesh -o $OUTPUT_FILES/NOISE_TOMOGRAPHY/2nd_$kernel.vtu

# at the end of this part, we obtain the 2nd contribution of the noise sensitivity kernel, stored as:
# $OUTPUT_FILES/NOISE_TOMOGRAPHY/2nd_$kernel.vtu

echo
echo `date`
echo "done"
