#!/bin/bash -eu

#############################################################################################
## USER PARAMETERS

# for 1st contribution, we inject adjoint source 1 at receiver 2
# (pay attention to "adj_sources_contribution1" & "X2.DB.MXZ.adj")
# (in 2nd contribution, we will be using "adj_sources_contribution2" & "X1.DB.MXZ.adj")
station_A=DB.X2
station_B=DB.X1

# fortran compiler (to compile adj_traveltime_filter.f90)
f90=gfortran    # or if preferred use: ifort

# specify which kernel we want to visualize
# since the Rayleigh wave is dominantly dependent on shear wave speed, we choose shear wave speed kernels
# you may visualize other kernels if you would like to
kernel="beta_kernel"

#############################################################################################



date
echo "running directory: `pwd`"
echo

#////////////////////////////// SIMULATION IS STARTING ////////////////////////////////////////////////////////////
# Note that one noise sensitivity kernel contains two contributions,
# both the 1st and the 2nd contributions may be obtained through THREE steps.
# Each contribution requires a distinct 'main' receiver, as shown in Tromp et al., 2010, GJI,
# each step requires slightly different Par_file, as documented in the Manual
#
# Please see the paper & Manual for more details on noise simulations
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

## first contribution
# in this part, we start noise simulations for 1st contribution of the noise sensitivity kernels
date
echo "*************************************"
echo "1. contribution..."
echo "*************************************"
echo
# the main receiver is receiver 1
cp -v NOISE_TOMOGRAPHY/irec_main_noise_contribution1  NOISE_TOMOGRAPHY/irec_main_noise

## step 1 of noise simulation
./run_single_step.sh 1
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

## step 2 of noise simulation
./run_single_step.sh 2
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

## calculating adjoint source
mkdir -p SEM
rm -f SEM/*.adj
cp -v OUTPUT_FILES/step_2/${station_A}.MXZ.semd  SEM/
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# note that "a.out" is compiled from "$f90 adj_traveltime_filter.f90"
# this program produces two traces --- adj_sources_contribution1 & adj_sources_contribution2
rm -f a.out
# copy and compile subroutine for adjoint source calculation
$f90 -o a.out adj_traveltime_filter.f90
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi
./a.out
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# setup adjoint traces
cp -v SEM/adj_sources_contribution1  SEM/${station_A}.MXZ.adj
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# zero out other components
awk '{print $1,0.0}' SEM/${station_A}.MXZ.adj > SEM/${station_A}.MXX.adj
awk '{print $1,0.0}' SEM/${station_A}.MXZ.adj > SEM/${station_A}.MXY.adj

# zero out other station B
awk '{print $1,0.0}' SEM/${station_A}.MXZ.adj > SEM/${station_B}.MXX.adj
awk '{print $1,0.0}' SEM/${station_A}.MXZ.adj > SEM/${station_B}.MXY.adj
awk '{print $1,0.0}' SEM/${station_A}.MXZ.adj > SEM/${station_B}.MXZ.adj

## step 3 of noise simulation
./run_single_step.sh 3
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# store kernels
rm -rf contribution_1
mv -v OUTPUT_FILES contribution_1
cp DATABASES_MPI/*kernel*  contribution_1/
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# visualization (refer to other examples, if you don't know the visualization process very well)
# note that "xcombine_vol_data_vtk" is compiled by "make all"
./bin/xcombine_vol_data_vtk 0 3 $kernel DATABASES_MPI/ contribution_1/ 1
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# at the end of this part, we obtain the 1st contribution of the noise sensitivity kernel, stored as:
# $OUTPUT_FILES/contribution_1/$kernel.vtu

echo

## second contribution
# in this part, we start noise simulations for 2nd contribution of the noise sensitivity kernels
date
echo "*************************************"
echo "2. contribution..."
echo "*************************************"
echo

# the main receiver is receiver 2
cp -v NOISE_TOMOGRAPHY/irec_main_noise_contribution2 NOISE_TOMOGRAPHY/irec_main_noise

## step 1 of noise simulation
./run_single_step.sh 1
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

## step 2 of noise simulation
./run_single_step.sh 2
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

## calculating adjoint source
# since it's 2nd contribution, we inject adjoint source 2 at receiver 1
# pay attention to "adj_sources_contribution2" & "X1.DB.MXZ.adj"
# in previous part, we have been using "adj_sources_contribution1" & "X2.DB.MXZ.adj" for the 1st contribution
rm -f SEM/*.adj
cp -v SEM/adj_sources_contribution2 SEM/${station_B}.MXZ.adj
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# zero out other components
awk '{print $1,0.0}' SEM/${station_B}.MXZ.adj > SEM/${station_B}.MXX.adj
awk '{print $1,0.0}' SEM/${station_B}.MXZ.adj > SEM/${station_B}.MXY.adj

# zero out other station
awk '{print $1,0.0}' SEM/${station_B}.MXZ.adj > SEM/${station_A}.MXX.adj
awk '{print $1,0.0}' SEM/${station_B}.MXZ.adj > SEM/${station_A}.MXY.adj
awk '{print $1,0.0}' SEM/${station_B}.MXZ.adj > SEM/${station_A}.MXZ.adj


## step 3 of noise simulation
./run_single_step.sh 3
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# store kernels
rm -rf contribution_2
mv -v OUTPUT_FILES contribution_2
cp DATABASES_MPI/*kernel*  contribution_2/
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# visualization (refer to other examples, if you don't know the visualization process very well)
# this program generates a file "$OUTPUT_FILES/$kernel.mesh"
./bin/xcombine_vol_data_vtk 0 3 $kernel DATABASES_MPI/ contribution_2/ 1
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# at the end of this part, we obtain the 2nd contribution of the noise sensitivity kernel, stored as:
# $OUTPUT_FILES/contribution_2/$kernel.vtu
mkdir -p OUTPUT_FILES
mv -v contribution_1 OUTPUT_FILES/
mv -v contribution_2 OUTPUT_FILES/

echo
echo `date`
echo "done"

