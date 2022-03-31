#!/bin/bash
#
# script runs mesher and solver (in serial)
# using this example setup
#
echo "running example: `date`"
currentdir=`pwd`

cd $currentdir

##############################################
## Acoustic benchmark

# Simulation type 1 == acoustic / 2 == elastic / 3 == coupled acoustic-elastic
SIM_TYPE=1

# perturbation model parameter rho/vp/vs (e.g. "rho" or "vp" or "rhovp")
perturb_param="vp"

# perturbation (should be small enough for approximating S(m - m0) ~ S(m) - S(m0)
perturb_percent=0.02

##############################################

echo
echo "setup:"
echo "  SIM_TYPE                : $SIM_TYPE     (1 == acoustic / 2 == elastic / 3 == coupled acoustic-elastic)"
echo "  perturbation parameter  : $perturb_param"
echo "  perturbation percent    : $perturb_percent"
echo

# gets Par_file parameters
# Get the number of processors
NPROC=`grep '^NPROC ' DATA/Par_file | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`
NSTEP=`grep '^NSTEP ' DATA/Par_file | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`
DT=`grep '^DT ' DATA/Par_file | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`
BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

echo "Par_file parameters:"
echo "  NPROC = $NPROC"
echo "  NSTEP = $NSTEP"
echo "  DT    = $DT"
echo

# create and compile all setup
do_setup=1

##############################################

if [ "$do_setup" == "1" ]; then
echo
echo "setting up example..."
echo

# cleans files
rm -rf DATA/*.bin

mkdir -p OUTPUT_FILES
rm -rf OUTPUT_FILES/*
mkdir -p OUTPUT_FILES/

mkdir -p SEM
rm -rf SEM/*
mkdir -p SEM/dat SEM/syn

mkdir -p MODELS
rm -rf MODELS/*

mkdir -p MODELS/initial_model MODELS/target_model

mkdir -p KERNELS
rm -rf KERNELS/*

rm -f adj_seismogram.py
ln -s ../adj_seismogram.py
rm -f model_add_Gaussian_perturbation.py
ln -s ../model_add_Gaussian_perturbation.py
rm -f model_update.py
ln -s ../model_update.py
rm -f kernel_evaluation_postprocessing.py
ln -s ../kernel_evaluation_postprocessing.py
rm -f helper_functions.py
ln -s ../helper_functions.py

rm -f change_simulation_type.pl
ln -s ../../../utils/change_simulation_type.pl

# SPECFEM3D root directory
ROOT=../../../

# re-compile executables for getting integration weights
cd $ROOT/
cp -v ./setup/constants.h ./setup/constants.h_backup
sed -i "s:SAVE_WEIGHTS.*=.*:SAVE_WEIGHTS = .true.:" ./setup/constants.h

# patch to write out acceleration, i.e. pressure, values in SU-format
# (- potential_dot_dot will be pressure for acoustic simulations, axd/ayd/azd will be pressure in 3-d vector format)
#cp ./src/specfem3D/write_seismograms.f90  ./src/specfem3D/write_seismograms.f90_backup
#patch ./src/specfem3D/write_seismograms.f90 < $currentdir/write_seismograms.patch

echo
echo "compiling binaries..."
echo

make -j4 >& $currentdir/make.log
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

cd $currentdir

# links executables
mkdir -p bin
cd bin/
rm -f ./*
ln -s ../$ROOT/bin/xmeshfem3D
ln -s ../$ROOT/bin/xgenerate_databases
ln -s ../$ROOT/bin/xspecfem3D
ln -s ../$ROOT/bin/xcombine_vol_data_vtk
cd ../

fi

########################### data ########################################

##
## data simulation
##

## forward simulation
echo
echo "running data forward simulation"
echo
./change_simulation_type.pl -f

# saving model files
sed -i "s:^MODEL .*=.*:MODEL = default:" ./DATA/Par_file
sed -i "s:^SAVE_MESH_FILES .*=.*:SAVE_MESH_FILES = .true.:" ./DATA/Par_file

./run_this_example.sh > output.log

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
mv -v output.log OUTPUT_FILES/

# backup copy
rm -rf OUTPUT_FILES.dat.forward
cp -rp OUTPUT_FILES OUTPUT_FILES.dat.forward

cp -v OUTPUT_FILES/*_SU SEM/dat/

# target model
cp -v $BASEMPIDIR/*rho.bin MODELS/target_model/
cp -v $BASEMPIDIR/*vp.bin  MODELS/target_model/
cp -v $BASEMPIDIR/*vs.bin  MODELS/target_model/
cp -v $BASEMPIDIR/*external_mesh.bin  MODELS/target_model/

# cleans OUTPUT_FILES, leaving model file
rm -f OUTPUT_FILES/*.* OUTPUT_FILES/*_SU

########################### model perturbation ################################

echo
echo "setting up perturbed model..."
echo "> ./model_add_Gaussian_perturbation.py $perturb_param $perturb_percent $NPROC "
echo
./model_add_Gaussian_perturbation.py $perturb_param $perturb_percent $NPROC

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# replaces model files with perturbed ones
for (( iproc = 0; iproc < $NPROC; iproc++ )); do
  rank=`printf "%06i\n" $iproc`
  cp -v $BASEMPIDIR/proc${rank}_rho_gaussian.bin $BASEMPIDIR/proc${rank}_rho.bin
  cp -v $BASEMPIDIR/proc${rank}_vp_gaussian.bin $BASEMPIDIR/proc${rank}_vp.bin
  cp -v $BASEMPIDIR/proc${rank}_vs_gaussian.bin $BASEMPIDIR/proc${rank}_vs.bin
  if [[ $? -ne 0 ]]; then exit 1; fi
done

########################### synthetics ################################

##
## synthetics simulation
##

echo
echo "running synthetics forward simulation (with saving forward wavefield)"
echo
./change_simulation_type.pl -F

# Par_file using GLL model
sed -i "s:^MODEL .*=.*:MODEL = gll:" ./DATA/Par_file
sed -i "s:^SAVE_MESH_FILES .*=.*:SAVE_MESH_FILES = .false.:" ./DATA/Par_file

./run_this_example.sh noclean > output.log

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
mv -v output.log OUTPUT_FILES/output.forward.log

# backup copy
rm -rf OUTPUT_FILES.syn.forward
cp -rp OUTPUT_FILES OUTPUT_FILES.syn.forward

# initial model
mv -v OUTPUT_FILES/*_SU SEM/syn/
cp -v $BASEMPIDIR/*rho.bin MODELS/initial_model/
cp -v $BASEMPIDIR/*vp.bin  MODELS/initial_model/
cp -v $BASEMPIDIR/*vs.bin  MODELS/initial_model/
cp -v $BASEMPIDIR/*external_mesh.bin  MODELS/initial_model/

########################### adj sources ################################
## adjoint sources
echo
echo "creating adjoint sources..."

# replaces model files with perturbed ones
for (( iproc = 0; iproc < $NPROC; iproc++ )); do
  rank=`printf "%i\n" $iproc`
  # pressure traces
  if [ -e OUTPUT_FILES.syn.forward/${rank}_p_SU ]; then
    syn=OUTPUT_FILES.syn.forward/${rank}_p_SU
    dat=OUTPUT_FILES.dat.forward/${rank}_p_SU
    echo "> ./adj_seismogram.py $syn $dat"
    echo
    # adjoint source f^adj = (s - d)
    ./adj_seismogram.py $syn $dat
    # checks exit code
    if [[ $? -ne 0 ]]; then exit 1; fi
  fi
done
########################### kernel ################################

## kernel simulation
echo
echo "running kernel simulation"
echo
./change_simulation_type.pl -b

# Par_file using GLL model
sed -i "s:^MODEL .*=.*:MODEL = gll:" ./DATA/Par_file
sed -i "s:^SAVE_MESH_FILES .*=.*:SAVE_MESH_FILES = .false.:" ./DATA/Par_file

# cleans OUTPUT_FILES, leaving model file
rm -f OUTPUT_FILES/*.* OUTPUT_FILES/*_SU

# In principle we do not need rerun meshing in the adjoint run.

./run_this_example.sh noclean > output.log
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
mv -v output.log OUTPUT_FILES/output.kernel.log

# backup
rm -rf OUTPUT_FILES.syn.adjoint
cp -rp OUTPUT_FILES OUTPUT_FILES.syn.adjoint

# kernels
cp -v OUTPUT_FILES/output.kernel.log KERNELS/
cp -v $BASEMPIDIR/*_kernel.* KERNELS/


########################### model update ################################

echo
echo "model update"
echo

# takes absolute value of percent
update_percent=$( sed "s/-//" <<< $perturb_percent )

# takes half of absolute perturbation for model update scaling length
update_percent=`echo "$update_percent" | awk '{print $0 * 0.5}'`

echo "> ./model_update.py $NPROC $SIM_TYPE $update_percent $perturb_param"
echo
./model_update.py $NPROC $SIM_TYPE $update_percent $perturb_param

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# replaces model files with perturbed ones
for (( iproc = 0; iproc < $NPROC; iproc++ )); do
  rank=`printf "%06i\n" $iproc`
  cp -v $BASEMPIDIR/proc${rank}_rho_new.bin $BASEMPIDIR/proc${rank}_rho.bin
  cp -v $BASEMPIDIR/proc${rank}_vp_new.bin $BASEMPIDIR/proc${rank}_vp.bin
  cp -v $BASEMPIDIR/proc${rank}_vs_new.bin $BASEMPIDIR/proc${rank}_vs.bin
  if [[ $? -ne 0 ]]; then exit 1; fi
done

########################### final forward ################################

## forward simulation
echo
echo "running forward simulation (updated model)"
echo
./change_simulation_type.pl -f

# cleans OUTPUT_FILES, leaving model file
rm -f OUTPUT_FILES/*.* OUTPUT_FILES/*_SU

# Par_file using GLL model
sed -i "s:^MODEL .*=.*:MODEL = gll:" ./DATA/Par_file
sed -i "s:^SAVE_MESH_FILES .*=.*:SAVE_MESH_FILES = .false.:" ./DATA/Par_file

# In principle we do not need rerun meshing in the adjoint run.

./run_this_example.sh noclean > output.log
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
mv -v output.log OUTPUT_FILES/output.log

# backup
rm -rf OUTPUT_FILES.syn.updated
cp -rp OUTPUT_FILES OUTPUT_FILES.syn.updated


########################### kernel ################################

echo
echo "postprocessing..."
echo "> ./kernel_evaluation_postprocessing.py $NSTEP $DT $NPROC $SIM_TYPE"
echo
./kernel_evaluation_postprocessing.py $NSTEP $DT $NPROC $SIM_TYPE

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# restore constants.h
if [ "$do_setup" == "1" ]; then
  # restore original files
  cd $ROOT
  mv -v ./setup/constants.h_backup ./setup/constants.h
fi

echo
echo "done: `date`"
echo
