#!/bin/bash -eu
echo "running example: `date`"
currentdir=`pwd`

# gets Par_file parameters
NSTEP=`grep ^NSTEP   ./DATA/Par_file | cut -d = -f 2 | sed 's/ //g'`
DT=`grep       ^DT   ./DATA/Par_file | cut -d = -f 2 | sed 's/ //g'`
NPROC=`grep ^NPROC   ./DATA/Par_file | cut -d = -f 2 `

echo "Par_file parameters:"
echo "  NSTEP = $NSTEP"
echo "  DT    = $DT"
echo "  NPROC = $NPROC"
echo

# perturbation (should be small enough for approximating S(m - m0) ~ S(m) - S(m0)
percent=0.02

echo "perturbation: $percent"
echo

# SPECFEM3D root directory
ROOT=../../../

## compiler
# intel
#FC="ifort -assume byterecl "
# gnu
FC="gfortran "

## MPI parallel run
MPIRUN="mpirun -np $NPROC "

# create and compile all setup
do_setup=1

##############################################

if [ "$do_setup" == "1" ]; then
echo
echo "setting up example..."
echo

# cleans output files
mkdir -p OUTPUT_FILES
rm -rf OUTPUT_FILES/*
mkdir -p OUTPUT_FILES/DATABASES_MPI

mkdir -p SEM
rm -rf SEM/*
mkdir -p SEM/dat SEM/syn

mkdir -p MODELS
rm -rf MODELS/*
mkdir -p MODELS/initial_model MODELS/target_model

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

make -j4 all >& $currentdir/make.log
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

fi # do_setup

echo
echo "compiling tools..."
echo

$FC ../random_model_generation.f90 -o ./bin/xrandom_model
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi
$FC ../adj_seismogram.f90 -o ./bin/xadj
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi
$FC ../postprocessing.f90 -o ./bin/xpostprocessing
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

########################### dat ########################################

sed -i "s:^SIMULATION_TYPE .*:SIMULATION_TYPE = 1:" ./DATA/Par_file
sed -i "s:^SAVE_FORWARD .*:SAVE_FORWARD = .false.:" ./DATA/Par_file

sed -i "s:^MODEL .*=.*:MODEL = default:" ./DATA/Par_file

echo "data simulation: xmeshfem3D ..."
$MPIRUN ./bin/xmeshfem3D
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo "data simulation: xgenerate_databases ..."
$MPIRUN ./bin/xgenerate_databases
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo "data simulation: xspecfem3D ..."
$MPIRUN ./bin/xspecfem3D
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# backup copy
echo
echo "backup output: OUTPUT_FILES.dat.forward"
rm -rf OUTPUT_FILES.dat.forward
cp -rp OUTPUT_FILES OUTPUT_FILES.dat.forward

echo "backup traces: SEM/dat"
mv -v OUTPUT_FILES/*SU SEM/dat/

# target model
echo
echo "backup model : MODELS/target_model"
cp -v OUTPUT_FILES/DATABASES_MPI/*rho.bin MODELS/target_model/
cp -v OUTPUT_FILES/DATABASES_MPI/*vp.bin  MODELS/target_model/
cp -v OUTPUT_FILES/DATABASES_MPI/*vs.bin  MODELS/target_model/

########################### syn ########################################
echo
echo "setting up perturbed model..."
echo "> ./bin/xrandom_model $percent $NPROC "
echo
./bin/xrandom_model $percent $NPROC
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

sed -i "s:^SIMULATION_TYPE .*:SIMULATION_TYPE = 1:" ./DATA/Par_file
sed -i "s:^SAVE_FORWARD .*:SAVE_FORWARD = .true.:" ./DATA/Par_file

sed -i "s:^MODEL .*=.*:MODEL = gll:" ./DATA/Par_file

echo "syn simulation: xgenerate_databases ..."
$MPIRUN ./bin/xgenerate_databases
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo "syn simulation: xspecfem3D ..."
$MPIRUN ./bin/xspecfem3D
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# backup
echo
echo "backup output: OUTPUT_FILES.syn.forward"
rm -rf OUTPUT_FILES.syn.forward
cp -rp OUTPUT_FILES OUTPUT_FILES.syn.forward

echo "backup traces: SEM/syn"
mv -v OUTPUT_FILES/*SU SEM/syn/

# initial model
echo
echo "backup model : MODELS/initial"
cp -v OUTPUT_FILES/DATABASES_MPI/*rho.bin MODELS/initial_model/
cp -v OUTPUT_FILES/DATABASES_MPI/*vp.bin  MODELS/initial_model/
cp -v OUTPUT_FILES/DATABASES_MPI/*vs.bin  MODELS/initial_model/

########################### adj sources ################################
echo
echo "creating adjoint sources..."
echo "> ./bin/xadj  $NSTEP $DT $NPROC acoustic"
echo
./bin/xadj $NSTEP $DT $NPROC acoustic
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

########################### adj ########################################

sed -i "s:^SIMULATION_TYPE .*:SIMULATION_TYPE = 3:" ./DATA/Par_file
sed -i "s:^SAVE_FORWARD .*:SAVE_FORWARD = .false.:" ./DATA/Par_file

echo "adj simulation: xspecfem3D ..."
$MPIRUN ./bin/xspecfem3D
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# backup
echo
echo "backup output: OUTPUT_FILES.syn.adjoint"
rm -rf OUTPUT_FILES.syn.adjoint
cp -rp OUTPUT_FILES OUTPUT_FILES.syn.adjoint

echo
echo "postprocessing..."
echo "> ./bin/xpostprocessing $NSTEP $DT $NPROC acoustic"
echo
./bin/xpostprocessing $NSTEP $DT $NPROC acoustic
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

if [ "$do_setup" == "1" ]; then
# restore original files
cd $ROOT
mv -v ./setup/constants.h_backup ./setup/constants.h
fi

echo
echo "done: `date`"
echo

