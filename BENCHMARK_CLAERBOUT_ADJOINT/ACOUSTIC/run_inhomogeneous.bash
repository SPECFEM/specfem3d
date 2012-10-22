#!/bin/bash -eu
echo
echo `date`
echo

NSTEP=`grep NSTEP   ./DATA/Par_file | cut -d = -f 2 | sed 's/ //g'`
DT=`grep       DT   ./DATA/Par_file | cut -d = -f 2 | sed 's/ //g'`
NPROC=`grep NPROC   ./DATA/Par_file | cut -d = -f 2 `

#daniel:
echo "Par_file parameters:"
echo "NSTEP=$NSTEP"
echo "DT=$DT"
echo "NPROC=$NPROC"
echo

#NSTEP=3000
#DT=0.001
#percent=0.0001

percent=0.02
echo "perturbation: $percent"
echo

#echo "NSTEP=$NSTEP"
#echo "DT=$DT"
#echo "NPROC=$NPROC"

bin="$PWD/bin"
DATA="$bin/../DATA"
OUTPUT_FILES="$bin/../OUTPUT_FILES"
models="$bin/../models"
SESAME="$bin/../../../../"


MPIRUN="mpirun -np $NPROC "
MPIFC="mpif90 -assume byterecl "

##############################################

#daniel
do_setup=$1

##############################################

if [ "$do_setup" == "" ]; then

echo
echo "setting up example..."
echo

rm -rf $OUTPUT_FILES $models $bin
mkdir -p $OUTPUT_FILES/DATABASES_MPI $OUTPUT_FILES $OUTPUT_FILES/SEM/dat $OUTPUT_FILES/SEM/syn
mkdir -p $models/initial_model $models/target_model $bin

cd $SESAME
cp ./src/shared/constants.h ./src/shared/constants.h_backup
cp ./src/specfem3D/save_adjoint_kernels.f90 ./src/specfem3D/save_adjoint_kernels.f90_backup
cp ./src/generate_databases/get_model.f90 ./src/generate_databases/get_model.f90_backup
cp ./src/specfem3D/write_seismograms.f90  ./src/specfem3D/write_seismograms.f90_backup

#daniel
#cp $bin/../constants.h ./src/shared/constants.h
sed -i "s:SU_FORMAT.*:SU_FORMAT = .true.:g" ./src/shared/constants.h
sed -i "s:FIX_UNDERFLOW_PROBLEM.*:FIX_UNDERFLOW_PROBLEM = .false.:g" ./src/shared/constants.h

#daniel
#cp $bin/../get_model_internal.f90 ./src/generate_databases/get_model.f90
sed -i "s:USE_EXTERNAL_FILES.*=.*:USE_EXTERNAL_FILES = .false.:" ./src/generate_databases/get_model.f90

#daniel
#cp $bin/../save_adjoint_kernels.f90 ./src/specfem3D/save_adjoint_kernels.f90
sed -i "s:SAVE_WEIGHTS.*=.*:SAVE_WEIGHTS = .true.:" ./src/specfem3D/save_adjoint_kernels.f90

#daniel
#cp $bin/../write_seismograms.f90 ./src/specfem3D/write_seismograms.f90
patch ./src/specfem3D/write_seismograms.f90 < $bin/../write_seismograms.patch

make >& $bin/../make.log

cp ./bin/xmeshfem3D            $bin/xmeshfem3D
cp ./bin/xgenerate_databases   $bin/xgenerate_databases_internal
cp ./bin/xspecfem3D            $bin/xspecfem3D

#daniel
make xcombine_vol_data >> $bin/../make.log
cp ./bin/xcombine_vol_data   $bin/

#daniel
#cp $bin/../get_model_external.f90 ./src/generate_databases/get_model.f90
sed -i "s:USE_EXTERNAL_FILES.*=.*:USE_EXTERNAL_FILES = .true.:" ./src/generate_databases/get_model.f90

make >& $bin/../make.log

cp ./bin/xgenerate_databases   $bin/xgenerate_databases_external

fi

cd $bin
$MPIFC $bin/../random_model_generation.f90 -o ./xrandom_model > $bin/../compile.log
$MPIFC $bin/../adj_seismogram.f90 -o ./xadj > $bin/../compile.log
$MPIFC $bin/../postprocessing.f90 -o ./xpostprocessing > $bin/../compile.log

########################### dat ########################################
FILE="$DATA/Par_file"
sed -e "s#^SIMULATION_TYPE.*#SIMULATION_TYPE = 1 #g"  < $FILE > ./tmp; mv ./tmp $FILE
sed -e "s#^SAVE_FORWARD.*#SAVE_FORWARD = .true. #g"   < $FILE > ./tmp; mv ./tmp $FILE
sed -e "s#^NSTEP.*#NSTEP = $NSTEP #g"                 < $FILE > ./tmp; mv ./tmp $FILE
sed -e "s#^DT.*#DT = $DT #g"                          < $FILE > ./tmp; mv ./tmp $FILE

cd $bin

echo "data simulation: $MPIRUN ./xmeshfem3D ..."
$MPIRUN ./xmeshfem3D
echo "data simulation: $MPIRUN ./xgenerate_databases ..."
$MPIRUN ./xgenerate_databases_internal
echo "data simulation: $MPIRUN ./xspecfem3D ..."
$MPIRUN ./xspecfem3D

#daniel
cp -rp $OUTPUT_FILES $OUTPUT_FILES.dat.forward

mv $OUTPUT_FILES/*SU $OUTPUT_FILES/SEM/dat/
cp $OUTPUT_FILES/DATABASES_MPI/*rho.bin $models/target_model/
cp $OUTPUT_FILES/DATABASES_MPI/*vp.bin  $models/target_model/
cp $OUTPUT_FILES/DATABASES_MPI/*vs.bin  $models/target_model/
########################### syn ########################################

$MPIRUN ./xrandom_model $percent

echo "syn simulation: $MPIRUN ./xgenerate_databases ..."
$MPIRUN ./xgenerate_databases_external
echo "syn simulation: $MPIRUN ./xspecfem3D ..."
$MPIRUN ./xspecfem3D

#daniel
cp -rp $OUTPUT_FILES $OUTPUT_FILES.syn.forward

mv $OUTPUT_FILES/*SU $OUTPUT_FILES/SEM/syn/
cp $OUTPUT_FILES/DATABASES_MPI/*rho.bin $models/initial_model/
cp $OUTPUT_FILES/DATABASES_MPI/*vp.bin  $models/initial_model/
cp $OUTPUT_FILES/DATABASES_MPI/*vs.bin  $models/initial_model/
########################### adj sources ################################

$MPIRUN ./xadj $NSTEP $DT

########################### adj ########################################
FILE="$DATA/Par_file"
sed -e "s#^SIMULATION_TYPE.*#SIMULATION_TYPE = 3 #g"  < $FILE > ./tmp; mv ./tmp $FILE
sed -e "s#^SAVE_FORWARD.*#SAVE_FORWARD = .false. #g"  < $FILE > ./tmp; mv ./tmp $FILE

echo "adj simulation: $MPIRUN ./xspecfem3D ..."
$MPIRUN ./xspecfem3D

#daniel
cp -rp $OUTPUT_FILES $OUTPUT_FILES.syn.adjoint

./xpostprocessing $NSTEP $DT $NPROC

if [ "$do_setup" == "" ]; then

cd $SESAME
mv ./src/generate_databases/get_model.f90_backup ./src/generate_databases/get_model.f90
mv ./src/specfem3D/save_adjoint_kernels.f90_backup ./src/specfem3D/save_adjoint_kernels.f90
mv ./src/shared/constants.h_backup ./src/shared/constants.h
mv ./src/specfem3D/write_seismograms.f90_backup ./src/specfem3D/write_seismograms.f90

#cd $bin/../
#rm -f compile.log make.log
#rm -rf $OUTPUT_FILES $models $bin

fi

echo
echo "done: `date`"
echo

