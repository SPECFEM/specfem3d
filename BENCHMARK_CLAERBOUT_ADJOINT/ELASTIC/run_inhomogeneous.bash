#!/bin/bash -eu

NSTEP=`grep NSTEP   ./in_data_files/Par_file | cut -d = -f 2 | sed 's/ //g'`
DT=`grep       DT   ./in_data_files/Par_file | cut -d = -f 2 | sed 's/ //g'`
NPROC=`grep NPROC   ./in_data_files/Par_file | cut -d = -f 2 `

NSTEP=3000
DT=0.001
percent=0.0001

echo "NSTEP=$NSTEP"
echo "DT=$DT"
echo "NPROC=$NPROC"

bin="$PWD/bin"
in_data_files="$bin/../in_data_files"
in_out_files="$bin/../in_out_files"
models="$bin/../models"
SESAME="$bin/../../../../"


MPIRUN="mpirun -np $NPROC "
MPIFC="mpif90 -assume byterecl "

rm -rf $in_out_files $models $bin
mkdir -p $in_out_files/DATABASES_MPI $in_out_files/OUTPUT_FILES $in_out_files/SEM/dat $in_out_files/SEM/syn
mkdir -p $models/initial_model $models/target_model $bin

cd $SESAME
cp ./src/shared/constants.h ./src/shared/constants.h_backup
cp ./src/specfem3D/save_adjoint_kernels.f90 ./src/specfem3D/save_adjoint_kernels.f90_backup
cp ./src/generate_databases/get_model.f90 ./src/generate_databases/get_model.f90_backup
cp $bin/../constants.h ./src/shared/constants.h
cp $bin/../get_model_internal.f90 ./src/generate_databases/get_model.f90
cp $bin/../save_adjoint_kernels.f90 ./src/specfem3D/save_adjoint_kernels.f90
make > $bin/../make.log

cp ./bin/xmeshfem3D            $bin/xmeshfem3D
cp ./bin/xgenerate_databases   $bin/xgenerate_databases_internal
cp ./bin/xspecfem3D            $bin/xspecfem3D

cp $bin/../get_model_external.f90 ./src/generate_databases/get_model.f90
make > $bin/../make.log
cp ./bin/xgenerate_databases   $bin/xgenerate_databases_external

cd $bin
$MPIFC $bin/../random_model_generation.f90 -o ./xrandom_model > $bin/../compile.log
$MPIFC $bin/../adj_seismogram.f90 -o ./xadj > $bin/../compile.log
$MPIFC $bin/../postprocessing.f90 -o ./xpostprocessing > $bin/../compile.log

########################### dat ########################################
FILE="$in_data_files/Par_file"
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

mv $in_out_files/OUTPUT_FILES/*SU $in_out_files/SEM/dat/
cp $in_out_files/DATABASES_MPI/*rho.bin $models/target_model/
cp $in_out_files/DATABASES_MPI/*vp.bin  $models/target_model/
cp $in_out_files/DATABASES_MPI/*vs.bin  $models/target_model/
########################### syn ########################################

$MPIRUN ./xrandom_model $percent

echo "syn simulation: $MPIRUN ./xgenerate_databases ..."
$MPIRUN ./xgenerate_databases_external
echo "syn simulation: $MPIRUN ./xspecfem3D ..."
$MPIRUN ./xspecfem3D

mv $in_out_files/OUTPUT_FILES/*SU $in_out_files/SEM/syn/
cp $in_out_files/DATABASES_MPI/*rho.bin $models/initial_model/
cp $in_out_files/DATABASES_MPI/*vp.bin  $models/initial_model/
cp $in_out_files/DATABASES_MPI/*vs.bin  $models/initial_model/
########################### adj sources ################################

$MPIRUN ./xadj $NSTEP $DT

########################### adj ########################################
FILE="$in_data_files/Par_file"
sed -e "s#^SIMULATION_TYPE.*#SIMULATION_TYPE = 3 #g"  < $FILE > ./tmp; mv ./tmp $FILE
sed -e "s#^SAVE_FORWARD.*#SAVE_FORWARD = .false. #g"  < $FILE > ./tmp; mv ./tmp $FILE

echo "adj simulation: $MPIRUN ./xspecfem3D ..."
$MPIRUN ./xspecfem3D

./xpostprocessing $NSTEP $DT $NPROC 

cd $SESAME
mv ./src/generate_databases/get_model.f90_backup ./src/generate_databases/get_model.f90
mv ./src/specfem3D/save_adjoint_kernels.f90_backup ./src/specfem3D/save_adjoint_kernels.f90
mv ./src/shared/constants.h_backup ./src/shared/constants.h

cd $bin/../
rm -f compile.log make.log
rm -rf $in_out_files $models $bin
