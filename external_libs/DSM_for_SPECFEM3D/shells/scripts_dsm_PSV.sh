#
#  scripts for computing tractions and velocity on chunk boundary
#
#
#
function setup_process ()
{
$BIN_DSM/xcreate_inputs_files<<EOF
parfile_for_benchmark
EOF

source params.in
source $SCRIPTS/scripts_specfem3D.sh
source $SCRIPTS/scripts_dsm.sh

}

function run_dsm_traction ()
{
#mkdir $DSM_tractions
cd $DSM_tractions
#echo "" >> $flog_file
#echo " compute expension coefs" >> $flog_file
#echo $(date) >> $flog_file
#make_dir_exp
#copy_input_files_exp
#compute_exp_coeff
#echo $(date) >> $flog_file
# faces
make_dir_faces
copy_input_files_faces

echo >> $flog_file
echo " FACE xmin" >> $flog_file
echo $(date) >> $flog_file
cd STXMIN
read_exp_vert inputIASP.infTra
echo $(date) >> $flog_file
fft_vert  inputIASP.infTra
echo $(date) >> $flog_file
change_format_vertical
echo $(date) >> $flog_file

echo >> $flog_file
echo "FACE xmax" >> $flog_file
echo $(date) >> $flog_file
cd ../STXMAX
read_exp_vert inputIASP.infTra
echo $(date) >> $flog_file
fft_vert inputIASP.infTra
echo $(date) >> $flog_file
change_format_vertical
echo $(date) >> $flog_file


echo >> $flog_file
echo "FACE ymin" >> $flog_file
echo $(date) >> $flog_file
cd ../STYMIN
read_exp_vert inputIASP.infTra
echo $(date) >> $flog_file
fft_vert  inputIASP.infTra
echo $(date) >> $flog_file
change_format_vertical
echo $(date) >> $flog_file

echo >> $flog_file
echo "FACE ymax" >> $flog_file
echo $(date) >> $flog_file
cd ../STYMAX
read_exp_vert inputIASP.infTra
echo $(date) >> $flog_file
fft_vert inputIASP.infTra
echo $(date) >> $flog_file
change_format_vertical
echo $(date) >> $flog_file

echo >> $flog_file
echo "FACE zmin" >> $flog_file
echo $(date) >> $flog_file
cd ../STZMIN
read_exp_zmin inputIASP.infTra
echo $(date) >> $flog_file
fft_zmin  inputIASP.infTra
echo $(date) >> $flog_file
change_format_zmin
echo $(date) >> $flog_file
cd ../..

# move output into specfem directory
move_output

}


# 1. Compute expansion coefficients
#
# NPROC : number of MPI processes
# BIN   : binary directory
# OTPION: mpirun options
#
function compute_exp_coeff ()
{
echo $NPROC
echo $OPTION
echo $BIN/xTraPSV_write_ceof_mpi_PSV
$MPIRUN $OPTION $BIN/xTraPSV_write_ceof_mpi_PSV
}


#
# clean expansion coefficients
#
function clean_exp_ceof ()
{
rm log/*
rm Displacement/*
rm Stress/*
rm ascii/*
}


#
# make directories to store expansion coefficients
#
function make_dir_exp ()
{
mkdir log
mkdir Displacement
mkdir Stress
mkdir ascii
}


#
# make directories for each face
#
function make_dir_faces ()
{
mkdir -p STXMIN/Displacement
mkdir -p STXMIN/Stress
mkdir -p STXMIN/log
mkdir -p STXMIN/out

mkdir -p STXMAX/Displacement
mkdir -p STXMAX/Stress
mkdir -p STXMAX/log
mkdir -p STXMAX/out

mkdir -p STYMIN/Displacement
mkdir -p STYMIN/Stress
mkdir -p STYMIN/log
mkdir -p STYMIN/out

mkdir -p STYMAX/Displacement
mkdir -p STYMAX/Stress
mkdir -p STYMAX/log
mkdir -p STYMAX/out

mkdir -p STZMIN/Displacement
mkdir -p STZMIN/Stress
mkdir -p STZMIN/log
mkdir -p STZMIN/out
}


#
# read expansion coefficients for vertical faces.
#
function read_exp_vert ()
{
$MPIRUN $OPTION $BIN/xTraPSV_MPI_read_ceof_vert_PSV $1
}


#
# read expansion coefficients for zmin
#
function read_exp_zmin ()
{
$MPIRUN $OPTION $BIN/xTraPSV_MPI_read_zmin_PSV $1
}


#
# fft for vertical faces
#
function fft_vert ()
{
$MPIRUN $OPTION $BIN/fft_face_vert_PSV $1
}


#
# fft for zmin
#
function fft_zmin ()
{
$MPIRUN $OPTION $BIN/TraFFT_MPI_face_zmin_PSV $1
}


#
# change format vertical
#
function change_format_vertical ()
{
$MPIRUN $OPTION $BIN/ChangeFormat
$MPIRUN $OPTION $BIN/ChangeFormat_disp
}

#
# change format zmin
#
function change_format_zmin ()
{
$MPIRUN $OPTION $BIN/ChangeFormat_zmin
$MPIRUN $OPTION $BIN/ChangeFormat_zmin_disp
}

#
# move outputs from $OUT (DSM_tractions) to $REP (Tract for SPECFEM)
#
function move_output ()
{
mkdir -p $REP

mv $OUT/STXMIN/Trac.bin $REP/tractxmin.bin
mv $OUT/STXMIN/Disp.bin $REP/velxmin.bin

mv $OUT/STXMAX/Trac.bin $REP/tractxmax.bin
mv $OUT/STXMAX/Disp.bin $REP/velxmax.bin

mv $OUT/STYMIN/Trac.bin $REP/tractymin.bin
mv $OUT/STYMIN/Disp.bin $REP/velymin.bin

mv $OUT/STYMAX/Trac.bin $REP/tractymax.bin
mv $OUT/STYMAX/Disp.bin $REP/velymax.bin

mv $OUT/STZMIN/Trac.bin $REP/tractzmin.bin
mv $OUT/STZMIN/Disp.bin $REP/velzmin.bin
}


#
# copy input file for expansion coeff calculation
#
function copy_input_files_exp ()
{
cp $IN_DSM/FrqsMpi.txt $DSM_tractions/.
cp $IN_DSM/inputIASP.infTra_for_coef $DSM_tractions/inputIASP.infTra
cp $IN_DSM/$MODELE_1D $DSM_tractions/.
cp $IN_DSM/st $DSM_tractions/.
cp $MESH/recdepth  $DSM_tractions/.
}


#
# copy  inptut files for each face
#
function copy_input_files_faces ()
{
STD=STXMIN
stf=stxmin
cp $IN_DSM/FrqsMpi.txt $STD/.
cp $IN_DSM/Double_para.txt $STD/.
cp $IN_DSM/$MODELE_1D $STD/.
cp $IN_DSM/inputIASP.infTra_${stf} $STD/inputIASP.infTra

cp $MESH/$stf $STD/.
cp $MESH/recdepth $STD/.
cp $MESH/OrigRepSpecfm $STD/.

STD=STXMAX
stf=stxmax
cp $IN_DSM/FrqsMpi.txt $STD/.
cp $IN_DSM/Double_para.txt $STD/.
cp $IN_DSM/$MODELE_1D $STD/.
cp $IN_DSM/inputIASP.infTra_${stf} $STD/inputIASP.infTra

cp $MESH/$stf $STD/.
cp $MESH/recdepth $STD/.
cp $MESH/OrigRepSpecfm $STD/.

STD=STYMIN
stf=stymin
cp $IN_DSM/FrqsMpi.txt $STD/.
cp $IN_DSM/Double_para.txt $STD/.
cp $IN_DSM/$MODELE_1D $STD/.
cp $IN_DSM/inputIASP.infTra_${stf} $STD/inputIASP.infTra

cp $MESH/$stf $STD/.
cp $MESH/recdepth $STD/.
cp $MESH/OrigRepSpecfm $STD/.

STD=STYMAX
stf=stymax
cp $IN_DSM/FrqsMpi.txt $STD/.
cp $IN_DSM/Double_para.txt $STD/.
cp $IN_DSM/$MODELE_1D $STD/.
cp $IN_DSM/inputIASP.infTra_${stf} $STD/inputIASP.infTra

cp $MESH/$stf $STD/.
cp $MESH/recdepth $STD/.
cp $MESH/OrigRepSpecfm $STD/.

STD=STZMIN
stf=stzmin
cp $IN_DSM/FrqsMpi.txt $STD/.
cp $IN_DSM/Double_para.txt $STD/.
cp $IN_DSM/$MODELE_1D $STD/.
cp $IN_DSM/inputIASP.infTra_${stf} $STD/inputIASP.infTra

cp $MESH/$stf $STD/.
cp $MESH/recdepth $STD/.
cp $MESH/OrigRepSpecfm $STD/.

}


#
##
function delete_directory_if_exist ()
{
if [ ! -d $1 ] ; then
 rm -fr $1
fi
}
#
##
function clean_and_make_dir_dsm ()
{
delete_directory_if_exist $OUT
delete_directory_if_exist $REP

mkdir -p $OUT
mkdir -p $REP
}




