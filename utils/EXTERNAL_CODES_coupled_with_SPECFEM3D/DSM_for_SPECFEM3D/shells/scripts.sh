#
#  scripts for computing tractions and velocity on chunk boundary
#
#
#
function run_dsm_traction ()
{
mkdir $DSM_tractions
cd $DSM_tractions

# expansion coeffs
make_dir_exp
copy_input_files_exp
compute_exp_coeff

# faces
make_dir_faces
copy_input_files_faces

cd STXMIN
read_exp_vert inputIASP.infTra
fft_vert  inputIASP.infTra
change_format_vertical

cd ../STXMAX
read_exp_vert inputIASP.infTra
fft_vert inputIASP.infTra
change_format_vertical

cd ../STYMIN
read_exp_vert inputIASP.infTra
fft_vert  inputIASP.infTra
change_format_vertical

cd ../STYMAX
read_exp_vert inputIASP.infTra
fft_vert inputIASP.infTra
change_format_vertical

cd ../STZMIN
read_exp_zmin inputIASP.infTra
fft_zmin  inputIASP.infTra
change_format_zmin

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
echo $BIN/TraPSV_write_ceof_mpi
mpirun -np $NPROC $OPTION $BIN/xTraPSV_write_ceof_mpi
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
mpirun -np $NPROC $OPTION $BIN/TraPSV_read_ceof_mpi_vert $1
}


#
# read expansion coefficients for zmin
#
function read_exp_zmin ()
{
mpirun -np $NPROC $OPTION $BIN/xTraPSV_MPI_read_zmin $1
}


#
# fft for vertical faces
#
function fft_vert ()
{
mpirun -np $NPROC $OPTION $BIN/fft_face_vert $1
}


#
# fft for zmin
#
function fft_zmin ()
{
mpirun -np $NPROC $OPTION $BIN/TraFFT_MPI_face_zmin $1
}


#
# change format vertical
#
function change_format_vertical ()
{
mpirun -np $NPROC $OPTION $BIN/ChangeFormat
mpirun -np $NPROC $OPTION $BIN/ChangeFormat_disp
}

#
# change format zmin
#
function change_format_zmin ()
{
mpirun -np $NPROC $OPTION $BIN/ChangeFormat_zmin
mpirun -np $NPROC $OPTION $BIN/ChangeFormat_zmin_disp
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
cp $IN_DSM/iasp91_dsm $DSM_tractions/.
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
cp $IN_DSM/iasp91_dsm $STD/.
cp $IN_DSM/inputIASP.infTra_${stf} $STD/inputIASP.infTra

cp $MESH/$stf $STD/.
cp $MESH/recdepth $STD/.
cp $MESH/OrigRepSpecfm $STD/.

STD=STXMAX
stf=stxmax
cp $IN_DSM/FrqsMpi.txt $STD/.
cp $IN_DSM/Double_para.txt $STD/.
cp $IN_DSM/iasp91_dsm $STD/.
cp $IN_DSM/inputIASP.infTra_${stf} $STD/inputIASP.infTra

cp $MESH/$stf $STD/.
cp $MESH/recdepth $STD/.
cp $MESH/OrigRepSpecfm $STD/.

STD=STYMIN
stf=stymin
cp $IN_DSM/FrqsMpi.txt $STD/.
cp $IN_DSM/Double_para.txt $STD/.
cp $IN_DSM/iasp91_dsm $STD/.
cp $IN_DSM/inputIASP.infTra_${stf} $STD/inputIASP.infTra

cp $MESH/$stf $STD/.
cp $MESH/recdepth $STD/.
cp $MESH/OrigRepSpecfm $STD/.

STD=STYMAX
stf=stymax
cp $IN_DSM/FrqsMpi.txt $STD/.
cp $IN_DSM/Double_para.txt $STD/.
cp $IN_DSM/iasp91_dsm $STD/.
cp $IN_DSM/inputIASP.infTra_${stf} $STD/inputIASP.infTra

cp $MESH/$stf $STD/.
cp $MESH/recdepth $STD/.
cp $MESH/OrigRepSpecfm $STD/.

STD=STZMIN
stf=stzmin
cp $IN_DSM/FrqsMpi.txt $STD/.
cp $IN_DSM/Double_para.txt $STD/.
cp $IN_DSM/iasp91_dsm $STD/.
cp $IN_DSM/inputIASP.infTra_${stf} $STD/inputIASP.infTra

cp $MESH/$stf $STD/.
cp $MESH/recdepth $STD/.
cp $MESH/OrigRepSpecfm $STD/.

}
