#
#
#
function run_pure_dsm ()
{
mkdir -p  $DSM_RUN
cd $DSM_RUN
echo "" >> $flog_file
echo " PART 2 compute expension coefs" >> $flog_file
echo $(date) >> $flog_file
make_dir_exp
copy_input_files_exp
compute_exp_coeff
echo $(date) >> $flog_file
make_dir_read
copy_input_files_read
echo >> $flog_file
echo " PART 3" >> $flog_file
echo $(date) >> $flog_file
cd ST
read_exp_vert inputIASP.infTra
echo $(date) >> $flog_file
fft_vert  inputIASP.infTra
echo $(date) >> $flog_file
change_format_vertical
echo $(date) >> $flog_file
read_dsm_out 6 1

}
function make_dir_read ()
{
mkdir -p ST/Displacement
mkdir -p ST/Stress
mkdir -p ST/log
mkdir -p ST/out
}
function make_dir_exp ()
{
mkdir log
mkdir Displacement
mkdir Stress
mkdir ascii
}
#
# copy input file for expansion coeff calculation
#
function copy_input_files_exp ()
{
cp $IN_DSM/FrqsMpi.txt $DSM_RUN/.
cp $IN_DSM/inputIASP.infTra_for_coef $DSM_RUN/inputIASP.infTra
cp $IN_DSM/iasp91_dsm $DSM_RUN/.
cp $IN_DSM/st $DSM_RUN/.
cp $MESH/recdepth  $DSM_RUN/.
}
function copy_input_files_read ()
{
STD=ST
stf=st
cp $IN_DSM/FrqsMpi.txt $STD/.
cp $IN_DSM/Double_para.txt $STD/.
cp $IN_DSM/iasp91_dsm $STD/.
cp $IN_DSM/inputIASP.infTra_${stf} $STD/inputIASP.infTra
cp $IN_DSM/$stf $STD/.

cp $MESH/recdepth $STD/.
cp $MESH/OrigRepSpecfm $STD/.
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
echo $BIN/xTraPSV_write_ceof_mpi
mpirun -np $NPROC $OPTION $BIN/xTraPSV_write_ceof_mpi
}
#
# read expansion coefficients for vertical faces.
#
function read_exp_vert ()
{
mpirun -np $NPROC $OPTION $BIN/xTraPSV_MPI_read_ceof_vert $1
}

#
# fft for vertical faces
#
function fft_vert ()
{
mpirun -np $NPROC $OPTION $BIN/fft_face_vert $1
}

#
# change format vertical
#
function change_format_vertical ()
{
mpirun -np $NPROC $OPTION $BIN/ChangeFormat
mpirun -np $NPROC $OPTION $BIN/ChangeFormat_disp
}

function read_dsm_out ()
# $1 numero station de st
# $2 numereo pronfondeur inverse de recdepth
{
$BIN/xLectureBinTrac<<EOF
Disp.bin
$1 $2
EOF
mv verifSismo sismo_$1_$2.txt
}

