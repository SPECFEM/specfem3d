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
function clean_and_make_dir ()
{
delete_directory_if_exist $MESH
delete_directory_if_exist OUTPUT_FILES
delete_directory_if_exist OUTPUT_FILES/DATABASES_MPI
delete_directory_if_exist DATA/DSM_tractions_for_specfem3D
delete_directory_if_exist bin

mkdir -p $MESH
mkdir -p OUTPUT_FILES/
mkdir -p OUTPUT_FILES/DATABASES_MPI/
mkdir -p DATA/DSM_tractions_for_specfem3D
mkdir bin/
}


function run_create_mesh ()
{
# fonction to create MESH for a chunk
# the output files for mesh are put in $MESH directory
#

current_dir=$(pwd)

cp ParFileMeshChunk $MESH/.
cp $IN_DSM/$MODELE_1D $MESH/.
cd $MESH
$BIN/xmesh_chunk_vm
cp $MESH/model_1D.in ../DATA/.
cd $current_dir

}


function run_create_specfem_databases ()
{

cp ParFileInterface bin/.

$BINSEM/xdecompose_mesh $NPROC $MESH OUTPUT_FILES/DATABASES_MPI/
mv Numglob2loc_elmn.txt $MESH/.

###cd bin
echo '!!!!!!!!!!!!!!!!!!!!! SCRPITS 1 !!!!!!!!!!!!!!!!'
pwd
#$MPIRUN $OPTION_SIMU /ccc/cont003/home/gen7165/wangyi/run/NEW_INVERSION_SYNTH/codes/SPECFEM3D_HYBRID/bin/xgenerate_databases
$MPIRUN $OPTION_SIMU $BINSEM/xgenerate_databases
###cd ..
}

function run_create_tractions_for_specfem ()
{
cp ParFileInterface bin/.

###cd bin
echo '!!!!!!!!!!!!!!!!!!!!! SCRPITS 2 !!!!!!!!!!!!!!!!'
pwd
$MPIRUN $OPTION_SIMU $BIN/xread_absorbing_interfaces > out_read.txt
###cd ..
cp out_read.txt bin/
}

function run_simu ()
{

###cd bin
echo '!!!!!!!!!!!!!!!!!!!!! SCRPITS 3 !!!!!!!!!!!!!!!!'
pwd
$MPIRUN $OPTION_SIMU $BINSEM/xspecfem3D > tmp_sem.out
###cd ..
cp out_read.txt bin/
echo '!!!!!!!!!!!!!!!!!!!!! SCRPITS 4 !!!!!!!!!!!!!!!!'
pwd
}

function create_movie ()
{

###cd bin

# $2 = IN (/scratch/vmonteil/BENCHMARK_COUPLAGE/chunk_15s/OUTPUT_FILES/DATABASES_MPI/)
# $3 = OUT (./movie)
# $4 istep
# $5 itmax
declare -i it itmax istep

PREFIX=$1 # (velocity_Z_it)
IN=$2 #/scratch/vmonteil/BENCHMARK_COUPLAGE/chunk_15s/OUTPUT_FILES/DATABASES_MPI/
OUT=$3 #./movie
istep=$4
itmax=$5

mkdir $OUT



it=$istep

while [ "$it" -le "$itmax" ] ; do

   if [ "$it" -lt 1000000 ]; then
     FICHIER=$PREFIX${it}
   fi;

   if [ "$it" -lt 100000 ]; then
    FICHIER=$PREFIX"0"${it}
   fi;

   if [ "$it" -lt 10000 ]; then
      FICHIER=$PREFIX"00"${it}
   fi;

   if [ "$it" -lt 1000 ]; then
     FICHIER=$PREFIX"000"${it}
   fi;

   if [ "$it" -lt 100 ]; then
      FICHIER=$PREFIX"0000"${it}
   fi;


   echo $FICHIER
   $BINSEM/xcombine_vol_data 0 $NPROC_MINUS_ONE $FICHIER $IN $OUT 0
   it="$it+$istep"

done;

tar -jcvf $OUT.tar.bz2 $OUT

cd ..

}
