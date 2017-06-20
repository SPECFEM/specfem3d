##

function delete_directory_if_exist ()
{
if [ ! -d $1 ] ; then
 rm -fr $1
fi
}

function clean_and_make_dir ()
{
delete_directory_if_exist $MESH
delete_directory_if_exist OUTPUT_FILES
delete_directory_if_exist OUTPUT_FILES/DATABASES_MPI
delete_directory_if_exist DATA/DSM_tractions_for_specfem3D

mkdir -p $MESH
mkdir -p OUTPUT_FILES/
mkdir -p OUTPUT_FILES/DATABASES_MPI/
mkdir -p DATA/DSM_tractions_for_specfem3D
}


function run_create_mesh ()
{
# fonction to create MESH for a chunk
# the output files for mesh are put in $MESH directory

current_dir=$(pwd)

cp ParFileMeshChunk $MESH/.
cp $IN_DSM/$MODELE_1D $MESH/.
###cd $MESH
###$BIN/xmesh_chunk_vm
cd $current_dir
$BINSEM/xmeshfem3D > Step1-create_mesh.out
cp $MESH/model_1D.in DATA/
}

function run_create_specfem_databases ()
{

$BINSEM/xdecompose_mesh $NPROC $MESH OUTPUT_FILES/DATABASES_MPI/
mv Numglob2loc_elmn.txt $MESH/.

pwd
$MPIRUN $OPTION_SIMU $BINSEM/xgenerate_databases > Step3-create_specfem3d_database.out
}

function run_create_tractions_for_specfem ()
{

pwd
$MPIRUN $OPTION_SIMU $BIN/xread_absorbing_interfaces > Step4-create_tractions_for_specfem3D_from_DSM.out
}

function run_simu ()
{
pwd
$MPIRUN $OPTION_SIMU $BINSEM/xspecfem3D > Step5-run_specfem3d_simulation.out

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
