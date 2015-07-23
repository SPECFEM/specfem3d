##

function delete_directory_if_exist ()
{
if [ ! -d $1 ] ; then
 rm -fr $1
fi
}

function run_create_dir_and_mesh ()
{
# fonction to create MESH for a chunk
# the output files for mesh are put in $MESH directory

delete_directory_if_exist OUTPUT_FILES/
mkdir -p OUTPUT_FILES/
mkdir -p OUTPUT_FILES/DATABASES_MPI/

delete_directory_if_exist MESH/
mkdir -p MESH/

cp ParFileMeshChunk MESH/.
cp $AxiSEM_FILE_1D_MODEL MESH/.

$SPECFEM3D_BINARY_PATH/xmeshfem3D > Step1-create_3D_chunk_mesh.out

cp MESH/model_1D.in DATA/
cp MESH/list_ggl_boundary_* ${HOME_AxiSEM}/AxiSEM_modif_for_coupling_with_specfem/SOLVER/

}

function run_create_partitioning_and_specfem_databases ()
{

$SPECFEM3D_BINARY_PATH/xdecompose_mesh $NPROC MESH/ OUTPUT_FILES/DATABASES_MPI/
mv Numglob2loc_elmn.txt MESH/.

pwd
$MPIRUN $OPTION_SIMU $SPECFEM3D_BINARY_PATH/xgenerate_databases > Step4_and_Step5-partitioning_and_specfem3d_databases.out
}

function run_specfem3d_simulation ()
{
pwd
$MPIRUN $OPTION_SIMU $SPECFEM3D_BINARY_PATH/xspecfem3D > Step8-run_specfem3d_simu.out

pwd
}

###!function create_movie ()
###!{
###!
###!###cd bin
###!
###!# $2 = IN (/scratch/vmonteil/BENCHMARK_COUPLAGE/chunk_15s/OUTPUT_FILES/DATABASES_MPI/)
###!# $3 = OUT (./movie)
###!# $4 istep
###!# $5 itmax
###!declare -i it itmax istep
###!
###!PREFIX=$1 # (velocity_Z_it)
###!IN=$2 #/scratch/vmonteil/BENCHMARK_COUPLAGE/chunk_15s/OUTPUT_FILES/DATABASES_MPI/
###!OUT=$3 #./movie
###!istep=$4
###!itmax=$5
###!
###!mkdir $OUT
###!
###!it=$istep
###!
###!while [ "$it" -le "$itmax" ] ; do
###!
###!if [ "$it" -lt 1000000 ]; then
###!     FICHIER=$PREFIX${it}
###!   fi;
###!
###!   if [ "$it" -lt 100000 ]; then
###!    FICHIER=$PREFIX"0"${it}
###!   fi;
###!
###!   if [ "$it" -lt 10000 ]; then
###!      FICHIER=$PREFIX"00"${it}
###!   fi;
###!
###!   if [ "$it" -lt 1000 ]; then
###!     FICHIER=$PREFIX"000"${it}
###!   fi;
###!
###!   if [ "$it" -lt 100 ]; then
###!      FICHIER=$PREFIX"0000"${it}
###!   fi;
###!
###!   echo $FICHIER
###!   $BINSEM/xcombine_vol_data 0 $NPROC_MINUS_ONE $FICHIER $IN $OUT 0
###!   it="$it+$istep"
###!
###!done;
###!
###!tar -jcvf $OUT.tar.bz2 $OUT
###!
###!cd ..
###!
###!}
