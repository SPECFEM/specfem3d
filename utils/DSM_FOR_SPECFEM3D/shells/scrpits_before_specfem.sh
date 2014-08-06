#
#
#
function create_out_directory_specfem ()
{
mkdir -p OUTPUT_FILES
mkdir -p OUTPUT_FILES/DATABASES_MPI
mkdir bin
cp ParFileInterface bin/.
}

#
#
#
function link_binary_specfem ()
{
cd $bin_local_specfem
ln -s $BIN_SPECFEM/xgenerate_databases_model_1D
ln -s $BIN_SPECFEM/xdecompose_mesh_SCOTCH
cd ..
}


#
#
#
function decompose_scotch ()
{
/bin/xdecompose_mesh_SCOTCH $NPROC $MESH OUTPUT_FILES/DATABASES_MPI/
}


#
#
#
function generate_database ()
{
cd $bin_local_specfem
mpirun -np $NPROC $OPTION ./xgenerate_databases_model_1D
cd ..
}
