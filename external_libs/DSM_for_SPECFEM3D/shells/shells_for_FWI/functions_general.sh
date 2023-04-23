#============================================================================
function create_directories ()
{
# WORKING_DIRECTORY=$(pwd)
# source ./global_parameters.in

# create directories

mkdir -p $WORKING_DIRECTORY/DATA
mkdir -p $WORKING_DIRECTORY/bin
mkdir -p $WORKING_DIRECTORY/OUTPUT_FILES/DATABASES_MPI
mkdir -p $WORKING_DIRECTORY/in_out_files/SEM
mkdir -p $WORKING_DIRECTORY/OUTPUT_FILES
mkdir -p $WORKING_DIRECTORY/OUTPUT_FILES/DATABASES_MPI_CURRENT

cp $PAR_FILE_DIRECTORY/* DATA/.

# create earthquake directories
isrc=1
while [ "$isrc" -le "$nsrc" ]; do
directory_to_create=$WORKING_DIRECTORY/in_out_files_${isrc}
mkdir -p $directory_to_create/DATABASES_MPI
mkdir -p $directory_to_create/OUTPUT_FILES/WF
isrc="$isrc+1"
done
}
#===========================================================================
function clean_directories ()
{
\rm -r bin DATA in_out_files
# create earthquake directories
isrc=1
while [ "$isrc" -le "$nsrc" ]; do
directory_to_remove=$WORKING_DIRECTORY/in_out_files_${isrc}
\rm -r  $directory_to_remove
isrc="$isrc+1"
done
}
#============================================================================
function copy_initial_database ()
{
cp $DATABASES_INIT_MODEL/* OUTPUT_FILES/DATABASES_MPI/.
}
#============================================================================
function copy_current_database ()
{
cp ./OUTPUT_FILES/DATABASES_MPI_CURRENT/* OUTPUT_FILES/DATABASES_MPI/.
}
#=============================================================================
function save_mv_current_database ()
{
mv ./OUTPUT_FILES/DATABASES_MPI/* ./OUTPUT_FILES/DATABASES_MPI_CURRENT/.
}
#=============================================================================
function save_cp_current_database ()
{
cp ./OUTPUT_FILES/DATABASES_MPI/* ./OUTPUT_FILES/DATABASES_MPI_CURRENT/.
}
#============================================================================
function initialise_process ()
{
copy_initial_database
cp $PAR_FILE_DIRECTORY/Par_for_projection_in_grid_tomo.par bin/.
rm log.linserach
rm bin/ajsutement.txt
rm bin/ajustement_t.txt
rm bin/ajustement_0.txt
rm bin/pas.txt
rm bin/tg.txt
rm bin/td.txt
rm bin/t_guess.txt
rm bin/derivee_cout.txt
}
#============================================================================
function continue_process ()
{
copy_current_database
cp $PAR_FILE_DIRECTORY/Par_for_projection_in_grid_tomo.par bin/.
rm bin/ajsutement.txt
rm bin/ajustement_t.txt
rm bin/ajustement_0.txt
rm bin/pas.txt
rm bin/tg.txt
rm bin/td.txt
rm bin/t_guess.txt
rm bin/derivee_cout.txt
}


