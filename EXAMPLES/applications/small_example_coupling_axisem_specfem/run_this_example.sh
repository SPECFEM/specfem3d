#!/bin/bash
#
# run script to couple SPECFEM3D and AxiSEM
#
#################################################################
## USER PARAMETERS

## choose example to run
# for box that reaches the free surface (default)
Param=Param_files

# for buried box
#Param=Param_files_for_buried_box

## choose AxiSEM resolution
# AxiSEM simulation minimum period
# default would be at 10 s (as SPECFEM3D box mesh goes down to 8 s)
# (choosing a larger period will run faster)
PERIOD_MIN=25

#################################################################

echo "running example: `date`"
currentdir=`pwd`

# sets up directory structure in current example directory
echo
echo "setting up example..."
echo

# directory with specfem3d
rootdir=../../../            # e.g., /mnt/Data1/vmont/GIT/specfem3d/

# checks if executables were compiled and available
if [ ! -e $rootdir/bin/xspecfem3D ]; then
  echo "Please compile first all binaries in the root directory, before running this example..."; echo
  exit 1
fi

# cleans output files
mkdir -p OUTPUT_FILES
rm -rf OUTPUT_FILES/*

# links executables
mkdir -p bin
cd bin/
rm -f x*
for f in ../$rootdir/bin/x*; do
  #echo "link executable: $f"
  ln -s $f
done
cd ../

#
#     ALL INPUTS FILES NEEDED TO RUN THIS SCRIPT ARE IN ./Param_files or ./Param_files_for_buried_box
#
#
#
#
#
echo
echo "  running setup for    : $Param"
echo "  using minimum period : $PERIOD_MIN (s)"
echo

# sets minimum period for AxiSEM simulation mesh
file=$Param/inputs_files_for_axisem/inparam_mesh
cp -v $file $file.org
sed -i "s:^DOMINANT_PERIOD .*:DOMINANT_PERIOD     $PERIOD_MIN:" $file
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# sets minimum period for AxiSEM simulation source
file=$Param/inputs_files_for_axisem/inparam_advanced
cp -v $file $file.org
sed -i "s:^SOURCE_PERIOD .*:SOURCE_PERIOD       $PERIOD_MIN:" $file
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# sets corresponding name for AxiSEM simulation
file=$Param/inputs_files_for_axisem/inparam_basic
cp -v $file $file.org
sed -i "s:^MESHNAME .*:MESHNAME           ak135_$PERIOD_MIN:" $file
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# sets minimum frequency for interpolation
#FREQ_MIN=10
#file=$Param/inputs_files_for_axisem/reformat.par
#cp -v $file $file.org
#sed -i "s:^10:$FREQ_MIN:" $file
# checks exit code
#if [[ $? -ne 0 ]]; then exit 1; fi

echo

########################    TRACTION DATABASES GENERATION   ########################################################

# directory containing axisem modified for specfem coupling
axisem_sources=$rootdir/external_libs/AxiSEM_for_SPECFEM3D/AxiSEM_modif_for_coupling_with_specfem

# directory containing utils for coupling axisem/specfem
axisem_utils_coupling=$rootdir/external_libs/AxiSEM_for_SPECFEM3D/UTILS_COUPLING_SpecFEM

# checks compilation of axisem tools
echo "checking AxiSEM MESHER compilation in directory: $axisem_sources ";
cd $axisem_sources/MESHER
make -j4 all
# checks exit code
if [[ $? -ne 0 ]]; then echo "AxiSEM MESHER compilation failed, please check in directory: $axisem_source/MESHER ...";exit 1; fi
cd $currentdir
echo

echo "checking AxiSEM SOLVER compilation in directory: $axisem_sources ";
cd $axisem_sources/SOLVER/
make -j4 all
# checks exit code
if [[ $? -ne 0 ]]; then echo "AxiSEM SOLVER compilation failed, please check in directory: $axisem_source/SOLVER ..."; exit 1; fi
cd UTILS/
make all
# checks exit code
if [[ $? -ne 0 ]]; then echo "AxiSEM SOLVER/UTILS compilation failed, please check in directory: $axisem_source/SOLVER/UTILS ..."; exit 1; fi
cd $currentdir
echo

echo "checking AxiSEM UTILS compilation in directory: $axisem_utils_coupling ...";
cd $axisem_utils_coupling
make all
# checks exit code
if [[ $? -ne 0 ]]; then echo "AxiSEM UTILS compilation failed, please check in directory: $axisem_utils_coupling ..."; exit 1; fi
cd $currentdir
echo

# link interpolation tools
cd bin/
ln -s ../$axisem_utils_coupling/xexpand_2D_3D
ln -s ../$axisem_utils_coupling/xreformat
cd ../
echo

# ------------- copy inputs files for specfem  ------------
#cp -r $Param/DATA/* DATA/
#cp -r $Param/MESH/* MESH/
# creates symbolic links
rm -f DATA MESH
ln -s $Param/DATA/
ln -s $Param/MESH/

# ------------ create traction directory ----
mkdir -p  DATA/AxiSEM_tractions/1
rm -rf DATA/AxiSEM_tractions/1/*

nproc_specfem=`grep ^NPROC DATA/Par_file_several_proc | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# ------ CREATE ALL DIRECTORIES
mkdir -p DATABASES_MPI OUTPUT_FILES run_axisem
rm -rf DATABASES_MPI/* OUTPUT_FILES/* run_axisem/*

# ------ copy axisem sources
cp -r $axisem_sources/* run_axisem/.

# copy input files defines by user
cp -v Param_files/inputs_files_for_axisem/inparam_mesh run_axisem/MESHER/.
cp -v Param_files/inputs_files_for_axisem/inparam_* run_axisem/SOLVER/.
cp -v Param_files/inputs_files_for_axisem/*.par run_axisem/SOLVER/.
echo


# --------------- CREATE INPUTS FOR SPECFEM -----------------
echo
echo "#######################################################"
echo " 1. step: inputs for coupling w/ SPECFEM"
echo "#######################################################"
echo `date`
echo

# run internal mesher (must be in serial mode in order to use scotch decomposer)
echo "running xmeshfem3D (in serial mode)..."
cp -v DATA/Par_file_one_proc DATA/Par_file         # copy Par_file for serial mode
echo


mpirun -np 1 ./bin/xmeshfem3D
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# run scotch decomposer (serial code)
echo "running xdecompose_mesh for partitioning into $nproc_specfem slices..."
cp -v DATA/Par_file_several_proc DATA/Par_file     # now copy Par_file for mpi parallel mode
echo

./bin/xdecompose_mesh $nproc_specfem  MESH/ DATABASES_MPI/ > OUTPUT_FILES/output_decompose_mesh.txt 2>&1
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# run generate database
echo "running xgenerate_databases..."
echo

mpirun -np $nproc_specfem ./bin/xgenerate_databases
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

cp -v Numglob2loc_elmn.txt MESH/.

# info output
echo
grep "Minimum period resolved" OUTPUT_FILES/output_generate_databases.txt
grep "Maximum suggested time step" OUTPUT_FILES/output_generate_databases.txt
echo

#--------------------------FOR FURTHER DEVELOPEMENTS------------------------------------
echo
echo "#######################################################"
echo " 2. step: normals"
echo "#######################################################"
echo `date`
echo
## info output (eventually for instaseis and get_rotation_matrix.py)
# get info abouts normals
./get_normal.sh $nproc_specfem
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi


# ------------------ RUNNING AXISEM ----------------------------
echo
echo "#######################################################"
echo " 3. step: AxiSEM"
echo "#######################################################"
echo `date`
echo

AXISEM_MESH_NAME=$(awk '$1=="MESHNAME" {print $2}' run_axisem/SOLVER/inparam_basic)

# mesher
echo "running AxiSEM mesher..."
echo

cd run_axisem/MESHER

./submit_called_from_matlab.csh > output_axisem_mesher.txt 2>&1
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi
echo

./movemesh_called_from_matlab.csh $AXISEM_MESH_NAME >> output_axisem_mesher.txt 2>&1
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi
echo

# solver
echo "running AxiSEM solver..."
echo
RUN_AXI_SOLVER=1

cd ../SOLVER

cp -v ../../MESH/list_ggl* .

./add_line_number_in_points_lists.sh
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi
echo

./submit_called_from_matlab.csh $RUN_AXI_SOLVER > output_axisem_solver.txt 2>&1
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi
echo

echo
echo "#######################################################"
echo " 4. step: reconstruct 3D wavefield"
echo "#######################################################"
echo `date`
echo

# reconstruct 3D wavefield on chunk edges from 2D axisem solution
echo "running 3D wavefield reconstruction xexpand_2D_3D..."
echo

cd $RUN_AXI_SOLVER

mpirun -np $nproc_specfem $currentdir/bin/xexpand_2D_3D > output_expand_2D_3D.txt 2>&1
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# interpolation of axisem solution for specfem time step
echo "running interpolation xreformat..."
echo

mpirun -np $nproc_specfem $currentdir/bin/xreformat > output_reformat.txt 2>&1
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# back to the launching directory
cd $currentdir


########################   END TRACTION DATABASES GENERATION   ########################################################


# ------ run specfem simulation
echo
echo "#######################################################"
echo " 5. step: coupled SPECFEM simulation"
echo "#######################################################"
echo `date`
echo

echo "running solver..."
echo

mpirun -np  $nproc_specfem ./bin/xspecfem3D
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi


#------- make vtk snapshots
echo
echo "#######################################################"
echo " 6. step: VTK snapshots"
echo "#######################################################"
echo `date`
echo

NTSTEP_BETWEEN_FRAMES=`grep ^NTSTEP_BETWEEN_FRAMES DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
NSTEP=`grep ^NSTEP DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

echo "creating snapshot files..."
echo

# log file
echo "snapshots:" > OUTPUT_FILES/output_snapshots.txt
echo "" >> OUTPUT_FILES/output_snapshots.txt

for IT in `seq $NTSTEP_BETWEEN_FRAMES $NTSTEP_BETWEEN_FRAMES $NSTEP`;
do
  ./create_one_snapshot.sh $IT >> OUTPUT_FILES/output_snapshots.txt 2>&1
  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi
done

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
echo `date`


