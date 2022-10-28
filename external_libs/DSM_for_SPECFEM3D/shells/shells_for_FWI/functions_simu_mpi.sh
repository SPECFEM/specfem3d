#=======================================================================================
function forward_simu ()
{

#
# arguement $1 is name of output directory where the results are
# saved
#

isrc=1
while [ "$isrc" -le "$nsrc" ]; do

  # forward simu for earthquake isrc
  cp ./DATA/Par_file_${isrc} ./DATA/Par_file
  cd bin
  $MPIRUN $OPTION_MPI $HYBRID_BINNARY/xspecfem3D > tmp.out
  cd ..

  # copy and save outputs
  mkdir -p in_out_files_${isrc}/OUTPUT_FILES/$1
  mv  OUTPUT_FILES/* in_out_files_${isrc}/OUTPUT_FILES/$1/.
  cp in_out_files_${isrc}/OUTPUT_FILES/$1/*semd in_out_files_${isrc}/OUTPUT_FILES/.

  # save direct field to prepare adjoint simu
  mv OUTPUT_FILES/DATABASES_MPI/*save_forward_arrays.bin $TRACTION/$EARTHQUAKE${isrc}/.
  mv OUTPUT_FILES/DATABASES_MPI/*_absorb_field.bin $TRACTION/$EARTHQUAKE${isrc}/.

  isrc="$isrc+1"

done
}
#=====================================================================================
function adjoint_simu ()
{
isrc=1
while [ "$isrc" -le "$nsrc" ]; do

  # copy adjoint source
  cp in_out_files_${isrc}/OUTPUT_FILES/WF/*.adj ./in_out_files/SEM/.

  # forward simu for earthquake isrc
  cp ./DATA/Par_file_adj_${isrc} ./DATA/Par_file
  cd bin
  $MPIRUN $OPTION_MPI $SEM_BINNARY/xspecfem3D > tmp.out
  cd ..

  # save kernel
  mv ./OUTPUT_FILES/DATABASES_MPI/*kernel.bin in_out_files_${isrc}/OUTPUT_FILES/WF/.
  mv ./OUTPUT_FILES/DATABASES_MPI/*vp.bin in_out_files_${isrc}/OUTPUT_FILES/WF/.
  mv ./OUTPUT_FILES/DATABASES_MPI/*vs.bin in_out_files_${isrc}/OUTPUT_FILES/WF/.

  isrc="$isrc+1"

done
}
#======================================================================================
function project_grad_in_tomo_grid ()
{

#
# argument $1 : numero d'iteration
#
isrc=1
PREFIX0="alpha_kernel"
cd bin

echo $nsrc > path_file_for_gradient_0_alpha.par
while [ "$isrc" -le "$nsrc" ]; do

  IN=../in_out_files_${isrc}/OUTPUT_FILES/WF/
  OUT=../in_out_files_${isrc}/OUTPUT_FILES/WF/
  echo $OUT >> path_file_for_gradient_0_alpha.par

  $MPIRUN $OPTION_MPI $HYBRID_BINNARY/xproject_tomo_grid 0 $SLICE $PREFIX0 $IN $OUT 0

  mv data_tomo.bin $OUT/grad_alpha_$1.bin

  isrc="$isrc+1"

done


echo $nsrc > path_file_for_gradient_0_beta.par
isrc=1
PREFIX0="beta_kernel"
while [ "$isrc" -le "$nsrc" ]; do

  IN=../in_out_files_${isrc}/OUTPUT_FILES/WF/
  OUT=../in_out_files_${isrc}/OUTPUT_FILES/WF/
  echo $OUT >> path_file_for_gradient_0_beta.par

  $MPIRUN $OPTION_MPI $HYBRID_BINNARY/xproject_tomo_grid 0 $SLICE $PREFIX0 $IN $OUT 0

  mv data_tomo.bin $OUT/grad_beta_$1.bin

  isrc="$isrc+1"

done

cd ..

}
#======================================================================================
function project_model_in_tomo_grid ()
{

#
# argument $1 : numero d'iteration
#

isrc=1
PREFIX0="vp"
cd bin
while [ "$isrc" -le "$nsrc" ]; do

  IN=../in_out_files_${isrc}/OUTPUT_FILES/WF/
  OUT=../in_out_files_${isrc}/OUTPUT_FILES/WF/
  # argument 1 at the end when model projection (0 for gradient)
  $MPIRUN $OPTION_MPI $HYBRID_BINNARY/xproject_tomo_grid 0 $SLICE $PREFIX0 $IN $OUT 1

  mv data_tomo.bin $OUT/mod_alpha_$1.bin

  isrc="$isrc+1"

done


isrc=1
PREFIX0="vs"
while [ "$isrc" -le "$nsrc" ]; do

  IN=../in_out_files_${isrc}/OUTPUT_FILES/WF/
  OUT=../in_out_files_${isrc}/OUTPUT_FILES/WF/
  # argument 1 at the end when model projection (0 for gradient)
  $MPIRUN $OPTION_MPI $HYBRID_BINNARY/xproject_tomo_grid 0 $SLICE $PREFIX0 $IN $OUT 1

  mv data_tomo.bin $OUT/mod_beta_$1.bin

  isrc="$isrc+1"

done

cd ..

}

#======================================================================================
