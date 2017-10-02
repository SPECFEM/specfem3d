# save OUTPUT_FILES/DATABASES_MPI/ from prevoius frequency group run
mkdir -p OUTPUT_FILES/$1
mv OUTPUT_FILES/DATABASES_MPI/* OUTPUT_FILES/$1/.
cp output_*txt  OUTPUT_FILES/$1/.

# get the number of processors
declare -i NPROC
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
NPROC="$NPROC-1"

# $2 : iteration
ITER=$(printf "%06d" $2)

for IP in `seq 0 1 $NPROC`;
do
  IPROC=$(printf "%06d" $IP)
  cp OUTPUT_FILES/$1/proc$IPROC'_model_vp_cr_'$ITER.bin OUTPUT_FILES/DATABASES_MPI/proc$IPROC'_model_vp_input.bin'
  cp OUTPUT_FILES/$1/proc$IPROC'_model_vs_cr_'$ITER.bin OUTPUT_FILES/DATABASES_MPI/proc$IPROC'_model_vs_input.bin'
  cp DATABASES_MPI/proc$IPROC'_rho'.bin OUTPUT_FILES/DATABASES_MPI/proc$IPROC'_model_rh_input.bin'
done
