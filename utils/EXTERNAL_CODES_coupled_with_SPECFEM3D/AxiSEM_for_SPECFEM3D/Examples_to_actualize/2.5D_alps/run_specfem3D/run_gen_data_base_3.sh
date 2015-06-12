source ./set_path_and_params.sh


cp MESH/model_1D.in DATA/.  # important copy the file produced by mesher
mpirun -np $NPROC $SEMBIN/xgenerate_databases
