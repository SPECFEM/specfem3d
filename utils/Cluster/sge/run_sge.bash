#!/bin/bash
# use the normal queue unless otherwise directed
# specific queue:
#     pe.q with 130 cpus,
#     all.q with 49 cpus
# can be added to qsub command at the end -p pe.q

rm -f OUTPUT_FILES/*

d=`date`
echo "Starting compilation $d"
make clean
make xgenerate_databases
make xspecfem3D
d=`date`
echo "Finished compilation $d"

echo "Submitting job"

# mesher
qsub go_mesher_sge.bash

# solver
#qsub go_solver_sge.bash
