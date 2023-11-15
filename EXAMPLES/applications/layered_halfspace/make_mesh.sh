#!/bin/bash
#
#
# Temporary instructions
#
# 1. set path to GEOCUBIT base directory (../../../CUBIT_GEOCUBIT/), for example:
#    export PYTHONPATH=$PYTHONPATH:/import/c/w/tape/3D/SPECFEM3D/CUBIT_GEOCUBIT
#    export PATH=$PATH:/import/c/w/tape/3D/SPECFEM3D/CUBIT_GEOCUBIT
#
#    check path:
#    which GEOCUBIT.py
#
# 2. run this script to generate mesh
#    ./make_mesh.sh
#

# checks if your GEOCUBIT.py is already in your path (which GEOCUBIT.py)
if ! [ -x "$(command -v GEOCUBIT.py)" ]; then
geocubit=../../../CUBIT_GEOCUBIT/GEOCUBIT.py
else
geocubit=GEOCUBIT.py
fi

# meshing
echo
echo "$geocubit --build_volume --mesh --cfg=layered_halfspace_tripling.cfg"
echo
$geocubit --build_volume --mesh --cfg=layered_halfspace_tripling.cfg


echo
echo "$geocubit --meshfiles=MESH_GEOCUBIT/mesh_vol_0.e --export2SPECFEM3D --SEMoutput=MESH"
echo

$geocubit --meshfiles=MESH_GEOCUBIT/mesh_vol_0.e --export2SPECFEM3D --SEMoutput=MESH

echo
echo

cp -v MESH-default/nummaterial_velocity_file.reference MESH/nummaterial_velocity_file


