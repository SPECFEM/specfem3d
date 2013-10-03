#!/bin/bash

GEOCUBIT.py --build_volume --mesh --cfg=homogeneous_halfspace.cfg
GEOCUBIT.py --collect --meshfiles=MESH_GEOCUBIT/mesh_vol_0.e --export2SPECFEM3D --SEMoutput=MESH
cp MESH-default/nummaterial_velocity_file MESH/