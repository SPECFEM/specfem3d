#!/bin/bash

../../CUBIT_GEOCUBIT/GEOCUBIT.py --build_volume --mesh --cfg=homogeneous_halfspace.cfg
../../CUBIT_GEOCUBIT/GEOCUBIT.py --collect --meshfiles=MESH_GEOCUBIT/mesh_vol_0.e --export2SPECFEM3D --SEMoutput=MESH
cp MESH-default/nummaterial_velocity_file MESH/
