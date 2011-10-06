#!/usr/bin/env python

import cubit
import boundary_definition
import cubit2specfem3d

import os
import sys

cubit.cmd('reset')
cubit.cmd('brick x 134000 y 134000 z 60000')
cubit.cmd('volume 1 move x 67000 y 67000 z -30000')


# Meshing the volumes
elementsize = 3750.0

cubit.cmd('volume 1 size '+str(elementsize))
cubit.cmd('mesh volume 1')


#### End of meshing

###### This is boundary_definition.py of GEOCUBIT
#..... which extracts the bounding faces and defines them into blocks
boundary_definition.entities=['face']
boundary_definition.define_bc(boundary_definition.entities,parallel=True)

#### Define material properties for the volume ################
cubit.cmd('#### DEFINE MATERIAL PROPERTIES #######################')


cubit.cmd('block 1 name "poroelastic 1" ')     # porous material region
cubit.cmd('block 1 attribute count 4')
cubit.cmd('block 1 attribute index 1 1')      # flag for material: 1 for 1. material
cubit.cmd('block 1 attribute index 2 3677 ')  # vp-dummy
cubit.cmd('block 1 attribute index 3 2378 ')  # vs-dummy
cubit.cmd('block 1 attribute index 4 1908 ')  # rho-dummy

cubit.cmd('export mesh "top.e" dimension 3 overwrite')
cubit.cmd('save as "meshing.cub" overwrite')

#### Export to SESAME format using cubit2specfem3d.py of GEOCUBIT

os.system('mkdir -p MESH')
cubit2specfem3d.export2SESAME('MESH')

# all files needed by SCOTCH are now in directory MESH
