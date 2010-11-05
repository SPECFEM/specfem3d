#!/usr/bin/env python

import cubit
import boundary_definition
import cubit2specfem3d

import os
import sys

# two volumes separating whole 134000x134000x60000 block in middle vertically
cubit.cmd('reset')
cubit.cmd('brick x 67000 y 134000 z 60000')
cubit.cmd('volume 1 move x 33500 y 67000 z -30000')
cubit.cmd('brick x 67000 y 134000 z 60000')
cubit.cmd('volume 2 move x 100500 y 67000 z -30000')
cubit.cmd('merge all')

# Meshing the volumes
elementsize = 3000.0

cubit.cmd('volume 1 size '+str(elementsize))
cubit.cmd('volume 2 size '+str(elementsize))
cubit.cmd('mesh volume 1 2')


#### End of meshing

###### This is boundary_definition.py of GEOCUBIT
#..... which extracts the bounding faces and defines them into blocks
boundary_definition.entities=['face']
boundary_definition.define_bc(boundary_definition.entities,parallel=True)

#### Define material properties for the 3 volumes ################
cubit.cmd('#### DEFINE MATERIAL PROPERTIES #######################')
cubit.cmd('block 1 name "acoustic 1" ')        # material region
cubit.cmd('block 1 attribute count 4')
cubit.cmd('block 1 attribute index 1 1')      # flag for material: 1 for 1. material
cubit.cmd('block 1 attribute index 2 3000')   # vp
cubit.cmd('block 1 attribute index 3 0')   # vs
cubit.cmd('block 1 attribute index 4 2300')   # rho

cubit.cmd('block 2 name "elastic 1" ')        # material region
cubit.cmd('block 2 attribute count 6')
cubit.cmd('block 2 attribute index 1 2')      # flag for material: 2 for 2. material
cubit.cmd('block 2 attribute index 2 3000')   # vp
cubit.cmd('block 2 attribute index 3 0.1')   # vs
cubit.cmd('block 2 attribute index 4 2300')   # rho
cubit.cmd('block 2 attribute index 5 9000.0')     # Q_mu
cubit.cmd('block 2 attribute index 6 0 ')      # anisotropy_flag


cubit.cmd('export mesh "top.e" dimension 3 overwrite')
cubit.cmd('save as "meshing.cub" overwrite')

#### Export to SESAME format using cubit2specfem3d.py of GEOCUBIT

os.system('mkdir -p MESH')
cubit2specfem3d.export2SESAME('MESH')

# all files needed by SCOTCH are now in directory MESH
