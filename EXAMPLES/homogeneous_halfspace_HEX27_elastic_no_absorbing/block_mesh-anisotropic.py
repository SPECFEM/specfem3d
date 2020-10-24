#!/usr/bin/env python
from __future__ import print_function

import cubit
cubit.init([""])
try:
	from geocubitlib import boundary_definition
	from geocubitlib import cubit2specfem3d
except:
    import boundary_definition
	import cubit2specfem3d

import os
import sys

# two volumes separating 134000x134000x60000 block horizontally
cubit.cmd('reset')
cubit.cmd('brick x 67000 y 134000 z 60000')
cubit.cmd('volume 1 move x 33500 y 67000 z -30000')
cubit.cmd('brick x 67000 y 134000 z 60000')
cubit.cmd('volume 2 move x 100500 y 67000 z -30000')
cubit.cmd('merge all')

# Meshing the volumes
elementsize = 3750.0

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
cubit.cmd('block 1 name "elastic" ')        # elastic material region
cubit.cmd('block 1 attribute count 6')
cubit.cmd('block 1 attribute index 1 1')      # flag for material: 1 for 1. material
cubit.cmd('block 1 attribute index 2 2800')   # vp
cubit.cmd('block 1 attribute index 3 1500')   # vs
cubit.cmd('block 1 attribute index 4 2300')   # rho
cubit.cmd('block 1 attribute index 5 9000.0')  # Qmu
cubit.cmd('block 1 attribute index 6 1 ')      # anisotropy_flag

cubit.cmd('block 2 name "elastic" ')        # elastic material region
cubit.cmd('block 2 attribute count 6')
cubit.cmd('block 2 attribute index 1 1')      # flag for material: 1 for 1. material
cubit.cmd('block 2 attribute index 2 2800')   # vp
cubit.cmd('block 2 attribute index 3 1500')   # vs
cubit.cmd('block 2 attribute index 4 2300')   # rho
cubit.cmd('block 2 attribute index 5 9000.0')  # Q_mu
cubit.cmd('block 2 attribute index 6 0 ')      # anisotropy_flag


# HEX27 DEFINITION: to be done AFTER the defintion of ALL the blocks
cubit.cmd('block 1 2 element type HEX27')
cubit.cmd('block 1001 element type SHELL9')
cubit.cmd('block 1002 element type SHELL9')
cubit.cmd('block 1003 element type SHELL9')
cubit.cmd('block 1004 element type SHELL9')
cubit.cmd('block 1005 element type SHELL9')
cubit.cmd('block 1006 element type SHELL9')


cubit.cmd('export mesh "top.e" dimension 3 overwrite')
cubit.cmd('save as "meshing.cub" overwrite')

#### Export to SPECFEM3D format using cubit2specfem3d.py of GEOCUBIT

os.system('mkdir -p MESH')
cubit2specfem3d.export2SPECFEM3D('MESH',hex27=True)

# all files needed by SCOTCH are now in directory MESH
