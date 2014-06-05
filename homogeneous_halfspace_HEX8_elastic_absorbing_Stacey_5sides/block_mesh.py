#!python
#!/usr/bin/env python

# block_mesh.py is a script that generates mesh specific to homogenous halfspace example i.e., a uniform mesh of 134 km x 134 km x 60 km with an element size 3.75 km.
# It is not applicable to other examples.

import cubit
import boundary_definition
import cubit2specfem3d

import os
import sys

cubit.cmd('reset')
cubit.cmd('brick x 134000 y 134000 z 60000')
cubit.cmd('volume 1 move x 67000 y 67000 z -30000')

elementsize = 3750.0

cubit.cmd('volume 1 size '+str(elementsize))
cubit.cmd('mesh volume 1')

boundary_definition.entities=['face']
boundary_definition.define_bc(boundary_definition.entities,parallel=True)

cubit.cmd('block 4 name "elastic 1" ') # elastic material region
cubit.cmd('block 4 attribute count 7')
cubit.cmd('block 4 attribute index 1 1') # flag for material: 1 for 1. material
cubit.cmd('block 4 attribute index 2 2300') # rho
cubit.cmd('block 4 attribute index 3 2800') # vp
cubit.cmd('block 4 attribute index 4 1500') # vs
cubit.cmd('block 4 attribute index 5 9000.0') # Qkappa
cubit.cmd('block 4 attribute index 6 9000.0') # Qmu
cubit.cmd('block 4 attribute index 7 0 ') # anisotropy_flag

os.system('mkdir -p MESH')
cubit2specfem3d.export2SPECFEM3D('MESH')

