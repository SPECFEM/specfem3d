#!/usr/bin/env python

# "create_mesh.py" is a script that generates mesh specific to homogenous halfspace example
# i.e., a uniform mesh of 134 km x 134 km x 60 km with an element size 3.75 km.
# It is not applicable to other examples.

import cubit

import os
import sys

# Creating the volumes
cubit.cmd('reset')

# single volume
cubit.cmd('brick x 134000 y 134000 z 60000')
cubit.cmd('volume 1 move x 67000 y 67000 z -30000')

# two merged volumes
#cubit.cmd('brick x 67000 y 134000 z 60000')
#cubit.cmd('volume 1 move x 33500 y 67000 z -30000')
#cubit.cmd('brick x 67000 y 134000 z 60000')
#cubit.cmd('volume 2 move x 100500 y 67000 z -30000')
#cubit.cmd('merge all')

# Meshing the volumes
elementsize = 3750.0

cubit.cmd('volume all size '+str(elementsize))
cubit.cmd('mesh volume all')

# End of meshing

#
# GEOCUBIT
#
# adds path to geocubit (if not setup yet)
sys.path.append('../../CUBIT_GEOCUBIT/')

print "path: "
print sys.path
print ""

from geocubitlib import boundary_definition,cubit2specfem3d

# bounding faces
boundary_definition.entities=['face']
boundary_definition.define_bc(boundary_definition.entities,parallel=True)

# sets the id of the volume block
# (volume block starts at id 4)
id_block = 4
print "cubit block:"
print "  volume block id = " + str(id_block)
print ""

# Define material properties
print "#### DEFINE MATERIAL PROPERTIES #######################"

# elastic material
cubit.cmd('block '+str(id_block)+' name "elastic 1" ')        # elastic material region
cubit.cmd('block '+str(id_block)+' attribute count 7')
cubit.cmd('block '+str(id_block)+' attribute index 1 1')      # flag for material: 1 for 1. material
cubit.cmd('block '+str(id_block)+' attribute index 2 2800')   # vp
cubit.cmd('block '+str(id_block)+' attribute index 3 1500')   # vs
cubit.cmd('block '+str(id_block)+' attribute index 4 2300')   # rho
cubit.cmd('block '+str(id_block)+' attribute index 5 9999.0')  # Qkappa
cubit.cmd('block '+str(id_block)+' attribute index 6 9000.0')  # Qmu
cubit.cmd('block '+str(id_block)+' attribute index 7 0')      # anisotropy_flag

# acoustic material
#cubit.cmd('block '+str(id_block)+' name "acoustic 1" ')       # acoustic material region
#cubit.cmd('block '+str(id_block)+' attribute count 4')
#cubit.cmd('block '+str(id_block)+' attribute index 1 1  ')     # material 1
#cubit.cmd('block '+str(id_block)+' attribute index 2 1480 ')  # vp
#cubit.cmd('block '+str(id_block)+' attribute index 3 0 ')      # vs
#cubit.cmd('block '+str(id_block)+' attribute index 4 1028 ')  # rho (ocean salt water density:

# creates MESH/ directory for file output
os.system('mkdir -p MESH')
cubit.cmd('export mesh "MESH/top.e" dimension 3 overwrite')
cubit.cmd('save as "MESH/meshing.cub" overwrite')

# Export to SPECFEM3D format

#note: exportlib-commands will overwrite the above material properties
# exportlib.collect(outdir='MESH/')
# exportlib.e2SEM(outdir='MESH/')

cubit2specfem3d.export2SPECFEM3D('MESH/')


# all files needed by xdecompose_mesh are now in directory MESH/
