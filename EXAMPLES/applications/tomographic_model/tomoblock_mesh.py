#!/usr/bin/env python
from __future__ import print_function

import os
import sys

import cubit
cubit.init([""])

# Creating the volumes
cubit.cmd('reset')
cubit.cmd('brick x 134000 y 134000 z 60000')
cubit.cmd('volume 1 move x 67000 y 67000 z -30000')


# Meshing the volumes
elementsize = 3750.0

cubit.cmd('volume 1 size '+str(elementsize))
cubit.cmd('mesh volume 1')


#
# GEOCUBIT
#
# adds path to geocubit (if not setup yet)
sys.path.append('../../../CUBIT_GEOCUBIT/')

print("path: ")
print(sys.path)
print("")

try:
    from geocubitlib import boundary_definition
    from geocubitlib import cubit2specfem3d
except:
    import boundary_definition
    import cubit2specfem3d

# bounding faces
boundary_definition.entities=['face']
boundary_definition.define_bc(boundary_definition.entities,parallel=True)

# sets the id of the volume block
# (volume block starts at id 4)
id_block = 4
print("cubit block:")
print("  volume block id = " + str(id_block))
print("")

# Define material properties
print("#### DEFINE MATERIAL PROPERTIES #######################")

# elastic model
cubit.cmd('block '+str(id_block)+' name "elastic tomography_model.xyz 1" ')        # elastic material region
cubit.cmd('block '+str(id_block)+' attribute count 2')
cubit.cmd('block '+str(id_block)+' attribute index 1 -1')      # flag for material: -1 for 1. undefined material
cubit.cmd('block '+str(id_block)+' attribute index 2 2')      # flag for tomographic model

# acoustic model
#cubit.cmd('block '+str(id_block)+' name "acoustic tomographic 1" ')       # acoustic material region
#cubit.cmd('block '+str(id_block)+' attribute count 2')
#cubit.cmd('block '+str(id_block)+' attribute index 1 -1')     # material 1
#cubit.cmd('block '+str(id_block)+' attribute index 2 2')      # tomographic model flag

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
