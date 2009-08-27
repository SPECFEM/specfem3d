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
elementsize = 10000.0
cubit.cmd('volume 1 size '+str(elementsize))
cubit.cmd('mesh volume 1')


#### End of meshing 

###### This is boundary_definition.py of GEOCUBIT
#..... which extracts the bounding faces and defines them into blocks
boundary_definition.entities=['face']
boundary_definition.define_bc(boundary_definition.entities,parallel=True)

#### Define material properties for the 3 volumes ################
cubit.cmd('#### DEFINE MATERIAL PROPERTIES #######################')
cubit.cmd('block 1 attribute count 5')
cubit.cmd('block 1 attribute index 1 1') # flag 1 for 1. material properties
cubit.cmd('block 1 attribute index 2 2800')
cubit.cmd('block 1 attribute index 3 1500')
cubit.cmd('block 1 attribute index 4 2300')
cubit.cmd('block 1 attribute index 5 100000.')

cubit.cmd('export mesh "top.e" dimension 3 overwrite')
cubit.cmd('save as "meshing.cub" overwrite')

#### Export to SESAME format using cubit2specfem3d.py of GEOCUBIT

os.system('mkdir -p MESH')
cubit2specfem3d.export2SESAME('MESH') 

# all files needed by SCOTCH are now in directory MESH
