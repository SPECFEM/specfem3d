#!/usr/bin/env python
# Surendra Nadh Somala, Caltech 2012
from __future__ import print_function

import os
import sys

import cubit
try:
    #cubit.init([""])
    cubit.init(["-noecho","-nojournal"])
except:
    pass

# gets version string
cubit_version = cubit.get_version()
print("version: ",cubit_version)

# extracts major number
v = cubit_version.split('.')
cubit_version_major = int(v[0])
print("major version number: ",cubit_version_major)

#
# GEOCUBIT
#
sys.path.append('../../../../CUBIT_GEOCUBIT')

# in case importing menu fails due to import utilities errors to find,
# this will add the geocubitlib/ folder to the sys.path:
import geocubitlib
sys.path.append(geocubitlib.__path__[0])

print("path: ")
print(sys.path)
print("")

#import cubit2specfem3d
from geocubitlib import absorbing_boundary
from geocubitlib import save_fault_nodes_elements
from geocubitlib import cubit2specfem3d

# current work directory
cubit.cmd('pwd')

# Creating the volumes
cubit.cmd('reset')

# avoids assigning empty blocks
cubit.cmd('set duplicate block elements on')

# runs journal file
cubit.cmd('playback "TPV102.jou" ')

xmin = [9,16]
xmax = [11,13]
ymin = [3]
ymax = [5]
zmax = [8,15]
zmin = [10,14]
entities=['face']

# bounding faces
absorbing_boundary.define_boundaries(entities,xmin,xmax,ymin,ymax,zmin,zmax)

# Define material properties
print("#### DEFINE MATERIAL PROPERTIES #######################")

cubit.cmd('block 1 name "elastic 1" ')         # elastic material
cubit.cmd('block 1 attribute count 6')
cubit.cmd('block 1 attribute index 1 1')       # flag for material
cubit.cmd('block 1 attribute index 2 6000')    # vp
cubit.cmd('block 1 attribute index 3 3464')    # vs
cubit.cmd('block 1 attribute index 4 2670')    # rho
cubit.cmd('block 1 attribute index 5 13')      # Qmu
cubit.cmd('block 1 attribute index 6 0')        # anisotropy_flag

cubit.cmd('block 2 name "elastic 2" ')
cubit.cmd('block 2 attribute count 6')
cubit.cmd('block 2 attribute index 1 1')
cubit.cmd('block 2 attribute index 2 6000')
cubit.cmd('block 2 attribute index 3 3464')
cubit.cmd('block 2 attribute index 4 2670')
cubit.cmd('block 2 attribute index 5 13')
cubit.cmd('block 2 attribute index 6 0')

#### Export to SPECFEM3D format using cubit2specfem3d.py of GEOCUBIT
os.system('mkdir -p MESH')

print("")
print("exporting to SPECFEM3D-format:")
print("")
# Export to SPECFEM3D format
cubit2specfem3d.export2SPECFEM3D('MESH')

## fault surfaces
Au = [8]    # A_up
Ad = [3]    # A_down

# create fault mesh files
faultA = save_fault_nodes_elements.fault_input(1,Au,Ad)

# backup cubit
cubit.cmd('export mesh "MESH/top.e" dimension 3 overwrite')
cubit.cmd('save as "MESH/mesh.cub" overwrite')
