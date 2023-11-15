#!/usr/bin/env python
###########################################################################
#
# example of a waterlayer on top, compacted poroelastic sediment layers and a final elastic layer below
#
# (similar to the compacted sedimentary layer model from Morency et al. 2008, section 14.1)
#
###########################################################################
from __future__ import print_function

import os
import sys

# checks path for modules
found_lib = False
for path in sys.path:
    if "geocubitlib" in path:
        found_lib = True
        break
if not found_lib:
    sys.path.append('../../../CUBIT_GEOCUBIT/geocubitlib')
    sys.path.append('../../../CUBIT_GEOCUBIT')
print("path:")
for path in sys.path: print("  ",path)
print("")

## choose your element mesh size
elementsize = 300.0

# default directories
outputdir='MESH'
os.system('mkdir -p '+ outputdir)

import cubit
try:
    #cubit.init([""])
    cubit.init(["-noecho","-nojournal"])
except:
    pass

version = cubit.get_version()
version_major = int(version.split(".")[0])
version_minor = int(version.split(".")[1])
print("cubit version: ",version)

cubit.cmd('reset')

# 4-layered model
#
if 1 == 1:
    # way 1: merges bricks to create volumes
    # layer 1
    cubit.cmd('brick x 4800 y 4800 z 2400')
    cubit.cmd('volume 1 move x 2400 y 2400 z -1200')  # between [0,-2400]
    # layer 2
    cubit.cmd('brick x 4800 y 4800 z 900')
    cubit.cmd('volume 2 move x 2400 y 2400 z -2850')  # between [-2400,-3300]
    # layer 3
    cubit.cmd('brick x 4800 y 4800 z 900')
    cubit.cmd('volume 3 move x 2400 y 2400 z -3750')  # between [-3300,-4200]
    # layer 4
    cubit.cmd('brick x 4800 y 4800 z 600')
    cubit.cmd('volume 4 move x 2400 y 2400 z -4500')  # between [-4200,-4800]
else:
    # way 2: uses surfaces and webcuts to create volumes
    # main volume
    cubit.cmd('brick x 4800 y 4800 z 4800')
    cubit.cmd('volume 1 move x 2400 y 2400 z -2400')
    # creates surfaces for webcut
    # 1. surface
    cubit.cmd('create vertex 0 0 -2400')
    cubit.cmd('create vertex 4800 0 -2400')
    cubit.cmd('create vertex 4800 4800 -2400')
    cubit.cmd('create vertex 0 4800 -2400')
    cubit.cmd('create surface vertex 9 10 11 12')
    # 2. surface
    cubit.cmd('create vertex 0 0 -3300')
    cubit.cmd('create vertex 4800 0 -3300')
    cubit.cmd('create vertex 4800 4800 -3300')
    cubit.cmd('create vertex 0 4800 -3300')
    cubit.cmd('create surface vertex 13 14 15 16')
    # 3. surface
    cubit.cmd('create vertex 0 0 -4200')
    cubit.cmd('create vertex 4800 0 -4200')
    cubit.cmd('create vertex 4800 4800 -4200')
    cubit.cmd('create vertex 0 4800 -4200')
    cubit.cmd('create surface vertex 17 18 19 20')
    # creates additional volumes with webcuts
    cubit.cmd('webcut volume 1 with sheet surface 7')
    cubit.cmd('webcut volume 5 with sheet surface 8')
    cubit.cmd('webcut volume 6 with sheet surface 9')
    # deletes sheet bodies
    cubit.cmd('delete volume 2')
    cubit.cmd('delete volume 3')
    cubit.cmd('delete volume 4')

# merges surfaces
cubit.cmd('merge all')
cubit.cmd('imprint all')

# resets numbering
cubit.cmd('compress all')

# Meshing the volumes
cubit.cmd('volume all size '+str(elementsize))
cubit.cmd('mesh volume all')

#### End of meshing

## directly from geocubitlib/ folder
import boundary_definition
import cubit2specfem3d
## new:
#from geocubitlib import boundary_definition
#from geocubitlib import cubit2specfem3d

###### This is boundary_definition.py of GEOCUBIT
#..... which extracts the bounding faces and defines them into blocks
boundary_definition.entities=['face']
boundary_definition.define_bc(boundary_definition.entities,parallel=True)

#### Define material properties for the 3 volumes ################
cubit.cmd('#### DEFINE MATERIAL PROPERTIES #######################')
cubit.cmd('block 1 name "acoustic 1" ')          # acoustic material region
cubit.cmd('block 1 attribute count 4')
cubit.cmd('block 1 attribute index 1 1  ')       # material 1
cubit.cmd('block 1 attribute index 2 1500 ')     # vp
cubit.cmd('block 1 attribute index 3 0 ')        # vs
cubit.cmd('block 1 attribute index 4 1020 ')     # rho

cubit.cmd('block 2 name "poroelastic 1" ')       # poroelastic material region (will be fully defined in nummaterial_poroelastic_file)
cubit.cmd('block 2 attribute count 7')
cubit.cmd('block 2 attribute index 1 2  ')       # material 2
cubit.cmd('block 2 attribute index 2 3200 ')     # vp
cubit.cmd('block 2 attribute index 3 1800 ')     # vs
cubit.cmd('block 2 attribute index 4 3000 ')     # rho
cubit.cmd('block 2 attribute index 5 9999.0')    # Q_kappa
cubit.cmd('block 2 attribute index 6 9999.0')    # Q_mu
cubit.cmd('block 2 attribute index 7 0 ')        # anisotropy_flag

cubit.cmd('block 3 name "poroelastic 2" ')       # poroelastic material region (will be fully defined in nummaterial_poroelastic_file)
cubit.cmd('block 3 attribute count 7')
cubit.cmd('block 3 attribute index 1 3  ')       # material 3
cubit.cmd('block 3 attribute index 2 3200 ')     # vp
cubit.cmd('block 3 attribute index 3 1800 ')     # vs
cubit.cmd('block 3 attribute index 4 3000 ')     # rho
cubit.cmd('block 3 attribute index 5 9999.0')    # Q_kappa
cubit.cmd('block 3 attribute index 6 9999.0')    # Q_mu
cubit.cmd('block 3 attribute index 7 0 ')        # anisotropy_flag

cubit.cmd('block 4 name "elastic 1" ')           # elastic material region
cubit.cmd('block 4 attribute count 7')
cubit.cmd('block 4 attribute index 1 4  ')       # material 4
cubit.cmd('block 4 attribute index 2 3399 ')     # vp
cubit.cmd('block 4 attribute index 3 1963 ')     # vs
cubit.cmd('block 4 attribute index 4 3200 ')     # rho
cubit.cmd('block 4 attribute index 5 9999.0')    # Q_kappa
cubit.cmd('block 4 attribute index 6 9999.0')    # Q_mu
cubit.cmd('block 4 attribute index 7 0 ')        # anisotropy_flag

# backup cubit mesh file
cubit.cmd('save as "' + outputdir + '/cubit_mesh.cub" overwrite')

#### Export to SPECFEM3D format using cubit2specfem3d.py of GEOCUBIT
cubit2specfem3d.export2SPECFEM3D(outputdir)

# screen shot
# (could crash version < 16.4
if version_major >= 16 and version_minor >= 4:
    cubit.cmd('view top')
    # needs to have a window -> only works when running this script from within Trelis/Cubit
    cubit.cmd('hardcopy "' + outputdir + '/cubit_mesh.png" png')

# all files needed by SCOTCH are now in directory MESH




