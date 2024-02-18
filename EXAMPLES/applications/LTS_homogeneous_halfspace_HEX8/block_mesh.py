#!/usr/bin/env python
from __future__ import print_function

import os
import sys

# default directories
SEMoutput='MESH'
os.system('mkdir -p '+ SEMoutput)

import cubit
try:
    #cubit.init([""])
    cubit.init(["-noecho","-nojournal"])
except:
    pass

version = cubit.get_version()
version_major = int(version.split(".")[0])
version_minor = int(version.split(".")[1])
print("cubit version: ",version," - major: ",version_major," minor: ",version_minor)

cubit.cmd('reset')
cubit.cmd('brick x 134000 y 134000 z 60000')
cubit.cmd('volume 1 move x 67000 y 67000 z -30000')

# Meshing the volumes
elementsize = 15000.0 #3750.0

cubit.cmd('volume 1 size '+str(elementsize))
cubit.cmd('mesh volume 1')

# refines single elements at topography surface
# element at top surface in the middle of the mesh
cubit.cmd('refine hex 161 numsplit 1 bias 1.0 depth 1')

#cubit.cmd('draw volume all'
#cubit.cmd('pause')

# refines again single element at topography surface
# element at top surface in the middle of the mesh
if version_major >= 2023:
    # cubit versions >= 2023.x
    cubit.cmd('refine hex 422 numsplit 1 bias 1.0 depth 1')
else:
    # older cubit/trelis
    cubit.cmd('refine hex 489 numsplit 1 bias 1.0 depth 1')

cubit.cmd('draw volume all')



#### End of meshing

# adds path to scripts (if not setup yet)
sys.path.append('../../../CUBIT_GEOCUBIT/geocubitlib')

## obsolete:
#import boundary_definition
#import cubit2specfem3d
## new:
from geocubitlib import boundary_definition
from geocubitlib import cubit2specfem3d


###### This is boundary_definition.py of GEOCUBIT
#..... which extracts the bounding faces and defines them into blocks
boundary_definition.entities=['face']
boundary_definition.define_bc(boundary_definition.entities,parallel=True)


#### Define material properties for the 3 volumes ################
cubit.cmd('#### DEFINE MATERIAL PROPERTIES #######################')

cubit.cmd('block 1 name "elastic 1" ')        # elastic material region
cubit.cmd('block 1 attribute count 7')
cubit.cmd('block 1 attribute index 1 1')      # flag for material: 1 for 1. material
cubit.cmd('block 1 attribute index 2 2800')   # vp
cubit.cmd('block 1 attribute index 3 1500')   # vs
cubit.cmd('block 1 attribute index 4 2300')   # rho
cubit.cmd('block 1 attribute index 5 9999.0') # Qkappa
cubit.cmd('block 1 attribute index 6 9999.0') # Qmu
cubit.cmd('block 1 attribute index 7 0 ')     # anisotropy_flag

cubit.cmd('export mesh "' + SEMoutput + '/top.e" dimension 3 overwrite')
cubit.cmd('save as "' + SEMoutput + '/top.cub" overwrite')

#### Export to SPECFEM3D format using cubit2specfem3d.py of GEOCUBIT

cubit2specfem3d.export2SPECFEM3D(SEMoutput)

# screen shot
# (could crash version < 16.4)
if version_major >= 2023:
    cubit.cmd('view iso')
    cubit.cmd('hardcopy "' + SEMoutput + '/block_mesh.png" png')
elif version_major >= 16 and version_minor >= 4:
    cubit.cmd('view top')
    cubit.cmd('hardcopy "' + SEMoutput + '/block_mesh.png" png')

# all files needed by SCOTCH are now in directory MESH
