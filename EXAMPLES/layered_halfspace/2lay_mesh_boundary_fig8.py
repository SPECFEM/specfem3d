#!/usr/bin/env python

###########################################################################
#### TNM: This is the mesh generation, adapted from a journal file
####      specific to the settings of Komatitsch and Tromp 1999, Fig.8
####      Aug 2009
###########################################################################

import cubit
import boundary_definition
import cubit2specfem3d 

import os
import sys

cubit.cmd('reset')
cubit.cmd('brick x 134000 y 134000 z 60000')

# This seems to conflict with boundary_definition.py
# ....which leaves the model space at e.g. x=[-67,67] km
cubit.cmd('volume 1 move x 67000 y 67000 z -30000')

# create vertices for discontinuity
cubit.cmd('split curve 9  distance 3000')
cubit.cmd('split curve 10  distance 3000')
cubit.cmd('split curve 11  distance 3000')
cubit.cmd('split curve 12  distance 3000')

# create surface for interface
cubit.cmd('create surface vertex 9 10 12 11')

cubit.cmd('section volume 1 with surface 7 keep normal')
cubit.cmd('section volume 1 with surface 7 reverse')

# create vertices for auxiliary interface to allow for refinement
cubit.cmd('split curve 29  distance 9000')
cubit.cmd('split curve 31  distance 9000')
cubit.cmd('split curve 32  distance 9000')
cubit.cmd('split curve 36  distance 9000')

# create surface for buffer interface to refine BELOW the discontinuity
cubit.cmd('create surface vertex 25 26 28 27')

cubit.cmd('section volume 3 with surface 19 keep normal')
cubit.cmd('section volume 3 with surface 19 reverse')

cubit.cmd('delete volume 2 4')

cubit.cmd('merge all')
cubit.cmd('imprint all')

# Meshing the volumes
cubit.cmd('volume 3 size 3589.2')
cubit.cmd('mesh volume 3')

cubit.cmd('refine surface 8 numsplit 1 bias 1.0 depth 1')

cubit.cmd('volume 1 size 1196.4')
cubit.cmd('mesh volume 1')

cubit.cmd('volume 5 size 4785.71')
cubit.cmd('mesh volume 5')

#### End of meshing 

###### This is boundary_definition.py of GEOCUBIT
#..... which extracts the bounding faces and defines them into blocks
boundary_definition.entities=['face']
boundary_definition.define_bc(boundary_definition.entities,parallel=True)

#### Define material properties for the 3 volumes ################
cubit.cmd('#### DEFINE MATERIAL PROPERTIES #######################')
cubit.cmd('block 1 name "elastic" ')        # elastic material region
cubit.cmd('block 1 attribute count 6')
cubit.cmd('block 1 attribute index 1 1  ')      # volume 1
cubit.cmd('block 1 attribute index 2 2800 ')   # vp 
cubit.cmd('block 1 attribute index 3 1500 ')   # vs 
cubit.cmd('block 1 attribute index 4 2300 ')   # rho 
cubit.cmd('block 1 attribute index 5 6 ')       # Q_flag
cubit.cmd('block 1 attribute index 6 0 ')     # anisotropy_flag

cubit.cmd('block 2 name "elastic" ')        # elastic material region
cubit.cmd('block 2 attribute count 6')
cubit.cmd('block 2 attribute index 1 2  ')      # volume 2
cubit.cmd('block 2 attribute index 2 7500 ')
cubit.cmd('block 2 attribute index 3 4300 ')
cubit.cmd('block 2 attribute index 4 3200 ')
cubit.cmd('block 2 attribute index 5 6 ')
cubit.cmd('block 2 attribute index 6 0 ')     # anisotropy_flag

cubit.cmd('block 3 name "elastic" ')        # elastic material region
cubit.cmd('block 3 attribute count 6')
cubit.cmd('block 3 attribute index 1 3  ')      # same properties as for volume 2
cubit.cmd('block 3 attribute index 2 7500 ')
cubit.cmd('block 3 attribute index 3 4300 ')
cubit.cmd('block 3 attribute index 4 3200 ')
cubit.cmd('block 3 attribute index 5 6 ')
cubit.cmd('block 3 attribute index 6 0 ')     # anisotropy_flag

cubit.cmd('export mesh "top.e" dimension 3 overwrite')
cubit.cmd('save as "meshing.cub" overwrite')

#### Export to SESAME format using cubit2specfem3d.py of GEOCUBIT

os.system('mkdir -p MESH')
cubit2specfem3d.export2SESAME('MESH') 

# all files needed by SCOTCH are now in directory MESH




