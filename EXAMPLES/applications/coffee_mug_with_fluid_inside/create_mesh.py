#!/usr/bin/env python

from __future__ import print_function

import os
import sys

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

# allows for elements in multiple block definitions (topo, absorbing faces,..)
cubit.cmd('set duplicate block elements on')

# mesh
#cubit.cmd('open "coffee_cup_with_handle_CUBIT_mesh.cub"')

## block definitions
#
# works with coffee_cup_with_handle_CUBIT_mesh.cub
cubit.cmd('delete block all')

# acoustic volume
if version_major <= 15:
    # older cubit versions
    cubit.cmd('block 1 hex in volume 2')
else:
    # version 16
    cubit.cmd('block 1 add hex in volume 2')
cubit.cmd('block 1 name "acoustic 1" ')        # acoustic material region
cubit.cmd('block 1 attribute count 4')
cubit.cmd('block 1 attribute index 1 1  ')     # material 1
cubit.cmd('block 1 attribute index 2 1400 ')   # vp
cubit.cmd('block 1 attribute index 3 0 ')      # vs
cubit.cmd('block 1 attribute index 4 1020 ')   # rho (water density)

# elastic volumes
if version_major <= 15:
    # older cubit versions
    cubit.cmd('block 2 hex in volume 1 3 4 5')
else:
    # version 16
    cubit.cmd('block 2 add hex in volume 1 3 4 5')
cubit.cmd('block 2 name "elastic 1" ')         # elastic material region
cubit.cmd('block 2 attribute count 7')
cubit.cmd('block 2 attribute index 1 2  ')     # material 2
cubit.cmd('block 2 attribute index 2 2300 ')   # vp
cubit.cmd('block 2 attribute index 3 1500 ')   # vs
cubit.cmd('block 2 attribute index 4 2800 ')   # rho
cubit.cmd('block 2 attribute index 5 9999.0')  # Q_kappa
cubit.cmd('block 2 attribute index 6 9999.0')  # Q_mu
cubit.cmd('block 2 attribute index 7 0 ')      # anisotropy_flag

# custom free surface (block id 1000)
if version_major <= 15:
    # older cubit versions
    cubit.cmd('block 1000 face in surface 1 2 4 6 7 9 11 12')
else:
    # version 16
    cubit.cmd('block 1000 add face in surface 1 2 4 6 7 9 11 12')
cubit.cmd('block 1000 name "free_or_absorbing_surface_file_zmax" ')

# saves updated mesh file
cubit.cmd('export mesh "cubit_mesh.e" dimension 3 overwrite')
cubit.cmd('save as "cubit_mesh.cub" overwrite')





