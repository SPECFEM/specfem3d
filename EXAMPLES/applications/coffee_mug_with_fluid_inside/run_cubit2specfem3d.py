#!python
#!/usr/bin/env python

import os
import sys

# checks for path for modules
found_lib = False
for path in sys.path:
    if "geocubitlib" in path:
        found_lib = True
        break
if not found_lib: sys.path.append('../../../CUBIT_GEOCUBIT/geocubitlib')
print("path:")
for path in sys.path: print("  ",path)
print("")


import cubit
import boundary_definition
import cubit2specfem3d

# gets version string
cubit_version = cubit.get_version()
print("CUBIT/Trelis version: ",cubit_version)
print("")

###### This is boundary_definition.py of GEOCUBIT
#..... which extracts the bounding faces and defines them into blocks
#import importlib
#importlib.(boundary_definition)
#boundary_definition.entities=['face']
#boundary_definition.define_bc(boundary_definition.entities,parallel=True)


#### Export to SESAME format using cubit2specfem3d.py of GEOCUBIT
os.system('mkdir -p MESH')

import importlib
importlib.reload(cubit2specfem3d)

cubit2specfem3d.export2SPECFEM3D('MESH')

# all files needed by SCOTCH are now in directory MESH

