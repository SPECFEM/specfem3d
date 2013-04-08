#!python
#!/usr/bin/env python

import cubit
import boundary_definition
import cubit2specfem3d 

import os
import sys


###### This is boundary_definition.py of GEOCUBIT
#..... which extracts the bounding faces and defines them into blocks
#reload(boundary_definition)
#boundary_definition.entities=['face']
#boundary_definition.define_bc(boundary_definition.entities,parallel=True)


#### Export to SESAME format using cubit2specfem3d.py of GEOCUBIT
os.system('mkdir -p MESH')

reload(cubit2specfem3d)
cubit2specfem3d.export2SESAME('MESH') 

# all files needed by SCOTCH are now in directory MESH

