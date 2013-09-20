#!/usr/bin/env python

import boundary_definition
import cubit2specfem3d 
import os
import sys

# put the path to your installation of the "bin" directory of the CUBIT package here
sys.path.append("/home/komatits/bin/cubit/bin")
import cubit

# let us also add the local GEOCUBIT library
sys.path.append("./geocubitlib")

# define the name of the CUBIT file to convert to SPECFEM format in the line below
#cubit.cmd('open "large_test_cpml.cub"')
cubit.cmd('open "test_cmpl_2layers.cub"')

###### This is boundary_definition.py of GEOCUBIT
#..... which extracts the bounding faces and defines them into blocks
#reload(boundary_definition)

# for CUBIT version 14 and higher (now called TRELIS) you may have to change 'face' to 'Face' in the line below
boundary_definition.entities=['face']

# this command can run in parallel or not
#boundary_definition.define_bc(boundary_definition.entities,parallel=True)
boundary_definition.define_bc(boundary_definition.entities,parallel=False)

#### Export to SPECFEM3D format using cubit2specfem3d.py of GEOCUBIT

os.system('mkdir -p MESH')
#reload(cubit2specfem3d)
cubit2specfem3d.export2SPECFEM3D('MESH') 

# all files needed by SCOTCH are now in directory MESH...

