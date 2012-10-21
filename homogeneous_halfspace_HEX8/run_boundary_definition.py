#!python
#!/usr/bin/env python

import cubit
import boundary_definition
import cubit2specfem3d 

import os
import sys


###### This is boundary_definition.py of GEOCUBIT
#..... which extracts the bounding faces and defines them into blocks
reload(boundary_definition)
boundary_definition.entities=['face']
boundary_definition.define_bc(boundary_definition.entities,parallel=True)


