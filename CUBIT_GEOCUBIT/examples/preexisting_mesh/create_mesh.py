#!/usr/bin/env python



import cubit
try:
    cubit.init([""])
except:
    pass

f='./cubit_mesh_files/mesh.e'
o='./specfem3d_mesh_files'


cubit.cmd('import mesh "'+f+'" block all  no_geom ')




from geocubitlib import  exportlib

exportlib.e2SEM(outdir=o)
