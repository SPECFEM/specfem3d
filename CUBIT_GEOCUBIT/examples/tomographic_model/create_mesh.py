#!/usr/bin/env python



import cubit
try:
    cubit.init([""])
except:
    pass

cfg='./tomographic_model//tomographic_model.cfg'
o='./tomographic_model/specfem3d_mesh_files'
o1='./tomographic_model/cubit_mesh_files'

from geocubitlib import volumes, mesh_volume, exportlib
volumes.volumes(cfg)
mesh_volume.mesh(cfg)
exportlib.collect(outdir=o1)
exportlib.e2SEM(outdir=o,listblock=[4],listflag=[-1])