#!/usr/bin/env python
import os
import cubit
try:
    cubit.init([""])
except:
    pass

cfg='./tomographic_model.cfg'
o='./MESH'
o1='./cubit_mesh_files'

from geocubitlib import volumes, mesh_volume, exportlib
volumes.volumes(cfg)
mesh_volume.mesh(cfg)

os.system('mkdir -p MESH/')
exportlib.collect(outdir=o1)
exportlib.e2SEM(outdir=o,listblock=[4],listflag=[-1])
