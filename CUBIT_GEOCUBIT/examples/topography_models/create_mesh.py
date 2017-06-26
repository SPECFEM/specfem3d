#!/usr/bin/env python



import cubit
try:
    cubit.init([""])
except:
    pass

cfg='fullpath_myparameterfile.cfg'
o='fullpath_specfem3d_meshfiles_outputdir'
o1='fullpath_cubit_meshfiles_outputdir'

from geocubitlib import volumes, mesh_volume, exportlib
volumes.volumes(cfg)
mesh_volume.mesh(cfg)
exportlib.collect(outdir=o1)
#exportlib.e2SEM(outdir=o,listblock=[4],listflag=[-1])
exportlib.e2SEM(outdir=o)
