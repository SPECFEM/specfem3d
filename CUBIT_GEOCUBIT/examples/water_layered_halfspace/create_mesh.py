#!/usr/bin/env python
import os
import cubit
try:
    cubit.init([""])
except:
    pass



cubit.cmd('reset')
cubit.cmd('brick x 67000 y 134000 z 60000')
cubit.cmd('volume 1 move x 33500 y 67000 z -30000')
cubit.cmd('brick x 67000 y 134000 z 60000')
cubit.cmd('volume 2 move x 100500 y 67000 z -30000')
cubit.cmd('merge all')

# Meshing the volumes
elementsize = 3000.0

cubit.cmd('volume 1 size '+str(elementsize))
cubit.cmd('volume 2 size '+str(elementsize))
cubit.cmd('mesh volume 1 2')


from geocubitlib import boundary_definition,exportlib

boundary_definition.define_bc(parallel=True)

os.system('mkdir -p MESH/')
exportlib.collect(outdir='MESH/')
exportlib.e2SEM(outdir='MESH/')
