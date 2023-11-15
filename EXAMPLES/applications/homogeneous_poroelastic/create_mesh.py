#!/usr/bin/env python
#
# "create_mesh.py" is a script that generates mesh specific to homogenous halfspace example
# i.e., a uniform mesh of 134 km x 134 km x 60 km with an element size 3.75 km.
# It is not applicable to other examples.
from __future__ import print_function

import os
import sys

# to run this script from command line, python must have its PATH environment set such that it
# includes the path to CUBIT/Trelis (cubit.py).
#
# you can also explicitly set it here e.g. like:
#sys.path.append('/opt/Trelis-15.0/bin/')

try:
    import cubit
except ImportError:
    print("Error: Importing cubit as python module failed")
    print("could not import cubit, please check your PYTHONPATH settings...")
    print("")
    print("current path: ")
    print(sys.path)
    print("")
    print("try to include path to directory which includes file cubit.py, e.g. /opt/Trelis-15.0/bin/")
    print("")
    sys.exit("Import cubit failed")

print(sys.path)

cubit.init([""])

# gets version string
cubit_version = cubit.get_version()
print("version: ",cubit_version)

# extracts major number
v = cubit_version.split('.')
cubit_version_major = int(v[0])
print("major version number: ",cubit_version_major)

# current work directory
cubit.cmd('pwd')

# Creating the volumes
cubit.cmd('reset')

# single volume
cubit.cmd('brick x 134000 y 134000 z 60000')
cubit.cmd('volume 1 move x 67000 y 67000 z -30000')

# two merged volumes
#cubit.cmd('brick x 67000 y 134000 z 60000')
#cubit.cmd('volume 1 move x 33500 y 67000 z -30000')
#cubit.cmd('brick x 67000 y 134000 z 60000')
#cubit.cmd('volume 2 move x 100500 y 67000 z -30000')
#cubit.cmd('merge all')

# Meshing the volumes
elementsize = 3750.0

cubit.cmd('volume all size '+str(elementsize))
cubit.cmd('mesh volume all')

# End of meshing

#
# GEOCUBIT
#
# adds path to geocubit (if not setup yet)
sys.path.append('../../../CUBIT_GEOCUBIT/')

print("path: ")
print(sys.path)
print("")

# avoids assigning empty blocks
cubit.cmd('set duplicate block elements on')

# creates MESH/ directory for file output
os.system('mkdir -p MESH')

# conversion to specfem-format
# use_explicit will explicitly assign material properties as block attributes
use_explicit=1

if use_explicit == 1:
    from geocubitlib import boundary_definition
    # bounding faces
    boundary_definition.entities=['face']
    boundary_definition.define_bc(boundary_definition.entities,parallel=True)
    from geocubitlib import cubit2specfem3d
    print("")
    print("material properties: assigned as block attributes")
    print("")
    # sets the id of the volume block
    # (volume block starts at id 4)
    id_block = 1
    print("cubit block:")
    print("  volume block id = " + str(id_block))
    print("")
    # Define material properties
    print("#### DEFINE MATERIAL PROPERTIES #######################")
    # elastic material
    cubit.cmd('block '+str(id_block)+' name "poroelastic 1" ')      # poroelastic material region
    cubit.cmd('block '+str(id_block)+' attribute count 7')
    cubit.cmd('block '+str(id_block)+' attribute index 1 1')        # flag for material: 1 for 1. material
    cubit.cmd('block '+str(id_block)+' attribute index 2 3371.0')   # vp
    cubit.cmd('block '+str(id_block)+' attribute index 3 2128.0')   # vs
    cubit.cmd('block '+str(id_block)+' attribute index 4 2473.0')   # rho
    cubit.cmd('block '+str(id_block)+' attribute index 5 9999.0')   # Qkappa
    cubit.cmd('block '+str(id_block)+' attribute index 6 9999.0')   # Qmu
    cubit.cmd('block '+str(id_block)+' attribute index 7 0')        # anisotropy_flag
    print("")
    print("exporting to SPECFEM3D-format:")
    print("")
    # Export to SPECFEM3D format
    cubit2specfem3d.export2SPECFEM3D('MESH/')
    # backup cubit
    cubit.cmd('export mesh "MESH/top.e" dimension 3 overwrite')
    cubit.cmd('save as "MESH/meshing.cub" overwrite')
else:
    from geocubitlib import exportlib
    print("")
    print("exporting to SPECFEM3D-format:")
    print("")
    # Export to SPECFEM3D format
    # note: exportlib-commands will overwrite material properties
    exportlib.define_blocks(outdir='MESH/',save_cubfile=True,outfilename='top')
    exportlib.e2SEM(outdir='MESH/')
    # Define material properties
    print("#### DEFINE MATERIAL PROPERTIES #######################")
    # elastic material
    material_cfg=[{'material region':'3','id_block':'1','vp':'3371.0','vs':'2128.0','rho':'2473.0','Qkappa':'9999.0','Qmu':'9999.0','anisotropy_flag':'0'}]
    # modifies material file
    nummaterial_velocity_file='MESH/nummaterial_velocity_file'
    f=open(nummaterial_velocity_file,'w')
    for block in material_cfg:
        print(block)
        s=block['material region']+' '
        s=s+block['id_block']+' '
        s=s+block['rho']+' '
        s=s+block['vp']+' '
        s=s+block['vs']+' '
        s=s+block['Qkappa']+' '
        s=s+block['Qmu']+' '
        s=s+block['anisotropy_flag']
        f.write(s+'\n')
    f.close()

# backup cubit

# all files needed by xdecompose_mesh are now in directory MESH/
