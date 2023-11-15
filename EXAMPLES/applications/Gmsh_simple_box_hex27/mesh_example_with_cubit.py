#!/usr/bin/env python2.7
#
# creates a simple mesh using cubit
#
# note: cubit/trelis still requires python 2.x
#
from __future__ import print_function

import os
import sys

#################################################################
## Model Parameters

# model parameter setup
domain_id = 1     # 1==acoustic/ 2==elastic

# model size
size = 1000.0

# meshing the volume with element size
elementsize = 100.0

#################################################################

# to run this script from command line, python must have its PATH environment set such that it
# includes the path to CUBIT/Trelis (cubit.py).
#
# you can also explicitly set it here e.g. like:
#sys.path.append('/opt/Trelis-15.0/bin/')

# version info
python_major_version = sys.version_info[0]
python_minor_version = sys.version_info[1]
print("Python version: ","{}.{}".format(python_major_version,python_minor_version))
print("")

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

cubit.init(["-noecho","-nojournal"])
print("")

# gets version string
cubit_version = cubit.get_version()
print("Cubit version: ",cubit_version)

# extracts major number
v = cubit_version.split('.')
cubit_version_major = int(v[0])
print("major version number: ",cubit_version_major)
print("")

# current work directory
cubit.cmd('pwd')

# Creating the volumes
cubit.cmd('reset')

# single volume
cubit.cmd('brick x '+str(size)+' y '+str(size)+' z '+str(size))
cubit.cmd('volume 1 move x '+str(size/2)+' y '+str(size/2)+' z '+str(-size/2))

cubit.cmd('volume all size '+str(elementsize))
cubit.cmd('mesh volume all')

# End of meshing

#
# GEOCUBIT
#
# adds path to geocubit (if not setup yet)
try:
    import geocubitlib
except:
    sys.path.append('../../../CUBIT_GEOCUBIT/')
    print("path: ")
    print(sys.path)
    print("")
    # test again
    try:
        import geocubitlib
    except:
        print("Error could not load module geocubitlib, please check your path...")
        sys.exit(1)

# avoids assigning empty blocks
cubit.cmd('set duplicate block elements on')

# creates MESH/ directory for file output
os.system('mkdir -p MESH')

# conversion to specfem-format
# use_explicit will explicitly assign material properties as block attributes
use_explicit = 1

if use_explicit == 1:
    from geocubitlib import boundary_definition
    # sets the id of the volume block
    # (volume block starts at id 4)
    id_block = 1
    print("cubit block:")
    print("  volume block id = " + str(id_block))
    print("")
    # for hex27 elements
    if 1 == 1:
        cubit.cmd('block '+str(id_block)+' vol 1 ')
        cubit.cmd('block '+str(id_block)+' element type hex27')
        cubit.cmd('reset block '+str(id_block))
    # bounding faces
    boundary_definition.entities=['face']
    boundary_definition.define_bc(boundary_definition.entities,parallel=True)
    from geocubitlib import cubit2specfem3d
    print("")
    print("material properties: assigned as block attributes")
    print("")
    # Define material properties
    print("#### DEFINE MATERIAL PROPERTIES #######################")
    if domain_id == 1:
        # acoustic material
        cubit.cmd('block '+str(id_block)+' name "acoustic 1" ')        # acoustic material region
        cubit.cmd('block '+str(id_block)+' attribute count 7')
        cubit.cmd('block '+str(id_block)+' attribute index 1 1  ')     # material 1
        cubit.cmd('block '+str(id_block)+' attribute index 2 1480 ')   # vp (1480 == ocean)
        cubit.cmd('block '+str(id_block)+' attribute index 3 0 ')      # vs
        cubit.cmd('block '+str(id_block)+' attribute index 4 1028 ')   # rho (1028 == ocean salt water density)
        cubit.cmd('block '+str(id_block)+' attribute index 5 400.0')   # Qkappa
        cubit.cmd('block '+str(id_block)+' attribute index 6 9999.0')  # Qmu
        cubit.cmd('block '+str(id_block)+' attribute index 7 0')       # anisotropy_flag
    elif domain_id == 2:
        # elastic material
        cubit.cmd('block '+str(id_block)+' name "elastic 1" ')         # elastic material region
        cubit.cmd('block '+str(id_block)+' attribute count 7')
        cubit.cmd('block '+str(id_block)+' attribute index 1 1')       # flag for material: 1 for 1. material
        cubit.cmd('block '+str(id_block)+' attribute index 2 2800')    # vp
        cubit.cmd('block '+str(id_block)+' attribute index 3 1500')    # vs
        cubit.cmd('block '+str(id_block)+' attribute index 4 2300')    # rho
        cubit.cmd('block '+str(id_block)+' attribute index 5 400.0')  # Qkappa
        cubit.cmd('block '+str(id_block)+' attribute index 6 200.0')  # Qmu
        cubit.cmd('block '+str(id_block)+' attribute index 7 0')       # anisotropy_flag
    else:
        print("domain setup not supported yet, exiting")
        sys.exit(1)


    print("")
    print("exporting to SPECFEM3D-format:")
    print("")
    # Export to SPECFEM3D format
    # for hex27 elements
    if 1 == 1:
        cubit2specfem3d.export2SPECFEM3D('MESH/',hex27=True)
    else:
        cubit2specfem3d.export2SPECFEM3D('MESH/')
    # backup cubit
    cubit.cmd('save as "MESH/meshing.cub" overwrite')
else:
    from geocubitlib import exportlib
    print("")
    print("exporting to SPECFEM3D-format:")
    print("")
    # Export to SPECFEM3D format
    # note: exportlib-commands will overwrite material properties
    # for hex27 elements
    if 1 == 1:
        exportlib.define_blocks(outdir='MESH/',save_cubfile=True,outfilename='top',hex27=True)
        exportlib.e2SEM(outdir='MESH/',hex27=True)
    else:
        exportlib.define_blocks(outdir='MESH/',save_cubfile=True,outfilename='top')
        exportlib.e2SEM(outdir='MESH/')
    # Define material properties
    print("#### DEFINE MATERIAL PROPERTIES #######################")
    if domain_id == 1:
        # acoustic material
        material_cfg=[{'material region':'1','id_block':'1','vp':'1480','vs':'0','rho':'1028','Qkappa':'400.0','Qmu':'9999.0','anisotropy_flag':'0'}]
    elif domain_id == 2:
        # elastic material
        material_cfg=[{'material region':'2','id_block':'1','vp':'2800','vs':'1500','rho':'2300','Qkappa':'400.0','Qmu':'200.0','anisotropy_flag':'0'}]
    else:
        print("domain setup not supported yet, exiting")
        sys.exit(1)

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
print("all mesh files in directory: MESH/")
print("")
print("done")
print("")

