#!/usr/bin/env python
#############################################################
#
# script uses STL surface formats
#
# note: this script seems only to work with CUBIT version 12.2
#      
# ( try out script mesh_mount.py when using a different CUBIT version)
#
#############################################################
import cubit
import boundary_definition
import cubit2specfem3d

import os
import sys
import os.path

# time stamp
print "#" + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())

# working directory
cwd = os.getcwd()
print "#current working directory: " + str(cwd)
if cwd[len(cwd)-14:len(cwd)] != "Mount_StHelens":
  print ""
  print "#please run this script from example directory: SPECFEM3D/example/Mount_StHelens/"
  print ""

cubit.cmd('version')
cubit.cmd('reset')

os.system('rm -f topo_brick.stl topo_vol.stl topo_vol2.stl')

print "running meshing script..."
print ""
print "note: this script uses topography surface in STL format"
print "         meshing will take around 2 min"
print ""

# note: this is a workaround to use STL file formats rather than ACIS formats.
#          for our purpose to create a simple mesh, this STL formats are faster 

# uses developer commands
cubit.cmd('set developer commands on')
cubit.cmd('set import mesh tolerance 100')


#############################################################
#
# 1. step: creates temporary brick volume in STL format
#
#############################################################

# creates temporary brick volume and exports in STL format
# new brick volume (depth will become 1/2 * 20,000 = 10,000 m)
cubit.cmd('brick x 15000 y 22000 z 20000')
# moves volume to UTM coordinates of topography surface
cubit.cmd('volume 1 move x 561738. y 5116370. z 0 ')
# saves as STL body
cubit.cmd('export stl ascii "topo_brick.stl" overwrite')
# clear 
cubit.cmd('reset')
#checks if new file available
if not os.path.exists("topo_brick.stl"):
  print ""
  print "error creating new STL file topo_brick.stl, please check manually..."
  print ""
  cubit.cmd('pause')

#############################################################
#
# 2. step: imports topography surface and brick volume in STL format
#
#############################################################

# topography surface
if os.path.exists("topo.stl"):
  print "opening existing topography surface"
  # previously run, just reopen the cubit file
  cubit.cmd('import stl "topo.stl" merge stitch')
else:
  # topo surface doesn't exist yet, this creates it:
  print "reading in topography surface"
  # reads in topography points and creates sheet surface
  execfile("./read_topo.py")
  # clear 
  cubit.cmd('reset')
  # now reopen the cubit file
  cubit.cmd('import stl "topo.stl" merge stitch')

# re-opens brick in STL format
cubit.cmd('import stl "topo_brick.stl" merge stitch')
# imprints topography surfaces into brick volume, creates connected surfaces
# note: this is a develop feature
cubit.cmd('imprint all')


cubit.cmd('version')

# cubit 12.2 specific
# exports as STL only surfaces which will create new volume
cubit.cmd('export stl ascii "topo_vol.stl" surface 9 5 6 8 11 13 overwrite')

# clear 
cubit.cmd('reset')
#checks if new file available
if not os.path.exists("topo_vol.stl"):
  print ""
  print "error creating new STL file topo_vol.stl, please check manually..."
  print ""
  cubit.cmd('pause')

#############################################################
#
# 3. step: manipulate STL file to create a single volume
#
#############################################################

os.system('awk \'BEGIN{print \"solid Body_1\";}{if($0 !~ /solid/) print $0;}END{print \"endsolid Body_1\";}\' topo_vol.stl > topo_vol2.stl')
#checks if new file available
if not os.path.exists("topo_vol2.stl"):
  print ""
  print "error creating new STL file topo_vol2.stl, please check manually..."
  print ""
  cubit.cmd('pause')
  
#############################################################
#
# 4. step: import STL volume and create mesh 
#
#############################################################

# re-opens new volume in STL format
cubit.cmd('import stl "topo_vol2.stl" merge stitch')

# Meshing the volumes
elementsize = 2000.0
cubit.cmd('volume all size '+str(elementsize))
# sets meshing type
# uses a sweep algorithm in vertical (Z) direction
cubit.cmd('volume all scheme sweep Vector 0 0 1')
# initial coarse mesh
cubit.cmd('mesh volume all')
# draw/update mesh lines for visualization
# this will draw also the tripling layer mesh lines
cubit.cmd('draw volume all')
# optional smoothing to improve mesh quality (takes up to 5 min)
cubit.cmd('volume all smooth scheme condition number beta 1.0 cpu 5')
cubit.cmd('smooth volume all')
# refines global mesh
cubit.cmd('refine volume 1 numsplit 1')
cubit.cmd('draw volume all')
# optional smoothing
cubit.cmd('smooth volume all')
# refines elements at topography surface
cubit.cmd('refine surface 2 numsplit 1 bias 1.0 depth 1')
cubit.cmd('draw volume all')
# optional smoothing to improve mesh quality (takes up all 10 min if beta is 2.0)
cubit.cmd('volume all smooth scheme condition number beta 2.3 cpu 10')
cubit.cmd('smooth volume all')
cubit.cmd('draw volume all')

#### End of meshing

###### This is boundary_definition.py of GEOCUBIT
#..... which extracts the bounding faces and defines them into blocks
boundary_definition.entities=['face']
boundary_definition.define_bc(boundary_definition.entities,parallel=True)

#### Define material properties for the 3 volumes ################
cubit.cmd('#### DEFINE MATERIAL PROPERTIES #######################')


cubit.cmd('block 1 name "elastic 1" ')        # elastic material region
cubit.cmd('block 1 attribute count 6')
cubit.cmd('block 1 attribute index 1 1')      # flag for material: 1 for 1. material
cubit.cmd('block 1 attribute index 2 2800')   # vp
cubit.cmd('block 1 attribute index 3 1500')   # vs
cubit.cmd('block 1 attribute index 4 2300')   # rho
cubit.cmd('block 1 attribute index 5 150.0')  # Qmu
cubit.cmd('block 1 attribute index 6 0 ')      # anisotropy_flag

# optional saves backups
cubit.cmd('export mesh "top.e" dimension 3 overwrite')
cubit.cmd('save as "meshing.cub" overwrite')
# cleanup
os.system('rm -f topo_brick.stl topo_vol.stl topo_vol2.stl')

#### Export to SESAME format using cubit2specfem3d.py of GEOCUBIT

os.system('mkdir -p MESH')
cubit2specfem3d.export2SESAME('MESH')

# all files needed by SCOTCH are now in directory MESH

# time stamp
print "#" + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
