#!/usr/bin/env python
#############################################################
#
# script uses ACIS surface formats
#
# note: this script seems to work with CUBIT version > 12.2
#          meshing takes about 15 minutes (without refinement)
#
#
#############################################################
import cubit
import boundary_definition
import cubit2specfem3d

import os
import sys
import os.path
import time

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

print "running meshing script..."
print ""
print "note: this script uses topography surface in ACIS format"
print "         meshing will take around 15 min"
print ""

# uses developer commands
cubit.cmd('set developer commands on')
cubit.cmd('set import mesh tolerance 1')

#############################################################
#
# 0. step: loading topography surface
#
#############################################################
print "#" + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print "#loading topo surface..."

# topography surface
if os.path.exists("topo.cub"):
  print "opening existing topography surface"
  # topography surface
  # previously run, just reopen the cubit file
  cubit.cmd('open "topo.cub"')
else:
  # topo surface doesn't exist yet, this creates it:
  print "reading in topography surface"
  # reads in topography points and creates sheet surface
  execfile("./read_topo.py")
  # clear 
  cubit.cmd('reset')
  # now reopen the cubit file
  cubit.cmd('open "topo.cub"')

# healing the surface...
cubit.cmd('Auto_clean volume 1 small_surfaces small_curve_size 10')
cubit.cmd('regularize volume 1')

#############################################################
#
# 1. step: creates temporary brick volume
#
#############################################################
print "#" + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print "#creating brick..."

# creates temporary brick volume 
# new brick volume (depth will become 1/2 * 20,000 = 10,000 m)
cubit.cmd('brick x 15000 y 22000 z 20000')

# moves volume to UTM coordinates of topography surface
cubit.cmd('volume 2 move x 561738. y 5116370. z 0 ')

# temporary backup
cubit.cmd('save as "topo_1.cub" overwrite')

#############################################################
#
# 2. step: creates volume with topography surface
#
#############################################################
print "#" + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print "#creating volume with topography..."

print "#" + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print "#imprinting volume, this will take around 1 min, please be patience..."
cubit.cmd('merge all')
cubit.cmd('imprint all')

# exports only surfaces which will create single volume
cubit.cmd('export acis "topo_2.acis" surface 3 10 12 14 15 9 ascii overwrite')

# backup
cubit.cmd('save as "topo_2.cub" overwrite')

#############################################################
#
# 3. step: manipulate ACIS file to create a single volume
#
#############################################################
# checks if new file available
if not os.path.exists("topo_2.acis"):
  print ""
  print "error creating new volume, please check manually..."
  print ""
  cubit.cmd('pause')
# clears workspace
cubit.cmd('reset')

# imports surfaces and merges to single volume
cubit.cmd('import acis "topo_2.acis" ascii merge_globally')

print "#" + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print "#creating new volume, this will take another 2 min..."
cubit.cmd('create volume surface all heal')

# backup
cubit.cmd('save as "topo_3.cub" overwrite')

#############################################################
#
# 4. step: create mesh 
#
#############################################################
print "#" + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print "#initial meshing..."
print "#(will take around 7 min)"

# optional: refining mesh at surface
#   
# note: refining using ACIS surface format takes a long time ... (up to 3 hours)
DO_TOPO_REFINEMENT = False

# Meshing the volumes
if DO_TOPO_REFINEMENT == False:
  elementsize = 500.0
  cubit.cmd('volume all size '+str(elementsize))
  # note: we will mesh first the topography surface, then sweep down the mesh
  # topography surface
  #cubit.cmd('control skew surface 12')
  cubit.cmd('surface 12 submap smooth off')
  cubit.cmd('surface 12 scheme submap')
  cubit.cmd('mesh surface 12')
  # propagates mesh down for whole volume
  cubit.cmd('volume 1  redistribute nodes off')
  cubit.cmd('volume 1 scheme sweep source surface 12 target surface 7')
  cubit.cmd('mesh volume all')
  # draw/update mesh lines for visualization
  # this will draw also the tripling layer mesh lines in case
  cubit.cmd('draw volume all')
else:
  # optional surface refinement
  # time stamp
  print "#" + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
  print "#refining surface mesh..."
  print "#(will take around 3 hours)"
  # starts with a crude mesh
  elementsize = 2000.0
  cubit.cmd('volume all size '+str(elementsize))
  # sets meshing type
  # explicitly sets scheme for topography surface
  cubit.cmd('surface 12 submap smooth off')
  cubit.cmd('surface 12 scheme submap')  
  # uses a sweep algorithm in vertical (Z) direction
  cubit.cmd('volume all scheme sweep Vector 0 0 1')
  # initial coarse mesh
  cubit.cmd('mesh volume all')  
  # optional smoothing to improve mesh quality (takes up to 5 min)
  cubit.cmd('volume all smooth scheme condition number beta 1.0 cpu 5')
  cubit.cmd('smooth volume all')
  # refines global mesh
  cubit.cmd('refine volume 1 numsplit 1')
  cubit.cmd('draw volume all')
  # optional smoothing
  cubit.cmd('smooth volume all')
  # refines elements at topography surface
  cubit.cmd('refine surface 12 numsplit 1 bias 1.0 depth 1')
  cubit.cmd('draw volume all')
  # optional smoothing to improve mesh quality (takes up to 10 min)
  cubit.cmd('volume all smooth scheme condition number beta 2.3 cpu 10')
  cubit.cmd('smooth volume all')
  # displays final mesh
  cubit.cmd('draw volume all')


# time stamp
print "#" + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print "#done meshing..."

# backup
cubit.cmd('save as "topo_4.cub" overwrite')

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

#### Export to SESAME format using cubit2specfem3d.py of GEOCUBIT

os.system('mkdir -p MESH')
cubit2specfem3d.export2SESAME('MESH')

# all files needed by SCOTCH are now in directory MESH

# time stamp
print "#" + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
