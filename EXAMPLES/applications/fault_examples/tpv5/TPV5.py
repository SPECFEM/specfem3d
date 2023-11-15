#!/usr/bin/env python
from __future__ import print_function

import math
import os
import sys

import cubit
try:
    #cubit.init([""])
    cubit.init(["-noecho","-nojournal"])
except:
    pass

# gets version string
cubit_version = cubit.get_version()
print("version: ",cubit_version)

# extracts major number
v = cubit_version.split('.')
cubit_version_major = int(v[0])
print("major version number: ",cubit_version_major)

#
# GEOCUBIT
#
# adds path to geocubit (if not setup yet)
sys.path.append('../../../../CUBIT_GEOCUBIT/')

# in case importing menu fails due to import utilities errors to find,
# this will add the geocubitlib/ folder to the sys.path:
import geocubitlib
sys.path.append(geocubitlib.__path__[0])

print("path: ")
print(sys.path)
print("")

from geocubitlib import absorbing_boundary
from geocubitlib import save_fault_nodes_elements
from geocubitlib import cubit2specfem3d

km = 1000
z_surf = 0*km

####  initializing coordinates x,y,z
x = []     # fault
y = []
z = []

xbulk = [] # bulk
ybulk = []
zbulk = []

xbulk.append(-21*km)   #x1
xbulk.append(21*km)    #x2
xbulk.append(21*km)    #x3
xbulk.append(-21*km)   #x4

ybulk.append(-21*km)  #y1
ybulk.append(-21*km)  #y2
ybulk.append(21*km)   #y3
ybulk.append(21*km)   #y4

zbulk = [z_surf]*4

x.append(-9*km) #x5
x.append(0*km)   #x6
x.append(9*km)  #x7
x.append(0*km)   #x8

y.append(0.0)       #y5
y.append(0.1)    #y6
y.append(0.0)       #y7
y.append(-0.1)   #y8

z = [z_surf]*4

# current work directory
cubit.cmd('pwd')

# Creating the volumes
cubit.cmd('reset')

####################  bulk ###########################################
for i in range(len(xbulk)):
   vert="create vertex x "+str(xbulk[i])+" y "+str(ybulk[i])+" z "+str(zbulk[i])
   cubit.cmd(vert)

################  Loading fault points profile#############################
for i in range(len(x)):
  vert="create vertex x "+str(x[i])+" y "+str(y[i])+" z "+str(z[i])
  cubit.cmd(vert)

################ creating fault domains #################################
bulk1="create curve vertex 1 2"
bulk2="create curve vertex 2 3"
bulk3="create curve vertex 3 4"
bulk4="create curve vertex 4 1"

fault_up="create curve spline vertex 5 6 7"
fault_down="create curve spline vertex 5 8 7"

cubit.cmd(bulk1)
cubit.cmd(bulk2)
cubit.cmd(bulk3)
cubit.cmd(bulk4)

cubit.cmd(fault_up)
cubit.cmd(fault_down)

surface="create surface curve 1 2 3 4 5 6"
cubit.cmd(surface)

cubit.cmd("sweep surface 1  vector 0  0 -1 distance "+str(21*km))
cubit.cmd("sweep curve 5 vector 0 0 -1 distance "+str(21*km))
cubit.cmd("sweep curve 6 vector 0 0 -1 distance "+str(21*km))

#####################################################

elementsize = 1000

cubit.cmd("imprint all")
cubit.cmd("merge all")
cubit.cmd("surface 1 size "+str(elementsize))
cubit.cmd("volume 1 size "+str(elementsize))
cubit.cmd("surface 1 scheme pave")
cubit.cmd("mesh surface 1")
cubit.cmd("mesh volume 1")
cubit.cmd("unmerge surface 2 3")

# clean up sheet bodies
cubit.cmd("delete volume 2 3")

os.system('mkdir -p MESH')

########### Fault elements and nodes ###############
# fault surfaces (up/down)
Au = [2]
Ad = [3]

faultA = save_fault_nodes_elements.fault_input(1,Au,Ad)

#  FOR THE BULK (Seismic wave propagation part for SPECFEM3D)

###### This is boundary_definition.py of GEOCUBIT
#..... which extracts the bounding faces and defines them into blocks
#entities=['face'] # this is a deprecated boundary definition function
#boundary_definition.define_bc(entities,parallel=True)

entities=['face']
absorbing_boundary.define_parallel_bc(entities) # in absorbing_boundary.py

#### Define material properties for the 2 volumes ################
cubit.cmd('#### DEFINE MATERIAL PROPERTIES #######################')

# Material properties in concordance with tpv5 benchmark.

cubit.cmd('block 1 name "elastic 1" ')        # material region
cubit.cmd('block 1 attribute count 6')
cubit.cmd('block 1 attribute index 1 1')      # flag for fault side 1
cubit.cmd('block 1 attribute index 2 6000')   # vp
cubit.cmd('block 1 attribute index 3 3464')    # vs
cubit.cmd('block 1 attribute index 4 2670')   # rho
cubit.cmd('block 1 attribute index 5 13')     # q flag (see constants.h: iattenuation_ ... )
cubit.cmd('block 1 attribute index 6 0')     # q flag (see constants.h: iattenuation_ ... )

#### Export to SPECFEM3D format using cubit2specfem3d.py of GEOCUBIT

cubit2specfem3d.export2SPECFEM3D('MESH')

# all files needed by SCOTCH are now in directory MESH

