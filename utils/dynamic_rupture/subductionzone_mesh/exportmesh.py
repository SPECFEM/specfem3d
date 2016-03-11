#!/usr/bin/env python
import cubit
import cubit2specfem3d 

import os
import sys
from save_fault_nodes_elements import *
from absorbing_boundary import *

def export_block(nb,vp,vs,rho,count=6,Q=13):
    cubit.cmd('block {0}  name "elastic {0}" '.format(nb))        # material region  
    cubit.cmd('block {0} attribute count {1}'.format(nb,count)) 
    cubit.cmd('block {0} attribute index 1 1'.format(nb))      # flag for fault side 1 
    cubit.cmd('block {0} attribute index 2 {1}'.format(nb,vp))   # vp 
    cubit.cmd('block {0} attribute index 3 {1}'.format(nb,vs))    # vs 
    cubit.cmd('block {0} attribute index 4 {1}'.format(nb,rho))   # rho 
    cubit.cmd('block {0} attribute index 5 {1}'.format(nb,Q))     # Q flag (see constants.h: #IATTENUATION_ ... ) 

def define_block_hex27(i):
    cubit.cmd('block {0}  vol {0} '.format(i))
    cubit.cmd('block {0} element type hex27'.format(i))
    cubit.cmd('reset block {0}'.format(i))
  


cubit.init([''])
#cubit.cmd('open "tpv29.cub"')
cubit.cmd('open "slab_rotate_refine.cub"')
cubit.cmd('vol all scale 1000')
########### Fault elements and nodes ###############
os.system('mkdir -p MESH') 
cubit.cmd('unmerge surf 3')
cubit.cmd('set node constraint on')
for iblock in range(1,11,1):
    define_block_hex27(iblock)
cubit.cmd('set node constraint off')
cubit.cmd('node in surf 3 move 0 0 -0.001')
cubit.cmd('node in surf 39 move 0 0 0.001')
cubit.cmd('compress all')

Au = [39]
Ad = [3]

#  FOR THE BULK (Seismic wave propagation part for SESAME)

###### This is boundary_definition.py of GEOCUBIT 
#..... which extracts the bounding faces and defines them into blocks 
xmin = [24,31]
xmax = [20]
ymin = [34,38]
ymax = [27,36]
zmax = [35,26,10,6,33,37,28,17]
zmin = [16]
entities=['face']
define_boundaries(entities,xmin,xmax,ymin,ymax,zmin,zmax)


#### Define material properties for the 2 volumes ################ 
cubit.cmd('#### DEFINE MATERIAL PROPERTIES #######################') 
for iblock in range(1,11,1):
    export_block(nb = iblock,vp=6000,vs=3464,rho=2670)
 
# Material properties in concordance with tpv5 benchmark. 
#cubit.cmd('block 1 name "elastic 1" ')        # material region  
#cubit.cmd('block 1 attribute count 6') 
#cubit.cmd('block 1 attribute index 1 1')      # flag for fault side 1 
#cubit.cmd('block 1 attribute index 2 6000')   # vp 
#cubit.cmd('block 1 attribute index 3 3464')    # vs 
#cubit.cmd('block 1 attribute index 4 2670')   # rho 
#cubit.cmd('block 1 attribute index 5 13')     # Q flag (see constants.h: #IATTENUATION_ ... ) 
#
#cubit.cmd('block 2 name "elastic 2" ')        # material region  
#cubit.cmd('block 2 attribute count 6') 
#cubit.cmd('block 2 attribute index 1 1')      # flag for fault domain 2 
#cubit.cmd('block 2 attribute index 2 6000')   # vp 
#cubit.cmd('block 2 attribute index 3 3464')    # vs 
#cubit.cmd('block 2 attribute index 4 2670')   # rho 
#cubit.cmd('block 2 attribute index 5 13')     # Q flag (see constants.h: #IATTENUATION_ ... )

 # Q flag (see constants.h: #IATTENUATION_ ... )
cubit2specfem3d.export2SPECFEM3D('MESH',hex27=True)  

faultA = fault_input(1,Au,Ad)


# all files needed by SCOTCH are now in directory MESH 

