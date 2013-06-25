#!python
#!/usr/bin/env python
# Surendra Nadh Somala, Caltech 2012

import cubit
import cubit2specfem3d 

import os
import sys
from save_fault_nodes_elements import *
from absorbing_boundary import *

cubit.cmd('playback "TPV16.jou" ') 

os.system('mkdir -p MESH') 

Au = [8]   # A_up
Ad = [13]  # A_down
faultA = fault_input(1,Au,Ad)

xmin = [4]
xmax = [6]
ymin = [17,10]
ymax = [12,15]
zmax = [9,14]
zmin = [11,16]
entities=['face']
define_boundaries(entities,xmin,xmax,ymin,ymax,zmin,zmax)

cubit.cmd('block 1 name "elastic 1" ') 
cubit.cmd('block 1 attribute count 5')
cubit.cmd('block 1 attribute index 1 1') 
cubit.cmd('block 1 attribute index 2 6000')
cubit.cmd('block 1 attribute index 3 3464')
cubit.cmd('block 1 attribute index 4 2670')
cubit.cmd('block 1 attribute index 5 13')  

cubit.cmd('block 2 name "elastic 2" ')    
cubit.cmd('block 2 attribute count 5') 
cubit.cmd('block 2 attribute index 1 1') 
cubit.cmd('block 2 attribute index 2 6000')
cubit.cmd('block 2 attribute index 3 3464')
cubit.cmd('block 2 attribute index 4 2670')
cubit.cmd('block 2 attribute index 5 13')  
 
#### Export to SESAME format using cubit2specfem3d.py of GEOCUBIT 
 
cubit2specfem3d.export2SESAME('MESH')  
