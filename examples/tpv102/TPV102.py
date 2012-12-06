#!python
#!/usr/bin/env python
# Surendra Nadh Somala, Caltech 2012

import cubit
import cubit2specfem3d 

import os
import sys
from save_fault_nodes_elements import *
from absorbing_boundary import *

cubit.cmd('playback "TPV102.jou" ') 

os.system('mkdir -p MESH') 


xmin = [9,16]
xmax = [11,13]
ymin = [3]
ymax = [5]
zmax = [8,15]
zmin = [10,14]
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

Au = [8]   # A_up
Ad = [3]  # A_down
faultA = fault_input(1,Au,Ad)
