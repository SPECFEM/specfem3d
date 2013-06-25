#!/usr/bin/env python

import cubit
import boundary_definition
import cubit2specfem3d 
import math
import os
import sys
import numarray
from save_fault_nodes_elements import *

cubit.cmd('reset')

km = 1000
y_vert=0*km
h  = 0.1#(h=0.1)
h1 = 0.1#(h1=0.1)
# dipping angle 1
d = math.pi/180
dip1=7.5*d
# dipping angle 2
dip2=15*d
# L , D and W (bulk)
L=230*km   #230*km
D=50*km    #50*km
W=160*km   #160*km

############# Domain size ############
Dy = 80*km+W
## redefine origin (input data)###

####  initializing coordinates x,y,z
x=[]     # fault
y=[]
z=[]

xbulk=[] # bulk
ybulk=[]
zbulk=[]

# Bulk points ###
xbulk.append(-10*km-L)    #x1
xbulk.append(250*km+L)    #x2  
xbulk.append(250*km+L)    #x3
xbulk.append(-10*km-L)    #x4

zbulk.append(0*km)   #y1
zbulk.append(0*km)   #y2   
zbulk.append(-50*km-D) #y3
zbulk.append(-50*km-D) #y4

ybulk=[y_vert]*4

### CRACKS ##########
x.append(0*km)                      #x5  
x.append(4*km/math.tan(dip1))       #x6  = x11 (triple joint)
x.append(x[1]+8*km/math.tan(dip2))  #x7
x.append(x[2])                      #x8  = x7 (kink point)
x.append(90*km-h1)                  #x9
x.append(90*km+h1)                  #x10 = x9+h1 (surface point, fault C)
x.append(x[1])                      #x11
x.append(227*km-h1)                 #x12
x.append(227*km+h1)                 #x13 = x12+h1 (surface point, fault A)

z.append(-30*km)     #z5
z.append(-26*km+h1)  #z6  (triple join - up)
z.append(-18*km+h1)  #z7  (kink point - up)
z.append(-18*km-h1)  #z8  (kink point - down)
z.append(0)          #z9
z.append(0)          #z10 = z9 
z.append(-26*km-h1)  #z11 (triple join - down)
z.append(0)          #z12
z.append(0)          #z13

y=[y_vert]*9
#####

####################  bulk ###########################################
for i in range(len(xbulk)): 
   vert="create vertex x "+str(xbulk[i])+" y "+str(ybulk[i])+" z "+str(zbulk[i]) 
   cubit.cmd(vert) 

################  Loading fault points profile#############################
for i in range(len(x)):
  vert="create vertex x "+str(x[i])+" y "+str(y[i])+" z "+str(z[i])
  cubit.cmd(vert)

################ creating fault domains #################################

bulk1="create curve vertex 1 9"   #c1
bulk2="create curve vertex 10 12" #c2
bulk3="create curve vertex 13 2"  #c3   
bulk4="create curve vertex 2 3"   #c4
bulk5="create curve vertex 3 4"   #c5
bulk6="create curve vertex 4 1"   #c6

cubit.cmd(bulk1)
cubit.cmd(bulk2)
cubit.cmd(bulk3)
cubit.cmd(bulk4)
cubit.cmd(bulk5)
cubit.cmd(bulk6)

fault_up_A1="create curve spline vertex 5 6"       #c7
fault_up_A2="create curve spline vertex 6 12"      #c8

fault_down_A1="create curve spline vertex 5 11"    #c9
fault_down_A2="create curve spline vertex 11 13"   #c10 

fault_up_BC1="create curve vertex 9 7"     #c11
fault_up_BC2="create curve vertex 7 6"     #c12
fault_down_BC1="create curve vertex 10 8"  #c13
fault_down_BC2="create curve vertex 8 6"   #c14
 
cubit.cmd(fault_up_A1) 
cubit.cmd(fault_up_A2) 
cubit.cmd(fault_down_A1) 
cubit.cmd(fault_down_A2) 

cubit.cmd(fault_up_BC1) 
cubit.cmd(fault_up_BC2) 
cubit.cmd(fault_down_BC1) 
cubit.cmd(fault_down_BC2) 

surface="create surface curve 1 11 12 7 9 10 3 4 5 6"
cubit.cmd(surface)

surface="create surface curve 2 13 14 8"
cubit.cmd(surface)

cubit.cmd("sweep surface 1 vector 0 1 0 distance "+str(Dy)) 
cubit.cmd("sweep surface 2 vector 0 1 0 distance "+str(Dy)) 

###  fault crack (Not necessary here) ###
# FAULT A
cubit.cmd('curve 8 7 9 10 merge off')
# FAULT B
cubit.cmd('curve 11 12  merge off')
# FAULT C
cubit.cmd('curve 13 14 merge off')
#cubit.cmd('merge tolerance 1e-3')

#####################################################
elementsize = 2000  #(2500)

# IMPRINTING 
cubit.cmd("imprint all")
# MERGING
cubit.cmd("merge all")

# Meshing coincide fault_A upper boundaries .
cubit.cmd('curve 8 10 7 9 size 2000')
cubit.cmd('curve 8 10 7 9 scheme equal')
cubit.cmd('mesh curve 8 10 7 9')
cubit.cmd("surface 13 18 size "+str(elementsize))
cubit.cmd("volume 1 2 3 size "+str(elementsize))
cubit.cmd("surface 13 18 scheme pave")
cubit.cmd("mesh surface 13 18 ")
cubit.cmd("mesh volume 1 2 ")
#cubit.cmd("unmerge surface 2 3")

########## loading cracks #######
#SAVING FAULT NODES AND ELEMENTS.
os.system('mkdir -p MESH')
########## FAULT A ##############################################################
Au = [8,15] #face_up 
Ad = [6,7]  #face_down 
faultA = fault_input(1,Au,Ad)

########## FAULT BC ##############################################################
BCu = [9,10]   #face_up 
BCd = [14,17]  #face_down 
faultBC = fault_input(2,BCu,BCd)

### Exporting the mesh to cubit.
boundary_definition.entities=['face'] 
boundary_definition.define_bc(boundary_definition.entities,parallel=True) 
 
#### Define material properties for the 2 volumes ################ 
cubit.cmd('#### DEFINE MATERIAL PROPERTIES #######################') 
 
# Material properties in concordance with tpv5 benchmark. 
 
cubit.cmd('block 1 name "elastic 1" ')        # material region  
cubit.cmd('block 1 attribute count 5') 
cubit.cmd('block 1 attribute index 1 1')      # flag for fault domain 1 
cubit.cmd('block 1 attribute index 2 5477.2')   # vp 
cubit.cmd('block 1 attribute index 3 3162.3')    # vs 
cubit.cmd('block 1 attribute index 4 3000')   # rho 
cubit.cmd('block 1 attribute index 5 13')     # Q flag (see constants.h: IATTENUATION_ ... )

# Material properties in concordance with tpv5 benchmark. 
 
cubit.cmd('block 2 name "elastic 2" ')        # material region  
cubit.cmd('block 2 attribute count 5') 
cubit.cmd('block 2 attribute index 1 1')      # flag for fault domain 2 
cubit.cmd('block 2 attribute index 2 5477.2')   # vp 
cubit.cmd('block 2 attribute index 3 3162.3')    # vs 
cubit.cmd('block 2 attribute index 4 3000')   # rho 
cubit.cmd('block 2 attribute index 5 13')     # Q flag (see constants.h: IATTENUATION_ ... )

#### Export to SESAME format using cubit2specfem3d.py of GEOCUBIT 
 
cubit2specfem3d.export2SESAME('MESH')  
 
# all files needed by SCOTCH are now in directory MESH 
