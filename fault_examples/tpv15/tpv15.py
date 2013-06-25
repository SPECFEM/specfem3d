#!/usr/bin/env python
# Percy. Script for TPV14-TPV15 SCEC benchmarks. 

import cubit
import cubit2specfem3d 
import math
import os
import sys
from save_fault_nodes_elements import *
from absorbing_boundary import *

cubit.cmd('reset')

km = 1000
z_surf=0*km
R = math.pi/180
h = 0.1 
q = math.sqrt(3)
####  initializing coordinates x,y,z
x=[]     # fault
y=[]
z=[]

xbulk=[] # bulk
ybulk=[]
zbulk=[]

xbulk.append(-18*km)   #x1
xbulk.append(12*km)    #x2
xbulk.append(12*km)    #x3
xbulk.append(-18*km)   #x4

ybulk.append(-8*km)  #y1
ybulk.append(-8*km)  #y2
ybulk.append(8*km)   #y3
ybulk.append(8*km)   #y4

zbulk=[z_surf]*4

### Main Fault #######################

x.append(-16*km) #x5
x.append(-50/math.tan(R*30))   #x6
x.append(12*km)  #x7
x.append(-50/math.tan(R*30))   #x8

y.append(0.0)    #y5
y.append(h)      #y6
y.append(0.0)    #y7
y.append(-h)     #y8

### Branch Fault ######################

x.append(x[3])                        #x9 = x8  (triple joint)
x.append((12*km/q)*math.cos(R*30))         #x10
x.append((24*km/q)*math.cos(R*30))        #x11
x.append((12*km/q)*math.cos(R*30))         #x12

y.append(y[3])                         #y9 = y8 (triple joint)
y.append(-50-(12*km/q)*math.sin(R*30)+h)   #y10
y.append(-50-(24*km/q)*math.sin(R*30))    #y11
y.append(-50-(12*km/q)*math.sin(R*30)-h)   #y12

z=[z_surf]*12

####################  bulk ###########################################
for i in range(len(xbulk)): 
   vert="create vertex x "+str(xbulk[i])+" y "+str(ybulk[i])+" z "+str(zbulk[i]) 
   cubit.cmd(vert) 

################  Loading fault points profile#############################
for i in range(len(x)):
  vert="create vertex x "+str(x[i])+" y "+str(y[i])+" z "+str(z[i])
  cubit.cmd(vert)


################ creating fault domains #################################
bulk1="create curve vertex 1 2"   # c1
bulk2="create curve vertex 2 11"  # c2
bulk3="create curve vertex 11 7"  # c3
bulk4="create curve vertex 7 3"   # c4
bulk5="create curve vertex 3 4"   # c5
bulk6="create curve vertex 4 1"   # c6

cubit.cmd(bulk1)
cubit.cmd(bulk2)
cubit.cmd(bulk3)
cubit.cmd(bulk4)
cubit.cmd(bulk5)
cubit.cmd(bulk6)

#### Main Fault ############################## 
fault_up_r="create curve spline vertex 5 6 "    #c7
fault_up_l="create curve spline vertex 6 7"     #c8

fault_down_r="create curve spline vertex 5 8"   #c9 
fault_down_l="create curve spline vertex 8 7"   #c10


cubit.cmd(fault_up_r) 
cubit.cmd(fault_up_l) 
cubit.cmd(fault_down_r)
cubit.cmd(fault_down_l)
 

#### Branch Fault ############################# 
fault_up="create curve spline vertex 9 10 11"    #c11 
fault_down="create curve spline vertex 9 12 11"  #c12

cubit.cmd(fault_up) 
cubit.cmd(fault_down) 
###############################################

surface1="create surface curve 2 1 6 5 4 7 8 9 12"
cubit.cmd(surface1)
surface2="create surface curve 3 10 11"
cubit.cmd(surface2)

cubit.cmd("sweep surface 1  vector 0  0 -1 distance "+str(30*km)) 
cubit.cmd("sweep surface 2  vector 0  0 -1 distance "+str(30*km)) 


######## Chuncks ########################################
L= 42*km

xchunk = []
ychunk = []
zchunk = []

xchunk.append(-L-6*km)   #x1..
xchunk.append(L)    #x
xchunk.append(L)    #x3
xchunk.append(-L-6*km)   #x4

ychunk.append(-L)  #y1
ychunk.append(-L)  #y2
ychunk.append(L)   #y3
ychunk.append(L)   #y4

zchunk=[z_surf]*4


####################  chunck ###########################################
for i in range(len(xchunk)): 
   vert="create vertex x "+str(xchunk[i])+" y "+str(ychunk[i])+" z "+str(zchunk[i]) 
   cubit.cmd(vert) 

########### creating chuncks ############################################

chunk1="create curve vertex 39 40"   # c1
chunk2="create curve vertex 40 41"   # c2
chunk3="create curve vertex 41 42"   # c3
chunk4="create curve vertex 42 39"   # c4

cubit.cmd(chunk1)
cubit.cmd(chunk2)
cubit.cmd(chunk3)
cubit.cmd(chunk4)

surface1="create surface curve 37 38 39 40"
cubit.cmd(surface1)

cubit.cmd("sweep surface 17  vector 0  0 -1 distance "+str(30*km)) 

cubit.cmd('subtract volume 1 2 from 3 keep') 
cubit.cmd('delete volume 3') 

#######  lateral faces of chuncks  #### 
cubit.cmd('create surface skin curve 15 90')   
cubit.cmd('create surface skin curve 25 92')   
cubit.cmd('create surface skin curve 27 93')   
cubit.cmd('create surface skin curve 29 91')   

### webcutting 4 chuncks ####
cubit.cmd("create body loft surface 45 48")
cubit.cmd("create body loft surface 48 47")
cubit.cmd("create body loft surface 47 46")
cubit.cmd("create body loft surface 46 45")

#####################################################
cubit.cmd('delete volume 4')

cubit.cmd('imprint all') 
cubit.cmd('merge all') 

## meshing inner solids ###
h_size = 300

## forcing aristas to have 300 m ##

cubit.cmd("surface 1 2 size "+str(h_size))
cubit.cmd("volume 1 2 size "+str(h_size))
cubit.cmd("surface 1 2 scheme pave")
cubit.cmd("mesh surface 1 2")
cubit.cmd("mesh volume 1 2")

bia=1.03
inter= 40

## Meshing Aristas #####
cubit.cmd('curve 95 97 99 101 103 105 107 109 interval '+str(inter))
cubit.cmd('curve 95 97 99 101 103 105 107 109 scheme bias factor '+str(bia))
cubit.cmd('mesh curve 95 97 99 101 103 105 107 109')

### meshing chuncks #######
### CHUNCK 1
cubit.cmd('volume 12 scheme Sweep  source surface 8 13 3 target surface 71 rotate off')
cubit.cmd('volume 12 sweep smooth Auto')
cubit.cmd('mesh volume 12')
#### CHUNCK 2
cubit.cmd('volume 9 scheme Sweep  source surface 11 target surface 53 rotate off')
cubit.cmd('volume 9 sweep smooth Auto')
cubit.cmd('mesh volume 9')
#### CHUNCK 3
cubit.cmd('volume 10 scheme Sweep  source surface 10 target surface 59 rotate off')
cubit.cmd('volume 10 sweep smooth Auto')
cubit.cmd('mesh volume 10')
#### CHUNCK 4
cubit.cmd('volume 11 scheme Sweep  source surface 9 target surface 65 rotate off')
cubit.cmd('volume 11 sweep smooth Auto')
cubit.cmd('mesh volume 11')

########### Fault elements and nodes ###############

### Main Fault ####################################################################
os.system('mkdir -p MESH') 

Au = [6,7]   # A_up
Ad = [5,14]  # A_down

faultA = fault_input(1,Au,Ad)

####### Branch Fault##############################################################################
Bu = [4]   # B_up
Bd = [15]  # B_down

faultB = fault_input(2,Bu,Bd)

##  FOR THE BULK (Seismic wave propagation part for SESAME)

####### This is boundary_definition.py of GEOCUBIT 
##..... which extracts the bounding faces and defines them into blocks 
entities=['face'] 
define_parallel_bc(entities) 
 
#### Define material properties for the 2 volumes ################ 
cubit.cmd('#### DEFINE MATERIAL PROPERTIES #######################') 
 
# Material properties in concordance with tpv14-15 benchmark. 
 
cubit.cmd('block 1 name "elastic 1" ')        # material region  
cubit.cmd('block 1 attribute count 5') 
cubit.cmd('block 1 attribute index 1 1')      # flag for fault side 1 
cubit.cmd('block 1 attribute index 2 6000')   # vp 
cubit.cmd('block 1 attribute index 3 3464')    # vs 
cubit.cmd('block 1 attribute index 4 2670')   # rho 
cubit.cmd('block 1 attribute index 5 13')     # Q flag (see constants.h: IATTENUATION_ ... )

cubit.cmd('block 2 name "elastic 2" ')        # material region  
cubit.cmd('block 2 attribute count 5') 
cubit.cmd('block 2 attribute index 1 1')      # flag for fault side 1 
cubit.cmd('block 2 attribute index 2 6000')   # vp 
cubit.cmd('block 2 attribute index 3 3464')    # vs 
cubit.cmd('block 2 attribute index 4 2670')   # rho 
cubit.cmd('block 2 attribute index 5 13')     # Q flag (see constants.h: IATTENUATION_ ... ) 

# Material properties in concordance with tpv14-15 benchmark chuncks. 
 
 
cubit.cmd('block 3 name "elastic 3" ')        # material region  
cubit.cmd('block 3 attribute count 5') 
cubit.cmd('block 3 attribute index 1 1')      # flag for fault side 1 
cubit.cmd('block 3 attribute index 2 6000')   # vp 
cubit.cmd('block 3 attribute index 3 3464')    # vs 
cubit.cmd('block 3 attribute index 4 2670')   # rho 
cubit.cmd('block 3 attribute index 5 13')     # Q flag (see constants.h: IATTENUATION_ ... ) 

cubit.cmd('block 4 name "elastic 4" ')        # material region  
cubit.cmd('block 4 attribute count 5') 
cubit.cmd('block 4 attribute index 1 1')      # flag for fault side 1 
cubit.cmd('block 4 attribute index 2 6000')   # vp 
cubit.cmd('block 4 attribute index 3 3464')    # vs 
cubit.cmd('block 4 attribute index 4 2670')   # rho 
cubit.cmd('block 4 attribute index 5 13')     # Q flag (see constants.h: IATTENUATION_ ... ) 

cubit.cmd('block 5 name "elastic 5" ')        # material region  
cubit.cmd('block 5 attribute count 5') 
cubit.cmd('block 5 attribute index 1 1')      # flag for fault side 1 
cubit.cmd('block 5 attribute index 2 6000')   # vp 
cubit.cmd('block 5 attribute index 3 3464')    # vs 
cubit.cmd('block 5 attribute index 4 2670')   # rho 
cubit.cmd('block 5 attribute index 5 13')     # Q flag (see constants.h: IATTENUATION_ ... ) 

cubit.cmd('block 6 name "elastic 6" ')        # material region  
cubit.cmd('block 6 attribute count 5') 
cubit.cmd('block 6 attribute index 1 1')      # flag for fault side 1 
cubit.cmd('block 6 attribute index 2 6000')   # vp 
cubit.cmd('block 6 attribute index 3 3464')    # vs 
cubit.cmd('block 6 attribute index 4 2670')   # rho 
cubit.cmd('block 6 attribute index 5 13')     # Q flag (see constants.h: IATTENUATION_ ... ) 

#### Export to SESAME format using cubit2specfem3d.py of GEOCUBIT 
 
cubit2specfem3d.export2SESAME('MESH')  
 
# all files needed by SCOTCH are now in directory MESH 

