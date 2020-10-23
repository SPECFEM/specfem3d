#!/usr/bin/env python
# Planar - Dipping Fault.
# P. Galvez, ETH-Zurich and J-P Ampuero (Caltech). August 26, 2014.
# K. Bai (Caltech). October 2016. Added fault kink (change of dip angle) and re-scaling trick.
from __future__ import print_function

import sys

# by STEP3
def export_block(nb,vp,vs,rho,count=6,Q=9999):
    cubit.cmd('block {0}  name "elastic {0}" '.format(nb))        # material region
    cubit.cmd('block {0} attribute count {1}'.format(nb,count))
    cubit.cmd('block {0} attribute index 1 1'.format(nb))      # flag for fault side 1
    cubit.cmd('block {0} attribute index 2 {1}'.format(nb,vp))   # vp
    cubit.cmd('block {0} attribute index 3 {1}'.format(nb,vs))    # vs
    cubit.cmd('block {0} attribute index 4 {1}'.format(nb,rho))   # rho
    cubit.cmd('block {0} attribute index 5 {1}'.format(nb,Q))     # Q flag (see constants.h: #IATTENUATION_ ... )

def define_bc_topo2(entities,self):
     # Temporal : Variable zmin should be obtained automatically.
     xmin = self.xmin
     xmax = self.xmax
     ymin = self.ymin
     ymax = self.ymax
     bottom = self.bottom
     topo = self.topo
     v_list,name_list=define_block()
     build_block(v_list,name_list)
     print(entities)
     for entity in entities:
         print("##entity: "+str(entity))
         build_block_side(xmin,entity+'_abs_xmin',obj=entity,id_0=1003)
         build_block_side(xmax,entity+'_abs_xmax',obj=entity,id_0=1005)
         build_block_side(ymin,entity+'_abs_ymin',obj=entity,id_0=1004)
         build_block_side(ymax,entity+'_abs_ymax',obj=entity,id_0=1006)
         build_block_side(bottom,entity+'_abs_bottom',obj=entity,id_0=1002)
         build_block_side(topo,entity+'_topo',obj=entity,id_0=1001)

# In ACIS, the third-party solid modeler used by Trelis,
# the smallest representable quantity is resabs
# and the largest one is resabs/resnor.
# The default values are resabs=1e-6 and resnor=1e-10.
# We have found that when the range of size in a model is out of ACIS' range
# Trelis takes very long to generate the mesh or fails.
# Corey Ernst (Trelis) recommends "modeling in ACIS'
# sweet spot [resabs, resabs/resnor], then scaling the model up at the end after
# all the solid modeling operations are done".
# The function get_global_scale computes the re-scaling factor such that the
# re-scaled range of sizes in the mesh, [min_size, max_size]*scale_factor,
# is centered on the ACIS range [resabs, resabs/resnor] in a logarithmic sense.
# It also gives a warning when the range of scales in the model
# exceeds ACIS' default range even after re-scaling.
def get_global_scale(Xmax,Xmin,Ymax,Ymin,Zmin,h_size):
    max_size = max((Xmax-Xmin),(Ymax-Ymin),(-Zmin))
    min_size = h_size
    resabs = 1e-6  # default smallest representable length in ACIS
    resnor = 1e-10 # default largest representable length in ACIS
    scale_factor = math.sqrt((resabs * resabs/resnor)/(max_size*min_size))
    if max_size/min_size > 1./resnor:
        exit('Trelis cannot handle the large dynamic range of the model.')
    return scale_factor




def main(parameter):
    print(parameter)
    cubit.init('[]')
    cubit.cmd('reset')
    Radian = math.pi/180
    L = float(parameter[0])        # Fault length
    alpha1 = float(parameter[3])   # Fault dip angle shallow
    alpha2 = float(parameter[4])   # Fault dip angle deep
    kink_depth = float(parameter[5])
    W1 = (kink_depth)/math.sin(alpha1*Radian)
    W2 = float(parameter[1])       # Fault width
    h_size = float(parameter[2])   # Element size on the fault domain
    small_opening = h_size/10000.0
    ratio = float(parameter[6])
    w1cos = W1*math.cos(alpha1*Radian)
    w2cos = W2*math.cos(alpha2*Radian)
    w1sin = W1*math.sin(alpha1*Radian)
    w2sin = W2*math.sin(alpha2*Radian)
    wsin = w1sin + w2sin
    wcos = w1cos + w2cos
    ########################################
    Xmax = (1+ratio) * wcos
    Xmin = -ratio * wcos
    Ymax = L/2*(1+ratio)
    Ymin = -Ymax
    Zmin = -wsin * (1+ratio)
    sf = get_global_scale(Xmax,Xmin,Ymax,Ymin,Zmin,small_opening)



    #### defining vertices ####
    x = np.array([0 ,0   ,0  ,0    ,w1cos   ,w1cos   ,wcos   ,Xmax   ,Xmax , Xmin   ,Xmin  ,Xmin  ,wcos   ,Xmax ,Xmax ,Xmin, wcos, wcos])*sf
    y = np.array([0  ,Ymax,0  ,Ymin ,Ymin   ,Ymax    ,Ymax    ,Ymax    ,Ymax    ,Ymax    ,Ymin   ,Ymin   ,Ymin   ,Ymin   ,Ymin   ,Ymax , Ymin, Ymax])*sf
    z = np.array([0  ,0  ,0  ,0    ,-w1sin  ,-w1sin  ,Zmin  ,Zmin    ,0      ,0      ,0      ,Zmin  ,Zmin  ,Zmin ,0      ,Zmin   , -wsin, -wsin])*sf

    ### creating vertices ####
    for i in range(len(x)):
      vert="create vertex x "+str(x[i])+" y "+str(y[i])+" z "+str(z[i])
      cubit.cmd(vert)

    cubit.cmd("vertex all move "+str(-0.5*wcos*sf)+" 0 0")
    ##### fault
    cubit.cmd("create curve spline vertex 4 1 2")
    cubit.cmd("create curve spline vertex 4 3 2")

    ### CURVES

    c=[[5,6],[4,5],[2,6],[13,7],[14,8],[15,9],[12,16],[11,10],[12,13],[11,4],[10,2],[16,7],[17,5],[17,13],[7,18],[18,6],[13,14],[7,8],[14,15],[8,9],[4,15],[2,9],[11,12],[10,16],[17,18]]

    for i in range(len(c)):
      curv="create curve vertex "+str(c[i][0])+ " "+ str(c[i][1])
      cubit.cmd(curv)

    # surfaces
    su=[[1,4,3,5],[18, 3, 15, 27],[16, 6, 17, 27],[19, 6, 20, 7],[21, 8, 22, 7],[23, 8, 24, 2],[13,10,12,1],[9,26,10,25],[6,14,9,11]]
    for i in range(len(su)):
       ss="create surface curve "+str(su[i][0])+" "+str(su[i][1])+" "+str(su[i][2])+" "+str(su[i][3])
       cubit.cmd(ss)

    s = [[14,26,13,5,18,17],[12,4,15,16,11,25],[4, 15, 16, 19, 21, 23],[5, 18, 17, 20, 22, 24]]
    for i in range(len(s)):
      ss="create surface curve "+str(s[i][0])+ " "+ str(s[i][1]) + " "+str(s[i][2])+ " "+ str(s[i][3])+" "+ str(s[i][4]) + " "+str(s[i][5])
      cubit.cmd(ss)

    cubit.cmd("create vol surface 1 2 3 7 8 9 10 11 ")
    cubit.cmd("create surface from surface 1")
    #this is because one surface cannot be used to create 2 vol which is ridiculous
    cubit.cmd("create surface from surface 2")
    cubit.cmd("create surface from surface 3")

    cubit.cmd("create vol surface 14 15 16 4 5 6 12 13 ")
    cubit.cmd("compress all")


    ##### MESHING fault ####
    cubit.cmd("imprint all")
    cubit.cmd("merge all")
    cubit.cmd("vol all size "+str(h_size*sf))
    cubit.cmd("mesh volume all")

    os.system('mkdir -p MESH')
    cubit.cmd('set node constraint off')
    cubit.cmd('unmerge surf 1 2')
    cubit.cmd('group "upp" add node in surf 14 15')
    cubit.cmd('group "downn" add node in surf 1 2')
    cubit.cmd('group "upp" remove node with y_coord< '+str(-L/2*sf))
    cubit.cmd('group "upp" remove node with y_coord> '+str(L/2*sf))
    cubit.cmd('group "downn" remove node with y_coord< '+str(-L/2*sf))
    cubit.cmd('group "downn" remove node with y_coord> '+str(L/2*sf))
    cubit.cmd('nodeset 100 group upp')
    cubit.cmd('nodeset 200 group downn')
    cubit.cmd('nodeset 200  move 0 0 '+str(-0.5*small_opening*sf))
    cubit.cmd('nodeset 100  move 0 0 '+str(0.5*small_opening*sf))
    cubit.cmd("vol all scale "+ str(1./sf * 1000)) #scale back, and make unit meter.
    cubit.cmd('delete vertex 1')


    fd = [1,2]  # fault face_down
    fu = [14,15]  # fault face_up


    ##  FOR THE BULK (Seismic wave propagation part for SESAME)
    ####### defining absorbing-boundary surface

    xmin   = [8]
    xmax   = [5]
    ymin   = [11,12]
    ymax   = [10,13]
    zmin   = [4,9]
    zmax   = [6,7] # Free surface.
    entities=['face']

    define_boundaries(entities,xmin,xmax,ymin,ymax,zmin,zmax)



    ##### USER: define material properties ################
    cubit.cmd('#### DEFINE MATERIAL PROPERTIES #######################')
    print('we are here')

    #for iblock in range(1,11,1):
    for iblock in range(1,3,1):
        export_block(nb = iblock,vp=6000,vs=3464,rho=2670)
    cubit.cmd('save as "subduction_kink.cub" overwrite')

    #### Export to SESAME format using cubit2specfem3d.py of GEOCUBIT

    cubit2specfem3d.export2SPECFEM3D('MESH',hex27=False)
    # all files needed by SCOTCH are now in directory MESH
    fault = fault_input(1,fu,fd)

def usage():
    print("Generate mesh for a planar dipping fault with a kink (change of dip angle)")
    print("usage: python dipping_fault_planar_dip10_kink.py L W H A1 A2 ZK DA")
    print("  L  = fault length")
    print("  W  = along-dip width of the fault segment below the kink")
    print("  H  = element size")
    print("  A1 = shallow dip angle")
    print("  A2 = deep dip angle")
    print("  ZK = kink depth")
    print("  DA = distance from fault edges to absorbing boundaries divided by fault length")
    print("  All input lengths in km, angles in degree. Output mesh is in m.")
    print("  To change material properties, edit the script in section commented as 'USER'")

if __name__ == '__main__':
    if(len(sys.argv) != 8):
        usage()
    else:
        import os
        sys.path.append('/opt/Trelis-15.1/bin')

        import cubit
        import numpy as np
        import cubit2specfem3d
        import exportlib
        import math
        import os
        import sys
        from save_fault_nodes_elements import *
        from absorbing_boundary import *
        from boundary_definition import *

        main(sys.argv[1:])

