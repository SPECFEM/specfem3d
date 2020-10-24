#!/usr/bin/python2.7
#from matplotlib.mlab import griddata
#this script generates topology and slab interface in cubit
from __future__ import print_function

import cubit
import os
import sys
import math
import numpy as np
#import matplotlib.pyplot as plt

cubit.init([''])
cubit.cmd('reset')


#-- BEGIN user settings ----------------------

FaultFileName = './kur_slab1.0_clip.xyz'  # fault surface file: longitude, latitude, depth (in km, negative)
#FaultFileName = './alu_slab1.0_clip.xyz'
#FaultFileName = './aluslab.xyz'

Latmin = 33          # regional box limits
Latmax = 44
Lonmin = 136
Lonmax = 150
Lat2dis = 100.0      # latitude to km conversion factor
                     # The true latitude to distance conversion ratio should be 111.195 km = 1 deg.
                     # We will reflect that at the end of exportmesh.py by scaling up the model
                     # by a factor of 1.1195. The reason we don't do that here is that it would
                     # change the whole mesh generation script.
Lon2dis = 76.0       # longitude to km conversion factor
zcutBottom = 100.0   # bottom depth of the fault in km, preliminary value (the final depth is set in exportmesh.py)
zcutTop = 15.0        # steepen the fault surface above this depth in km
                     # to avoid elements with small angles at the trench
                     # Use the matlab script demo_subduction_smoothing.m to explore this feature
rotate = -15         # set this to minus average strike. Approximately aligns the trench with the Y axis to facilitate meshing
# dimensions (in km) of the box containing the fault:
xLim = 600           # along the X axis of the rotated frame (it must fit inside the lat-lon box)
yLim = 800           # along the Y axis of the rotated frame (it must fit inside the lat-lon box)
zLim = 250           # depth (it must be deeper than zcutBottom)
radius = 1000.0      # radius of semi-spherical absorbing boundary in km
Meshsize = 4.0       # mesh size in km
refine_slab = False  # refine the mesh near the subduction interface. Can be very slow
Mesh = False         # set to true to trigger meshing
Plotsquare=False

#-- END user settings ----------------------


xc = 0.5*(Lonmin+Lonmax)
yc = 0.5*(Latmin+Latmax)

# rotate function doesn't seem to work
#def rotate(X,deg):
#    rad = math.radians(deg)
#    Y=X
#    Y[:,0] = X[:,0]* math.cos(rad) + X[:,1] * math.sin(rad)
#    Y[:,1] = - X[:,0]* math.sin(rad) + X[:,1]*math.cos(rad)
#    return Y

# this function modifies the slab geometry near the surface
# by making the fault dip angle steeper above depth "z"
# to avoid elements with small angles at the trench.
# See demo_subduction_smoothing.m
def surf(z, minz):
#    c = 0.5/zcutTop
#    f = -math.log(math.exp(-c*min(z - minz, -0.1))-1)/c + minz
    if(z > -zcutTop):
        f = z + 2.0 * (z + zcutTop);
    else:
        f = z
    return max(f, -zcutBottom)


def import_elev_data():
    cmd = 'awk \'{{if($1>={0} && $1<={1} && $2>={2} && $2<={3}) print }}\' ./etopo2.xyz > ./subduct_slab.xyz'.format(Lonmin,Lonmax,Latmin,Latmax)
    print('extract data from the global topography file etopo2.xyz')
    os.system(cmd)
    print('done!')
    filename = './subduct_slab.xyz'
    data = np.loadtxt(filename)
    MASK = np.bitwise_not(data[:,2]==-9999)
    MASK = np.bitwise_and(MASK,data[:,1]>=Latmin)
    MASK = np.bitwise_and(MASK,data[:,1]<=Latmax)
    MASK = np.bitwise_and(MASK,data[:,0]>=Lonmin)
    MASK = np.bitwise_and(MASK,data[:,0]<=Lonmax)
    data = data[MASK,:]
    return data


def import_slab_data(filename):
    data = np.loadtxt(filename)
    data = data[np.bitwise_not(np.isnan(data[:,2])),:]
    data = data[np.bitwise_and(np.bitwise_and(data[:,1]>=Latmin, data[:,1]<=Latmax),data[:,2]>-zcutBottom*2.5)]
    return data


def generating_contour(X,Y,Z):
    xi = np.linspace(min(X),max(X),10)
    yi = np.linspace(min(Y),max(Y),10)
    Xm,Ym = np.meshgrid(xi,yi)
    Zm = griddata(X,Y,Z,xi,yi)
    return Xm,Ym,Zm

data = import_slab_data(FaultFileName)

#plt.scatter(data[:,0],data[:,1],1.0,c=data[:,2],edgecolors='none')

if Plotsquare:
    XB = [393.3,186.3,-393.3,-186.3,393.3]
    YB = [308.7,-464.0,-308.7,464.0,308.7]
    XB = np.array(XB)/Lon2dis+xc
    YB = np.array(YB)/Lat2dis+yc
    print(XB)
    print(YB)
    #plt.plot(XB,YB)
#Xm,Ym,Zm = generating_contour(data[:,0],data[:,1],data[:,2])
#plt.contourf(Xm,Ym,Zm)
    #plt.colorbar()
    #plt.show()
np.savetxt('outputslab.txt',data,delimiter=',')

#exit()

N,D = data.shape
X = range(0,N,5) # down sampling
data = data[X,:]
#xc = np.mean(data[:,0])
#yc = np.mean(data[:,1])
xc = 0.5*(Lonmax+Lonmin)
yc = 0.5*(Latmax+Latmin)
data[:,0]=data[:,0]-xc
data[:,1]=data[:,1]-yc
#data[:,0:2] = rotate(data[:,0:2],45)
N,D = data.shape
mindepth = max(data[:,2])
# falt surface: create one spline per latitude
print(N)
start = 1
n_curve = 0
n = 0
for ii in range(0,N):
    vert = "create vertex x "+str(data[ii,0]*Lon2dis)+" y "+str(data[ii,1]*Lat2dis)+" z "+str(surf(data[ii,2], mindepth))
    cubit.cmd(vert)
    if(ii<N-1):
        if(data[ii,1]!=data[ii+1,1]):
            end = ii+1
            mindepth = data[ii,2]
            if(n%10 == 0):
                cc = "create curve spline vertex "+str(start)+" to "+str(end)+" delete"
                cubit.cmd(cc)
                n_curve = n_curve + 1
            n = n+1
            start = ii+2
n_curve_slab = n_curve
print('total curve %d'%n_curve)

#cubit.cmd('delete vertex  all')
#cubit.cmd('create surface skin curve 1 to {0}'.format(n_curve))
#for ii in range(1,n_curve,20):#
#cubit.cmd('create surface skin curve %d to %d'%(ii,ii+20))
data = import_elev_data()
if Plotsquare:
    #plt.scatter(data[:,0],data[:,1],1.0,c=data[:,2],edgecolors='none')
    XB = [393.3,186.3,-393.3,-186.3,393.3]
    YB = [308.7,-464.0,-308.7,464.0,308.7]
    XB = np.array(XB)/Lon2dis+xc
    YB = np.array(YB)/Lat2dis+yc
    #plt.plot(XB,YB)
#Xm,Ym,Zm = generating_contour(data[:,0],data[:,1],data[:,2])
#plt.contourf(Xm,Ym,Zm)
    #plt.colorbar()
    #plt.show()
np.savetxt('outputelev.txt',data,delimiter=',')


N2,D = data.shape
X = range(0,N2,10) # down sampling
data = data[X,:]
data[:,0]=data[:,0]-xc
data[:,1]=data[:,1]-yc
#data[:,0:2] = rotate(data[:,0:2],45)
N2,D = data.shape

# topography surface: create one spline per latitude
print(N2)
#exit()
start = N+1
start_curve = n_curve+1
for ii in range(0,N2):
    vert = "create vertex x "+str(data[ii,0]*Lon2dis)+" y "+str(data[ii,1]*Lat2dis)+" z "+str(data[ii,2]/1000.0)
    # topography file: longitude, latitude and depth
    # depths are negative and in meters
    cubit.cmd(vert)
    if(ii<N2-1):
        if(data[ii,1]!=data[ii+1,1]):
            end = N+ii+1
            if(n_curve%1 == 0):
                cc = "create curve spline vertex "+str(start)+" to "+str(end)+" delete"
                cubit.cmd(cc)
                n_curve = n_curve+1
            start = N+ii+2

print('total curve %d'%n_curve)

# create fault and topography surfaces
cubit.cmd('delete vertex  all')
cubit.cmd('create surface skin curve {0} to {1}'.format(1,n_curve_slab))
cubit.cmd('create surface skin curve {0} to {1}'.format(n_curve_slab+1,n_curve))
#cubit.cmd('create surface skin curve {0} to {1}'.format(start_curve,n_curve))
#for ii in range(1,n_curve,20):#
#cubit.cmd('create surface skin curve %d to %d'%(ii,ii+20))
cubit.cmd("delete curve all")

cubit.cmd('set node constraint on')

# create a box and cut it using fault and topography smooth surfaces
cubit.cmd('create brick x {0} y {1} z {2}'.format(xLim,yLim,zLim*2))
# rotate to align the trench (on average) with the Y axis
cubit.cmd('rotate vol 3 about 0 0 0 direction 0 0 1 angle {0}'.format(rotate))
cubit.cmd("save as 'slab_rotate.cub' overwrite")
cubit.cmd('webcut vol 3 with sheet body 2')
cubit.cmd('delete vol 4')
cubit.cmd('webcut vol 3 with sheet body 1')
cubit.cmd('delete vol 1 2')
cubit.cmd('compress all')
cubit.cmd("merge all")

if True:

    cubit.cmd("create sphere radius {0}".format(radius))
    cubit.cmd('project surface 1 onto surf 13 imprint keepcurve keepbody')
    cubit.cmd('project surface 2 onto surf 13 imprint keepcurve keepbody')
    cubit.cmd('project surface 4 onto surf 13 imprint keepcurve keepbody')
    cubit.cmd('project surface 5 onto surf 13 imprint keepcurve keepbody')
    cubit.cmd('project surface 7 onto surf 13 imprint keepcurve keepbody')
    cubit.cmd('project surface 9 onto surf 13 imprint keepcurve keepbody')
    cubit.cmd('project surface 10 onto surf 13 imprint keepcurve keepbody')
    cubit.cmd('project surface 12 onto surf 13 imprint keepcurve keepbody')
    loft_surf = np.array([[1,15],[2,16],[4,18],[5,20],[7,23],[9,25],[10,27],[12,29]])
    cubit.cmd("save as 'slab_rotate_before_loft.cub' overwrite")

    exit()
#--- Everything below this line is currently ignored. The rest of the meshing is completed by two other scripts.

    for ii in range(0,8):
        cubit.cmd('create volume loft surface {0} {1}'.format(loft_surf[ii,0],loft_surf[ii,1]))


if Mesh:
    cubit.cmd("merge all")
    cubit.cmd("vol all size {0}".format(Meshsize))
    cubit.cmd("surf 5 12 scheme pave")
    cubit.cmd("mesh surf 5 12")
    cubit.cmd("vol 1 scheme sweep source 5 target 7")
    cubit.cmd("vol 2 scheme sweep source 12 target 10")
    cubit.cmd("mesh vol 1")
    cubit.cmd("mesh vol 2")
    if refine_slab:
        cubit.cmd("refine surf 3 depth 3")
cubit.cmd("save as 'slab_rotate.cub' overwrite")
