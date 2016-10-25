#!/usr/bin/python2.7
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.mlab import griddata
#this script generates topology and slab interface in cubit
import os
import sys
import cubit
import os
import sys
import math
 

cubit.init([''])
cubit.cmd('reset')
Latmin = 33
Latmax = 44
Lonmin = 136
Lonmax = 150
xc = 0.5*(Lonmin+Lonmax)
yc = 0.5*(Latmin+Latmax)
Lat2dis = 111.195 # 100km = 1deg in latitude
Lon2dis = Lat2dis * math.cos(math.radians(0.5*(Latmin+Latmax)))
Meshsize = 4.0 # mesh size set to 4.0 km
radius = 1000.0
cuttingdepth = -100.0
refine_slab = False # refine mesh near the subduction interface , can be very slow
Mesh = False # set to true to trigger meshing 
rotate =-15
Plotsquare=False
# rotate function doesn't seem to work
#def rotate(X,deg):
#    rad = math.radians(deg)
#    Y=X
#    Y[:,0] = X[:,0]* math.cos(rad) + X[:,1] * math.sin(rad)
#    Y[:,1] = - X[:,0]* math.sin(rad) + X[:,1]*math.cos(rad)
#    return Y

def surf(z):#this function modifies the slab geometry at surface
    c = 5.0e-2
    f = -math.log(math.exp(-c*z)-1)/c
    f = max(f,cuttingdepth)
    return f


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


def import_slab_data():
    filename = './alu_slab1.0_clip.xyz'
    filename = './kur_slab1.0_clip.xyz'
#    filename = './aluslab.xyz'
    data = np.loadtxt(filename)
    data = data[np.bitwise_not(np.isnan(data[:,2])),:]
    data = data[np.bitwise_and(np.bitwise_and(data[:,1]>=Latmin, data[:,1]<=Latmax),data[:,2]>-250)] # select data between latitude 45 and 50 degree and depth above 250km 
    return data


def generating_contour(X,Y,Z):
    xi = np.linspace(min(X),max(X),10)
    yi = np.linspace(min(Y),max(Y),10)
    Xm,Ym = np.meshgrid(xi,yi)
    Zm = griddata(X,Y,Z,xi,yi)
    return Xm,Ym,Zm

data = import_slab_data()

plt.scatter(data[:,0],data[:,1],1.0,c=data[:,2],edgecolors='none')

if Plotsquare:
    XB = [393.3,186.3,-393.3,-186.3,393.3]
    YB = [308.7,-464.0,-308.7,464.0,308.7]
    XB = np.array(XB)/Lon2dis+xc
    YB = np.array(YB)/Lat2dis+yc
    print(XB)
    print(YB)
    plt.plot(XB,YB)
#Xm,Ym,Zm = generating_contour(data[:,0],data[:,1],data[:,2])
#plt.contourf(Xm,Ym,Zm)
    plt.colorbar()
    plt.show()
np.savetxt('outputslab.txt',data,delimiter=',')

#exit()

N,D = data.shape
X = range(0,N,5) # down sampling 
data = data[X,:]
xc = np.mean(data[:,0])
yc = np.mean(data[:,1])
xc = 0.5*(Lonmax+Lonmin)
yc = 0.5*(Latmax+Latmin)
data[:,0]=data[:,0]-xc
data[:,1]=data[:,1]-yc
#data[:,0:2] = rotate(data[:,0:2],45)
N,D = data.shape

print(N)
start = 1
n_curve = 0
n = 0
for ii in range(0,N):
    vert = "create vertex x "+str(data[ii,0]*Lon2dis)+" y "+str(data[ii,1]*Lat2dis)+" z "+str(surf(data[ii,2]))
    cubit.cmd(vert)
    if(ii<N-1):
        if(data[ii,1]!=data[ii+1,1]):
            end = ii+1
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
    plt.scatter(data[:,0],data[:,1],1.0,c=data[:,2],edgecolors='none')
    XB = [393.3,186.3,-393.3,-186.3,393.3]
    YB = [308.7,-464.0,-308.7,464.0,308.7]
    XB = np.array(XB)/Lon2dis+xc
    YB = np.array(YB)/Lat2dis+yc
    plt.plot(XB,YB)
#Xm,Ym,Zm = generating_contour(data[:,0],data[:,1],data[:,2])
#plt.contourf(Xm,Ym,Zm)
    plt.colorbar()
    plt.show()
np.savetxt('outputelev.txt',data,delimiter=',')


N2,D = data.shape
X = range(0,N2,10) # down sampling 
data = data[X,:]
data[:,0]=data[:,0]-xc
data[:,1]=data[:,1]-yc
#data[:,0:2] = rotate(data[:,0:2],45)
N2,D = data.shape

print(N2)
#exit()
start = N+1
start_curve = n_curve+1
for ii in range(0,N2):
    vert = "create vertex x "+str(data[ii,0]*Lon2dis)+" y "+str(data[ii,1]*Lat2dis)+" z "+str(data[ii,2]/1000.0) # the topography profile longitude latitude and depth.depth are minus signed and of unit meters. A crude conversion of 1 degree=100 km , should be improved later. 
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

cubit.cmd('delete vertex  all')
cubit.cmd('create surface skin curve {0} to {1}'.format(1,n_curve_slab))
cubit.cmd('create surface skin curve {0} to {1}'.format(n_curve_slab+1,n_curve))
#cubit.cmd('create surface skin curve {0} to {1}'.format(start_curve,n_curve))
#for ii in range(1,n_curve,20):#
#cubit.cmd('create surface skin curve %d to %d'%(ii,ii+20))
cubit.cmd("delete curve all")
cubit.cmd('set node constraint on')
cubit.cmd('create brick x 600 y 800 z 400')
cubit.cmd("save as 'slab_rotate.cub' overwrite")

cubit.cmd('vol 3  move 0 0 -50')
cubit.cmd('rotate vol 3 about 0 0 0 direction 0 0 1 angle {0}'.format(rotate))

cubit.cmd("save as 'slab_rotate.cub' overwrite")
cubit.cmd('webcut vol 3 with sheet body 2')
cubit.cmd('delete vol 4')
cubit.cmd('webcut vol 3 with sheet body 1')
cubit.cmd('delete vol 1 2')
cubit.cmd('compress all')
#cubit.cmd("vol all move -100 -350 0")
cubit.cmd("merge all")
if True:

# make wrapping sphere not working 
    cubit.cmd("create sphere radius {0}".format(radius))
#cubit.cmd("section volume 3 with zplane offset 0 reverse")

#!python
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
