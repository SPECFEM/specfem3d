#############################################################################
# local_volume.py                                                    
# this file is part of GEOCUBIT                                             #
#                                                                           #
# Created by Emanuele Casarotti                                             #
# Copyright (c) 2008 Istituto Nazionale di Geofisica e Vulcanologia         #
#                                                                           #
#############################################################################
#                                                                           #
# GEOCUBIT is free software: you can redistribute it and/or modify          #
# it under the terms of the GNU General Public License as published by      #
# the Free Software Foundation, either version 3 of the License, or         #
# (at your option) any later version.                                       #
#                                                                           #
# GEOCUBIT is distributed in the hope that it will be useful,               #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
# GNU General Public License for more details.                              #
#                                                                           #
# You should have received a copy of the GNU General Public License         #
# along with GEOCUBIT.  If not, see <http://www.gnu.org/licenses/>.         #
#                                                                           #
#############################################################################
try:
    import start as start
    cubit                   = start.start_cubit()
except:
    try:
        import cubit
    except:
        print 'error importing cubit, check if cubit is installed'
        pass

numpy                       = start.start_numpy()

def check_orientation(grdfileNAME):

    try:
         grdfile = open(grdfileNAME, 'r')
         print 'reading ',grdfileNAME
    except:
         txt='check_orintation ->error reading: '+  str( grdfile )
         raise Exception(txt)
    diff=1
    txt=grdfile.readline()
    x0,y0,z=map(float,txt.split())
    while diff>0:
        try:
            txt=grdfile.readline()
        except:
            break
        x,y,z=map(float,txt.split())
        diff=x-x0
        x0=x
    diff=y-y0
    if diff>0:
        orientation= 'SOUTH2NORTH'
    else:
        orientation= 'NORTH2SOUTH'
    grdfile.close()
    return orientation
    
def read_irregular_surf(filename):
    "read irregular grid"
    try:
        xyz = numpy.loadtxt(filename)
    except:
        txt='error reading '+filename
        raise Exception(txt)
    gridpoints = xyz[:,0:2]
    z = xyz[:,2]
    return gridpoints,z
    
def get_interpolated_elevation(point,gridpoints,z,k=1):
    """for x0 and y0 return the interpolated z point of a irregular x,y,z grid
    K=1 nearest
    K>1 number of points used in a inverse distance weighted interpolation
    
    point=(x0,y0)
    gridpoints=numpy.array([[x1,y1],[x2,y2],...)
    
    """
    dist=numpy.sum((point-gridpoints)**2,axis=1)
    zindex=dist.argsort()[:k]
    kmindist=dist[zindex]
    w=1/kmindist
    w/=w.sum()
    zi = numpy.dot(w.T, z[zindex])
    return zi

def create_grid(xmin,xmax,ymin,ymax,xstep,ystep):
    """create regular grid with xmin,xmax by xstep and  ymin,ymax by ystep"""
    x,y=numpy.mgrid[xmin:xmax+xstep/2.:xstep,ymin:ymax+ystep/2.:ystep] #this includes the bounds
    gridpoints = numpy.vstack([x.ravel(), y.ravel()]).T
    return x,y,gridpoints
    

def process_surfacefiles(iproc,nx,ny,nstep,grdfile,unit,lat_orientation):
        from utilities import geo2utm
        numpy                       = start.start_numpy()
        elev=numpy.zeros([nx,ny],float)
        coordx=numpy.zeros([nx,ny],float)
        coordy=numpy.zeros([nx,ny],float)
        icoord=0
        
        
        lat_orientation=check_orientation(grdfile)
        
        try:
             grdfile = open(grdfile, 'r')
             #print 'reading ',grdfile
        except:
             txt='error reading: '+  str( grdfile )
             raise Exception(txt)
        
        
        if lat_orientation is 'SOUTH2NORTH':
            rangey=range(0,ny)
        else:
            rangey=range(ny-1,-1,-1)
            lat_orientation='NORTH2SOUTH'
        print lat_orientation
        for iy in rangey:
            for ix in range(0,nx):
                txt=grdfile.readline()
                try:
                    if len(txt) != 0:
                        x,y,z=map(float,txt.split())
                        if iy%nstep == 0 and ix%nstep == 0:
                            icoord=icoord+1
                            x_current,y_current=geo2utm(x,y,unit)
                            jx=min(nx-1,ix/nstep)
                            jy=min(ny-1,iy/nstep)
                            coordx[jx,jy]=x_current
                            coordy[jx,jy]=y_current
                            elev[jx,jy]=z      
                except:
                    print 'error reading point ',iy*nx+ix,txt, grdfile.name, ' proc '
                    raise NameError, 'error reading point'

        
        if  (nx)*(ny) != icoord: 
            if iproc == 0: print 'error in the surface file '+grdfile.name
            if iproc == 0: print 'x points ' +str(nx)+ ' y points ' +str(ny)+ ' tot points '+str((nx)*(ny)) 
            if iproc == 0: print 'points read in '+grdfile.name+': '+str(icoord)
            raise NameError
        
        grdfile.close()
        
        return coordx,coordy,elev







def process_irregular_surfacefiles(iproc,nx,ny,xmin,xmax,ymin,ymax,xstep,ystep,grdfile,unit_surf,lat_orientation):
    gridpoints,z=read_irregular_surf(grdfile)
    coordx,coordy,points=create_grid(xmin,xmax,ymin,ymax,xstep,ystep)

    elev = numpy.empty([len(points)])
    for i in xrange(len(points)):
        elev[i] = get_interpolated_elevation(points[i],gridpoints,z,k=4)
        
    coordx.shape=(nx,ny)
    coordy.shape=(nx,ny)
    elev.shape=(nx,ny)
    
    return coordx,coordy,elev


def read_grid(filename=None):
    import sys
    import start as start
    mpiflag,iproc,numproc,mpi   = start.start_mpi()
    #
    numpy                       = start.start_numpy()
    cfg                         = start.start_cfg(filename=filename)
    from utilities import geo2utm
    
    #if cfg.irregulargridded_surf==True then cfg.nx and cfg.ny are the desired number of point along the axis....
    if cfg.nx and cfg.ny:
        nx=cfg.nx
        ny=cfg.ny
        if cfg.nstep:
            nx=min(cfg.nx,int(cfg.nx/cfg.nstep)+1)
            ny=min(cfg.ny,int(cfg.ny/cfg.nstep)+1)
            nstep=cfg.nstep
        else:
            nstep=1
    else:
        try:
            xstep=cfg.step
            ystep=cfg.step
        except:
            xstep=cfg.xstep
            ystep=cfg.ystep
        nx= int((cfg.longitude_max-cfg.longitude_min)/xstep)+1
        ny= int((cfg.latitude_max-cfg.latitude_min)/ystep)+1
        nstep=1
    #
    
    if cfg.irregulargridded_surf:
        xt,xstep=numpy.linspace(cfg.xmin, cfg.xmax, num=nx, retstep=True)
        yt,ystep=numpy.linspace(cfg.ymin, cfg.ymax, num=ny, retstep=True)
    
    elev=numpy.zeros([nx,ny,cfg.nz],float)
    #
    if  cfg.bottomflat: 
        elev[:,:,0] = cfg.depth_bottom
        bottomsurface=1
    else:
        bottomsurface=0
        
    for inz in range(bottomsurface,cfg.nz-1):
        grdfilename=cfg.filename[inz-bottomsurface]

        if cfg.irregulargridded_surf:
            coordx,coordy,elev_1=process_irregular_surfacefiles(iproc,nx,ny,cfg.xmin,cfg.xmax,cfg.ymin,cfg.ymax,xstep,ystep,grdfilename)
        else:
            coordx,coordy,elev_1=process_surfacefiles(iproc,nx,ny,nstep,grdfilename,cfg.unit,cfg.lat_orientation)
        elev[:,:,inz]=elev_1[:,:]
        #
    
    inz=cfg.nz-1 #last surface
    if cfg.sea:
        elev[:,:,inz]=elev[:,:,inz-1]
    else:
        #try:
        grdfile = cfg.filename[inz-bottomsurface]
        print 'reading ',cfg.filename[inz-bottomsurface]
        if cfg.irregulargridded_surf:
            coordx,coordy,elev_1=process_irregular_surfacefiles(iproc,nx,ny,cfg.xmin,cfg.xmax,cfg.ymin,cfg.ymax,xstep,ystep,grdfile)
        else:
            coordx,coordy,elev_1=process_surfacefiles(iproc,nx,ny,nstep,grdfile,cfg.unit,cfg.lat_orientation)
        elev[:,:,inz]=elev_1[:,:]
        #except:
        #     txt='error reading: '+  str( cfg.filename[inz-bottomsurface] )
        #    raise NameError, txt
        
        
        if cfg.subduction:
          print 'subduction'
          top=elev[:,:,inz]
          slab=elev[:,:,inz-1]
          subcrit=numpy.abs(top-slab)<cfg.subduction_thres
          top[subcrit]=slab[subcrit]+cfg.subduction_thres
          print len(top[subcrit])
          elev[:,:,inz]=top
    return coordx,coordy,elev,nx,ny

