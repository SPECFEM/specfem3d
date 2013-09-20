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

def read_grid(filename=None):
    import sys
    import start as start
    mpiflag,iproc,numproc,mpi   = start.start_mpi()
    #
    numpy                       = start.start_numpy()
    cfg                         = start.start_cfg(filename=filename)
    from utilities import geo2utm
    
    #     
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
    elev=numpy.zeros([nx,ny,cfg.nz],float)
    coordx=numpy.zeros([nx,ny],float)
    coordy=numpy.zeros([nx,ny],float)
    #
    if  cfg.bottomflat: 
        elev[:,:,0] = cfg.depth_bottom
        bottomsurface=1
    else:
        bottomsurface=0
            #
    for inz in range(bottomsurface,cfg.nz):
        try:
             grdfile = open(cfg.filename[inz-bottomsurface], 'r')
             print 'reading ',cfg.filename[inz-bottomsurface]
        except:
             txt='error reading: '+  str( cfg.filename[inz-bottomsurface] )
             raise NameError, txt
        #
        icoord=0
        for iy in range(0,ny):
            for ix in range(0,nx):
                txt=grdfile.readline()
                try:
                    if len(txt) != 0:
                        x,y,z=map(float,txt.split())
                        if iy%nstep == 0 and ix%nstep == 0:
                            icoord=icoord+1
                            x_current,y_current=geo2utm(x,y,cfg.unit)
                            jx=min(nx-1,ix/nstep)
                            jy=min(ny-1,iy/nstep)
                            coordx[jx,jy]=x_current
                            coordy[jx,jy]=y_current
                            elev[jx,jy,inz]=z      
                except:
                    print 'error reading point ',iy*cfg.nx+ix,txt, cfg.filename[inz-bottomsurface], ' proc ',iproc
                    raise NameError, 'error reading point'
                    #
        if  (nx)*(ny) != icoord: 
            if iproc == 0: print 'error in the surface file '+cfg.filename[inz-bottomsurface]
            if iproc == 0: print 'x points ' +str(nx)+ ' y points ' +str(ny)+ ' tot points '+str((nx)*(ny)) 
            if iproc == 0: print 'points read in '+cfg.filename[inz-bottomsurface]+': '+str(icoord)
            raise NameError
            
        #if iproc == 0: print 'end of reading grd ascii file '+cfg.filename[inz-bottomsurface]+' '+str(icoord)+ ' points'
        grdfile.close()
    
    
    return coordx,coordy,elev,nx,ny
    

def extract_volume(xmin,ymin,xmax,ymax,coordx,coordy,elev,nx,ny,filename=None):
    import sys
    import start as start
    #
    mpiflag,iproc,numproc,mpi   = start.start_mpi()
    #
    numpy                       = start.start_numpy()
    cfg                         = start.start_cfg(filename=filename)             
    
    from utilities import geo2utm
    #
    rxstep=coordx[1,0]-coordx[0,0]
    rystep=coordy[0,1]-coordy[0,0]
    
    nxmin_cpu=min(0,int((x0-cfg.xmin)/rxstep)+1-10)
    nymin_cpu=min(0,int((y0-cfg.ymin)/rxstep)+1-10)
    nxmax_cpu=min(nx,int((x0-cfg.xmin)/rystep)+1+10)
    nymax_cpu=min(ny,int((y0-cfg.ymin)/rystep)+1+10)
    #
    #
    icurve=0
    isurf=0
    ivertex=0
    #
    #create vertex
    last_surface=cubit.get_last_id('surface')
    for inz in range(0,cfg.nz):
        if  cfg.bottomflat and inz == 0: #bottom layer
                    
                    x_current,y_current=geo2utm(coordx[nxmin_cpu,nymin_cpu],coordy[nxmin_cpu,nymin_cpu],cfg.unit)
                    cubitcommand= 'create vertex '+ str( x_current )+ ' ' + str( y_current) +' '+ str( cfg.depth_bottom )
                    cubit.cmd(cubitcommand)
                    #
                    x_current,y_current=geo2utm(coordx[nxmin_cpu,nymax_cpu],coordy[nxmin_cpu,nymax_cpu],cfg.unit)
                    cubitcommand= 'create vertex '+ str( x_current )+ ' ' + str( y_current) +' '+ str( cfg.depth_bottom )
                    cubit.cmd(cubitcommand)                                                                              
                    #
                    x_current,y_current=geo2utm(coordx[nxmax_cpu,nymax_cpu],coordy[nxmax_cpu,nymax_cpu],cfg.unit)
                    cubitcommand= 'create vertex '+ str( x_current )+ ' ' + str( y_current) +' '+ str( cfg.depth_bottom )
                    cubit.cmd(cubitcommand)
                    #
                    x_current,y_current=geo2utm(coordx[nxmax_cpu,nymin_cpu],coordy[nxmax_cpu,nymin_cpu],cfg.unit)
                    cubitcommand= 'create vertex '+ str( x_current )+ ' ' + str( y_current) +' '+ str( cfg.depth_bottom )
                    cubit.cmd(cubitcommand)
                    #
                    cubitcommand= 'create surface vertex 1 2 3 4'
                    cubit.cmd(cubitcommand)
                    #
                    isurf = isurf + 1
                    
        else:
                vertex=[]
                
                for iy in range(nymin_cpu,nymax_cpu+1):
                    ivx=0
                    for ix in range(nxmin_cpu,nxmax_cpu+1):
                        zvertex=elev[ix,iy,inz]
                        x_current,y_current=geo2utm(coordx[ix,iy],coordy[ix,iy],cfg.unit)
                        #
                        vertex.append(' Position '+ str( x_current ) +' '+ str( y_current )+' '+ str( zvertex ) )
                #
                print iproc, 'vertex created....'
                n=max(nx,ny)
                uline=[]
                vline=[]
                iv=0
                
                cubit.cmd("set info off")
                cubit.cmd("set echo off")
                cubit.cmd("set journal off")
                
                for iy in range(0,nymax_cpu-nymin_cpu+1):
                    positionx=''
                    for ix in range(0,nxmax_cpu-nxmin_cpu+1):
                        positionx=positionx+vertex[iv]
                        iv=iv+1
                    command='create curve spline '+positionx
                    cubit.cmd(command)
                    uline.append( cubit.get_last_id("curve") )
                for ix in range(0,nxmax_cpu-nxmin_cpu+1):
                    positiony=''
                    for iy in range(0,nymax_cpu-nymin_cpu+1):
                        positiony=positiony+vertex[ix+iy*(nxmax_cpu-nxmin_cpu+1)]
                    command='create curve spline '+positiony
                    cubit.cmd(command)
                    vline.append( cubit.get_last_id("curve") )
                #
                cubit.cmd("set info "+cfg.cubit_info)
                cubit.cmd("set echo "+cfg.echo_info)
                cubit.cmd("set journal "+cfg.jou_info)
                #
                #
                print iproc,'line created....'
                umax=max(uline)
                umin=min(uline)
                vmax=max(vline)
                vmin=min(vline)
                cubitcommand= 'create surface net u curve '+ str( umin )+' to '+str( umax )+ ' v curve '+ str( vmin )+ ' to '+str( vmax )+' heal'
                cubit.cmd(cubitcommand)
                command = "del curve all"
                cubit.cmd(command)
                isurf=isurf+1
                #
                #
        cubitcommand= 'del vertex all'
        cubit.cmd(cubitcommand)
        #cubit_error_stop(iproc,cubitcommand,ner)
    cubitcommand= 'del curve all'
    cubit.cmd(cubitcommand)
    #
    last_surface_2=cubit.get_last_id('surface')
    #
    for inz in range(1,cfg.nz):
        #!cubit cmd
        cubitcommand= 'create volume loft surface '+ str( inz+1 )+' '+str( inz )
        cubit.cmd(cubitcommand)
        #cubit_error_stop(iproc,cubitcommand,ner)
        isurf=isurf+6
    cubitcommand= 'del surface '+str(last_surface+1)+' to '+ str( last_surface_2 )
    cubit.cmd(cubitcommand)
    #cubit_error_stop(iproc,cubitcommand,ner)
    #
    #        
    cubit.cmd("set info "+cfg.cubit_info)
    cubit.cmd("set echo "+cfg.echo_info)
    cubit.cmd("set journal "+cfg.jou_info)
    command = "compress all"
    cubit.cmd(command)
    