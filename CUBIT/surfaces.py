#############################################################################
# surfaces.py                                                    
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

def surfaces(filename=None):
    """creating the surfaces defined in the parameter files
       #
       nsurf = number of surfaces
       surf_type = list of stype of method for the creation of the surface
                   -- regulare_grid (u and v lines)
                   -- skin (u lines)       
    """
    #
    import start as start
    cfg                     = start.start_cfg(filename=filename)
    #
    #
    for isurface in range(0,cfg.nsurf):
        surf_type=cfg.surf_type[isurface]
        if surf_type == 'regular_grid':
            surface_regular_grid(isurface,cfgname=filename)
        elif surf_type == 'skin':
            surface_skin(isurface,cfgname=filename)

def surface_regular_grid(isurface=0,cfgname=None):
    """
    create an acis surface from a regular lon/lat/z grid
    """
    import sys,os
    from math import sqrt
    from utilities import geo2utm
    import start as start
    #
    #
    cfg                     = start.start_cfg(cfgname)
    numpy                   = start.start_numpy()
    #
    def create_line_u(ind,n,step,data,unit):
        last_curve_store=cubit.get_last_id("curve")
        command='create curve spline '
        for i in range(0,n):
            if i%step == 0:
                lon,lat,z=data[i+ind][0],data[i+ind][1],data[i+ind][2]
                x,y=geo2utm(lon,lat,unit)
                txt=' Position ' +   str(x)  +' '+  str(y) +' '+  str(z)
                command=command+txt
                #print command
        cubit.silent_cmd(command)
        last_curve=cubit.get_last_id("curve")
        if last_curve != last_curve_store:
            return last_curve
        else:
            return 0
    
    def create_line_v(ind,n,n2,step,data,unit):
        last_curve_store=cubit.get_last_id("curve")
        command='create curve spline '
        for i in range(0,n):
            if i%step == 0:
                lon,lat,z=data[n2*i+ind][0],data[n2*i+ind][1],data[n2*i+ind][2]
                x,y=geo2utm(lon,lat,unit)
                txt=' Position ' +   str(x)  +' '+  str(y) +' '+  str(z)
                command=command+txt
                #print command
        cubit.silent_cmd(command)
        last_curve=cubit.get_last_id("curve")
        if last_curve != last_curve_store:
            return last_curve
        else:
            return 0
    #
    #
    cubit.cmd("reset")
    #
    position=True
    #
    #
    nu= cfg.num_x[isurface]
    nv= cfg.num_y[isurface]
    ustep= cfg.xstep[isurface]
    vstep= cfg.ystep[isurface]
    exag=1.
    unit=cfg.unit2[isurface]
    #
    #
    data=numpy.loadtxt(cfg.surface_name[isurface])
    if len(data) > 100:
        command = "set echo off"
        cubit.cmd(command)
        command = "set journal off"
        cubit.cmd(command)
    #
    u_curve=[]
    v_curve=[]
    #
    for iv in range(0,nv):
        if iv%vstep == 0.:
            u=create_line_u(iv*(nu),nu,ustep,data,unit)
            u_curve.append(u)
    for iu in range(0,nu):
        if iu%ustep == 0.:
            v=create_line_v(iu,nv,nu,ustep,data,unit)
            v_curve.append(v)
    #
    umax=max(u_curve)
    umin=min(u_curve)
    vmax=max(v_curve)
    vmin=min(v_curve)
    cubitcommand= 'create surface net u curve '+ str( umin )+' to '+str( umax )+ ' v curve '+ str( vmin )+ ' to '+str( vmax )+' heal'
    cubit.cmd(cubitcommand)
    command = "del curve all"
    cubit.cmd(command)
    suff=cfg.surface_name[isurface].split('/')
    command = "save as '"+cfg.working_dir+"/surf_"+suff[-1]+".cub' overwrite"
    cubit.cmd(command)
    #
    #
    #        
    cubit.cmd("set info "+cfg.cubit_info)
    cubit.cmd("set echo "+cfg.echo_info)
    cubit.cmd("set journal "+cfg.jou_info)

def surface_skin(isurface=0,cfgname=None):
    """
    create an acis surface interpolating no-intersecting lines
    """
    import sys,os
    from math import sqrt
    from utilities import geo2utm
    import start as start
    #
    #
    cubit                   = start.start_cubit()
    cfg                     = start.start_cfg(cfgname)
    #
    def define_next_line(directionx,directiony,n,data):
        ndata=len(data)
        command=''
        ind=n
        try:
            record=data[ind]
        except:
            return False,False
        try:
            x,y,z=map(float,record.split())
        except:
            return False,False
        txt=' Position ' +   record
        command=command+txt
        x_store,y_store,z_store = x,y,z
        icount=1
        while True:
                    ind+=1
                    if ind >= ndata: return ind,command
                    record=data[ind]
                    try:
                        x,y,z=map(float,record.split())
                    except:
                        return ind,command
                    dx,dy = x-x_store,y-y_store        
                    if  directionx == 0 and dy/abs(dy) * directiony >= 0:
                        txt=' Position ' +   record
                        command=command+txt
                        icount+=1
                        x_store,y_store,z_store = x,y,z
                    elif  directiony == 0 and dx/abs(dx) == directionx :
                        txt=' Position ' +   record
                        command=command+txt
                        icount+=1
                        x_store,y_store,z_store = x,y,z
                    else:
                        if icount==1:
                           x,y,z=x_store+1e-4*directionx,y_store+1e-4*directiony,z_store
                           txt=' Position ' +str(x)+ ' '+str(y)+ ' '+str(z)
                           command=command+txt
                        return ind,command    
    def create_line(position):
        if position:
            last_curve_store=cubit.get_last_id("curve")
            command='create curve spline '+position
            cubit.silent_cmd(command)
            last_curve=cubit.get_last_id("curve")
            if last_curve != last_curve_store:
                return last_curve
            else:
                return False
        else:
            return False
                        
    command = "reset"
    cubit.cmd(command)
    #
    position=True
    #
    try:
         grdfile = open(cfg.surface_name[isurface], 'r')
    except:
         raise NameError, 'No such file or directory: '+  str( cfg.surface_name[isurface] )
    #
    directionx=cfg.directionx[isurface]
    directiony=cfg.directiony[isurface]
    step=cfg.step[isurface]
    position=True
    curveskin=[]
    count_line=0
    data=grdfile.read().split('\n')
    ndata=len(data)
    n=0
    #
    #
    command = "set echo off"
    cubit.cmd(command)
    command = "set journal off"
    cubit.cmd(command)
    command = "set info off"
    cubit.cmd(command)         
    #
    while position:
          index,position=define_next_line(directionx,directiony,n,data)
          if n%step == 0:
              curve=create_line(position)
              if curve: curveskin.append(curve)
          elif n%step != 0 and not position:
              curve=create_line(position)
              if curve: curveskin.append(curve)
          n=index
    umax=max(curveskin)
    umin=min(curveskin)
    print 'create surface skin curve '+ str( umin )+' to '+str( umax )
    cubitcommand= 'create surface skin curve '+ str( umin )+' to '+str( umax )
    cubit.cmd(cubitcommand)
    command = "del curve all"
    cubit.cmd(command)
    last_surface=cubit.get_last_id("surface")
    command = "regularize surf "+str(last_surface)
    cubit.cmd(command)
    #
    suff=cfg.surface_name[isurface].split('/')
    command = "save as '"+cfg.working_dir+"/surf_"+suff[-1]+".cub' overwrite"
    cubit.cmd(command)
    #
    #
    #        
    cubit.cmd("set info "+cfg.cubit_info)
    cubit.cmd("set echo "+cfg.echo_info)
    cubit.cmd("set journal "+cfg.jou_info)
    
    
def plane(cfg):
    import sys,os
    from math import sqrt
    from utilities import geo2utm
    import start as start
    #
    #
    cubit                   = start.start_cubit()
    #
    #
    command = "reset"
    cubit.cmd(command)
    #
    #
    for p in [cfg.x1,cfg.x2,cfg.x3,cfg.x4]:
        x_current,y_current=geo2utm(p[0],p[1],cfg.unit)
        cubitcommand= 'create vertex '+ str( x_current )+ ' ' + str( y_current) +' '+ str( p[2] )
        cubit.cmd(cubitcommand)
    #
    cubitcommand= 'create surface vertex 1 2 3 4'
    cubit.cmd(cubitcommand)
    command = "del vertex all"
    cubit.cmd(command)
    command = "save as 'plane.cub' overwrite"
    cubit.cmd(command)
    #
    #
    #        
    cubit.cmd("set info "+cfg.cubit_info)
    cubit.cmd("set echo "+cfg.echo_info)
    cubit.cmd("set journal "+cfg.jou_info)