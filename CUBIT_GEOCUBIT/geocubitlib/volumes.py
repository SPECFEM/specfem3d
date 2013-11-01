#############################################################################
# volumes.py                                                    
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

    
def volumes(filename=None):
    """create the volumes"""
    import start as start
    print'volume'
    cfg                     = start.start_cfg(filename=filename)
    print cfg
    #
    if cfg.volume_type == 'layercake_volume_ascii_regulargrid_regularmap':
            layercake_volume_ascii_regulargrid_mpiregularmap(filename=filename)
    elif cfg.volume_type == 'layercake_volume_fromacis_mpiregularmap':
            layercake_volume_fromacis_mpiregularmap(filename=filename)
    elif cfg.volume_type == 'verticalsandwich_volume_ascii_regulargrid_mpiregularmap':
            layercake_volume_ascii_regulargrid_mpiregularmap(filename=filename,verticalsandwich=True)
    
def ordering_surfaces(list_surfaces):
    list_z=[]
    for s in list_surfaces:
        _,_,z=cubit.get_center_point("surface",s)
        list_z.append(z)
    ord_list_surfaces=[s for s,z in sorted(zip(list_surfaces,list_z),key=lambda x: (x[1]))]
    return ord_list_surfaces


def onlyvolumes():
    list_surfaces=cubit.parse_cubit_list("surface","all")
    ord_list_surfaces=ordering_surfaces(list_surfaces)
    for s1,s2 in zip(ord_list_surfaces[:-1],ord_list_surfaces[1:]):
        create_volume(s1,s2,method=None)
    cubitcommand= 'del surface all'
    cubit.cmd(cubitcommand)
    list_vol=cubit.parse_cubit_list("volume","all")
    if len(list_vol) > 1:     
        cubitcommand= 'imprint volume all'
        cubit.cmd(cubitcommand)
        cubitcommand= 'merge all'
        cubit.cmd(cubitcommand)


def surfaces(filename=None):
    """create the volumes"""
    import start as start
    print'volume'
    cfg                     = start.start_cfg(filename=filename)
    print cfg
    if cfg.volume_type == 'layercake_volume_ascii_regulargrid_regularmap':
            layercake_volume_ascii_regulargrid_mpiregularmap(filename=filename,onlysurface=True)
    elif cfg.volume_type == 'layercake_volume_fromacis_mpiregularmap':
            layercake_volume_fromacis_mpiregularmap(filename=filename,onlysurface=True)
    elif cfg.volume_type == 'verticalsandwich_volume_ascii_regulargrid_mpiregularmap':
            layercake_volume_ascii_regulargrid_mpiregularmap(filename=filename,verticalsandwich=True,onlysurface=True)
    
    
    
    
    
    


def layercake_volume_ascii_regulargrid_mpiregularmap(filename=None,verticalsandwich=False,onlysurface=False):
    import sys
    import start as start
    #
    mpiflag,iproc,numproc,mpi   = start.start_mpi()
    #
    numpy                       = start.start_numpy()
    cfg                         = start.start_cfg(filename=filename)                       
    
    from utilities import geo2utm, savegeometry,savesurf,cubit_command_check
    
    from math import sqrt
    #
    try:
        mpi.barrier()
    except:
        pass
    #
    #
    command = "comment '"+"PROC: "+str(iproc)+"/"+str(numproc)+" '"
    cubit_command_check(iproc,command,stop=True)
    if verticalsandwich: cubit.cmd("comment 'Starting Vertical Sandwich'")
    #
    #get icpuy,icpux values
    if mpiflag:
        icpux = iproc % cfg.nproc_xi
        icpuy = int(iproc / cfg.nproc_xi)
    else:
        icpuy=int(cfg.id_proc/cfg.nproc_xi)
        icpux=cfg.id_proc%cfg.nproc_xi
    
    #
    if  cfg.geometry_format == 'ascii':
        #for the original surfaces
        #number of points in the files that describe the topography
        import local_volume
        if cfg.localdir_is_globaldir:
            if iproc == 0 or not mpiflag:
                coordx_0,coordy_0,elev_0,nx_0,ny_0=local_volume.read_grid(filename)
                print 'end of reading grd files '+str(nx_0*ny_0)+ ' points'
            else:
                pass
            if iproc == 0 or not mpiflag:
                coordx=mpi.bcast(coordx_0)
            else:
                coordx=mpi.bcast()
            if iproc == 0 or not mpiflag:
                coordy=mpi.bcast(coordy_0)
            else:
                coordy=mpi.bcast()
            if iproc == 0 or not mpiflag:
                elev=mpi.bcast(elev_0)
            else:
                elev=mpi.bcast()
            if iproc == 0 or not mpiflag:
                nx=mpi.bcast(nx_0)
            else:
                nx=mpi.bcast()       
            if iproc == 0 or not mpiflag:
                ny=mpi.bcast(ny_0)
            else:
                ny=mpi.bcast()
        else:
            coordx,coordy,elev,nx,ny=local_volume.read_grid(filename)
        print str(iproc)+ ' end of receving grd files '
        nx_segment=int(nx/cfg.nproc_xi)+1
        ny_segment=int(ny/cfg.nproc_eta)+1
        
    elif cfg.geometry_format=='regmesh': # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        if cfg.depth_bottom != cfg.zdepth[0]:
            if iproc == 0: print 'the bottom of the block is at different depth than depth[0] in the configuration file'
        nx= cfg.nproc_xi+1
        ny= cfg.nproc_eta+1
        nx_segment=2
        ny_segment=2
        #if iproc == 0: print nx,ny,cfg.cpux,cfg.cpuy
        xp=(cfg.xmax-cfg.xmin)/float((nx-1))
        yp=(cfg.ymax-cfg.ymin)/float((ny-1))
        #
        elev=numpy.zeros([nx,ny,cfg.nz],float)
        coordx=numpy.zeros([nx,ny],float)
        coordy=numpy.zeros([nx,ny],float)
        #
        #
        xlength=(cfg.xmax-cfg.xmin)/float(cfg.nproc_xi) #length of x slide for chunk
        ylength=(cfg.ymax-cfg.ymin)/float(cfg.nproc_eta) #length of y slide for chunk
        nelem_chunk_x=1    
        nelem_chunk_y=1
        ivxtot=nelem_chunk_x+1
        ivytot=nelem_chunk_y+1 
        xstep=xlength #distance between vertex on x
        ystep=ylength
        for i in range(0,cfg.nz):
            elev[:,:,i] = cfg.zdepth[i]
        
        icoord=0
        for iy in range(0,ny):
            for ix in range(0,nx):
                icoord=icoord+1
                coordx[ix,iy]=cfg.xmin+xlength*(ix)
                coordy[ix,iy]=cfg.ymin+ylength*(iy)
        
        #print coordx,coordy,nx,ny
    #
    print 'end of building grid '+str(iproc)
    print 'number of point: ', len(coordx)*len(coordy)
    #
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #for each processor
    #
    nxmin_cpu=(nx_segment-1)*(icpux)
    nymin_cpu=(ny_segment-1)*(icpuy)
    nxmax_cpu=min(nx-1,(nx_segment-1)*(icpux+1))
    nymax_cpu=min(ny-1,(ny_segment-1)*(icpuy+1))
    #if iproc == 0:
    #    print nx_segment,ny_segment,nx,ny
    #    print icpux,icpuy,nxmin_cpu,nxmax_cpu
    #    print icpux,icpuy,nymin_cpu,nymax_cpu
    #    print coordx[0,0],coordx[nx-1,ny-1]
    #    print coordy[0,0],coordy[nx-1,ny-1]
    #
    #
    icurve=0
    isurf=0
    ivertex=0
    #
    #create vertex
    for inz in range(0,cfg.nz):
        if cfg.sea and inz==cfg.nz-1: #sea layer
            sealevel=True
            bathymetry=False
        elif cfg.sea and inz==cfg.nz-2: #bathymetry layer
            sealevel=False
            bathymetry=True
        else:
            sealevel=False
            bathymetry=False
        print sealevel,bathymetry
        
        if  cfg.bottomflat and inz == 0: #bottom layer
                #
                if cfg.geometry_format == 'ascii':
                    lv=cubit.get_last_id("vertex")     
                    
                    x_current,y_current=(coordx[nxmin_cpu,nymin_cpu],coordy[nxmin_cpu,nymin_cpu])
                    cubitcommand= 'create vertex '+ str( x_current )+ ' ' + str( y_current) +' '+ str( cfg.depth_bottom )
                    cubit.cmd(cubitcommand)
                    #
                    x_current,y_current=(coordx[nxmin_cpu,nymax_cpu],coordy[nxmin_cpu,nymax_cpu])
                    cubitcommand= 'create vertex '+ str( x_current )+ ' ' + str( y_current) +' '+ str( cfg.depth_bottom )
                    cubit.cmd(cubitcommand)                                                                              
                    #
                    x_current,y_current=(coordx[nxmax_cpu,nymax_cpu],coordy[nxmax_cpu,nymax_cpu])
                    cubitcommand= 'create vertex '+ str( x_current )+ ' ' + str( y_current) +' '+ str( cfg.depth_bottom )
                    cubit.cmd(cubitcommand)
                    #
                    x_current,y_current=(coordx[nxmax_cpu,nymin_cpu],coordy[nxmax_cpu,nymin_cpu])
                    cubitcommand= 'create vertex '+ str( x_current )+ ' ' + str( y_current) +' '+ str( cfg.depth_bottom )
                    cubit.cmd(cubitcommand)
                    #
                    lv2=cubit.get_last_id("vertex")     
                    
                    cubitcommand= 'create surface vertex '+str(lv+1)+' to '+str(lv2)
                    cubit.cmd(cubitcommand)
                    #
                    isurf = isurf + 1
                else:
                    lv=cubit.get_last_id("vertex") 
                    x_current,y_current=geo2utm(coordx[nxmin_cpu,nymin_cpu],coordy[nxmin_cpu,nymin_cpu],'utm')
                    cubitcommand= 'create vertex '+ str( x_current )+ ' ' + str( y_current) +' '+ str( cfg.depth_bottom )
                    cubit.cmd(cubitcommand)
                    #
                    x_current,y_current=geo2utm(coordx[nxmin_cpu,nymax_cpu],coordy[nxmin_cpu,nymax_cpu],'utm')
                    cubitcommand= 'create vertex '+ str( x_current )+ ' ' + str( y_current) +' '+ str( cfg.depth_bottom )
                    cubit.cmd(cubitcommand)                                                                              
                    #
                    x_current,y_current=geo2utm(coordx[nxmax_cpu,nymax_cpu],coordy[nxmax_cpu,nymax_cpu],'utm')
                    cubitcommand= 'create vertex '+ str( x_current )+ ' ' + str( y_current) +' '+ str( cfg.depth_bottom )
                    cubit.cmd(cubitcommand)
                    #
                    x_current,y_current=geo2utm(coordx[nxmax_cpu,nymin_cpu],coordy[nxmax_cpu,nymin_cpu],'utm')
                    cubitcommand= 'create vertex '+ str( x_current )+ ' ' + str( y_current) +' '+ str( cfg.depth_bottom )
                    cubit.cmd(cubitcommand)
                    #
                    lv2=cubit.get_last_id("vertex") 
                    cubitcommand= 'create surface vertex '+str(lv+1)+' to '+str(lv2)
                    cubit.cmd(cubitcommand)
                    #
                    isurf = isurf + 1
        else:
            if cfg.geometry_format == 'regmesh':
                zvertex=cfg.zdepth[inz]
                lv=cubit.get_last_id("vertex")                        
                x_current,y_current=geo2utm(coordx[nxmin_cpu,nymin_cpu],coordy[nxmin_cpu,nymin_cpu],'utm')
                cubitcommand= 'create vertex '+ str( x_current )+ ' ' + str( y_current) +' '+ str( zvertex )
                cubit.cmd(cubitcommand)
                #
                x_current,y_current=geo2utm(coordx[nxmin_cpu,nymax_cpu],coordy[nxmin_cpu,nymax_cpu],'utm')
                cubitcommand= 'create vertex '+ str( x_current )+ ' ' + str( y_current) +' '+ str( zvertex )
                cubit.cmd(cubitcommand)                                                                              
                #
                x_current,y_current=geo2utm(coordx[nxmax_cpu,nymax_cpu],coordy[nxmax_cpu,nymax_cpu],'utm')
                cubitcommand= 'create vertex '+ str( x_current )+ ' ' + str( y_current) +' '+ str( zvertex )
                cubit.cmd(cubitcommand)
                #
                x_current,y_current=geo2utm(coordx[nxmax_cpu,nymin_cpu],coordy[nxmax_cpu,nymin_cpu],'utm')
                cubitcommand= 'create vertex '+ str( x_current )+ ' ' + str( y_current) +' '+ str( zvertex )
                cubit.cmd(cubitcommand)
                #
                cubitcommand= 'create surface vertex '+str(lv+1)+' '+str(lv+2)+' '+str(lv+3)+' '+str(lv+4)
                cubit.cmd(cubitcommand)
                #
                isurf = isurf + 1
            elif cfg.geometry_format == 'ascii':
                
                vertex=[]
                
                for iy in range(nymin_cpu,nymax_cpu+1):
                    ivx=0
                    for ix in range(nxmin_cpu,nxmax_cpu+1):
                        zvertex=elev[ix,iy,inz]
                        #zvertex=adjust_sea_layers(zvertex,sealevel,bathymetry,cfg)
                        x_current,y_current=(coordx[ix,iy],coordy[ix,iy])
                        #
                        vertex.append(' Position '+ str( x_current ) +' '+ str( y_current )+' '+ str( zvertex ) )
                #
                print 'proc',iproc, 'vertex list created....',len(vertex)
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
                    #print command
                    uline.append( cubit.get_last_id("curve") )
                for ix in range(0,nxmax_cpu-nxmin_cpu+1):
                    positiony=''
                    for iy in range(0,nymax_cpu-nymin_cpu+1):
                        positiony=positiony+vertex[ix+iy*(nxmax_cpu-nxmin_cpu+1)]
                    command='create curve spline '+positiony
                    cubit.cmd(command)
                    #print command
                    vline.append( cubit.get_last_id("curve") )
                #
                cubit.cmd("set info "+cfg.cubit_info)
                cubit.cmd("set echo "+cfg.echo_info)
                cubit.cmd("set journal "+cfg.jou_info)
                #
                #
                print 'proc',iproc, 'lines created....',len(uline),'*',len(vline)
                umax=max(uline)
                umin=min(uline)
                vmax=max(vline)
                vmin=min(vline)
                ner=cubit.get_error_count()
                cubitcommand= 'create surface net u curve '+ str( umin )+' to '+str( umax )+ ' v curve '+ str( vmin )+ ' to '+str( vmax )+' heal'
                cubit.cmd(cubitcommand)
                ner2=cubit.get_error_count()
                if ner == ner2: 
                    command = "del curve all"
                    cubit.cmd(command)
                    isurf=isurf+1
                #
            else:
                raise NameError, 'error, check geometry_format, it should be ascii or regmesh'   #
                #
        cubitcommand= 'del vertex all'
        cubit.cmd(cubitcommand)
    if cfg.save_surface_cubit:
        savegeometry(iproc=iproc,surf=True,filename=filename)
    #
    #
    #!create volume
    if not onlysurface:
        if cfg.nz == 1:
            nsurface=2
        else:
            nsurface=cfg.nz
        for inz in range(1,nsurface):
            ner=cubit.get_error_count()
            create_volume(inz,inz+1,method=cfg.volumecreation_method)
            ner2=cubit.get_error_count()
        if ner == ner2 and not cfg.debug_geometry:
            #cubitcommand= 'del surface 1 to '+ str( cfg.nz )
            cubitcommand= 'del surface all'
            cubit.cmd(cubitcommand)
            list_vol=cubit.parse_cubit_list("volume","all")
            if len(list_vol) > 1:     
                cubitcommand= 'imprint volume all'
                cubit.cmd(cubitcommand)
                cubitcommand= 'merge all'
                cubit.cmd(cubitcommand)
            #ner=cubit.get_error_count()
            #cubitcommand= 'composite create curve in vol all'
            #cubit.cmd(cubitcommand)
    savegeometry(iproc,filename=filename)
    #if cfg.geological_imprint:
    #    curvesname=[cfg.outlinebasin_curve,cfg.transition_curve,cfg.faulttrace_curve]
    #    outdir=cfg.working_dir
    #    imprint_topography_with_geological_outline(curvesname,outdir)
    #
    #        
    cubit.cmd("set info "+cfg.cubit_info)
    cubit.cmd("set echo "+cfg.echo_info)
    cubit.cmd("set journal "+cfg.jou_info)
    
    
def layercake_volume_fromacis_mpiregularmap(filename=None):
    import sys
    import start as start
    #
    mpiflag,iproc,numproc,mpi   = start.start_mpi()
    #
    cfg                         = start.start_cfg(filename=filename)                       
    #
    from utilities import geo2utm, savegeometry
    #
    from math import sqrt
    #
    try:
        mpi.barrier()
    except:
        pass
    #
    #
    command = "comment '"+"PROC: "+str(iproc)+"/"+str(numproc)+" '"
    cubit.cmd(command)
    #
    #get the limit of the volume considering the cpu
    def xwebcut(x):
        command='create planar surface with plane xplane offset '+str(x)
        cubit.cmd(command)
        last_surface=cubit.get_last_id("surface")
        command="webcut volume all tool volume in surf "+str(last_surface)
        cubit.cmd(command)
        command="del surf "+str(last_surface)
        cubit.cmd(command)
        
    def ywebcut(x):
        command='create planar surface with plane yplane offset '+str(x)
        cubit.cmd(command)
        last_surface=cubit.get_last_id("surface")
        command="webcut volume all tool volume in surf "+str(last_surface)
        cubit.cmd(command)
        command="del surf "+str(last_surface)
        cubit.cmd(command)
        
    def translate2zero():
        ss=cubit.parse_cubit_list('surface','all')
        box = cubit.get_total_bounding_box("surface", ss)
        xmin=box[0]
        ymin=box[3]
        cubit.cmd('move surface all x '+str(-1*xmin)+' y '+str(-1*ymin))
        return xmin,ymin
        
    def translate2original(xmin,ymin):
        cubit.cmd('move surface all x '+str(xmin)+' y '+str(ymin))
        
    if mpiflag:
        icpux = iproc % cfg.nproc_xi
        icpuy = int(iproc / cfg.nproc_xi)
    else:
        icpuy=int(cfg.id_proc/cfg.nproc_xi)
        icpux=cfg.id_proc%cfg.nproc_xi
    #
    ner=cubit.get_error_count()
    #
    icurve=0
    isurf=0
    ivertex=0
    #
    xlength=(cfg.xmax-cfg.xmin)/float(cfg.cpux) #length of x slide for chunk
    ylength=(cfg.ymax-cfg.ymin)/float(cfg.cpuy) #length of y slide for chunk
    xmin_cpu=cfg.xmin+(xlength*(icpux))
    ymin_cpu=cfg.ymin+(ylength*(icpuy))
    xmax_cpu=xmin_cpu+xlength
    ymax_cpu=ymin_cpu+ylength
    #
    #importing the surfaces
    for inz in range(cfg.nz-2,-2,-1):
        if cfg.bottomflat and inz==-1:
            command = "create planar surface with plane zplane offset "+str(cfg.depth_bottom)
            cubit.cmd(command)
        else:
            command = "import cubit '"+cfg.filename[inz]+"'"
            cubit.cmd(command)
            
            
    #translate
    xmin,ymin=translate2zero()
    print 'translate ...', -xmin,-ymin
    xmin_cpu=xmin_cpu-xmin
    ymin_cpu=ymin_cpu-ymin
    xmax_cpu=xmax_cpu-xmin
    ymax_cpu=ymax_cpu-ymin
    
    ss=cubit.parse_cubit_list('surface','all')
    box = cubit.get_total_bounding_box("surface", ss)
    print 'dimension... ', box
    #cutting the surfaces
    xwebcut(xmin_cpu)
    xwebcut(xmax_cpu)
    ywebcut(ymin_cpu)
    ywebcut(ymax_cpu)
    #
    list_surface_all=cubit.parse_cubit_list("surface","all")
    #condisidering only the surfaces inside the boundaries
    dict_surf={}
    for isurf in list_surface_all:
        p=cubit.get_center_point("surface",isurf)
        if p[0] < xmin_cpu or p[0] > xmax_cpu or p[1] > ymax_cpu or p[1] < ymin_cpu:
            command = "del surf "+str(isurf)
            cubit.cmd(command)
        else:
            dict_surf[str(isurf)]=p[2]
    z=dict_surf.values()
    z.sort()
    list_surf=[]
    for val in z:
        isurf=[k for k, v in dict_surf.iteritems() if v == val][0]
        list_surf.append(int(isurf))
    #
    
    #lofting the volume
    for i,j in zip(list_surf,list_surf[1:]):
        ner=cubit.get_error_count()
        create_volume(i,j,method=cfg.volumecreation_method)
        #cubitcommand= 'create volume loft surface '+ str(i)+' '+str(j)
        #cubit.cmd(cubitcommand)
        ner2=cubit.get_error_count()
    #
    translate2original(xmin,ymin)
    
    
    if ner == ner2:
        cubitcommand= 'del surface all'
        cubit.cmd(cubitcommand)
        #
        #
        #cubitcommand= 'composite create curve in vol all'
        #cubit.cmd(cubitcommand)
        list_vol=cubit.parse_cubit_list("volume","all")
        if len(list_vol) > 1:     
            cubitcommand= 'imprint volume all'
            cubit.cmd(cubitcommand)
            #cubit_error_stop(iproc,cubitcommand,ner)
            #
            cubitcommand= 'merge all'
            cubit.cmd(cubitcommand)
    #
    savegeometry(iproc,filename=filename)
    #if cfg.geological_imprint:
    #    curvesname=[cfg.outlinebasin_curve,cfg.transition_curve,cfg.faulttrace_curve]
    #    outdir=cfg.working_dir
    #    imprint_topography_with_geological_outline(curvesname,outdir)


def hor_distance(c1,c2):
    p1=cubit.get_center_point("curve", c1)
    p2=cubit.get_center_point("curve", c2)
    d=(p1[0]-p2[0])**2+(p1[1]-p2[1])**2
    return d



def coupling_curve(lcurve1,lcurve2):
    """get the couple of curve that we  use to get the skin surface"""
    import operator
    couples=[]
    for c1 in lcurve1:
        distance=[]
        for c2 in lcurve2:
            d=hor_distance(c1,c2)
            distance.append(d)
        min_index, min_value = min(enumerate(distance), key=operator.itemgetter(1))
        couples.append((c1,lcurve2[min_index]))
    return couples

def create_volume(surf1,surf2,method='loft'):
    if method == 'loft':
        cmd='create volume loft surface '+str(surf1)+' '+str(surf2)
        cubit.cmd(cmd)
    else:
        lcurve1=cubit.get_relatives("surface",surf1,"curve")
        lcurve2=cubit.get_relatives("surface",surf2,"curve")
        couples=coupling_curve(lcurve1,lcurve2)
        is_start=cubit.get_last_id('surface')+1
        for cs in couples:
            cmd='create surface skin curve '+str(cs[1])+' '+str(cs[0])
            cubit.cmd(cmd)
        is_stop=cubit.get_last_id('surface')
        cmd="create volume surface "+str(surf1)+' '+str(surf2)+' '+str(is_start)+' to '+str(is_stop)+"  heal keep"
        cubit.cmd(cmd)












#def imprint_topography_with_geological_outline(curvesname,outdir='.'):
#    import sys,os
#    from sets import Set
#    #
#    from utilities import load_curves,project_curves,get_v_h_list
#    
#    list_vol=cubit.parse_cubit_list("volume","all")
#    surf_or,surf_vertical,list_curve_or,list_curve_vertical,bottom,top=get_v_h_list(list_vol)
#    
#    
#    
#    outlinebasin_curve=load_curves(curvesname[0])
#    transition_curve=load_curves(curvesname[1])
#    faulttrace_curve=load_curves(curvesname[2])
#    
#    curves=[]
#    if outlinebasin_curve: curves=curves+outlinebasin_curve
#    if transition_curve: curves=curves+transition_curve
#    if faulttrace_curve: curves=curves+faulttrace_curve
#    
#    if curves:
#            command='imprint tolerant surface '+str(top)+' with curve '+' '.join(str(x) for x in curves)+'  merge'
#            cubit.cmd(command)
#             
#    
#    command = "merge surf all"
#    cubit.cmd(command)
#    command = "compress vol all"
#    cubit.cmd(command)
#    command = "compress surf all"
#    cubit.cmd(command)
#    command = "save as '"+outdirs+"/"+"imprinted_vol_"+str(iproc)+".cub' overwrite"
#    cubit.cmd(command)    
#
#
