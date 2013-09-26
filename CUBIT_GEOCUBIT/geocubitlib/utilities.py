#############################################################################
# utilities.py                                                    
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


def get_cubit_version():
    v=cubit.get_version()
    try:
        v=float(v[0:4])
    except:
        v=float(v[0:2])
    return v
    

def snapshot(name=None,i=0,viewnumber=1):
    """
    it takes a snapshot of the figure, following the predefined view position.
    view 1: vector 1 1 1 z up
    """
    if name is None:
        name='snapshot_'+str(i)
    i=i+1
    if viewnumber == 1:
        command = "at 0"
        cubit.cmd(command)
        command = "from 1 1 1"
        cubit.cmd(command)
        command = "up 0 0 1"
        cubit.cmd(command)
    cubit.cmd('graphics autocenter on')
    cubit.cmd("zoom reset")
    command = "hardcopy '"+name+".png' png"
    cubit.cmd(command)
    return i

def cubit_command_check(iproc,command,stop=True):
    """
    Run a cubit command, checking if it performs correctly. 
    If the command fails, it writes the result on a file "error_[processor number]" and stop the meshing process
    
    iproc = process number
    command = cubit command
    stop = if command fails, stop the meshing process (Default: True)
    
    return status variable (0 ok, -1 fail)
    
    """
    er=cubit.get_error_count()
    cubit.cmd(command)
    ner=cubit.get_error_count()
    flag=0
    if ner > er:
        text='"Proc: '+str(iproc)+' ERROR '+str(command)+' number of error '+str(er)+'/'+str(ner)+'"'
        cubitcommand = 'comment '+text
        cubit.cmd(cubitcommand)
        f=open('error_'+str(iproc)+'.log','a')
        f.write("CUBIT ERROR: \n"+text)
        f.close()
        if stop: raise Exception("CUBIT ERROR: "+text)
        flag=-1
    return flag

    
    

def cubit_error_stop(iproc,command,ner):
    """obsolete"""
    er=cubit.get_error_count()
    if er > ner: 
       text='"Proc: '+str(iproc)+' ERROR '+str(command)+' number of error '+str(er)+'/'+str(ner)+'"'
       cubitcommand = 'comment '+text
       cubit.cmd(cubitcommand)
       raise NameError, text

def cubit_error_continue(iproc,command,n_er):
    """obsolete"""
    er=cubit.get_error_count()
    if  er >= n_er: 
        text='"Proc: '+str(iproc)+' ERROR continue '+str(command)+' number of error '+str(er)+' '+str(n_er)+'"'
        cubit.cmd('comment '+ text)
        print 'error: ',text
    return er
    
def savemesh(mpiflag,iproc=0,filename=None):
    import start as start
    cfg                         = start.start_cfg(filename=filename)
    mpiflag,iproc,numproc,mpi   = start.start_mpi()
    
    def runsave(meshfile,iproc,filename=None):
        import start as start
        cubit                   = start.start_cubit()
        cfg                         = start.start_cfg(filename=filename)
        flag=0
        ner=cubit.get_error_count()
        cubitcommand= 'save as "'+ cfg.output_dir+'/'+meshfile+'.cub'+ '" overwrite' 
        cubit.cmd(cubitcommand)
        ner2=cubit.get_error_count()
        if ner == ner2:
            cubitcommand= 'export mesh "'+ cfg.output_dir+'/'+meshfile+'.e'+ '" dimension 3 block all overwrite' 
            cubit.cmd(cubitcommand)
            ner2=cubit.get_error_count()                                                    
        if ner == ner2:
            flag=1
        return flag
    
    
    meshfile='mesh_vol_'+str(iproc)
    
    flagsaved=0
    infosave=(iproc,flagsaved)
    
    mpi.barrier()
    total_saved=mpi.allgather(flagsaved)
    if isinstance(total_saved,int): total_saved=[total_saved]
    
    ind=0
    saving=True
    while saving:
        if len(total_saved) != sum(total_saved):
            #
            if not flagsaved: 
                flagsaved=runsave(meshfile,iproc,filename=filename)
                if flagsaved:
                    infosave=(iproc,flagsaved)        
                    if numproc > 1:
                        f=open('mesh_saved'+str(iproc),'w')
                        f.close()
            mpi.barrier()
            total_saved=mpi.allgather(flagsaved)
            if isinstance(total_saved,int): total_saved=[total_saved]
            ind=ind+1
        else:
            saving=False
        if ind > len(total_saved)+10: saving=False
        print sum(total_saved),'/',len(total_saved),' saved'
    
    info_total_saved=mpi.allgather(infosave)
    if isinstance(info_total_saved,int): info_total_saved=[info_total_saved]
    
    if iproc==0:
        f=open('mesh_saving.log','w')
        f.write('\n'.join(str(x) for x in info_total_saved))                
        f.close()                           
            
    f=open(cfg.output_dir+'/'+'blocks_'+str(iproc).zfill(5),'w')
    blocks=cubit.get_block_id_list()
    
    for block in blocks:
        name=cubit.get_exodus_entity_name('block',block)
        element_count = cubit.get_exodus_element_count(block, "block")
        nattrib=cubit.get_block_attribute_count(block)
        attr=[cubit.get_block_attribute_value(block,x) for x in range(0,nattrib)]
        ty=cubit.get_block_element_type(block)
        f.write(str(block)+' ; '+name+' ; nattr '+str(nattrib)+' ; '+' '.join(str(x) for x in attr)+' ; '+ty+' '+str(element_count)+'\n')
    f.close()
    
    import quality_log
    f=open(cfg.output_dir+'/'+'quality_'+str(iproc).zfill(5),'w')
    max_skewness,min_length=quality_log.quality_log(f)
    f.close()
    
    
    count_hex=[cubit.get_hex_count()]
    count_node=[cubit.get_node_count()]
    max_skew=[(iproc,max_skewness)]
    min_l=[(iproc,min_length)]
    
    mpi.barrier()
    total_min_l=mpi.gather(min_l)
    total_hex=mpi.gather(count_hex)        
    total_node=mpi.gather(count_node)      
    total_max_skew=mpi.gather(max_skew)    
    
    
    mpi.barrier()                          
    if iproc == 0:
        min_total_min_l=min([ms[1] for ms in total_min_l])
        max_total_max_skew=max([ms[1] for ms in total_max_skew])
        sum_total_node=sum(total_node)
        sum_total_hex=sum(total_hex)
        
        totstat_file=open(cfg.output_dir+'/totstat.log','w')
        text='hex total number,node total number,max skew, min length\n'
        totstat_file.write(text)
        
        text=str(sum_total_hex)+' , '+str(sum_total_node)+' , '+str(max_total_max_skew)+' , '+str(min_total_min_l)+'\n'
        totstat_file.write(text)
        
        totstat_file.write(str(total_max_skew))    
        totstat_file.close()
    
    print 'meshing process end... proc ',iproc 

def importgeometry(geometryfile,iproc=0,filename=None):
    import start as start
    cfg                         = start.start_cfg(filename=filename)
    mpiflag,iproc,numproc,mpi   = start.start_mpi()
    
    if iproc == 0: print 'importing geometry....'
    a=['ok from '+str(iproc)]
    
    mpi.barrier()
    total_a=mpi.allgather(a)
    if iproc == 0: print total_a
    
    def runimport(geometryfile,iproc,filename=None):
        import start as start
        cubit                   = start.start_cubit()
        cfg                         = start.start_cfg(filename=filename)
        file1=cfg.output_dir+'/'+geometryfile
        cubitcommand= 'open "'+ file1+ '"  ' 
        cubit.cmd(cubitcommand)                
        
    if cfg.parallel_import:
        runimport(geometryfile,iproc,filename=filename)
    else:
        if iproc == 0:
            runimport(geometryfile,iproc,filename=filename)
            for i in range(1,mpi.size):
                mpi.send('import',i)
                msg,status=mpi.recv(i)
        else:
            msg,status=mpi.recv(0)
            runimport(geometryfile,iproc,filename=filename)
            mpi.send('ok'+str(iproc),0)

        
def savesurf(iproc=0):
    savegeometry(iproc=iproc,surf=True)
    
###################################################################################### BELOW OK
    
def load_curves(acis_filename):
    """
    load the curves from acis files
    """
    import os
    #
    #
    print acis_filename
    if acis_filename and os.path.exists(acis_filename):
        tmp_curve=cubit.get_last_id("curve")
        command = "import acis '"+acis_filename+"'"
        cubit.cmd(command)
        tmp_curve_after=cubit.get_last_id("curve")
        curves=' '.join(str(x) for x in range(tmp_curve+1,tmp_curve_after+1))
    elif not os.path.exists(acis_filename):
        print str(acis_filename)+' not found'
        curves=None
    return [curves]

def project_curves(curves,top_surface):
    """
    project curves on surface
    """
    if not isinstance(curves,list): curves=curves.split()
    tmpc=[]
    for curve in curves:
        command = "project curve "+str(curve)+" onto surface "+str(top_surface)
        cubit.cmd(command)
        tmp_curve_after=cubit.get_last_id("curve")
        tmpc.append(tmp_curve_after)
        command = "del curve "+str(curve)
        cubit.cmd(command)
    return tmpc
    
def geo2utm(lon,lat,unit,ellipsoid=23):
    """conversion geocoodinates from geographical to utm
    
    usage: x,y=geo2utm(lon,lat,unit,ellipsoid=23)
    
    dafault ellipsoid is 23 = WGS-84, 
        ellipsoid:
        1, "Airy"
        2, "Australian National"
        3, "Bessel 1841"
        4, "Bessel 1841 (Nambia] "
        5, "Clarke 1866"
        6, "Clarke 1880"
        7, "Everest"
        8, "Fischer 1960 (Mercury] "
        9, "Fischer 1968"
        10, "GRS 1967"
        11, "GRS 1980"
        12, "Helmert 1906"
        13, "Hough"
        14, "International"
        15, "Krassovsky"
        16, "Modified Airy"
        17, "Modified Everest"
        18, "Modified Fischer 1960"
        19, "South American 1969"
        20, "WGS 60"
        21, "WGS 66"
        22, "WGS-72"
        23, "WGS-84"
        
    unit:  'geo' if the coordinates of the model (lon,lat) are geographical
           'utm' if the coordinates of the model (lon,lat) are utm
           
    x,y: the function return the easting, northing utm coordinates 
    """
    import LatLongUTMconversion
    if unit == 'geo' :          
       (zone, x, y) = LatLongUTMconversion.LLtoUTM(ellipsoid, lat, lon)
    elif unit == 'utm' : 
       x=lon
       y=lat
    return x,y


def savegeometry(iproc=0,surf=False,filename=None):
    import start as start
    cfg                         = start.start_cfg(filename=filename)
    mpiflag,iproc,numproc,mpi   = start.start_mpi()
    
    def runsave(geometryfile,iproc,filename=None):
        import start as start
        cubit                   = start.start_cubit()
        cfg                         = start.start_cfg(filename=filename)
        flag=0
        ner=cubit.get_error_count()
        cubitcommand= 'save as "'+ cfg.output_dir+'/'+geometryfile+ '"  overwrite' 
        cubit.cmd(cubitcommand)                                                    
        ner2=cubit.get_error_count()                                             
        if ner == ner2:
            flag=1
        return flag
        
    if surf:
        geometryfile='surf_vol_'+str(iproc)+'.cub'
    else:
        geometryfile='geometry_vol_'+str(iproc)+'.cub'
        
    flagsaved=0
    infosave=(iproc,flagsaved)
    
    mpi.barrier()
    total_saved=mpi.allgather(flagsaved)
    if isinstance(total_saved,int): total_saved=[total_saved]
    
    ind=0
    saving=True
    while saving:
        if len(total_saved) != sum(total_saved):
            #
            if not flagsaved: 
                flagsaved=runsave(geometryfile,iproc,filename=filename)
                if flagsaved:
                    infosave=(iproc,flagsaved)        
                    if numproc > 1:
                        f=open('geometry_saved'+str(iproc),'w')
                        f.close()
            mpi.barrier()
            total_saved=mpi.allgather(flagsaved)
            if isinstance(total_saved,int): total_saved=[total_saved]
            ind=ind+1
        else:
            saving=False
        if ind > len(total_saved)+10: saving=False
        print sum(total_saved),'/',len(total_saved),' saved'
    
    info_total_saved=mpi.allgather(infosave)
    if isinstance(info_total_saved,int): info_total_saved=[info_total_saved]
    
    if iproc==0:
        f=open('geometry_saving.log','w')
        f.write('\n'.join(str(x) for x in info_total_saved))                
        f.close()


def get_v_h_list(vol_id_list,chktop=False):
    """return the lists of the cubit ID of vertical/horizontal surface and vertical/horizontal curves
    where v/h is defined by the distance of the z normal component from the axis direction
    the parameter cfg.tres is the threshold as for example if 
    -tres <= normal[2] <= tres 
    then the surface is vertical
    #
    usage: surf_or,surf_vertical,list_curve_or,list_curve_vertical,bottom,top = get_v_h_list(list_vol,chktop=False)
    """
    #
    tres=0.3
    
    try:
        nvol=len(vol_id_list)
    except:
        nvol=1
        vol_id_list=[vol_id_list]
    surf_vertical=[]
    surf_or=[]
    list_curve_vertical=[]
    list_curve_or=[]
    #
    #
    for id_vol in vol_id_list:
        lsurf=cubit.get_relatives("volume",id_vol,"surface")
        for k in lsurf:
            normal=cubit.get_surface_normal(k)
            center_point = cubit.get_center_point("surface", k)
            if  -1*tres <= normal[2] <= tres:
               surf_vertical.append(k)
               lcurve=cubit.get_relatives("surface",k,"curve")
               list_curve_vertical=list_curve_vertical+list(lcurve)                                                                                                          
            else:
                surf_or.append(k)
                lcurve=cubit.get_relatives("surface",k,"curve")
                list_curve_or=list_curve_or+list(lcurve)
    for x in list_curve_or:
        try:
            list_curve_vertical.remove(x)
        except:
            pass
            
    #find the top and the bottom surfaces
    k=surf_or[0]
    center_point = cubit.get_center_point("surface", k)[2]
    center_point_top=center_point
    center_point_bottom=center_point
    top=k
    bottom=k 
    for k in surf_or[1:]:
        center_point = cubit.get_center_point("surface", k)[2]
        if center_point > center_point_top:
            center_point_top=center_point
            top=k
        elif center_point < center_point_bottom:
            center_point_bottom=center_point
            bottom=k
    #check that a top surface exists
    #it assume that the z coord of the center point 
    if chktop:
        k=lsurf[0]
        vertical_centerpoint_top = cubit.get_center_point("surface", k)[2]
        vertical_zmax_box_top=cubit.get_bounding_box('surface',k)[7]
        normal_top=cubit.get_surface_normal(k)    
        top=k        
        for k in lsurf:
            vertical_centerpoint = cubit.get_center_point("surface", k)[2]
            vertical_zmax_box=cubit.get_bounding_box('surface',k)[7]
            normal=cubit.get_surface_normal(k)
            check=(vertical_centerpoint >= vertical_centerpoint_top) and (vertical_zmax_box >= vertical_zmax_box_top) and (normal >= normal_top)
            if check:
                top=k
        if top in surf_vertical:
            surf_vertical.remove(top)
        if top not in surf_or: 
            surf_or.append(top)
    #if more than one surf is on the top, I get all the surfaces that are in touch with top surface but not the vertical surfaces
    surftop=list(cubit.get_adjacent_surfaces("surface", top)) #top is included in the list
    for s in surf_vertical:
        try:                          
            surftop.remove(s)    
        except:                       
            pass                      
    top=surftop
    #check that all the surf are Horizontal or vertical
    surf_all=surf_vertical+surf_or
    if len(surf_all)!=len(lsurf):
        print 'not all the surf are horizontal or vertical, check the normals'
        print 'list of surfaces: ',surf_all
        print 'list of vertical surface',surf_vertical
        print 'list of horizontal surface',surf_or        

    bottom=[bottom]
    return surf_or,surf_vertical,list_curve_or,list_curve_vertical,bottom,top

def list2str(l):
    if not isinstance(l,list): l=list(l)
    return ' '.join(str(x) for x in l)
    
def highlight(ent,l):
    txt=list2str(l)
    txt='highlight '+ent+' '+txt
    cubit.cmd(txt)
