#!/usr/bin/env python
#############################################################################
# cpml.py                                                    
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
#
try:
    import start as start
    cubit                   = start.start_cubit()
except:
    try:
        import cubit
    except:
        print 'error importing cubit, check if cubit is installed'
        pass


def create_single_pml(n,coord,operator,limit,normal,distance,layers):
    idv1=cubit.get_last_id('volume')+1
    if normal=='0 0 -1':
        txt="create element extrude face in surface with %s %s %s direction %s distance %s layers %s"
        cmd=txt % (coord,operator,limit,normal,distance,layers)    
    else:
        txt="create element extrude face in (surf in vol %s) with %s %s %s direction %s distance %s layers %s"
        cmd=txt % (n,coord,operator,limit,normal,distance,layers)
    cubit.cmd(cmd)
    cmd='create mesh geometry Hex all except hex in vol all feature_angle 135.0'
    cubit.cmd(cmd)
    cmd='del hex all except hex in vol all'
    cubit.cmd(cmd)
    idv2=cubit.get_last_id('volume')
    idv=' '.join(str(x) for x in range(idv1,idv2+1))
    return idv

def create_volume_pml(*args,**keys):
    ymin=keys.get("ymin",False)
    xmin=keys.get("xmin",False)
    ymax=keys.get("ymax",False)
    xmax=keys.get("xmax",False)
    zmin=keys.get("zmin",False)
    zmax=keys.get("zmax",False)
    nvol=keys.get("vol",False)
    layers=int(keys.get("layers",2))
    size=float(keys.get("size",0))
    thres=keys.get("thres",1)
    distance=str(size*layers)
    layers=str(layers)
    #
    #
    for n in nvol:
        normal='1 0 0'
        coord='X_coord'
        operator='>'
        limit=str(xmax-thres)
        idv=create_single_pml(n,coord,operator,limit,normal,distance,layers)
        ##
        ##
        normal='-1 0 0'
        coord='X_coord'
        operator='<'
        limit=str(xmin+thres)
        idv2=create_single_pml(n,coord,operator,limit,normal,distance,layers)
        #
        #
        #
        normal='0 -1 0'
        coord='Y_coord'
        operator='<'
        limit=str(ymin+thres)
        idv=create_single_pml(n,coord,operator,limit,normal,distance,layers)
        #
        #
        normal='0 1 0'
        coord='Y_coord'
        operator='>'
        limit=str(ymax-thres)
        idv=create_single_pml(n,coord,operator,limit,normal,distance,layers)
        #
        #
    normal='0 0 -1'
    coord='Z_coord'
    operator='<'
    limit=str(zmin+thres)
    nvol='all'
    idv=create_single_pml(nvol,coord,operator,limit,normal,distance,layers)


def mesh_cpml(list_vol,remesh=True,refinement=None,top_surf=None,size=None):
    from utilities import list2str
    top_surf=list2str(top_surf)
    if remesh:
        cubit.cmd('reset vol all')
        cubit.cmf('set dev on')
        cubit.cmd('imprint vol all')
        cubit.cmd('merge vol all')
        cubit.cmd('vol all size '+str(size))
        cubit.cmd('mesh vol all')
        try:
            for refdepth in refinement:
                cubit.cmd('refine surf '+top_surf+' numsplit 1 bias 1 depth '+str(refdepth)) 
        except:
            print 'DEBUG: error in refinement cpml'
        xmin=xmin-size
        xmax=xmax+size
        ymin=ymin-size
        ymax=ymax+size
        zmin=zmin-size
        zmin=zmax-size
        txt="group 'vol_xmin' add vol in surf  with X_coord < "+str(xmin)
        cubit.cmd(txt)
        txt="group 'vol_xmax' add vol in surf  with X_coord > "+str(xmax)
        cubit.cmd(txt)
        txt="group 'vol_ymin' add vol in surf  with Y_coord < "+str(ymin)
        cubit.cmd(txt)
        txt="group 'vol_ymax' add vol in surf  with Y_coord > "+str(ymax)
        cubit.cmd(txt)
        txt="group 'vol_zmin' add vol in surf  with Z_coord < "+str(zmin)
        cubit.cmd(txt)


    
def collecting_cpml(ip,size=None,cpuxmin=0,cpuxmax=1,cpuymin=0,cpuymax=1,cpux=1,cpuy=1,cubfiles=False,decimate=False,layers=2):
    import glob
    import re
    from utilities import list2str
    #
    if not size:
        print 'cpml size must be specified'
        return
        
        
    boundary_dict={}
    ##
    try:
        from boundary_definition import check_bc, map_boundary, define_surf
    except:
        pass
    #
    xmin,xmax,ymin,ymax,listfull=map_boundary(cpuxmin,cpuxmax,cpuymin,cpuymax,cpux,cpuy)
    #
    if cubfiles:        
        nf,listip,filenames,cubflag=importing_cubfiles(cubfiles)
    else:
        nf=0
        filenames=[]
        ip=0
    #
    if nf > 0:
        for ip,filename in zip(listip,filenames):
            try:
                if ip in listfull:
                    if cubflag:
                        cubit.cmd('import cubit "'+filename+'"')
                    else:
                        cubit.cmd('import mesh geometry "'+filename+'" block all use nodeset sideset feature_angle 135.00 linear merge')
                    if decimate: cubit.cmd('refine volume all numsplit 1 bias 1.0 depth 1 ')
            except:
                cubit.cmd('import mesh geometry "'+filename+'" block all use nodeset sideset feature_angle 135.00 linear merge')
                if decimate: cubit.cmd('refine volume all numsplit 1 bias 1.0 depth 1 ')
                ip=0
        if decimate: cubit.cmd('export mesh "decimated_before_cmpl.e" dimension 3 block all overwrite')
    else:
        if decimate: 
            cubit.cmd('refine volume all numsplit 1 bias 1.0 depth 1 ')
            cubit.cmd('export mesh "decimated_before_cmpl.e" dimension 3 block all overwrite')

    
    #
    #
    #print boundary_dict
    block_list=cubit.get_block_id_list()
    for block in block_list:
        ty=cubit.get_block_element_type(block)
        if ty == 'HEX8':
            cubit.cmd('block '+str(block)+' name "vol'+str(block)+'"')
            
    list_vol=list(cubit.parse_cubit_list('volume','all'))
    create_pml(xmin=xmin,xmax=xmax,ymax=ymax,ymin=ymin,zmin=zmin,zmax=zmax,size=size,layers=layers,vol=list_vol)

