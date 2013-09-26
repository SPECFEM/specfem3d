#############################################################################
# boundary_definition.py                                                    
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

from utilities import list2str


def map_boundary(cpuxmin=0,cpuxmax=1,cpuymin=0,cpuymax=1,cpux=1,cpuy=1):
    ymin=[]
    xmax=[]
    xmin=[]
    ymax=[]
    listfull=[]
    #
    for ix in range(cpuxmin,cpuxmax):
        for iy in range(cpuymin,cpuymax):
            ip=iy*cpux+ix
            if ix == cpuxmin:
                xmin.append(ip)
            if ix == cpuxmax-1:
                xmax.append(ip)
            if iy == cpuymin:
                ymin.append(ip)
            if iy == cpuymax-1:
                ymax.append(ip)
                #
            listfull.append(ip)
    return xmin,xmax,ymin,ymax,listfull


def define_4side_lateral_surfaces():
    list_vol=cubit.parse_cubit_list("volume","all")
    surf_xmin=[]
    surf_ymin=[]
    surf_xmax=[]
    surf_ymax=[]
    for id_vol in list_vol:
        surf_vertical=[]
        xsurf=[]
        ysurf=[]
        tres=0.3
        lsurf=cubit.get_relatives("volume",id_vol,"surface")
        for k in lsurf:
            normal=cubit.get_surface_normal(k)
            center_point = cubit.get_center_point("surface", k)
            if normal[2] >= -1*tres and normal[2] <= tres:
                surf_vertical.append(k)
                xsurf.append(center_point[0])
                ysurf.append(center_point[1])
        surf_xmin.append(surf_vertical[xsurf.index(min(xsurf))])
        surf_ymin.append(surf_vertical[ysurf.index(min(ysurf))])
        surf_xmax.append(surf_vertical[xsurf.index(max(xsurf))])
        surf_ymax.append(surf_vertical[ysurf.index(max(ysurf))])
    return surf_xmin,surf_ymin,surf_xmax,surf_ymax




def lateral_boundary_are_absorbing(ip=0,cpuxmin=0,cpuxmax=1,cpuymin=0,cpuymax=1,cpux=1,cpuy=1):
    #
    xmin,ymin,xmax,ymax=define_4side_lateral_surfaces()
    ip_xmin,ip_xmax,ip_ymin,ip_ymax,listfull=map_boundary(cpuxmin,cpuxmax,cpuymin,cpuymax,cpux,cpuy)
    if  not isinstance(ip_xmin, list): ip_xmin=[ip_xmin]
    if  not isinstance(ip_ymin, list): ip_ymin=[ip_ymin]
    if  not isinstance(ip_xmax, list): ip_xmax=[ip_xmax]
    if  not isinstance(ip_ymax, list): ip_ymax=[ip_ymax]
    #
    abs_xmin=[]
    abs_ymin=[]
    abs_xmax=[]
    abs_ymax=[]
    #
    if ip in ip_xmin: 
        abs_xmin=xmin
        print 'proc ',ip,' is has absorbing boundary xmin'
    if ip in ip_ymin:     
        print 'proc ',ip,' is has absorbing boundary ymin'
        abs_ymin=ymin
    if ip in ip_xmax:     
        print 'proc ',ip,' is has absorbing boundary xmax'
        abs_xmax=xmax
    if ip in ip_ymax:     
        print 'proc ',ip,' is has absorbing boundary ymax'
        abs_ymax=ymax
    return abs_xmin,abs_xmax,abs_ymin,abs_ymax


def define_surf(ip=0,cpuxmin=0,cpuxmax=1,cpuymin=0,cpuymax=1,cpux=1,cpuy=1):
    """
    define the absorbing surfaces for a layered topological box where boundary are surfaces parallel to the axis.
    it returns absorbing_surf,absorbing_surf_xmin,absorbing_surf_xmax,absorbing_surf_ymin,absorbing_surf_ymax,absorbing_surf_bottom,topo_surf
    where
    absorbing_surf is the list of all the absorbing boundary surf
    absorbing_surf_xmin is the list of the absorbing boundary surfaces that correnspond to x=xmin
    ...
    absorbing_surf_bottom is the list of the absorbing boundary surfaces that correspond to z=zmin
    """
    from utilities import get_v_h_list
    #
    from sets import Set
    def product(*args, **kwds):
        # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
        # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
        # for compatibility with python2.5
        pools = map(tuple, args) * kwds.get('repeat', 1)
        result = [[]]
        for pool in pools:
            result = [x+[y] for x in result for y in pool]
        return result
    absorbing_surf=[]
    xmin=[]
    xmax=[]
    ymin=[]
    ymax=[]
    #
    top_surf=[]
    bottom_surf=[]
    list_vol=cubit.parse_cubit_list("volume","all")
    zmax_box=cubit.get_total_bounding_box("volume",list_vol)[7]
    zmin_box=cubit.get_total_bounding_box("volume",list_vol)[6] #it is the z_min of the box ... box= xmin,xmax,d,ymin,ymax,d,zmin...    
    xmin_box=cubit.get_total_bounding_box("volume",list_vol)[0]
    xmax_box=cubit.get_total_bounding_box("volume",list_vol)[1]
    ymin_box=cubit.get_total_bounding_box("volume",list_vol)[3]
    ymax_box=cubit.get_total_bounding_box("volume",list_vol)[4]
    list_surf=cubit.parse_cubit_list("surface","all")
    
    absorbing_surface_distance_tolerance=0.001
    topographic_surface_distance_tolerance=0.1
    topographic_surface_normal_tolerance=0.4
    
    lv=[]
    for k in list_surf:
        sbox=cubit.get_bounding_box('surface',k)
        if zmax_box == 0 and sbox[7] == 0:
             dzmax=0
        elif zmax_box == 0 or sbox[7] == 0:
            dzmax=abs(sbox[7] - zmax_box)
        else:
            dzmax=abs(sbox[7] - zmax_box)/max(abs(sbox[7]),abs(zmax_box))
        if zmin_box == 0 and sbox[6] == 0:
             dzmin=0
        elif zmin_box == 0 or sbox[6] == 0:
            dzmin=abs(sbox[6] - zmin_box)                                            
        else:                                                    
            dzmin=abs(sbox[6] - zmin_box)/max(abs(sbox[6]),abs(zmin_box))
        normal=cubit.get_surface_normal(k)
        zn=normal[2]
        if dzmax <= topographic_surface_distance_tolerance and zn > topographic_surface_normal_tolerance:
            top_surf.append(k)
            list_vertex=cubit.get_relatives('surface',k,'vertex')
            for v in list_vertex:
                valence=cubit.get_valence(v)
                if valence <= 4: #valence 3 is a corner, 4 is a vertex between 2 volumes, > 4 is a vertex not in the boundaries
                    lv.append(v)
        elif dzmin <= 0.001 and zn < -1+topographic_surface_normal_tolerance:
            bottom_surf.append(k)
    if len(top_surf) ==0: #assuming that one topo surface need to be selected
            _,_,_,_,_,top_surf=get_v_h_list(list_vol,chktop=False)
    lp=[]
    labelp=[]
    combs=product(lv,lv)
    for comb in combs:
        v1=comb[0]
        v2=comb[1]
        c=Set(cubit.get_relatives("vertex",v1,"curve")) & Set(cubit.get_relatives("vertex",v2,"curve"))
        if len(c) == 1:
            p=cubit.get_center_point("curve",list(c)[0])
            labelp.append(list(c)[0])
    labelps=Set(labelp)
    for c in labelps:
        p=cubit.get_center_point("curve",c)
        lp.append(p)
    
    for k in list_surf: 
        center_point = cubit.get_center_point("surface", k)
        for p in lp:
            try:
                if abs((center_point[0] - p[0])/p[0]) <= absorbing_surface_distance_tolerance and abs((center_point[1] - p[1])/p[1]) <= absorbing_surface_distance_tolerance:
                    absorbing_surf.append(k)
                    break
            except:
                if -1 <= center_point[0] <= 1 and -1 <= center_point[1] <= 1:
                    absorbing_surf.append(k)
                    break
    #
    four_side=True
    if four_side:
        xmintmp,ymintmp,xmaxtmp,ymaxtmp=define_4side_lateral_surfaces()
        xmin=list(Set(xmintmp)-Set(xmaxtmp))
        xmax=list(Set(xmaxtmp)-Set(xmintmp))
        ymin=list(Set(ymintmp)-Set(ymaxtmp))
        ymax=list(Set(ymaxtmp)-Set(ymintmp))
        abs_xmintmp,abs_xmaxtmp,abs_ymintmp,abs_ymaxtmp=lateral_boundary_are_absorbing(ip,cpuxmin,cpuxmax,cpuymin,cpuymax,cpux,cpuy)
        abs_xmin=list(Set(abs_xmintmp)-Set(abs_xmaxtmp))
        abs_xmax=list(Set(abs_xmaxtmp)-Set(abs_xmintmp))
        abs_ymin=list(Set(abs_ymintmp)-Set(abs_ymaxtmp))
        abs_ymax=list(Set(abs_ymaxtmp)-Set(abs_ymintmp))
    return absorbing_surf,abs_xmin,abs_xmax,abs_ymin,abs_ymax,top_surf,bottom_surf,xmin,ymin,xmax,ymax



def define_block():
    list_vol=cubit.parse_cubit_list("volume","all")
    init_n_vol=len(list_vol)
    list_name=map(lambda x: 'vol'+x,map(str,list_vol))
    return list_vol,list_name

def build_block(vol_list,name,id_0=1,top_surf=None,optionsea=False):
    #
    from sets import Set
    if optionsea:
        sea=optionsea['sea']
        seaup=optionsea['seaup']
        sealevel=optionsea['sealevel']
        seathres=optionsea['seathres']
    else:
        sea=False
        seaup=False
        sealevel=False
        seathres=False
    
    #
    block_list=cubit.get_block_id_list()
    if len(block_list) > 0:
        id_block=max(max(block_list),2)+id_0
    else:
        id_block=2+id_0
    for v,n in zip(vol_list,name):
        id_block+=1
        v_other=Set(vol_list)-Set([v])
        #command= 'block '+str(id_block)+' hex in node in vol '+str(v)+' except hex in vol '+str(list(v_other))
        if sea and v == vol_list[-1]:
            cubit.cmd('set duplicate block elements off')
            tsurf_string=" ".join(str(x) for x in top_surf)
            #sea
            command= 'block '+str(id_block)+' hex in node in surf '+tsurf_string+' with Z_coord < '+str(seathres)
            cubit.cmd(command)
            command = "block "+str(id_block)+" name 'sea"+n+"'"
            cubit.cmd(command)
            if not seaup:
                id_block+=1
                command= 'block '+str(id_block)+' hex in node in surf '+tsurf_string+' with (Z_coord > '+str(seathres)+' and Z_coord < '+str(sealevel)
                cubit.cmd(command)
                command = "block "+str(id_block)+" name 'shwater"+n+"'"
                cubit.cmd(command)
            id_block+=1
            command= 'block '+str(id_block)+' hex in node in surf '+tsurf_string+' with Z_coord >= '+str(sealevel)
            cubit.cmd(command)
            command = "block "+str(id_block)+" name 'continent"+n+"'"
            cubit.cmd(command)
        else:
            command= 'block '+str(id_block)+' hex in vol '+str(v)+' except hex in vol '+str(list(v_other))
            print command
            command = command.replace("["," ").replace("]"," ")
            cubit.cmd(command) 
            command = "block "+str(id_block)+" name '"+n+"'"
            cubit.cmd(command)



def build_block_side(surf_list,name,obj='surface',id_0=1):
    id_nodeset=cubit.get_next_nodeset_id()
    id_block=id_0
    if obj == 'hex':
        txt='hex in node in surface'
        txt1='block '+str(id_block)+ ' '+ txt +' '+str(list(surf_list))
        txt2="block "+str(id_block)+" name '"+name+"'"
        txt1=txt1.replace("["," ").replace("]"," ")
        cubit.cmd(txt1)
        cubit.cmd(txt2)
    elif obj == 'node':
        txt=obj+' in surface'
        txt1= 'nodeset '+str(id_nodeset)+ ' '+ txt +' '+str(list(surf_list))
        txt1 = txt1.replace("["," ").replace("]"," ")
        txt2 = "nodeset "+str(id_nodeset)+" name '"+name+"'"
        cubit.cmd(txt1)
        cubit.cmd(txt2)
    elif obj == 'face' or obj == 'edge':
        txt=obj+' in surface'
        txt1= 'block '+str(id_block)+ ' '+ txt +' '+str(list(surf_list))
        txt1 = txt1.replace("["," ").replace("]"," ")
        txt2 = "block "+str(id_block)+" name '"+name+"'"
        cubit.cmd(txt1)
        cubit.cmd(txt2)
    else:
        txt1=''
        txt2="block "+str(id_block)+" name "+name+"_notsupported (only hex,face,edge,node)"
        try:
            cubit.cmd('comment "'+txt1+'"')
            cubit.cmd('comment "'+txt2+'"')
        except:
            pass

def define_bc(*args,**keys):
    id_0=1
    #
    #
    parallel=keys.get('parallel',True)
    closed=keys.get('closed',False)
    ip=keys.get("iproc",0)
    cpuxmin=keys.get("cpuxmin",0)
    cpuymin=keys.get("cpuymin",0)
    cpux=keys.get("cpux",1)
    cpuy=keys.get("cpuy",1)
    cpuxmax=keys.get("cpuxmax",cpux)
    cpuymax=keys.get("cpuymax",cpuy)
    optionsea=keys.get("optionsea",False)
    #
    if parallel:
        absorbing_surf,abs_xmin,abs_xmax,abs_ymin,abs_ymax,top_surf,bottom_surf,xmin,ymin,xmax,ymax=define_surf(ip=ip,cpuxmin=cpuxmin,cpuxmax=cpuxmax,cpuymin=cpuymin,cpuymax=cpuymax,cpux=cpux,cpuy=cpuy)
        id_0=cubit.get_next_block_id()
        v_list,name_list=define_block()
        build_block(v_list,name_list,id_0,top_surf,optionsea=optionsea)
        #
    elif closed:
        surf=define_absorbing_surf_sphere()
        v_list,name_list=define_block()
        build_block(v_list,name_list,id_0)
        entities=args[0]
        id_side=1
        for entity in entities:
            build_block_side(surf,entity+'_closedvol',obj=entity,id_0=id_side)
            id_side=id_side+1
            
#########################################


def extract_bottom_curves(surf_xmin,surf_ymin,surf_xmax,surf_ymax):
    curve_xmin=[]
    curve_ymin=[]
    curve_xmax=[]
    curve_ymax=[]
    from sets import Set #UPGRADE.... the sets module is deprecated after python 2.6
    for s in surf_xmin:
        lcs=cubit.get_relatives("surface",s,"curve")
        for lc in lcs:
            curve_xmin.append(lc)
    for s in surf_xmax:
        lcs=cubit.get_relatives("surface",s,"curve")
        for lc in lcs:
            curve_xmax.append(lc)
    for s in surf_ymin:
        lcs=cubit.get_relatives("surface",s,"curve")
        for lc in lcs:
            curve_ymin.append(lc)
    for s in surf_ymax:
        lcs=cubit.get_relatives("surface",s,"curve")
        for lc in lcs:
            curve_ymax.append(lc)
    curve_xmin=list(Set(curve_xmin))
    curve_ymin=list(Set(curve_ymin))
    curve_xmax=list(Set(curve_xmax))
    curve_ymax=list(Set(curve_ymax))
    curve_bottom_xmin=select_bottom_curve(curve_xmin)
    curve_bottom_ymin=select_bottom_curve(curve_ymin)
    curve_bottom_xmax=select_bottom_curve(curve_xmax)
    curve_bottom_ymax=select_bottom_curve(curve_ymax)
    #
    return curve_bottom_xmin,curve_bottom_ymin,curve_bottom_xmax,curve_bottom_ymax

def select_bottom_curve(lc):
    z=[]
    for l in lc:
        center_point = cubit.get_center_point("curve", l)
        z.append(center_point[2])
    result=zip(z,lc)
    result.sort()
    return result[0][1]


def get_ordered_node_surf(lsurface,icurve):
    if not isinstance(lsurface,str): 
        lsurf=list2str(lsurface)
    #
    if not isinstance(icurve,str): 
        icurvestr=str(icurve)
    orient_nodes_surf=[]
    #
    #get the nodes on a surface, I don't use the method get_surface_nodes since it has different behavior in cubit12.2 and cubit13.2+
    k=cubit.get_id_from_name('sl')
    if k!=0:
        cubit.cmd('del group sl')
    else:
        print 'initializing group sl'
    cubit.cmd("group 'sl' add node in surf "+lsurf)
    group1 = cubit.get_id_from_name("sl")
    nodes_ls =list(cubit.get_group_nodes(group1))
    nnode=len(nodes_ls)
    #
    #get the nodes on curves
    orient=[]
    k=cubit.get_id_from_name('n1')
    if k!=0:
        cubit.cmd('del group n1')
    else:
        print 'initializing group n1'
    cubit.cmd("group 'n1' add node in curve "+icurvestr)
    x=cubit.get_bounding_box('curve', icurve)
    if x[2]>x[5]:
        idx=0
    else:
        idx=1
    group1 = cubit.get_id_from_name("n1")
    nodes1 = list(cubit.get_group_nodes(group1))
    for n in nodes1:
        v = cubit.get_nodal_coordinates(n)
        orient.append(v[idx])
    result=zip(orient,nodes1)
    result.sort()
    nodes2=[c[1] for c in result]
    for n in nodes2:
        try:
            nodes_ls.remove(n)
        except:
            pass             
    orient_nodes_surf=orient_nodes_surf+nodes2
    #
    while len(orient_nodes_surf) < nnode:
        cubit.cmd('del group n1')
        cubit.cmd("group 'n1' add node in edge in node "+str(nodes2).replace('[',' ').replace(']',' '))
        group1 = cubit.get_id_from_name("n1")
        nodes1 = list(cubit.get_group_nodes(group1))
        orient=[]
        nd=[]
        for n in nodes1:
            if n in nodes_ls:
                v = cubit.get_nodal_coordinates(n)
                orient.append(v[idx])
                nd.append(n)
        result=zip(orient,nd)
        result.sort()
        nodes2=[c[1] for c in result]
        for n in nodes2:
            try:
                nodes_ls.remove(n)
            except:
                pass
        orient_nodes_surf=orient_nodes_surf+nodes2
    #get the vertical curve
    curve_vertical=[]
    for s in lsurface:
        lcs=cubit.get_relatives("surface",s,"curve")
        for l in lcs:
            x=cubit.get_bounding_box('curve', l)
            length=[(x[2],1),(x[5],2),(x[8],3)]
            length.sort()
            if length[-1][1] == 3:
                curve_vertical.append(l)
    #
    kcurve=list2str(curve_vertical)
    k=cubit.get_id_from_name('curve_vertical')
    if k!=0:
        cubit.cmd('del group curve_vertical')
    else:
        print 'initializing group curve_vertical'
    cubit.cmd("group 'curve_vertical' add node in curve "+kcurve)
    group1 = cubit.get_id_from_name('curve_vertical')
    nodes_curve = list(cubit.get_group_nodes(group1))
    for n in nodes_curve:
        try:
            orient_nodes_surf.remove(n)
        except:
            pass
    #
    return nodes_curve,orient_nodes_surf



def check_bc(iproc,xmin,xmax,ymin,ymax,cpux,cpuy,cpuxmin,cpuxmax,cpuymin,cpuymax):
    """
    boundary=check_bc(iproc,xmin,xmax,ymin,ymax,cpux,cpuy,cpuxmin,cpuxmax,cpuymin,cpuymax)
    #
    set the boundary condition during the collecting phase and group the nodes of the vertical surface in groups for the merging phase
    iproc is the value of the processor 
    xmin,ymin,ymax,ymin are the list of iproc that have at least one absorbing boundary condition
    """
    #
    absorbing_surf,abs_xmin,abs_xmax,abs_ymin,abs_ymax,top_surf,bottom_surf,surf_xmin,surf_ymin,surf_xmax,surf_ymax=define_surf(ip=iproc,cpuxmin=cpuxmin,cpuxmax=cpuxmax,cpuymin=cpuymin,cpuymax=cpuymax,cpux=cpux,cpuy=cpuy)
    curve_bottom_xmin,curve_bottom_ymin,curve_bottom_xmax,curve_bottom_ymax=extract_bottom_curves(surf_xmin,surf_ymin,surf_xmax,surf_ymax)
    print absorbing_surf,abs_xmin,abs_xmax,abs_ymin,abs_ymax,top_surf,bottom_surf,surf_xmin,surf_ymin,surf_xmax,surf_ymax
    #
    #
    #
    from sets import Set #UPGRADE.... the sets module is deprecated after python 2.6
    
    cubit.cmd('set info off')
    cubit.cmd('set echo off')
    nodes_curve_ymax,orient_nodes_surf_ymax=get_ordered_node_surf(surf_ymax,curve_bottom_ymax)
    nodes_curve_xmax,orient_nodes_surf_xmax=get_ordered_node_surf(surf_xmax,curve_bottom_xmax)
    nodes_curve_ymin,orient_nodes_surf_ymin=get_ordered_node_surf(surf_ymin,curve_bottom_ymin)
    nodes_curve_xmin,orient_nodes_surf_xmin=get_ordered_node_surf(surf_xmin,curve_bottom_xmin)
    c_xminymin=Set(nodes_curve_xmin).intersection(nodes_curve_ymin)
    c_xminymax=Set(nodes_curve_xmin).intersection(nodes_curve_ymax)
    c_xmaxymin=Set(nodes_curve_xmax).intersection(nodes_curve_ymin)
    c_xmaxymax=Set(nodes_curve_xmax).intersection(nodes_curve_ymax)
    cubit.cmd('set info on')
    cubit.cmd('set echo on')
    
    
    #
    orient=[]
    nd=[]
    for n in c_xminymin:
        v = cubit.get_nodal_coordinates(n)
        orient.append(v[2])
        nd.append(n)
    result=zip(orient,nd)
    result.sort()
    c_xminymin=[c[1] for c in result]
    #
    orient=[]
    nd=[]    
    for n in c_xminymax:
        v = cubit.get_nodal_coordinates(n)
        orient.append(v[2])
        nd.append(n)
    result=zip(orient,nd)
    result.sort()
    c_xminymax=[c[1] for c in result]
    #
    orient=[]
    nd=[]
    for n in c_xmaxymin:
        v = cubit.get_nodal_coordinates(n)
        orient.append(v[2])
        nd.append(n)
    result=zip(orient,nd)
    result.sort()
    c_xmaxymin=[c[1] for c in result]
    #
    orient=[]
    nd=[]
    for n in c_xmaxymax:
        v = cubit.get_nodal_coordinates(n)
        orient.append(v[2])
        nd.append(n)
    result=zip(orient,nd)
    result.sort()
    c_xmaxymax=[c[1] for c in result]
    #
    boundary={}
    boundary['id']=iproc
    #
    boundary['nodes_surf_xmin']=orient_nodes_surf_xmin
    boundary['nodes_surf_xmax']=orient_nodes_surf_xmax
    boundary['nodes_surf_ymin']=orient_nodes_surf_ymin
    boundary['nodes_surf_ymax']=orient_nodes_surf_ymax
    #
    boundary['node_curve_xminymin']=c_xminymin
    boundary['node_curve_xminymax']=c_xminymax
    boundary['node_curve_xmaxymin']=c_xmaxymin
    boundary['node_curve_xmaxymax']=c_xmaxymax
    
    #print boundary['node_curve_xminymin']
    #print     boundary['node_curve_xminymax']
    #print     boundary['node_curve_xmaxymin']
    #print     boundary['node_curve_xmaxymax']
    entities=['face']
    #
    for entity in entities:
        if len(abs_xmin) != 0:
            refname=entity+'_abs_xmin'
            build_block_side(abs_xmin,refname,obj=entity,id_0=1003)
        #
        if len(abs_ymin) != 0:
            refname=entity+'_abs_ymin'
            build_block_side(abs_ymin,refname,obj=entity,id_0=1004)
        #
        if len(abs_xmax) != 0:
            refname=entity+'_abs_xmax'
            build_block_side(abs_xmax,refname,obj=entity,id_0=1005)
        #
        if len(abs_ymax) != 0:
            refname=entity+'_abs_ymax'
            build_block_side(abs_ymax,refname,obj=entity,id_0=1006)
        ##
        refname=entity+'_topo'
        block=3 #change here..... must be 1 @@@@@@@@@
        ty=None
        ty=cubit.get_block_element_type(block)
        if ty=='':
            pass
        elif ty == 'HEX8':
            pass
        else:
            cubit.cmd('del block '+str(block))
        build_block_side(top_surf,refname,obj=entity,id_0=1001)
        #
        refname=entity+'_bottom'
        block=4 #change here..... must be 2 @@@@@@@@
        ty=None
        ty=cubit.get_block_element_type(block)
        if ty=='':
            pass
        elif ty == 'HEX8':
            pass
        else:
            cubit.cmd('del block '+str(block))
        build_block_side(bottom_surf,refname,obj=entity,id_0=1002)
    #
    #
    return boundary