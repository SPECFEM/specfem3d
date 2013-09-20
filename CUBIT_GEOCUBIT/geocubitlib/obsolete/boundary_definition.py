#############################################################################
# boundary_definition.py                                                    #
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


def define_absorbing_surf():
    """
    define the absorbing surfaces for a layered topological box where boundary are surfaces parallel to the axis.
    it returns absorbing_surf,absorbing_surf_xmin,absorbing_surf_xmax,absorbing_surf_ymin,absorbing_surf_ymax,absorbing_surf_bottom,topo_surf
    where
    absorbing_surf is the list of all the absorbing boundary surf
    absorbing_surf_xmin is the list of the absorbing boundary surfaces that correnspond to x=xmin
    ...
    absorbing_surf_bottom is the list of the absorbing boundary surfaces that correspond to z=zmin
    """
    import initializing as initializing
    cubit                   = initializing.initializing_cubit()
    mpiflag,iproc,numproc,mpi   = initializing.initializing_mpi()
    #     
    absorbing_surf=[]
    absorbing_surf_xmin=[]
    absorbing_surf_xmax=[]
    absorbing_surf_ymin=[]
    absorbing_surf_ymax=[]
    absorbing_surf_bottom=[]
    top_surf=[]
    list_vol=cubit.parse_cubit_list("volume","all")
    init_n_vol=len(list_vol)
    cfg                         = initializing.initializing_cfg()
    
    zmax_box=cubit.get_total_bounding_box("volume",list_vol)[7]
    zmin_box=cubit.get_total_bounding_box("volume",list_vol)[6] #it is the z_min of the box ... box= xmin,xmax,d,ymin,ymax,d,zmin...    
    #        [                                                  ]
    #xmin_box=[cubit.get_total_bounding_box("volume",list_vol)[0]]
    #xmax_box=[cubit.get_total_bounding_box("volume",list_vol)[1]]
    #ymin_box=[cubit.get_total_bounding_box("volume",list_vol)[3]]
    #ymax_box=[cubit.get_total_bounding_box("volume",list_vol)[4]]
    
    #mpi.barrier()
    #
    #total_xmin=mpi.allgather(xmin_box)
    #xmin_box=min(total_xmin)
    #
    #total_xmax=mpi.allgather(xmax_box)
    #xmax_box=max(total_xmax)       
    #
    #total_ymin=mpi.allgather(ymin_box)
    #ymin_box=min(total_ymin)       
    #
    #total_ymax=mpi.allgather(ymax_box)
    #ymax_box=max(total_ymax)     
    #
    #mpi.barrier()                        
    #
    def distance(x0,y0,x1,y1,x2,y2):
        import math
        d=abs((x2-x1)*(y1-y0)-(x1-x0)*(y2-y1))/math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))
        return d
    
    x1=cfg.x1_box
    x2=cfg.x2_box
    x3=cfg.x3_box
    x4=cfg.x4_box
    y1=cfg.y1_box
    y2=cfg.y2_box
    y3=cfg.y3_box
    y4=cfg.y4_box
    #tres=cfg.tres_boundarydetection
    tres=cfg.tres
    list_surf=cubit.parse_cubit_list("surface","all")
    surface_hor=[]
    for k in list_surf:
        center_point = cubit.get_center_point("surface", k)
        x0=center_point[0]
        y0=center_point[1]
        normal=cubit.get_surface_normal(k)[2]
        vertical= normal >= -1*tres and normal <= tres
        if distance(x0,y0,x1,y1,x4,y4) < abs(x4-x1)*.2:
            absorbing_surf_xmin.append(k)
            absorbing_surf.append(k)
        elif distance(x0,y0,x2,y2,x3,y3) < abs(x3-x2)*.2:
            absorbing_surf_xmax.append(k)
            absorbing_surf.append(k)
        elif distance(x0,y0,x1,y1,x2,y2) < abs(y2-y1)*.2:
            absorbing_surf_ymin.append(k)
            absorbing_surf.append(k)
        elif distance(x0,y0,x3,y3,x4,y4) < abs(y3-y4)*.2:
            absorbing_surf_ymax.append(k)
            absorbing_surf.append(k)
        elif not vertical: 
            surface_hor.append(k)
    ztop=zmin_box
    zbottom=zmax_box
    for k in surface_hor:
        sbox=cubit.get_bounding_box('surface',k)
        zsurf_min=sbox[6]
        zsurf_max=sbox[7]
        if zsurf_max >= ztop:
            ktop=k
            ztop=zsurf_max
        if zsurf_min <= zbottom:
            kbottom=k
            zbottom=zsurf_min
    absorbing_surf_bottom.append(kbottom)
    absorbing_surf.append(kbottom)           
    top_surf.append(ktop)
                
    f=open('bc_ok_'+str(iproc),'w')
    txt=str(absorbing_surf)+str(absorbing_surf_xmin)+str(absorbing_surf_xmax)+str(absorbing_surf_ymin)+str(absorbing_surf_ymax)+str(absorbing_surf_bottom)+str(top_surf)
    f.write(txt)
    f.close()




    return absorbing_surf,absorbing_surf_xmin,absorbing_surf_xmax,absorbing_surf_ymin,absorbing_surf_ymax,absorbing_surf_bottom,top_surf

def define_absorbing_surf_nopar():
    """
    define the absorbing surfaces for a layered topological box where boundary surfaces are not parallel to the axis.
    it returns absorbing_surf,topo_surf
    where
    absorbing_surf is the list of all the absorbing boundary surf
    """
    #import initializing as initializing
    #cubit                   = initializing.initializing_cubit()
    #
    from sets import Set
    def product(*args, **kwds):
        # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
        # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
        pools = map(tuple, args) * kwds.get('repeat', 1)
        result = [[]]
        for pool in pools:
            result = [x+[y] for x in result for y in pool]
        return result
    absorbing_surf=[]
    absorbing_surf_xmin=[]
    absorbing_surf_xmax=[]
    absorbing_surf_ymin=[]
    absorbing_surf_ymax=[]
    absorbing_surf_bottom=[]
    top_surf=[]
    bottom_surf=[]
    list_vol=cubit.parse_cubit_list("volume","all")
    init_n_vol=len(list_vol)
    zmax_box=cubit.get_total_bounding_box("volume",list_vol)[7]
    zmin_box=cubit.get_total_bounding_box("volume",list_vol)[6] #it is the z_min of the box ... box= xmin,xmax,d,ymin,ymax,d,zmin...    
    xmin_box=cubit.get_total_bounding_box("volume",list_vol)[0]
    xmax_box=cubit.get_total_bounding_box("volume",list_vol)[1]
    ymin_box=cubit.get_total_bounding_box("volume",list_vol)[3]
    ymax_box=cubit.get_total_bounding_box("volume",list_vol)[4]
    list_surf=cubit.parse_cubit_list("surface","all")
    lv=[]
    for k in list_surf:
            sbox=cubit.get_bounding_box('surface',k)
            dzmax=abs((sbox[7] - zmax_box)/zmax_box)
            dzmin=abs((sbox[6] - zmin_box)/zmin_box)
            normal=cubit.get_surface_normal(k)
            zn=normal[2]
            if dzmax <= 0.1 and zn > 0.4:
                top_surf.append(k)
                list_vertex=cubit.get_relatives('surface',k,'vertex')
                for v in list_vertex:
                    valence=cubit.get_valence(v)
                    if valence <= 4: #valence 3 is a corner, 4 is a vertex between 2 volumes, > 4 is a vertex not in the boundaries
                        lv.append(v)
            elif dzmin <= 0.001 and zn < -0.7:
                bottom_surf.append(k)
                absorbing_surf.append(k)
    lp=[]
    combs=product(lv,lv)
    for comb in combs:
        v1=comb[0]
        v2=comb[1]
        c=Set(cubit.get_relatives("vertex",v1,"curve")) & Set(cubit.get_relatives("vertex",v2,"curve"))
        if len(c) == 1:
            p=cubit.get_center_point("curve",list(c)[0])
            lp.append(p)
    for k in list_surf: 
        center_point = cubit.get_center_point("surface", k)
        for p in lp:
            if abs((center_point[0] - p[0])/p[0]) <= 0.001 and abs((center_point[1] - p[1])/p[1]) <= 0.001:
             absorbing_surf.append(k)
             break
    return absorbing_surf,top_surf,bottom_surf

def define_absorbing_surf_sphere():
    import initializing as initializing
    cubit                   = initializing.initializing_cubit()
    #
    surf=[]
    list_surf=cubit.parse_cubit_list("surface","all")
    for s in list_surf:
       v=cubit.get_relatives('surface',s,'volume')
       if len(v) == 1:
           surf.append(s)
    return surf

def define_block():
    #try:
    #    import initializing as initializing
    #    cubit                   = initializing.initializing_cubit()
    #except:
    #    pass
    #
    list_vol=cubit.parse_cubit_list("volume","all")
    init_n_vol=len(list_vol)
    list_name=map(lambda x: 'vol'+x,map(str,list_vol))
    return list_vol,list_name

def build_block(vol_list,name,id_0):
    #try:
    #    import initializing as initializing
    #    cubit                   = initializing.initializing_cubit()
    #except:
    #    pass
    #
    from sets import Set
    #
    block_list=cubit.get_block_id_list()
    if len(block_list) > 0:
         id_block=max(max(block_list),2)+id_0
    else:
        id_block=1+id_0
    for v,n in zip(vol_list,name):
       id_block+=1
       v_other=Set(vol_list)-Set([v])
       #command= 'block '+str(id_block)+' hex in node in vol '+str(v)+' except hex in vol '+str(list(v_other))
       command= 'block '+str(id_block)+' hex in vol '+str(v)+' except hex in vol '+str(list(v_other))
       print command
       command = command.replace("["," ").replace("]"," ")
       cubit.cmd(command) 
       command = "block "+str(id_block)+" name '"+n+"'"
       cubit.cmd(command)

def build_block_side(surf_list,name,obj='surface',id_0=1):
    #try:
    #    import initializing as initializing
    #    cubit                   = initializing.initializing_cubit()
    #except:
    #    pass
    #
    id_nodeset=cubit.get_next_nodeset_id()
    #id_block=cubit.get_next_block_id()
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
    import initializing as initializing
    mpiflag,iproc,numproc,mpi   = initializing.initializing_mpi()
    cubit                   = initializing.initializing_cubit()
    
    parallel=keys.get('parallel',True)
    closed=keys.get('closed',False)
    #partitioning_type=keys.get('partitioning_type',None)
    id_0=1
    #if partitioning_type:
    #    if partitioning_type == 't':
    #        id_0=10
    #    elif partitioning_type == 'in':
    #        id_0=100
    #    elif partitioning_type == 'out':
    #        id_0=1000
    #    else:
    #        print 'check the partitioning definition..', partitioning_type
    #        return
    ##
    cubit.cmd('comment "def bc"')
    if not closed:
        if parallel:
            surf,xmin,xmax,ymin,ymax,bottom,topo=define_absorbing_surf()
            cubit.cmd('comment "'+str([surf,xmin,xmax,ymin,ymax,bottom,topo])+'"')
        else:
            surf,topo,bottom=define_absorbing_surf_nopar()
        v_list,name_list=define_block()
        build_block(v_list,name_list,id_0)
        entities=args[0]
        id_side=cubit.get_next_block_id()
        for entity in entities:
            build_block_side(topo,entity+'_topo',obj=entity,id_0=1) #topo ha block 1
            id_side=cubit.get_next_block_id()
            build_block_side(bottom,entity+'_abs_bottom',obj=entity,id_0=id_side)
            id_side=id_side+1
            build_block_side(surf,entity+'_abs',obj=entity,id_0=id_side)
            id_side=id_side+1
            if parallel: 
                build_block_side(xmin,entity+'_abs_xmin',obj=entity,id_0=id_side)
                id_side=id_side+1
                build_block_side(xmax,entity+'_abs_xmax',obj=entity,id_0=id_side)
                id_side=id_side+1
                build_block_side(ymin,entity+'_abs_ymin',obj=entity,id_0=id_side)
                id_side=id_side+1
                build_block_side(ymax,entity+'_abs_ymax',obj=entity,id_0=id_side)
                id_side=id_side+1
    else:
        surf=define_absorbing_surf_sphere()
        v_list,name_list=define_block()
        build_block(v_list,name_list,id_0)
        entities=args[0]
        id_side=1
        for entity in entities:
            build_block_side(surf,entity+'_closedvol',obj=entity,id_0=id_side)
            id_side=id_side+1


def list2str(l):
    if not isinstance(l,list): l=list(l)
    return ' '.join(str(x) for x in l)



def get_ordered_node_surf(lsurface,icurve):
    if not isinstance(lsurface,str): 
        lsurf=list2str(lsurface)
    #
    if not isinstance(icurve,str): 
        icurvestr=str(icurve)
    orient_nodes_surf=[]
    #
    cubit.cmd('del group sl')
    cubit.cmd("group 'sl' add node in surf "+lsurf)
    group1 = cubit.get_id_from_name("sl")
    nodes_ls =list(cubit.get_group_nodes(group1))
    cubit.cmd('del group sl')
    nnode=len(nodes_ls)
    #
    orient=[]
    cubit.cmd('del group n1')
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
    icurve=list2str(curve_vertical)
    orientx=[]
    orienty=[]
    orientz=[]
    cubit.cmd('del group curve_vertical')
    cubit.cmd("group 'curve_vertical' add node in curve "+icurve)
    group1 = cubit.get_id_from_name('curve_vertical')
    nodes_curve = list(cubit.get_group_nodes(group1))
    for n in nodes_curve:
        try:
             orient_nodes_surf.remove(n)
        except:
             pass
    #
    return nodes_curve,orient_nodes_surf



def select_bottom_curve(lc):
    z=[]
    for l in lc:
        center_point = cubit.get_center_point("curve", l)
        z.append(center_point[2])
    result=zip(z,lc)
    result.sort()
    print result
    return result[0][1]
    


def check_bc(iproc,xmin,xmax,ymin,ymax,cpux,cpuy):
    """
    set the boundary condition during the collectiong phase and group the nodes of the vertical surface in groups for the merging phase
    iproc is the value of the processor 
    xmin,ymin,ymax,ymin are the list of iproc that have at least one absorbing boundary condition
    """
    try:
        import initializing as initializing
        cubit                   = initializing.initializing_cubit()
    except:
        pass
    list_vol=cubit.parse_cubit_list("volume","all")
    surf_xmin=[]
    surf_ymin=[]
    surf_xmax=[]
    surf_ymax=[]
    if  not isinstance(xmin, list): xmin=[xmin]
    if  not isinstance(ymin, list): ymin=[ymin]
    if  not isinstance(xmax, list): xmax=[xmax]
    if  not isinstance(ymax, list): ymax=[ymax]
    for id_vol in list_vol:
        surf_vertical=[]
        xsurf=[]
        ysurf=[]
        tres=0.3
        lsurf=cubit.get_relatives("volume",id_vol,"surface")
        icheck=0
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
    print curve_bottom_xmin,curve_bottom_ymin,curve_bottom_xmax,curve_bottom_ymax
    #
    #
    nodes_curve_ymax,orient_nodes_surf_ymax=get_ordered_node_surf(surf_ymax,curve_bottom_ymax)
    nodes_curve_xmax,orient_nodes_surf_xmax=get_ordered_node_surf(surf_xmax,curve_bottom_xmax)
    nodes_curve_ymin,orient_nodes_surf_ymin=get_ordered_node_surf(surf_ymin,curve_bottom_ymin)
    nodes_curve_xmin,orient_nodes_surf_xmin=get_ordered_node_surf(surf_xmin,curve_bottom_xmin)
    c_xminymin=Set(nodes_curve_xmin).intersection(nodes_curve_ymin)
    c_xminymax=Set(nodes_curve_xmin).intersection(nodes_curve_ymax)
    c_xmaxymin=Set(nodes_curve_xmax).intersection(nodes_curve_ymin)
    c_xmaxymax=Set(nodes_curve_xmax).intersection(nodes_curve_ymax)
    
    for n in c_xminymin:
             v = cubit.get_nodal_coordinates(n)
             orient.append(v[2])
             nd.append(n)
    result=zip(orient,nd)
    result.sort()
    c_xminymin=[c[1] for c in result]
    
    for n in c_xminymax:
             v = cubit.get_nodal_coordinates(n)
             orient.append(v[2])
             nd.append(n)
    result=zip(orient,nd)
    result.sort()
    c_xminymax=[c[1] for c in result]
    
    for n in c_xmaxymin:
             v = cubit.get_nodal_coordinates(n)
             orient.append(v[2])
             nd.append(n)
    result=zip(orient,nd)
    result.sort()
    c_xmaxymin=[c[1] for c in result]
    
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
    #
    #
    block_list=cubit.get_block_id_list()
    entities=['face']
    entity='face'
    
    if len(block_list) == 0:
        boundary_definition.define_bc(entities,parallel=False)
    #
    if iproc in list(xmin):
        xminflag=False
        refname=entity+'_abs_xmin'
        if block_list != 0:
           for block in block_list:
                ty=cubit.get_block_element_type(block)
                if ty != 'HEX8':
                    name=cubit.get_exodus_entity_name('block',block)
                    if name == refname:
                        build_block_side(surf_xmin,refname,obj=entity,id_0=block)
                        #build_block_side(surf_xmin,entity+'_abs',obj=entity,id_0=block+1)
                        xminflag=True
        if not xminflag:
            block=cubit.get_next_block_id()
            build_block_side(surf_xmin,refname,obj=entity,id_0=block)
            #build_block_side(surf_xmin,entity+'_abs',obj=entity,id_0=block+1)
    #
    if iproc in list(xmax):
        xmaxflag=False
        refname=entity+'_abs_xmax'
        if block_list != 0:
           for block in block_list:
                ty=cubit.get_block_element_type(block)
                if ty != 'HEX8':
                    name=cubit.get_exodus_entity_name('block',block)
                    if name == refname:
                        build_block_side(surf_xmax,refname,obj=entity,id_0=block)
                        #build_block_side(surf_xmax,entity+'_abs',obj=entity,id_0=block+1)
                        xmaxflag=True
        if not xmaxflag:
            block=cubit.get_next_block_id()
            build_block_side(surf_xmax,refname,obj=entity,id_0=block)
            #build_block_side(surf_xmax,entity+'_abs',obj=entity,id_0=block)
            #
    if iproc in list(ymin):
        yminflag=False
        refname=entity+'_abs_ymin'
        if block_list != 0:
           for block in block_list:
                ty=cubit.get_block_element_type(block)
                if ty != 'HEX8':
                    name=cubit.get_exodus_entity_name('block',block)
                    if name == refname:
                        build_block_side(surf_ymin,refname,obj=entity,id_0=block)
                        #build_block_side(surf_ymin,entity+'_abs',obj=entity,id_0=block)
                        yminflag=True
        if not yminflag:
            block=cubit.get_next_block_id()
            build_block_side(surf_ymin,refname,obj=entity,id_0=block)
            #build_block_side(surf_ymin,entity+'_abs',obj=entity,id_0=block)
    #
    if iproc in list(ymax):
        ymaxflag=False
        refname=entity+'_abs_ymax'
        if block_list != 0:
           for block in block_list:
                ty=cubit.get_block_element_type(block)
                if ty != 'HEX8':
                    name=cubit.get_exodus_entity_name('block',block)
                    if name == refname:
                        build_block_side(surf_ymax,refname,obj=entity,id_0=block)
                        #build_block_side(surf_ymax,entity+'_abs',obj=entity,id_0=block)
                        ymaxflag=True
        if not ymaxflag:
            block=cubit.get_next_block_id()
            build_block_side(surf_ymax,refname,obj=entity,id_0=block)
            #build_block_side(surf_ymax,entity+'_abs',obj=entity,id_0=block)
    
    return boundary
    


    