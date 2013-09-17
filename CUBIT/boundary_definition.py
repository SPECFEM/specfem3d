try:
    import start as start
    cubit                   = start.start_cubit()
except:
    try:
        import cubit
    except:
        print 'error importing cubit, check if cubit is installed'
        pass

def list2str(l):
    if not isinstance(l,list): l=list(l)
    return ' '.join(str(x) for x in l)


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
    define the surfaces defining the boundaries of the volume
    """
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
        if dzmax <= 0.1 and zn > 0.4:
            top_surf.append(k)
            list_vertex=cubit.get_relatives('surface',k,'vertex')
            for v in list_vertex:
                valence=cubit.get_valence(v)
                if valence <= 4: #valence 3 is a corner, 4 is a vertex between 2 volumes, > 4 is a vertex not in the boundaries
                    lv.append(v)
        elif dzmin <= 0.001 and zn < -0.7:
            bottom_surf.append(k)
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
    #
    four_side=True
    if four_side:
        xmin,ymin,xmax,ymax=define_4side_lateral_surfaces()
        abs_xmin,abs_xmax,abs_ymin,abs_ymax=lateral_boundary_are_absorbing(ip,cpuxmin,cpuxmax,cpuymin,cpuymax,cpux,cpuy)

    return absorbing_surf,abs_xmin,abs_xmax,abs_ymin,abs_ymax,top_surf,bottom_surf,xmin,ymin,xmax,ymax



def define_block():
    list_vol=cubit.parse_cubit_list("volume","all")
    init_n_vol=len(list_vol)
    list_name=map(lambda x: 'vol'+x,map(str,list_vol))
    return list_vol,list_name

def build_block(vol_list,name,id_0=1):
    #
    from sets import Set
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
    #
    if parallel:
        absorbing_surf,abs_xmin,abs_xmax,abs_ymin,abs_ymax,top_surf,bottom_surf,xmin,ymin,xmax,ymax=define_surf(ip=ip,cpuxmin=cpuxmin,cpuxmax=cpuxmax,cpuymin=cpuymin,cpuymax=cpuymax,cpux=cpux,cpuy=cpuy)
        #
        #id_side=cubit.get_next_block_id()
        #try:
        #    entities=args[0]
        #except:
        #    entities=['face']
        #for entity in entities:
        #    build_block_side(top_surf,entity+'_topo',obj=entity,id_0=1) #topo is block 1 so id_0=1
        #    build_block_side(bottom_surf,entity+'_abs_bottom',obj=entity,id_0=2)
        #    build_block_side(xmin,entity+'_xmin',obj=entity,id_0=3)
        #    build_block_side(ymin,entity+'_ymin',obj=entity,id_0=4)
        #    build_block_side(xmax,entity+'_xmax',obj=entity,id_0=5)
        #    build_block_side(ymax,entity+'_ymax',obj=entity,id_0=6)
        #     
        id_0=cubit.get_next_block_id()
        v_list,name_list=define_block()
        build_block(v_list,name_list,id_0)
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
    cubit.cmd('del group sl')
    cubit.cmd("group 'sl' add node in surf "+lsurf)
    group1 = cubit.get_id_from_name("sl")
    nodes_ls =list(cubit.get_group_nodes(group1))
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
    #
    #
    #
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
        if ty != 'HEX8': cubit.cmd('del block '+str(block))
        build_block_side(top_surf,refname,obj=entity,id_0=1001)
        #
        refname=entity+'_bottom'
        block=4 #change here..... must be 2 @@@@@@@@
        ty=None
        ty=cubit.get_block_element_type(block)
        if ty != 'HEX8': cubit.cmd('del block '+str(block))
        build_block_side(bottom_surf,refname,obj=entity,id_0=1002)
    #
    #
    return boundary