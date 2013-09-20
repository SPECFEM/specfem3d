#!python
# Retrieving absorbing boundaries.
#    P. Galvez (ETH-Zurich, 10.09.2011):
#    This function is based on Emmanuele Cassarotti , boundary_definition.py routine. 
#    
#    It returns absorbing_surf,absorbing_surf_xmin,absorbing_surf_xmax,
#    absorbing_surf_ymin,absorbing_surf_ymax,absorbing_surf_bottom,topo_surf
#    where absorbing_surf is the list of all the absorbing boundary surf
#    absorbing_surf_xmin is the list of the absorbing boundary surfaces that correnspond to x=xmin
#    ...
#    absorbing_surf_bottom is the list of the absorbing boundary surfaces that correspond to z=zmin

class abs_surface:
   def __init__(self,xmin,xmax,ymin,ymax): 
       self.xmin = xmin
       self.xmax = xmax
       self.ymin = ymin
       self.ymax = ymax

class abs_surface_topo:
   def __init__(self,xmin,xmax,ymin,ymax,bottom,topo): 
       self.xmin = xmin
       self.xmax = xmax
       self.ymin = ymin
       self.ymax = ymax
       self.bottom = bottom
       self.topo = topo
 
# Emmanuele Cassarotti function for Parallel absorbing boundaries.
# WARNING : absorbing.surf deleted due to CUBIT 13.0 does not allow elements beloging to diferent blocks.

def define_parallel_absorbing_surf():
    """
    define the absorbing surfaces for a layered topological box where boundary are surfaces parallel to the axis.
    it returns absorbing_surf,absorbing_surf_xmin,absorbing_surf_xmax,absorbing_surf_ymin,absorbing_surf_ymax,absorbing_surf_bottom,topo_surf
    where
    absorbing_surf is the list of all the absorbing boundary surf
    absorbing_surf_xmin is the list of the absorbing boundary surfaces that correnspond to x=xmin
    ...
    absorbing_surf_bottom is the list of the absorbing boundary surfaces that correspond to z=zmin
    """
    try:
        cubit.cmd('comment')
    except:
        try:
            import cubit
            cubit.init([""])
        except:
            print 'error importing cubit'
            import sys
            sys.exit()
    absorbing_surf_xmin=[]
    absorbing_surf_xmax=[]
    absorbing_surf_ymin=[]
    absorbing_surf_ymax=[]
    absorbing_surf_bottom=[]
    top_surf=[]
    
    
    list_vol=cubit.parse_cubit_list("volume","all")
    init_n_vol=len(list_vol)
    zmax_box=cubit.get_total_bounding_box("volume",list_vol)[7]
    zmin_box=cubit.get_total_bounding_box("volume",list_vol)[6] #it is the z_min of the box ... box= xmin,xmax,d,ymin,ymax,d,zmin...    
    xmin_box=cubit.get_total_bounding_box("volume",list_vol)[0]
    xmax_box=cubit.get_total_bounding_box("volume",list_vol)[1]
    ymin_box=cubit.get_total_bounding_box("volume",list_vol)[3]
    ymax_box=cubit.get_total_bounding_box("volume",list_vol)[4]
    list_surf=cubit.parse_cubit_list("surface","all")
    print '##boundary box: '
    print '##  x min: ' + str(xmin_box)
    print '##  y min: ' + str(ymin_box)
    print '##  z min: ' + str(zmin_box)
    print '##  x max: ' + str(xmax_box)
    print '##  y max: ' + str(ymax_box)
    print '##  z max: ' + str(zmax_box)

    #box lengths
    x_len = abs( xmax_box - xmin_box)
    y_len = abs( ymax_box - ymin_box)
    z_len = abs( zmax_box - zmin_box)
    
    print '##boundary box: '
    print '##  x length: ' + str(x_len)
    print '##  y length: ' + str(y_len)
    print '##  z length: ' + str(z_len)
    
    # tolerance parameters 
    absorbing_surface_distance_tolerance=0.005
    topographic_surface_distance_tolerance=0.001
    topographic_surface_normal_tolerance=0.2
        
    for k in list_surf:
        center_point = cubit.get_center_point("surface", k)
        if abs((center_point[0] - xmin_box)/x_len) <= absorbing_surface_distance_tolerance:
             absorbing_surf_xmin.append(k)
        elif abs((center_point[0] - xmax_box)/x_len) <= absorbing_surface_distance_tolerance:
             absorbing_surf_xmax.append(k)
        elif abs((center_point[1] - ymin_box)/y_len) <= absorbing_surface_distance_tolerance:
             absorbing_surf_ymin.append(k)
        elif abs((center_point[1] - ymax_box)/y_len) <= absorbing_surface_distance_tolerance:
             absorbing_surf_ymax.append(k)
        elif abs((center_point[2] - zmin_box)/z_len) <= absorbing_surface_distance_tolerance:
             print 'center_point[2]' + str(center_point[2])
             print 'kz:' + str(k)
             absorbing_surf_bottom.append(k)
                       
        else:
            sbox=cubit.get_bounding_box('surface',k)
            dz=abs((sbox[7] - zmax_box)/z_len)
            normal=cubit.get_surface_normal(k)
            zn=normal[2]
            dn=abs(zn-1)
            if dz <= topographic_surface_distance_tolerance and dn < topographic_surface_normal_tolerance:
                top_surf.append(k)

    return absorbing_surf_xmin,absorbing_surf_xmax,absorbing_surf_ymin,absorbing_surf_ymax,absorbing_surf_bottom,top_surf

def define_top_bottom_absorbing_surf(zmin_box,zmax_box):
    """
      absorbing_surf_bottom is the list of the absorbing boundary surfaces that correspond to z=zmin
    """
    try:
        cubit.cmd('comment')
    except:
        try:
            import cubit
            cubit.init([""])
        except:
            print 'error importing cubit'
            import sys
            sys.exit()
    absorbing_surf_bottom=[]
    top_surf = []
    
    list_vol=cubit.parse_cubit_list("volume","all")
    init_n_vol=len(list_vol)
#   TO DO : Make zmin_box work properly.
#   zmax_box=cubit.get_total_bounding_box("volume",list_vol)[7]
#   zmin_box=cubit.get_total_bounding_box("volume",list_vol)[6] #it is the z_min of the box ... box= xmin,xmax,d,ymin,ymax,d,zmin...
    xmin_box=cubit.get_total_bounding_box("volume",list_vol)[0]
    xmax_box=cubit.get_total_bounding_box("volume",list_vol)[1]
    ymin_box=cubit.get_total_bounding_box("volume",list_vol)[3]
    ymax_box=cubit.get_total_bounding_box("volume",list_vol)[4]
    list_surf=cubit.parse_cubit_list("surface","all")
   
    print '##boundary box: '
    print '##  x min: ' + str(xmin_box)
    print '##  y min: ' + str(ymin_box)
    print '##  z min: ' + str(zmin_box)
    print '##  x max: ' + str(xmax_box)
    print '##  y max: ' + str(ymax_box)
    print '##  z max: ' + str(zmax_box)

    #box lengths
    x_len = abs( xmax_box - xmin_box)
    y_len = abs( ymax_box - ymin_box)
    z_len = abs( zmax_box - zmin_box)
    
    print '##boundary box: '
    print '##  x length: ' + str(x_len)
    print '##  y length: ' + str(y_len)
    print '##  z length: ' + str(z_len)
    
    # tolerance parameters 
    absorbing_surface_distance_tolerance=0.005
    topographic_surface_distance_tolerance=0.001
    topographic_surface_normal_tolerance=0.2

    for k in list_surf:
        center_point = cubit.get_center_point("surface", k)
        if abs((center_point[2] - zmin_box)/z_len) <= absorbing_surface_distance_tolerance:
             print 'center_point[2]' + str(center_point[2])
             print 'kz:' + str(k)
             absorbing_surf_bottom.append(k)
   
        else:
            sbox=cubit.get_bounding_box('surface',k)
            dz=abs((sbox[7] - zmax_box)/z_len)
            normal=cubit.get_surface_normal(k)
            zn=normal[2]
            dn=abs(zn-1)
            if dz <= topographic_surface_distance_tolerance and dn < topographic_surface_normal_tolerance:
                top_surf.append(k)
    
    return absorbing_surf_bottom,top_surf


def build_block(vol_list,name):
    from sets import Set
    try:
            cubit.cmd('comment')
    except:
            try:
                import cubit
                cubit.init([""])
            except:
                print 'error importing cubit'
                import sys
                sys.exit()
    block_list=cubit.get_block_id_list()
    if len(block_list) > 0:
        id_block=max(block_list)
    else:
        id_block=0
    for v,n in zip(vol_list,name):
       id_block+=1
       v_other=Set(vol_list)-Set([v])
       command= 'block '+str(id_block)+' hex in vol '+str(v)
       command = command.replace("["," ").replace("]"," ")
       cubit.cmd(command) 
       command = "block "+str(id_block)+" name '"+n+"'"
       cubit.cmd(command)


def define_block():
    """ 
     Renumbering number of volumes from 1 to NVOLUMES.
    """ 
    try:
            cubit.cmd('comment')
    except:
            try:
                import cubit
                cubit.init([""])
            except:
                print 'error importing cubit'
                import sys
                sys.exit()
    list_vol=cubit.parse_cubit_list("volume","all")
    init_n_vol=len(list_vol)
    list_name=map(lambda x: 'vol'+x,map(str,list_vol))
    return list_vol,list_name

 
def build_block_side(surf_list,name,obj='surface'):
    try:
            cubit.cmd('comment')
    except:
            try:
                import cubit
                cubit.init([""])
            except:
                print 'error importing cubit'
                import sys
                sys.exit()
    id_nodeset=cubit.get_next_nodeset_id()
    id_block=cubit.get_next_block_id()
    
    if obj == 'hex':
        txt='hex in node in surface'
        txt1='block '+str(id_block)+ ' '+ txt +' '+str(list(surf_list))
        txt2="block "+str(id_block)+" name '"+name+"'"
        txt1=txt1.replace("["," ").replace("]"," ")
    elif obj == 'node':
         txt=obj+' in surface'
         txt1= 'nodeset '+str(id_nodeset)+ ' '+ txt +' '+str(list(surf_list))
         txt1 = txt1.replace("["," ").replace("]"," ")
         txt2 = "nodeset "+str(id_nodeset)+" name '"+name+"'"
    elif obj == 'face' or obj == 'edge':
        txt=obj+' in surface'
        txt1= 'block '+str(id_block)+ ' '+ txt +' '+str(list(surf_list))
        txt1 = txt1.replace("["," ").replace("]"," ")
        txt2 = "block "+str(id_block)+" name '"+name+"'"
    else:
        txt1=''
        # do not execute: block id might be wrong
        print "##block "+str(id_block)+" name '"+name+"_notsupported (only hex,face,edge,node)'"
        txt2=''
   
 
    cubit.cmd(txt1)
    cubit.cmd(txt2)


def define_bc(entities,zmin,zmax,self):
     # Temporal : Variable zmin should be obtained automatically. 
     xmin = self.xmin
     xmax = self.xmax
     ymin = self.ymin
     ymax = self.ymax
     bottom,topo=define_top_bottom_absorbing_surf(zmin,zmax)
     v_list,name_list=define_block()
     build_block(v_list,name_list)
     print entities
     for entity in entities:
         print "##entity: "+str(entity)
         build_block_side(xmin,entity+'_abs_xmin',obj=entity)
         build_block_side(xmax,entity+'_abs_xmax',obj=entity)
         build_block_side(ymin,entity+'_abs_ymin',obj=entity)
         build_block_side(ymax,entity+'_abs_ymax',obj=entity)
         build_block_side(bottom,entity+'_abs_bottom',obj=entity)
         build_block_side(topo,entity+'_topo',obj=entity)

def define_parallel_bc(entities):
     xmax = []
     ymin = []
     ymax = []
     zmin = []
     zmax = []
     #Extracting parallel surfaces.
     xmin,xmax,ymin,ymax,bottom,topo=define_parallel_absorbing_surf()
     v_list,name_list=define_block()
     build_block(v_list,name_list)
     print entities
     for entity in entities:
         print "##entity: "+str(entity)
         build_block_side(xmin,entity+'_abs_xmin',obj=entity)
         build_block_side(xmax,entity+'_abs_xmax',obj=entity)
         build_block_side(ymin,entity+'_abs_ymin',obj=entity)
         build_block_side(ymax,entity+'_abs_ymax',obj=entity)
         build_block_side(bottom,entity+'_abs_bottom',obj=entity)
         build_block_side(topo,entity+'_topo',obj=entity)


def define_boundaries(entities,xmin,xmax,ymin,ymax,zmin,zmax):
     bottom=zmin
     topo=zmax
     v_list,name_list=define_block()
     build_block(v_list,name_list)
     print entities
     for entity in entities:
         print "##entity: "+str(entity)
         build_block_side(xmin,entity+'_abs_xmin',obj=entity)
         build_block_side(xmax,entity+'_abs_xmax',obj=entity)
         build_block_side(ymin,entity+'_abs_ymin',obj=entity)
         build_block_side(ymax,entity+'_abs_ymax',obj=entity)
         build_block_side(bottom,entity+'_abs_bottom',obj=entity)
         build_block_side(topo,entity+'_topo',obj=entity)

def define_bc_topo(entities,self):
     # Temporal : Variable zmin should be obtained automatically. 
     xmin = self.xmin
     xmax = self.xmax
     ymin = self.ymin
     ymax = self.ymax
     bottom = self.bottom
     topo = self.topo
     v_list,name_list=define_block()
     build_block(v_list,name_list)
     print entities
     for entity in entities:
         print "##entity: "+str(entity)
         build_block_side(xmin,entity+'_abs_xmin',obj=entity)
         build_block_side(xmax,entity+'_abs_xmax',obj=entity)
         build_block_side(ymin,entity+'_abs_ymin',obj=entity)
         build_block_side(ymax,entity+'_abs_ymax',obj=entity)
         build_block_side(bottom,entity+'_abs_bottom',obj=entity)
         build_block_side(topo,entity+'_topo',obj=entity)


