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
    absorbing_surf=[]
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
#    for k in list_surf:
#        center_point = cubit.get_center_point("surface", k)
#        if abs((center_point[0] - xmin_box)/xmin_box) <= 0.005:
#             absorbing_surf_xmin.append(k)
#             absorbing_surf.append(k)
#        elif abs((center_point[0] - xmax_box)/xmax_box) <= 0.005:
#             absorbing_surf_xmax.append(k)
#             absorbing_surf.append(k)
#        elif abs((center_point[1] - ymin_box)/ymin_box) <= 0.005:
#             absorbing_surf_ymin.append(k)
#             absorbing_surf.append(k)
#        elif abs((center_point[1] - ymax_box)/ymax_box) <= 0.005:
#             absorbing_surf_ymax.append(k)
#             absorbing_surf.append(k)
#        elif abs((center_point[2] - zmin_box)/zmin_box) <= 0.005:
#             absorbing_surf_bottom.append(k)
#             absorbing_surf.append(k)
#        else:
#            sbox=cubit.get_bounding_box('surface',k)
#            dz=abs((sbox[7] - zmax_box)/zmax_box)
#            normal=cubit.get_surface_normal(k)
#            zn=normal[2]
#            dn=abs(zn-1)
#            if dz <= 0.001 and dn < 0.2:
#                top_surf.append(k)

    #box lengths
    x_len = abs( xmax_box - xmin_box)
    y_len = abs( ymax_box - ymin_box)
    z_len = abs( zmax_box - zmin_box)

    print '##boundary box: '
    print '##  x length: ' + str(x_len)
    print '##  y length: ' + str(y_len)
    print '##  z length: ' + str(z_len)

    # debug
    print '##  xmin: ' + str(xmin_box)
    print '##  xmax: ' + str(xmax_box)
    print '##  ymin: ' + str(ymin_box)
    print '##  ymax: ' + str(ymax_box)
    print '##  zmin: ' + str(zmin_box)
    print '##  zmax: ' + str(zmax_box)

    ############################################
    ##
    ## tolerance parameters
    ##
    ## modified for surface topography
    ############################################
    absorbing_surface_distance_tolerance=0.1
    topographic_surface_distance_tolerance=0.1
    topographic_surface_normal_tolerance=0.3

    for k in list_surf:
        center_point = cubit.get_center_point("surface", k)
        
        #debug
        print '##surface: ' + str(k)
        print '## center point: ' + str(center_point)

        if abs((center_point[0] - xmin_box)/x_len) <= absorbing_surface_distance_tolerance:
          #debug 
          print '## xmin surface: ' + str(k)
          absorbing_surf_xmin.append(k)
          absorbing_surf.append(k)
        elif abs((center_point[0] - xmax_box)/x_len) <= absorbing_surface_distance_tolerance:
          #debug 
          print '## xmax surface: ' + str(k)
          absorbing_surf_xmax.append(k)
          absorbing_surf.append(k)
        elif abs((center_point[1] - ymin_box)/y_len) <= absorbing_surface_distance_tolerance:
          #debug 
          print '## ymin surface: ' + str(k)
          absorbing_surf_ymin.append(k)
          absorbing_surf.append(k)
        elif abs((center_point[1] - ymax_box)/y_len) <= absorbing_surface_distance_tolerance:
          #debug 
          print '## ymax surface: ' + str(k)
          absorbing_surf_ymax.append(k)
          absorbing_surf.append(k)
        elif abs((center_point[2] - zmin_box)/z_len) <= absorbing_surface_distance_tolerance:
          #debug 
          print '## bottom surface: ' + str(k)
          absorbing_surf_bottom.append(k)
          absorbing_surf.append(k)
        else:
          sbox=cubit.get_bounding_box('surface',k)
          dz=abs((sbox[7] - zmax_box)/z_len)
          normal=cubit.get_surface_normal(k)
          zn=normal[2]
          dn=abs(abs(zn)-1)
          #debug 
          #print '## surface element: ' + str(k)
          #print '## surface element: zn ' + str(zn)
          #print '## surface element: dn ' + str(dn)
          #print '## surface element: dz ' + str(dz)
          if dz <= topographic_surface_distance_tolerance and dn < topographic_surface_normal_tolerance:
            #debug 
            print '## topo surface: ' + str(k)          
            top_surf.append(k)
    
    return absorbing_surf,absorbing_surf_xmin,absorbing_surf_xmax,absorbing_surf_ymin,absorbing_surf_ymax,absorbing_surf_bottom,top_surf

def define_absorbing_surf_nopar():
    """
    define the absorbing surfaces for a layered topological box where boundary surfaces are not parallel to the axis.
    it returns absorbing_surf,topo_surf
    where
    absorbing_surf is the list of all the absorbing boundary surf
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
            if dzmax <= 0.001 and zn > 0.7:
                top_surf.append(k)
                list_vertex=cubit.get_relatives('surface',k,'vertex')
                for v in list_vertex:
                    valence=cubit.get_valence(v)
                    if valence <= 4: #valence 3 is a corner, 4 is a vertex between 2 volumes, > 4 is a vertex not in the boundaries
                        lv.append(v)
            elif dzmin <= 0.001 and zn < -0.7:
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
            if abs((center_point[0] - p[0])/p[0]) <= 0.005 and abs((center_point[1] - p[1])/p[1]) <= 0.005:
             absorbing_surf.append(k)
             break
    return absorbing_surf,top_surf

def define_absorbing_surf_sphere():
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
    surf=[]
    list_surf=cubit.parse_cubit_list("surface","all")
    for s in list_surf:
       v=cubit.get_relatives('surface',s,'volume')
       if len(v) == 1:
           surf.append(s)
    return surf

def define_block():
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
       #command= 'block '+str(id_block)+' hex in node in vol '+str(v)+' except hex in vol '+str(list(v_other))
       command= 'block '+str(id_block)+' hex in vol '+str(v)
       command = command.replace("["," ").replace("]"," ")
       cubit.cmd(command)
       command = "block "+str(id_block)+" name '"+n+"'"
       cubit.cmd(command)

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

    # creates command string
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

    # executes commands
    print "# command: " + txt1
    print "# command: " + txt2
    cubit.cmd(txt1)
    cubit.cmd(txt2)

def define_bc(*args,**keys):
    parallel=keys.get('parallel',True)
    closed=keys.get('closed',False)
    if not closed:
        print "##open region"

        # model with parallel sides (e.g. a block)
        if parallel:
            surf,xmin,xmax,ymin,ymax,bottom,topo=define_absorbing_surf()
        else:
            # arbitrary geometry
            surf,topo=define_absorbing_surf_nopar()

        v_list,name_list=define_block()
        build_block(v_list,name_list)
        entities=args[0]
        print entities
        for entity in entities:
            print "##entity: "+str(entity)

            # block for free surface (w/ topography)
            #print '## topo surface block: ' + str(topo)
            if len(topo) == 0:
              print ""
              print "no topo surface found, please create block face_topo manually..."
              print ""
            else:
              build_block_side(topo,entity+'_topo',obj=entity)

            # model has parallel sides (e.g. a block model )
            if parallel:
            
              # blocks for each side
              if len(xmin) == 0:
                print ""
                print "no abs_xmin surface found, please create block manually..."
                print ""
              else:
                build_block_side(xmin,entity+'_abs_xmin',obj=entity)
                
              # blocks for each side
              if len(xmax) == 0:
                print ""
                print "no abs_xmax surface found, please create block manually..."
                print ""
              else:
                build_block_side(xmax,entity+'_abs_xmax',obj=entity)

              # blocks for each side
              if len(ymin) == 0:
                print ""
                print "no abs_xmin surface found, please create block manually..."
                print ""
              else:
                build_block_side(ymin,entity+'_abs_ymin',obj=entity)

              # blocks for each side
              if len(ymax) == 0:
                print ""
                print "no abs_ymax surface found, please create block manually..."
                print ""
              else:
                build_block_side(ymax,entity+'_abs_ymax',obj=entity)

              # blocks for each side
              if len(bottom) == 0:
                print ""
                print "no abs_bottom surface found, please create block manually..."
                print ""
              else:
                build_block_side(bottom,entity+'_abs_bottom',obj=entity)

              # block for all sides together
              # NOTE:
              #    this might fail in some CUBIT versions, when elements are already
              #    assigned to other blocks
              try:
                build_block_side(surf,entity+'_abs',obj=entity)
              except:
                print "no combined surface with all sides created"

            else:
                # arbitrary geometry
                # puts all elements in single block
                build_block_side(surf,entity+'_abs',obj=entity)

    else:
        print "##closed region"

        # model without absorbing boundaries, only one surface, e.g. a sphere
        surf=define_absorbing_surf_sphere()

        v_list,name_list=define_block()
        build_block(v_list,name_list)

        entities=args[0]
        for entity in entities:
            # puts all elements in single block
            build_block_side(surf,entity+'_closedvol',obj=entity)





## calling example:

#entities=['face']
#define_bc(entities,parallel=True)
#define_bc(entities,parallel=False)
#define_bc(entities,parallel=False,closed=True)

## block material assigning example:

#block 1  attribute count 5
#block 2  attribute count 0
#block 2  name '/prova/interface1'
#block 2  attribute count 3
#block 3  attribute count 5
#block 1  attribute index 2 1500
#block 1  attribute count 2
#block 3  attribute index 1 2
#block 3  attribute index 2 5800
#block 3  attribute index 3 3900
#block 3  attribute index 4 1500
#block 3  attribute index 5 3.5
#block 1  name 'top'
#block 3  name 'bottom'
#block 2  attribute index 1 -1
#block 2  attribute index 2 2
#block 2  attribute count 4
#block 2  attribute index 2 1
#block 2  attribute index 3 2

# to create block manually:
#
# use a commands like:
#
# e.g. surface with topography
# block 2 face in surface 6 
# block 2 name "face_topo"
#
# e.g. all surface which are absorbing
# block 3 face in surface 1 2 3 4 5
# block 2 name "face_abs"


