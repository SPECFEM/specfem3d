#!python
#############################################################################
# cubit2specfem3d.py                                                           #
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
#for a complete definition of the format of the mesh in SPECFEM3D_SESAME check the manual (http:/):
#
#USAGE
#
#############################################################################
#PREREQUISITE
#The mesh must be prepared
#   automatically using the module boundary_definition (see boundary_definition.py for more information)
#or
#   manually following the convention:
#     - each material should have a block defined by:
#         material domain_flag (acoustic/elastic/poroelastic)name,flag of the material (integer),p velocity
#       (or the full description: name, flag, vp, vs, rho, Q ... if not present these last 3 parameters will be
#       interpolated by module mat_parameter)
#     - each mesh should have the block definition for the face on the free_surface (topography),
#       the name of this block must be 'face_topo' or you can change the default name in mesh.topo defined in profile.
#     - each mesh should have the block definition for the faces on the absorbing boundaries,
#       one block for each surface with x=Xmin,x=Xmax,y=Ymin,y=Ymax and z=bottom. The names of
#       the blocks should contain the strings "xmin,xmax,ymin,ymax,bottom"
#
#############################################################################
#RUN
#In a python script or in the cubit python tab call:
#
#           export2SESAME(path_exporting_mesh_SPECFEM3D_SESAME)
#
#the modele create a python class for the mesh: ex. profile=mesh()
#and it export the files of the mesh needed by the partitioner of SESAME
#
#############################################################################
#OUTPUT
#The default output are 11 ASCII files:
#__________________________________________________________________________________________
#mesh_name='mesh_file' -> the file that contains the connectity of the all mesh
#    format:
#        number of elements
#        id_elements id_node1 id_node2 id_node3 id_node4 id_node5 id_node6 id_node7 id_node8
#        .....
#
#__________________________________________________________________________________________
##nodecoord_name='nodes_coords_file' -> the file that contains the coordinates of the nodes of the all mesh
#    format:
#        number of nodes
#        id_node x_coordinate y_coordinate z_coordinate
#        .....
#
#__________________________________________________________________________________________
##material_name='materials_file' -> the file that contains the material flag of the elements
#    format:
#        id_element flag
#        .....
#
#__________________________________________________________________________________________
##nummaterial_name='nummaterial_velocity_file' -> table of the material properties
#    format:
#        flag rho vp vs 0 0 #full definition of the properties, flag > 0
#        .....
#        flag 'tomography' file_name #for interpolation with tomography
#        .....
#        flag 'interface' file_name flag_for_the_gll_below_the_interface
#        flag_for_the_gll_above_the_interface #for interpolation with interface
#__________________________________________________________________________________________
##absname='absorbing_surface_file' -> this file contains all the face in all the absorbing  boundaries
##absname_local='absorbing_surface_file'+'_xmin' -> this file contains all the face in the
#                                                                                    absorbing  boundary defined by x=Xmin
##absname_local='absorbing_surface_file'+'_xmax' -> this file contains all the face in the
#                                                                                    absorbing  boundary defined by x=Xmax
##absname_local='absorbing_surface_file'+'_ymin' -> this file contains all the face in the
#                                                                                    absorbing  boundary defined by y=Ymin
##absname_local='absorbing_surface_file'+'_ymax' -> this file contains all the face in the
#                                                                                    absorbing  boundary defined by y=Ymax
##absname_local='absorbing_surface_file'+'_bottom' -> this file contains all the face in the
#                                                                                     absorbing  boundary defined by z=bottom
#    format:
#        number of faces
#        id_(element containg the face) id_node1_face id_node2_face id_node3_face id_node4_face
#        ....
#
#__________________________________________________________________________________________
##freename='free_surface_file' -> file with the hex on the free surface (usually the topography)
#    format:
#        number of faces
#        id_(element containg the face) id_node1_face id_node2_face id_node3_face id_node4_face
#
#__________________________________________________________________________________________
# it is possible save only one (or more) file singularly: for example if you want only the nodecoord_file
# call the module mesh.nodescoord_write(full path name)
#
#############################################################################

import cubit

class mtools(object):
    """docstring for ciao"""
    def __init__(self,frequency,list_surf,list_vp):
        super(mtools, self).__init__()
        self.frequency = frequency
        self.list_surf = list_surf
        self.list_vp = list_vp
        self.ngll=5
        self.percent_gll=0.172
        self.point_wavelength=5
    def __repr__(self):
        txt='Meshing for frequency up to '+str(self.frequency)+'Hz\n'
        for surf,vp in zip(self.list_surf,self.list_vp):
            txt=txt+'surface '+str(surf)+', vp ='+str(vp)+'  -> size '+str(self.freq2meshsize(vp)[0])\
                                                      +' -> dt '+str(self.freq2meshsize(vp)[0])+'\n'
        return txt
    def freq2meshsize(self,vp):
        velocity=vp*.5
        self.size=(1/2.5)*velocity/self.frequency*(self.ngll-1)/self.point_wavelength
        self.dt=.4*self.size/vp*self.percent_gll
        return self.size,self.dt
    def mesh_it(self):
        for surf,vp in zip(self.list_surf,self.list_vp):
            command = "surface "+str(surf)+" size "+str(self.freq2meshsize(vp)[0])
            cubit.cmd(command)
            command = "surface "+str(surf)+ 'scheme pave'
            cubit.cmd(command)
            command = "mesh surf "+str(surf)
            cubit.cmd(command)

class block_tools:
    def __int__(self):
        pass
    def create_blocks(self,mesh_entity,list_entity=None,):
        if mesh_entity =='surface':
            txt=' face in surface '
        elif mesh_entity == 'curve':
            txt=' edge in curve '
        elif mesh_entity == 'group':
            txt=' face in group '
        if list_entity:
            if not isinstance(list_entity,list):
                list_entity=[list_entity]
        for entity in list_entity:
            iblock=cubit.get_next_block_id()
            command = "block "+str(iblock)+ txt +str(entity)
            cubit.cmd(command)
    def material_file(self,filename):
        matfile=open(filename,'w')
        material=[]
        for record in matfile:
            mat_name,vp_str=record.split()
            vp=float(vp_str)
            material.append([mat_name,vp])
        self.material=dict(material)
    def assign_block_material(self,id_block,mat_name,vp=None):
        try:
            material=self.material
        except:
            material=None
        cubit.cmd('block '+str(id_block)+' attribute count 2')
        cubit.cmd('block '+str(id_block)+'  attribute index 1 '+str(id_block))
        if material:
            if material.has_key(mat_name):
                cubit.cmd('block '+str(id_block)+'  attribute index 2 '+str(material[mat_name]))
                print 'block '+str(id_block)+' - material '+mat_name+' - vp '+str(material[mat_name])+' from database'
        elif vp:
            cubit.cmd('block '+str(id_block)+'  attribute index 2 '+str(vp))
            print 'block '+str(id_block)+' - material '+mat_name+' - vp '+str(vp)
        else:
            print 'assignment impossible: check if '+mat_name+' is in the database or specify vp'

class mesh_tools(block_tools):
    """Tools for the mesh
    #########
    dt,edge_dt,freq,edge_freq=seismic_resolution(edges,velocity,bins_d=None,bins_u=None,sidelist=None,ngll=5,np=8)
        Given the velocity of a list of edges, seismic_resolution provides the minimum Dt
        required for the stability condition (and the corrisponding edge).
        Furthermore, given the number of gll point in the element (ngll) and the number
        of GLL point for wavelength, it provide the maximum resolved frequency.
    #########
    length=edge_length(edge)
        return the length of a edge
    #########
    edge_min,length=edge_min_length(surface)
        given the cubit id of a surface, it return the edge with minimun length
    #########
    """
    def __int__(self):
        pass
    def seismic_resolution(self,edges,velocity,bins_d=None,bins_u=None,sidelist=None):
        """
        dt,edge_dt,freq,edge_freq=seismic_resolution(edges,velocity,bins_d=None,bins_u=None,sidelist=None,ngll=5,np=8)
            Given the velocity of a list of edges, seismic_resolution provides the minimum Dt
            required for the stability condition (and the corrisponding edge).
            Furthermore, given the number of gll point in the element (ngll) and the number
            of GLL point for wavelength, it provide the maximum resolved frequency.
        """
        ratiostore=1e10
        dtstore=1e10
        edgedtstore=-1
        edgeratiostore=-1
        for edge in edges:
            d=self.edge_length(edge)
            ratio=(1/2.5)*velocity/d*(self.ngll-1)/self.point_wavelength
            dt=.4*d/velocity*self.percent_gll
            if dt<dtstore:
               dtstore=dt
               edgedtstore=edge
            if ratio < ratiostore:
                ratiostore=ratio
                edgeratiostore=edge
            try:
                for bin_d,bin_u,side in zip(bins_d,bins_u,sidelist):
                    if ratio >= bin_d and ratio < bin_u:
                        command = "sideset "+str(side)+" edge "+str(edge)
                        cubit.cmd(command)
                        #print command
                        break
            except:
                pass
        return dtstore,edgedtstore,ratiostore,edgeratiostore
    def edge_length(self,edge):
        """
        length=edge_length(edge)
            return the length of a edge
        """
        from math import sqrt
        nodes=cubit.get_connectivity('Edge',edge)
        x0,y0,z0=cubit.get_nodal_coordinates(nodes[0])
        x1,y1,z1=cubit.get_nodal_coordinates(nodes[1])
        d=sqrt((x1-x0)**2+(y1-y0)**2+(z1-z0)**2)
        return d
    def edge_min_length(self,surface):
        """
        edge_min,length=edge_min_length(surface)
            given the cubit id of a surface, it return the edge with minimun length
        """
        from math import sqrt
        self.dmin=99999
        edge_store=0
        command = "group 'list_edge' add edge in surf "+str(surface)
        command = command.replace("["," ").replace("]"," ")
        #print command
        cubit.cmd(command)
        group=cubit.get_id_from_name("list_edge")
        edges=cubit.get_group_edges(group)
        command = "delete group "+ str(group)
        cubit.cmd(command)
        for edge in edges:
            d=self.edge_length(edge)
            if d<dmin:
                self.dmin=d
                edge_store=edge
        self.edgemin=edge_store
        return self.edgemin,self.dmin
    def normal_check(self,nodes,normal):
        tres=.2
        p0=cubit.get_nodal_coordinates(nodes[0])
        p1=cubit.get_nodal_coordinates(nodes[1])
        p2=cubit.get_nodal_coordinates(nodes[2])
        a=[p1[0]-p0[0],p1[1]-p0[1],p1[2]-p0[2]]
        b=[p2[0]-p1[0],p2[1]-p1[1],p2[2]-p1[2]]
        axb=[a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]]
        dot=0.0
        for i in (0,1,2):
            dot=dot+axb[i]*normal[i]
        if  dot > 0:
            return nodes
        elif dot < 0:
            return nodes[0],nodes[3],nodes[2],nodes[1]
        else:
            print 'error: surface normal, dot=0', axb,normal,dot,p0,p1,p2
    def mesh_analysis(self,frequency):
        from sets import Set
        cubit.cmd('set info off')
        cubit.cmd('set echo off')
        cubit.cmd('set journal off')
        bins_d=[0.0001]+range(0,int(frequency)+1)+[1000]
        bins_u=bins_d[1:]
        dt=[]
        ed_dt=[]
        r=[]
        ed_r=[]
        nstart=cubit.get_next_sideset_id()
        command = "del sideset all"
        cubit.cmd(command)
        for bin_d,bin_u in zip(bins_d,bins_u):
            nsideset=cubit.get_next_sideset_id()
            command='create sideset '+str(nsideset)
            cubit.cmd(command)
            command = "sideset "+str(nsideset)+ " name "+ "'ratio-["+str(bin_d)+"_"+str(bin_u)+"['"
            cubit.cmd(command)
        nend=cubit.get_next_sideset_id()
        sidelist=range(nstart,nend)
        for block in self.block_mat:
            name=cubit.get_exodus_entity_name('block',block)
            velocity=self.material[name][1]
            if velocity > 0:
                faces=cubit.get_block_faces(block)
                edges=[]
                for face in faces:
                    es=cubit.get_sub_elements("face", face, 1)
                    edges=edges+list(es)
                edges=Set(edges)
                dtstore,edgedtstore,ratiostore,edgeratiostore=self.seismic_resolution(edges,\
                                                              velocity,bins_d,bins_u,sidelist)
                dt.append(dtstore)
                ed_dt.append(edgedtstore)
                r.append(ratiostore)
                ed_r.append(edgeratiostore)
        self.ddt=zip(ed_dt,dt)
        self.dr=zip(ed_r,r)
        def sorter(x, y):
            return cmp(x[1],y[1])
        self.ddt.sort(sorter)
        self.dr.sort(sorter)
        print self.ddt,self.dr
        print 'Deltat minimum => edge:'+str(self.ddt[0][0])+' dt: '+str(self.ddt[0][1])
        print 'minimum frequency resolved => edge:'+str(self.dr[0][0])+' frequency: '+str(self.dr[0][1])
        return self.ddt[0],self.dr[0]

class mesh(object,mesh_tools):
    def __init__(self):
        super(mesh, self).__init__()
        self.mesh_name='mesh_file'
        self.nodecoord_name='nodes_coords_file'
        self.material_name='materials_file'
        self.nummaterial_name='nummaterial_velocity_file'
        self.absname='absorbing_surface_file'
        self.freename='free_surface_file'
        self.recname='STATIONS'
        self.face='QUAD4'
        self.face2='SHELL4'
        self.hex='HEX8'
        self.edge='BAR2'
        self.topo='face_topo'
        self.rec='receivers'
        self.ngll=5
        self.percent_gll=0.172
        self.point_wavelength=5
        self.block_definition()
        cubit.cmd('compress')
    def __repr__(self):
        pass
    def block_definition(self):
        block_flag=[]
        block_mat=[]
        block_bc=[]
        block_bc_flag=[]
        material={}
        bc={}
        blocks=cubit.get_block_id_list()
        for block in blocks:
            name=cubit.get_exodus_entity_name('block',block)
            type=cubit.get_block_element_type(block)
            print block,name,blocks,type,self.hex,self.face
            # block has hexahedral elements (HEX8)
            if type == self.hex:
                flag=None
                vel=None
                vs=None
                rho=None
                q=0
                ani=0
                # material domain id
                if name.find("acoustic") >= 0 :
                  imaterial = 1
                elif name.find("elastic") >= 0 :
                  imaterial = 2
                elif name.find("poroelastic") >= 0 :
                  imaterial = 3
                else :
                  imaterial = 0
                  print "block: ",name
                  print "  could not find appropriate material for this block..."
                  print ""
                  break

                nattrib=cubit.get_block_attribute_count(block)
                if nattrib != 0:
                    # material flag:
                    #   positive => material properties,
                    #   negative => interface/tomography domain
                    flag=int(cubit.get_block_attribute_value(block,0))
                    if flag > 0 and nattrib >= 2:
                      # vp
                      vel=cubit.get_block_attribute_value(block,1)
                      if nattrib >= 3:
                        # vs
                        vs=cubit.get_block_attribute_value(block,2)
                        if nattrib >= 4:
                          #density
                          rho=cubit.get_block_attribute_value(block,3)
                          if nattrib >= 5:
                            #Q_mu
                            q=cubit.get_block_attribute_value(block,4)
                            # for q to be valid: it must be positive
                            if q < 0 :
                              print 'error, q value invalid:', q
                              break
                            if nattrib == 6:
                              #anisotropy_flag
                              ani=cubit.get_block_attribute_value(block,5)
                    elif flag < 0:
                        # velocity model
                        vel=name
                        attrib=cubit.get_block_attribute_value(block,1)
                        if attrib == 1:
                            kind='interface'
                            flag_down=cubit.get_block_attribute_value(block,2)
                            flag_up=cubit.get_block_attribute_value(block,3)
                        elif attrib == 2:
                            kind='tomography'
                else:
                    flag=block
                    vel,vs,rho,q,ani=(name,0,0,0,0)
                block_flag.append(int(flag))
                block_mat.append(block)
                if flag > 0:
                    par=tuple([imaterial,flag,vel,vs,rho,q,ani])
                elif flag < 0:
                    if kind=='interface':
                        par=tuple([imaterial,flag,kind,name,flag_down,flag_up])
                    elif kind=='tomography':
                        par=tuple([imaterial,flag,kind,name])
                elif flag==0:
                    par=tuple([imaterial,flag,name])
                material[block]=par
            elif (type == self.face) or (type == self.face2) :
                # block has surface elements (QUAD4 or SHELL4)
                block_bc_flag.append(4)
                block_bc.append(block)
                bc[block]=4 #face has connectivity = 4
                if name == self.topo: topography_face=block
            else:
                # block elements differ from HEX8/QUAD4/SHELL4
                print '****************************************'
                print 'block not properly defined:'
                print '  name:',name
                print '  type:',type
                print
                print 'please check your block definitions!'
                print
                print 'only supported types are:'
                print '  HEX8  for volumes'
                print '  QUAD4 for surface'
                print '  SHELL4 for surface'
                print '****************************************'
                continue

        nsets=cubit.get_nodeset_id_list()
        if len(nsets) == 0: self.receivers=None
        for nset in nsets:
            name=cubit.get_exodus_entity_name('nodeset',nset)
            if name == self.rec:
                self.receivers=nset
            else:
                print 'nodeset '+name+' not defined'
                self.receivers=None
        try:
            self.block_mat=block_mat
            self.block_flag=block_flag
            self.block_bc=block_bc
            self.block_bc_flag=block_bc_flag
            self.material=material
            self.bc=bc
            self.topography=topography_face
        except:
            print '****************************************'
            print 'sorry, no blocks or blocks not properly defined'
            print block_mat
            print block_flag
            print block_bc
            print block_bc_flag
            print material
            print bc
            print topography
            print '****************************************'
    def mat_parameter(self,properties):
        #note: material property acoustic/elastic/poroelastic are defined by the block's name
        print "#material properties:"
        print properties
        imaterial=properties[0]
        flag=properties[1]
        if flag > 0:
            vel=properties[2]
            if properties[3] is None and type(vel) != str:
                # velocity model scales with given vp value
                if vel >= 30:
                    m2km=1000.
                else:
                    m2km=1.
                vp=vel/m2km
                rho=(1.6612*vp-0.472*vp**2+0.0671*vp**3-0.0043*vp**4+0.000106*vp**4)*m2km
                txt='%1i %3i %20f %20f %20f %1i %1i\n' % (properties[0],properties[1],rho,vel,vel/(3**.5),0,0)
            elif type(vel) != str:
                # velocity model given as vp,vs,rho,..
                #format nummaterials file: #material_domain_id #material_id #rho #vp #vs #Q_mu #anisotropy_flag
                txt='%1i %3i %20f %20f %20f %20f %2i\n' % (properties[0],properties[1],properties[4], \
                         properties[2],properties[3],properties[5],properties[6])
            else:
                txt='%1i %3i %s \n' % (properties[0],properties[1],properties[2])
        elif flag < 0:
            if properties[2] == 'tomography':
                txt='%1i %3i %s %s\n' % (properties[0],properties[1],properties[2],properties[3])
            elif properties[2] == 'interface':
                txt='%1i %3i %s %s %1i %1i\n' % (properties[0],properties[1],properties[2],properties[3],\
                                            properties[4],properties[5])
        return txt
    def nummaterial_write(self,nummaterial_name):
        print 'Writing '+nummaterial_name+'.....'
        nummaterial=open(nummaterial_name,'w')
        for block in self.block_mat:
            #name=cubit.get_exodus_entity_name('block',block)
            nummaterial.write(self.mat_parameter(self.material[block]))
        nummaterial.close()
        print 'Ok'
    def mesh_write(self,mesh_name):
        meshfile=open(mesh_name,'w')
        print 'Writing '+mesh_name+'.....'
        num_elems=cubit.get_hex_count()
        print '  number of elements:',str(num_elems)
        meshfile.write(str(num_elems)+'\n')
        num_write=0
        for block,flag in zip(self.block_mat,self.block_flag):
            #print block,flag
            hexes=cubit.get_block_hexes(block)
            #print len(hexes)
            for hexa in hexes:
                #print hexa
                nodes=cubit.get_connectivity('Hex',hexa)
                #nodes=self.jac_check(nodes) #is it valid for 3D? TODO
                txt=('%10i ')% hexa
                txt=txt+('%10i %10i %10i %10i %10i %10i %10i %10i\n')% nodes[:]
                meshfile.write(txt)
        meshfile.close()
        print 'Ok'
    def material_write(self,mat_name):
        mat=open(mat_name,'w')
        print 'Writing '+mat_name+'.....'
        for block,flag in zip(self.block_mat,self.block_flag):
                hexes=cubit.get_block_hexes(block)
                for hexa in hexes:
                    mat.write(('%10i %10i\n') % (hexa,flag))
        mat.close()
        print 'Ok'
    def nodescoord_write(self,nodecoord_name):
        nodecoord=open(nodecoord_name,'w')
        print 'Writing '+nodecoord_name+'.....'
        node_list=cubit.parse_cubit_list('node','all')
        num_nodes=len(node_list)
        print '  number of nodes:',str(num_nodes)
        nodecoord.write('%10i\n' % num_nodes)
        #
        for node in node_list:
            x,y,z=cubit.get_nodal_coordinates(node)
            txt=('%10i %20f %20f %20f\n') % (node,x,y,z)
            nodecoord.write(txt)
        nodecoord.close()
        print 'Ok'
    def free_write(self,freename=None):
        # free surface
        cubit.cmd('set info off')
        cubit.cmd('set echo off')
        cubit.cmd('set journal off')
        from sets import Set
        normal=(0,0,1)
        if not freename: freename=self.freename
        # writes free surface file
        print 'Writing '+freename+'.....'
        freehex=open(freename,'w')
        # searches block definition with name face_topo
        for block,flag in zip(self.block_bc,self.block_bc_flag):
            if block == self.topography:
                name=cubit.get_exodus_entity_name('block',block)
                print '  block name:',name,'id:',block
                quads_all=cubit.get_block_faces(block)
                print '  number of faces = ',len(quads_all)
                dic_quads_all=dict(zip(quads_all,quads_all))
                freehex.write('%10i\n' % len(quads_all))
                list_hex=cubit.parse_cubit_list('hex','all')
                for h in list_hex:
                    faces=cubit.get_sub_elements('hex',h,2)
                    for f in faces:
                        if dic_quads_all.has_key(f):
                            #print f
                            nodes=cubit.get_connectivity('Face',f)
                            nodes_ok=self.normal_check(nodes,normal)
                            txt='%10i %10i %10i %10i %10i\n' % (h,nodes_ok[0],\
                                         nodes_ok[1],nodes_ok[2],nodes_ok[3])
                            freehex.write(txt)
        freehex.close()
        print 'Ok'
        cubit.cmd('set info on')
        cubit.cmd('set echo on')
    def abs_write(self,absname=None):
        # absorbing boundaries
        import re
        cubit.cmd('set info off')
        cubit.cmd('set echo off')
        cubit.cmd('set journal off')
        from sets import Set
        if not absname: absname=self.absname
        #
        # loops through all block definitions
        list_hex=cubit.parse_cubit_list('hex','all')
        for block,flag in zip(self.block_bc,self.block_bc_flag):
            if block != self.topography:
                name=cubit.get_exodus_entity_name('block',block)
                print name,block
                absflag=False
                if re.search('xmin',name):
                    filename=absname+'_xmin'
                    normal=(-1,0,0)
                elif re.search('xmax',name):
                    filename=absname+'_xmax'
                    normal=(1,0,0)
                elif re.search('ymin',name):
                    filename=absname+'_ymin'
                    normal=(0,-1,0)
                elif re.search('ymax',name):
                    filename=absname+'_ymax'
                    normal=(0,1,0)
                elif re.search('bottom',name):
                    filename=absname+'_bottom'
                    normal=(0,0,-1)
                elif re.search('abs',name):
                    print "  ...face_abs - not used so far..."
                    continue
                else:
                    continue
                # opens file
                print 'Writing '+filename+'.....'
                abshex_local=open(filename,'w')
                # gets face elements
                quads_all=cubit.get_block_faces(block)
                dic_quads_all=dict(zip(quads_all,quads_all))
                print '  number of faces = ',len(quads_all)
                abshex_local.write('%10i\n' % len(quads_all))
                #command = "group 'list_hex' add hex in face "+str(quads_all)
                #command = command.replace("["," ").replace("]"," ").replace("("," ").replace(")"," ")
                #cubit.cmd(command)
                #group=cubit.get_id_from_name("list_hex")
                #list_hex=cubit.get_group_hexes(group)
                #command = "delete group "+ str(group)
                #cubit.cmd(command)
                for h in list_hex:
                    faces=cubit.get_sub_elements('hex',h,2)
                    for f in faces:
                        if dic_quads_all.has_key(f):
                            nodes=cubit.get_connectivity('Face',f)
                            if not absflag:
                                # checks with specified normal
                                nodes_ok=self.normal_check(nodes,normal)
                                txt='%10i %10i %10i %10i %10i\n' % (h,nodes_ok[0],\
                                             nodes_ok[1],nodes_ok[2],nodes_ok[3])
                            else:
                                txt='%10i %10i %10i %10i %10i\n' % (h,nodes[0],\
                                             nodes[1],nodes[2],nodes[3])
                            abshex_local.write(txt)
                # closes file
                abshex_local.close()
        print 'Ok'
        cubit.cmd('set info on')
        cubit.cmd('set echo on')
    def surface_write(self,pathdir=None):
        # optional surfaces, e.g. moho_surface
        # should be created like e.g.:
        #  > block 10 face in surface 2
        #  > block 10 name 'moho_surface'
        import re
        from sets import Set
        for block in self.block_bc :
            if block != self.topography:
                name=cubit.get_exodus_entity_name('block',block)
                # skips block names like face_abs**, face_topo**
                if re.search('abs',name):
                  continue
                elif re.search('topo',name):
                  continue
                elif re.search('surface',name):
                  filename=pathdir+name+'_file'
                else:
                  continue
                # gets face elements
                print '  surface block name: ',name,'id: ',block
                quads_all=cubit.get_block_faces(block)
                print '  face = ',len(quads_all)
                if len(quads_all) == 0 :
                  continue
                # writes out surface infos to file
                print 'Writing '+filename+'.....'
                surfhex_local=open(filename,'w')
                dic_quads_all=dict(zip(quads_all,quads_all))
                # writes number of surface elements
                surfhex_local.write('%10i\n' % len(quads_all))
                # writes out element node ids
                list_hex=cubit.parse_cubit_list('hex','all')
                for h in list_hex:
                    faces=cubit.get_sub_elements('hex',h,2)
                    for f in faces:
                        if dic_quads_all.has_key(f):
                            nodes=cubit.get_connectivity('Face',f)
                            txt='%10i %10i %10i %10i %10i\n' % (h,nodes[0],\
                                             nodes[1],nodes[2],nodes[3])
                            surfhex_local.write(txt)
                # closes file
                surfhex_local.close()
        print 'Ok'
    def rec_write(self,recname):
        print 'Writing '+self.recname+'.....'
        recfile=open(self.recname,'w')
        nodes=cubit.get_nodeset_nodes(self.receivers)
        for i,n in enumerate(nodes):
            x,y,z=cubit.get_nodal_coordinates(n)
            recfile.write('ST%i XX %20f %20f 0.0 0.0 \n' % (i,x,z))
        recfile.close()
        print 'Ok'
    def write(self,path=''):
        cubit.cmd('set info off')
        cubit.cmd('set echo off')
        cubit.cmd('set journal off')
        if len(path) != 0:
            if path[-1] != '/': path=path+'/'
        # mesh file
        self.mesh_write(path+self.mesh_name)
        # mesh material
        self.material_write(path+self.material_name)
        # mesh coordinates
        self.nodescoord_write(path+self.nodecoord_name)
        # material definitions
        self.nummaterial_write(path+self.nummaterial_name)
        # free surface: face_top
        self.free_write(path+self.freename)
        # absorbing surfaces: abs_***
        self.abs_write(path+self.absname)
        # any other surfaces: ***surface***
        self.surface_write(path)
        # receivers
        if self.receivers: self.rec_write(path+self.recname)
        cubit.cmd('set info on')
        cubit.cmd('set echo on')

def export2SESAME(path_exporting_mesh_SPECFEM3D_SESAME):
    cubit.cmd('set info on')
    cubit.cmd('set echo on')
    sem_mesh=mesh()
    sem_mesh.write(path=path_exporting_mesh_SPECFEM3D_SESAME)


if __name__ == '__main__':
    path='MESH/'
    export2SESAME(path)

# call by:
# import cubit2specfem3d
# cubit2specfem3d.export2SESAME('MESH')
