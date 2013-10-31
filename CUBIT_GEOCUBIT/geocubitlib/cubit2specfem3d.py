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
#for a complete definition of the format of the mesh in SPECFEM3D check the manual (http://www.geodynamics.org/cig/software/specfem3d):
#
#USAGE 
#
#############################################################################
#PREREQUISITE
#The mesh must be prepared 
#   automatically using the module boundary_definition (see boundary_definition.py for more information)
#or 
#   manually following the convention:
#     - each material should have a block defined by material domain_flag (acoustic/elastic/poroelastic) name,flag of the material (integer),p velocity 
#       (or the full description: name, flag, vp, vs, rho, Q ... if not present these last 3 parameters will be interpolated by module mat_parameter)
#     - each mesh should have the block definition for the face on the free_surface (topography), 
#       the name of this block must be 'face_topo' or you can change the default name in mesh.topo defined in profile.
#     - each mesh should have the block definition for the faces on the absorbing boundaries, 
#       one block for each surface with x=Xmin,x=Xmax,y=Ymin,y=Ymax and z=bottom. The names of the blocks should contain the strings "xmin,xmax,ymin,ymax,bottom"
#
#############################################################################
#RUN
#In a python script or in the cubit python tab call:
#           
#           export2SPECFEM3D(path_exporting_mesh_SPECFEM3D) 
#
#the module creates a python class for the mesh: ex. profile=mesh()
#and it export the files of the mesh needed by the partitioner of SPECFEM3D
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
#        #material_domain_id #material_id #rho #vp #vs #Q_mu #anisotropy
#        .....
#        #material_domain_id 'tomography' file_name #for interpolation with tomography
#        .....
#        #material_domain_id 'interface' file_name flag_for_the_gll_below_the_interface flag_for_the_gll_above_the_interface #for interpolation with interface
#__________________________________________________________________________________________        
##absname='absorbing_surface_file' -> this file contains all the face in all the absorbing  boundaries
##absname_local='absorbing_surface_file'+'_xmin' -> this file contains all the face in the absorbing  boundary defined by x=Xmin
##absname_local='absorbing_surface_file'+'_xmax' -> this file contains all the face in the absorbing  boundary defined by x=Xmax
##absname_local='absorbing_surface_file'+'_ymin' -> this file contains all the face in the absorbing  boundary defined by y=Ymin
##absname_local='absorbing_surface_file'+'_ymax' -> this file contains all the face in the absorbing  boundary defined by y=Ymax
##absname_local='absorbing_surface_file'+'_bottom' -> this file contains all the face in the absorbing  boundary defined by z=bottom
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
#__________________________________________________________________________________________
##surface='*_surface_file' -> file with the hex on any surface (define by the word 'surface' in the name of the block, ex: moho_surface)
# optional surfaces, e.g. moho_surface
# should be created like e.g.:
#  > block 10 face in surface 2
#  > block 10 name 'moho_surface'
#
#
#    format:
#        number of faces
#        id_(element containg the face) id_node1_face id_node2_face id_node3_face id_node4_face
#
#__________________________________________________________________________________________
# it is possible save only one (or more) file singularly: for example if you want only the nodecoord_file call the module mesh.nodescoord_write(full path name)
#
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

from utilities import get_cubit_version

class mtools(object):
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
            txt=txt+'surface '+str(surf)+', vp ='+str(vp)+'  -> size '+str(self.freq2meshsize(vp)[0])+' -> dt '+str(self.freq2meshsize(vp)[0])+'\n' 
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

class block_tools():
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
        Given the velocity of a list of edges, seismic_resolution provides the minimum Dt required for the stability condition (and the corrisponding edge).
        Furthermore, given the number of gll point in the element (ngll) and the number of GLL point for wavelength, it provide the maximum resolved frequency.
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
            Given the velocity of a list of edges, seismic_resolution provides the minimum Dt required for the stability condition (and the corrisponding edge).
            Furthermore, given the number of gll point in the element (ngll) and the number of GLL point for wavelength, it provide the maximum resolved frequency.
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
                dtstore,edgedtstore,ratiostore,edgeratiostore=self.seismic_resolution(edges,velocity,bins_d,bins_u,sidelist)
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
    def __init__(self,hex27=False,cpml=False,cpml_size=False,top_absorbing=False):
        super(mesh, self).__init__()
        self.mesh_name='mesh_file'
        self.nodecoord_name='nodes_coords_file'
        self.material_name='materials_file'
        self.nummaterial_name='nummaterial_velocity_file'
        self.absname='absorbing_surface_file'
        self.cpmlname='absorbing_cpml_file'
        self.freename='free_or_absorbing_surface_file_zmax'
        self.recname='STATIONS'
        version_cubit=get_cubit_version()
        if version_cubit >= 12:
            self.face='SHELL4'
        else:
            self.face='QUAD4'
        self.hex='HEX'
        if version_cubit <= 13:
            if hex27: 
                print "ATTENTION **********************\n\nCubit <= 12.2 doesn't support HEX27\nassuming HEX8 .....\n\n"
            self.hex27=False
        else:
            self.hex27=hex27
        self.edge='BAR2'
        self.topo='face_topo'
        self.topography=None
        self.free=None
        self.freetxt='free'
        self.rec='receivers'
        self.cpml=cpml
        if cpml:
            if cpml_size:
                self.size=cpml_size
            else:
                print 'please specify cmpl size if you want to use cpml'
        self.top_absorbing=top_absorbing
        if hex27: cubit.cmd('block all except block 1001 1002 1003 1004 1005 1006 element type hex27')
        self.block_definition()
        self.ngll=5
        self.percent_gll=0.172
        self.point_wavelength=5
        self.xmin=False
        self.ymin=False
        self.zmin=False
        self.xmax=False
        self.ymax=False
        self.zmax=False
        cubit.cmd('compress all')
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
            ty=cubit.get_block_element_type(block)
            #print block,blocks,ty,self.hex,self.face
            if self.hex in ty:
                nattrib=cubit.get_block_attribute_count(block)
                flag=None
                vel=None
                vs=None
                rho=None
                q=0
                ani=0
                # material domain id
                if "acoustic" in name :
                  imaterial = 1
                elif "elastic"  in name:
                  imaterial = 2
                elif "poroelastic"  in name:
                  imaterial = 3
                else :
                  imaterial = 0
                #
                if nattrib > 1:
                    # material flag:
                    #   positive => material properties,
                    #   negative => interface/tomography domain
                    flag=int(cubit.get_block_attribute_value(block,0))
                    if 0< flag and nattrib >= 2:
                        vel=cubit.get_block_attribute_value(block,1)
                        if nattrib >= 3:
                            vs=cubit.get_block_attribute_value(block,2)
                            if nattrib >= 4:
                                rho=cubit.get_block_attribute_value(block,3)
                                if nattrib >= 5:
                                    qk=cubit.get_block_attribute_value(block,4)
                                    if nattrib >= 6:
                                        qmu=cubit.get_block_attribute_value(block,4)
                                    # for q to be valid: it must be positive
                                    if qk < 0 or qmu<0:
                                      print 'error, q value invalid:', qk,qmu
                                      break                                                   
                                    if nattrib == 7:
                                        ani=cubit.get_block_attribute_value(block,5)
                    elif flag < 0:
                        vel=name
                        attrib=cubit.get_block_attribute_value(block,1)
                        if attrib == 1: 
                            kind='interface'
                            flag_down=cubit.get_block_attribute_value(block,2)
                            flag_up=cubit.get_block_attribute_value(block,3)
                        elif attrib == 2:
                            kind='tomography'
                elif  nattrib == 1:
                    flag=cubit.get_block_attribute_value(block,0)
                    #print 'only 1 attribute ', name,block,flag
                    vel,vs,rho,qk,qmu,ani=(0,0,0,9999.,9999.,0)
                else:
                    flag=block
                    vel,vs,rho,qk,qmu,ani=(name,0,0,9999.,9999.,0)
                block_flag.append(int(flag))
                block_mat.append(block)
                if (flag > 0) and nattrib != 1:
                    par=tuple([imaterial,flag,vel,vs,rho,qk,qmu,ani])
                elif flag < 0 and nattrib != 1:
                    if kind=='interface':
                        par=tuple([imaterial,flag,kind,name,flag_down,flag_up])
                    elif kind=='tomography':
                        par=tuple([imaterial,flag,kind,name])
                elif flag==0 or nattrib == 1:
                    par=tuple([imaterial,flag,name])
                material[block]=par
            elif ty == self.face:
                block_bc_flag.append(4)
                block_bc.append(block)
                bc[block]=4 #face has connectivity = 4
                if name == self.topo or block == 1001: 
                    self.topography=block
                if self.freetxt in name:
                    self.free=block
            elif ty == 'SPHERE':
                pass
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
                print '  HEX/HEX8  for volumes'
                print '  QUAD4 for surface'
                print '  SHELL4 for surface'
                print '****************************************'
                continue
                return None, None,None,None,None,None,None,None
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
            print 'HEX Blocks:'
            for m,f in zip(self.block_mat,self.block_flag):
                print 'block ',m,'material flag ',f
            print 'Absorbing Boundary Conditions:'
            for m,f in zip(self.block_bc,self.block_bc_flag):
                print  'bc ',m,'bc flag ',f
            print 'Topography (free surface)'
            print self.topography
            print 'Free surface'
            print self.free
        except:
            print '****************************************'
            print 'sorry, no blocks or blocks not properly defined'
            print block_mat
            print block_flag
            print block_bc
            print block_bc_flag
            print material
            print bc
            print '****************************************'
            
            
            
            
            
            
    def get_hex_connectivity(self,ind):
        if self.hex27:
                cubit.silent_cmd('group "nh" add Node in hex '+str(ind))
                group1 = cubit.get_id_from_name("nh")
                result=cubit.get_group_nodes(group1)
                if len(result) != 27: print 'error hex27'
                cubit.cmd('del group '+str(group1))
        else:
            result=cubit.get_connectivity('hex',ind)
        return result
    #
    def get_face_connectivity(self,ind):
        if self.hex27:
                cubit.silent_cmd('group "nf" add Node in face '+str(ind))
                group1 = cubit.get_id_from_name("nf")
                result=cubit.get_group_nodes(group1)
                cubit.cmd('del group '+str(group1))
        else:
            result=cubit.get_connectivity('face',ind)
        return result        
    
    
    def mat_parameter(self,properties): 
        #print properties
        #format nummaterials file: #material_domain_id #material_id #rho #vp #vs #Q_kappa #Q_mu #anisotropy_flag
        imaterial=properties[0]
        flag=properties[1]
        #print 'prop',flag
        if flag > 0:
            vel=properties[2]
            if properties[2] is None and type(vel) != str:
                if vel >= 30:
                    m2km=1000.
                else:
                    m2km=1.
                vp=vel/m2km
                rho=(1.6612*vp-0.472*vp**2+0.0671*vp**3-0.0043*vp**4+0.000106*vp**4)*m2km
                txt='%1i %3i %20f %20f %20f %1i %1i\n' % (properties[0],properties[1],rho,vel,vel/(3**.5),0,0)     
            elif type(vel) != str and vel != 0.:
                try: 
                    qmu=properties[6]
                except:
                    qmu=9999.
                try: 
                    qk=properties[5]
                except:
                    qk=9999.
                try:
                    ani=properties[7]
                except:
                    ani=0.
                #print properties[0],properties[3],properties[1],properties[2],q,ani
                txt='%1i %3i %20f %20f %20f %20f %20f\n' % (properties[0],properties[1],properties[4],properties[2],properties[3],qk,qmu,ani)
            elif type(vel) != str and vel != 0.:
                helpstring="#material_domain_id #material_id #rho #vp #vs #Q_kappa #Q_mu #anisotropy"
                txt='%1i %3i %s \n' % (properties[0],properties[1],helpstring)
            else:
                helpstring=" -->       sintax: #material_domain_id #material_id #rho #vp #vs #Q_kappa #Q_mu #anisotropy"
                txt='%1i %3i %s %s\n' % (properties[0],properties[1],properties[2],helpstring)
        elif flag < 0:
            if properties[2] == 'tomography':
                txt='%1i %3i %s %s\n' % (properties[0],properties[1],properties[2],properties[3])
            elif properties[2] == 'interface':
                txt='%1i %3i %s %s %1i %1i\n' % (properties[0],properties[1],properties[2],properties[3],properties[4],properties[5])
            else:
                helpstring=" -->       sintax: #material_domain_id 'tomography' #file_name "
                txt='%1i %3i %s %s \n' % (properties[0],properties[1],properties[2],helpstring)
                #
        #print txt
        return txt
    def nummaterial_write(self,nummaterial_name,placeholder=True):
        print 'Writing '+nummaterial_name+'.....'
        nummaterial=open(nummaterial_name,'w')
        for block in self.block_mat:
            #name=cubit.get_exodus_entity_name('block',block)
            nummaterial.write(str(self.mat_parameter(self.material[block])))
        if placeholder:
            txt='''



! note: format of nummaterial_velocity_file must be


! #(1)material_domain_id #(2)material_id  #(3)rho  #(4)vp   #(5)vs   #(6)Q_kappa   #(7)Q_mu  #(8)anisotropy_flag
!
! where
!     material_domain_id : 1=acoustic / 2=elastic 
!     material_id        : POSITIVE integer identifier corresponding to the identifier of material block
!     rho                : density
!     vp                 : P-velocity
!     vs                 : S-velocity
!     Q_kappa            : 9999 = no Q_kappa attenuation
!     Q_mu               : 9999 = no Q_mu attenuation
!     anisotropy_flag    : 0=no anisotropy/ 1,2,.. check with implementation in aniso_model.f90            
!
!example:
!2   1 2300 2800 1500 9999.0 9999.0 0

!or

! #(1)material_domain_id #(2)material_id  tomography elastic  #(3)tomography_filename #(4)positive_unique_number
!         
! where
!     material_domain_id : 1=acoustic / 2=elastic 
!     material_id        : NEGATIVE integer identifier corresponding to the identifier of material block
!     tomography_filename: filename of the tomography file       
!     positive_unique_number: a positive unique identifier   
!
!example:
!2  -1 tomography elastic tomo.xyz 1


'''
            nummaterial.write(txt)
        nummaterial.close()
    
    def create_hexnode_string(self,hexa):
        nodes=self.get_hex_connectivity(hexa)
        #nodes=self.jac_check(nodes) #is it valid for 3D? TODO
        if self.hex27:
            ordered_nodes=[hexa]+list(nodes[:20])+[nodes[21]]+[nodes[25]]+[nodes[24]]+[nodes[26]]+[nodes[23]]+[nodes[22]]+[nodes[20]]
            txt=' '.join(str(x) for x in ordered_nodes)
            txt=txt+'\n'
            #txt=('%10i %10i %10i %10i %10i %10i %10i %10i ')% nodes[:8] #first 8 nodes following specfem3d numbering convenction..
            #txt=txt+('%10i %10i %10i %10i %10i %10i %10i %10i ')% nodes[8:16] #middle 12 nodes following specfem3d numbering convenction..
            #txt=txt+('%10i %10i %10i %10i ')% nodes[16:20]
            #txt=txt+('%10i %10i %10i %10i %10i %10i ')% (nodes[21], nodes[25], nodes[24], nodes[26], nodes[23], nodes[22])
            #txt=txt+('%10i\n ')% nodes[20] #center volume
        else:
            txt=str(hexa)+' '+' '.join(str(x) for x in nodes)
            txt=txt+'\n'
            #txt=('%10i %10i %10i %10i %10i %10i %10i %10i\n')% nodes[:]
        return txt
        
    def create_facenode_string(self,hexa,face,normal=None,cknormal=True):
        nodes=self.get_face_connectivity(face)
        if cknormal:
            nodes_ok=self.normal_check(nodes[0:4],normal)
            if self.hex27: nodes_ok2=self.normal_check(nodes[4:8],normal)
        else:
            nodes_ok=nodes[0:4]
            if self.hex27: nodes_ok2=nodes[4:8]
        #
        if self.hex27:
            ordered_nodes=[hexa]+list(nodes_ok)+list(nodes_ok2)+[nodes[8]]
            txt=' '.join(str(x) for x in ordered_nodes)
            txt=txt+'\n'
            #txt=('%10i %10i %10i %10i %10i ') % (hexa,nodes_ok[0],nodes_ok[1],nodes_ok[2],nodes_ok[3]) #first 4 nodes following specfem3d numbering convenction..
            #txt=txt+('%10i %10i %10i %10i ')% (nodes_ok2[0],nodes_ok2[1],nodes_ok2[2],nodes_ok2[3]) #middle 4 nodes following specfem3d numbering convenction..
            #txt=txt+('%10i\n')% nodes[8]
        else:
            txt=str(hexa)+' '+' '.join(str(x) for x in nodes_ok)
            txt=txt+'\n'
            #txt=('%10i %10i %10i %10i %10i\n') % (hexa,nodes_ok[0],nodes_ok[1],nodes_ok[2],nodes_ok[3])
        return txt
    
    
    def mesh_write(self,mesh_name):
        meshfile=open(mesh_name,'w')
        print 'Writing '+mesh_name+'.....'
        num_elems=cubit.get_hex_count()
        print ' total number of elements:',str(num_elems)
        meshfile.write(str(num_elems)+'\n')
        for block,flag in zip(self.block_mat,self.block_flag):
            hexes=cubit.get_block_hexes(block)
            print 'block ',block,' hexes ',len(hexes)
            for hexa in hexes:
                txt=self.create_hexnode_string(hexa)
                meshfile.write(txt)
        meshfile.close()
    def material_write(self,mat_name):
        mat=open(mat_name,'w')
        print 'Writing '+mat_name+'.....'
        for block,flag in zip(self.block_mat,self.block_flag):
                print 'block ',block,'flag ',flag
                hexes=cubit.get_block_hexes(block)
                for hexa in hexes:
                    mat.write(('%10i %10i\n') % (hexa,flag))
        mat.close()
    def get_extreme(self,c,cmin,cmax):
        if not cmin and not cmax:
            cmin=c
            cmax=c
        else:
            if c<cmin: cmin=c
            if c>cmax: cmax=c
        return cmin,cmax
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
            self.xmin,self.xmax=self.get_extreme(x,self.xmin,self.xmax)
            self.ymin,self.ymax=self.get_extreme(y,self.ymin,self.ymax)
            self.zmin,self.zmax=self.get_extreme(z,self.zmin,self.zmax)
            txt=('%10i %20f %20f %20f\n') % (node,x,y,z)
            nodecoord.write(txt)
        nodecoord.close()
    def free_write(self,freename=None):
        cubit.cmd('set info off')
        cubit.cmd('set echo off')
        cubit.cmd('set journal off')
        from sets import Set
        normal=(0,0,1)
        if not freename: freename=self.freename
        freehex=open(freename,'w')
        print 'Writing '+freename+'.....'
        #
        #
        for block,flag in zip(self.block_bc,self.block_bc_flag):
            if block == self.topography:
                name=cubit.get_exodus_entity_name('block',block)
                print 'free surface (topography) block name:',name,'id:',block
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
                            txt=self.create_facenode_string(h,f,normal,cknormal=True)
                            freehex.write(txt)
                freehex.close()
            elif block == self.free: 
                name=cubit.get_exodus_entity_name('block',block)
                print 'free surface block name:',name,'id:',block
                quads_all=cubit.get_block_faces(block)
                print '  number of faces = ',len(quads_all)
                dic_quads_all=dict(zip(quads_all,quads_all))
                freehex.write('%10i\n' % len(quads_all))
                list_hex=cubit.parse_cubit_list('hex','all')
                for h in list_hex:
                    faces=cubit.get_sub_elements('hex',h,2)
                    for f in faces:
                        if dic_quads_all.has_key(f):
                            txt=self.create_facenode_string(h,f,normal,cknormal=False)
                            freehex.write(txt)
                freehex.close()
        cubit.cmd('set info on')
        cubit.cmd('set echo on')
    def check_cmpl_size(self,case='x'):
        if case=='x':
            vmaxtmp=self.xmax
            vmintmp=self.xmin
        elif case=='y':
            vmaxtmp=self.ymax
            vmintmp=self.ymin
        elif case=='z':
            vmaxtmp=self.zmax
            vmintmp=self.zmin
            
        if self.size > .3*(vmaxtmp-vmintmp):
            print 'please select the size of cpml less than 30% of the '+case+' size of the volume'
            print vmaxtmp-vmintmp,.3*(vmaxtmp-vmintmp)
            print 'cmpl set to false, no '+self.cpmlname+' file will be created'
            return False,False
        else:
            vmin=vmintmp+self.size
            vmax=vmaxtmp-self.size
        return vmin,vmax
    def select_cpml(self):
        xmin,xmax=self.check_cmpl_size(case='x')
        ymin,ymax=self.check_cmpl_size(case='y')
        zmin,zmax=self.check_cmpl_size(case='z')
        #
        if xmin is False or xmax is False or ymin is False or ymax is False or zmin is False or zmax is False:
            return False
        else:
            txt="group 'hxmin' add hex  with X_coord < "+str(xmin)
            cubit.cmd(txt)        
            txt="group 'hxmax' add hex  with X_coord > "+str(xmax)
            cubit.cmd(txt)        
            txt="group 'hymin' add hex  with Y_coord < "+str(ymin)
            cubit.cmd(txt)         
            txt="group 'hymax' add hex  with Y_coord > "+str(ymax)
            cubit.cmd(txt)        
            txt="group 'hzmin' add hex  with Z_coord < "+str(zmin)
            cubit.cmd(txt)       
            txt="group 'hzmax' add hex  with Z_coord > "+str(zmax)
            cubit.cmd(txt)
            from sets import Set
            group1 = cubit.get_id_from_name("hxmin")
            cpml_xmin =Set(list(cubit.get_group_hexes(group1)))
            group1 = cubit.get_id_from_name("hymin")
            cpml_ymin =Set(list(cubit.get_group_hexes(group1)))
            group1 = cubit.get_id_from_name("hxmax")
            cpml_xmax =Set(list(cubit.get_group_hexes(group1)))
            group1 = cubit.get_id_from_name("hymax")
            cpml_ymax =Set(list(cubit.get_group_hexes(group1)))
            group1 = cubit.get_id_from_name("hzmin")
            cpml_zmin =Set(list(cubit.get_group_hexes(group1)))
            if self.top_absorbing:
                group1 = cubit.get_id_from_name("hzmax")
                cpml_zmax =Set(list(cubit.get_group_hexes(group1)))
            else:
                cpml_zmax =Set([])
            cpml_all=cpml_ymin | cpml_ymax | cpml_xmin | cpml_xmax | cpml_zmin | cpml_zmax
            cpml_x=cpml_all-cpml_zmin-cpml_ymin-cpml_ymax-cpml_zmax
            cpml_y=cpml_all-cpml_zmin-cpml_xmin-cpml_xmax-cpml_zmax
            cpml_xy=cpml_all-cpml_zmin-cpml_y-cpml_x-cpml_zmax
            cpml_z=cpml_all-cpml_xmin-cpml_ymin-cpml_ymax-cpml_xmax
            cpml_xz=cpml_zmin-cpml_ymin-cpml_ymax-cpml_z
            cpml_yz=cpml_zmin-cpml_xmin-cpml_xmax-cpml_z
            cpml_xyz=cpml_zmin-cpml_xz-cpml_yz-cpml_z
            txt=' '.join(str(h) for h in cpml_x)
            cubit.cmd("group 'x_cpml' add hex "+txt)
            txt=' '.join(str(h) for h in cpml_y)
            cubit.cmd("group 'y_cpml' add hex "+txt)
            txt=' '.join(str(h) for h in cpml_z)
            cubit.cmd("group 'z_cpml' add hex "+txt)
            txt=' '.join(str(h) for h in cpml_xy)
            cubit.cmd("group 'xy_cpml' add hex "+txt)
            txt=' '.join(str(h) for h in cpml_xz)
            cubit.cmd("group 'xz_cpml' add hex "+txt)
            txt=' '.join(str(h) for h in cpml_yz)
            cubit.cmd("group 'yz_cpml' add hex "+txt)            
            txt=' '.join(str(h) for h in cpml_xyz)
            cubit.cmd("group 'xyz_cpml' add hex "+txt)
            return cpml_x,cpml_y,cpml_z,cpml_xy,cpml_xz,cpml_yz,cpml_xyz
        
        
        
        
        
        
        
    
    
    def abs_write(self,absname=None):
        import re
        cubit.cmd('set info off')
        cubit.cmd('set echo off')
        cubit.cmd('set journal off')
        from sets import Set
        if not absname: absname=self.absname
        
        if self.cpml:
            if not absname: absname=self.cpmlname
            print 'Writing cpml'+absname+'.....'
            list_cpml=self.select_cpml()
            if list_cpml is False:
                print 'error writing cpml files'
                return
            else:
                abshex_cpml=open(absname,'w')
                hexcount=sum(map(len,list_cpml))
                abshex_cpml.write(('%10i\n') % (hexcount))
                for icpml,lcpml in enumerate(list_cpml):
                    for hexa in lcpml:
                        abshex_cpml.write(('%10i %10i\n') % (hexa,icpml))
            
            
        stacey_absorb=True
        if stacey_absorb:      
            #
            #
            if not absname: absname=self.absname
            list_hex=cubit.parse_cubit_list('hex','all')
            for block,flag in zip(self.block_bc,self.block_bc_flag):
                if block != self.topography:
                    name=cubit.get_exodus_entity_name('block',block)
                    print '  block name:',name,'id:',block
                    cknormal=True
                    abshex_local=False
                    if re.search('xmin',name):
                        print 'xmin'
                        abshex_local=open(absname+'_xmin','w')
                        normal=(-1,0,0)
                    elif re.search('xmax',name):
                        print "xmax"
                        abshex_local=open(absname+'_xmax','w')
                        normal=(1,0,0)
                    elif re.search('ymin',name):
                        print "ymin"
                        abshex_local=open(absname+'_ymin','w')
                        normal=(0,-1,0)
                    elif re.search('ymax',name):
                        print "ymax"
                        abshex_local=open(absname+'_ymax','w')
                        normal=(0,1,0)
                    elif re.search('bottom',name):
                        print "bottom"
                        abshex_local=open(absname+'_bottom','w')
                        normal=(0,0,-1)
                    elif re.search('abs',name):
                        print "abs all - experimental, check the output"
                        cknormal=False
                        abshex_local=open(absname,'w')
                    else:
                        if block == 1003:
                            print 'xmin'
                            abshex_local=open(absname+'_xmin','w')
                            normal=(-1,0,0)
                        elif block == 1004:
                            print "ymin"
                            abshex_local=open(absname+'_ymin','w')
                            normal=(0,-1,0)
                        elif block == 1005:
                            print "xmax"
                            abshex_local=open(absname+'_xmax','w')
                            normal=(1,0,0)
                        elif block == 1006:
                            print "ymax"
                            abshex_local=open(absname+'_ymax','w')
                            normal=(0,1,0)
                        elif block == 1002:
                            print "bottom"
                            abshex_local=open(absname+'_bottom','w')
                            normal=(0,0,-1)
                        elif block == 1000:
                            print "custumized"
                            abshex_local=open(absname,'w')
                            cknormal=False
                            normal=None
                    #
                    #
                    if abshex_local:
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
                                    txt=self.create_facenode_string(h,f,normal=normal,cknormal=cknormal)
                                    abshex_local.write(txt)
                        abshex_local.close()   
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
                            txt=self.create_facenode_string(h,f,cknormal=False)
                            surfhex_local.write(txt)
                # closes file
                surfhex_local.close()
    def rec_write(self,recname):
        print 'Writing '+self.recname+'.....'
        recfile=open(self.recname,'w')
        nodes=cubit.get_nodeset_nodes(self.receivers)
        for i,n in enumerate(nodes):
            x,y,z=cubit.get_nodal_coordinates(n)
            recfile.write('ST%i XX %20f %20f 0.0 0.0 \n' % (i,x,z))
        recfile.close()
    def write(self,path=''):
        cubit.cmd('set info off')
        cubit.cmd('set echo off')
        cubit.cmd('set journal off')
        cubit.cmd('compress all')
        if len(path) != 0:
            if path[-1] != '/': path=path+'/'
        self.mesh_write(path+self.mesh_name)
        self.material_write(path+self.material_name)
        self.nodescoord_write(path+self.nodecoord_name)
        self.free_write(path+self.freename)
        if self.cpml:
            self.abs_write(path+self.cpmlname)
        else:
            self.abs_write(path+self.absname)
        self.nummaterial_write(path+self.nummaterial_name)
        # any other surfaces: ***surface***
        self.surface_write(path)
        if self.receivers: self.rec_write(path+self.recname)
        cubit.cmd('set info on')
        cubit.cmd('set echo on')

def export2SPECFEM3D(path_exporting_mesh_SPECFEM3D='.',hex27=False,cpml=False,cpml_size=False,top_absorbing=False):
    sem_mesh=mesh(hex27,cpml,cpml_size,top_absorbing)
    #sem_mesh.block_definition()
    #print sem_mesh.block_mat
    #print sem_mesh.block_flag
    #
    sem_mesh.write(path=path_exporting_mesh_SPECFEM3D)
    print 'END SPECFEM3D exporting process......'
    if cpml:
        cmd='save as "cpml.cub" overwrite'
        cubit.cmd(cmd)



