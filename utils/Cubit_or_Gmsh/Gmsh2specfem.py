#!/usr/bin/env python
#
# uses meshio:
#   https://github.com/nschloe/meshio
#
# install with:
#   pip install -U meshio
#
# this script reads in a Gmsh file (in Gmsh-ASCII format) and outputs SPECFEM-readable mesh files
#
# Boundary surfaces must be defined prior when generating the mesh as Physical Surfaces in Gmsh.
# Surfaces must have names:
#    - "top"
#    - "bottom"
#    - "xmin"
#    - "xmax"
#    - "ymin"
#    - "ymax"
#
from __future__ import print_function

import sys
import os
import math

import numpy as np

try:
    import meshio
except:
    print("importing module meshio failed, please make sure to install it via: pip install -U meshio")
    sys.exit(1)

# globals
cell_type = ''
quad_type = ''
meshio_major_version = 0

#--------------------------------------------------------------------------------------------------

def read_mesh_file_msh(file):
    """
    reads in Gmsh file
    """
    global cell_type,quad_type,meshio_major_version

    print("")
    print("reading file: ",file)
    print("")

    # reads *.msh file (Gmsh Ascii-format)
    mesh = meshio.read(file)

    # version changes:
    #   meshio 3.3.1: cells and cell_data are returned as dictionary to mesh object
    #   meshio >= 4.0.0: cells and cell_data are returns as lists to mesh object,
    #                    instead use cells_dict and cell_data_dict
    #
    # we would have something like:
    #
    # if meshio_major_version >= 4:
    #     points, cells, point_data, cell_data, field_data = (mesh.points, mesh.cells_dict,
    #                                                         mesh.point_data, mesh.cell_data_dict,
    #                                                         mesh.field_data)
    # else:
    #     points, cells, point_data, cell_data, field_data = (mesh.points, mesh.cells,
    #                                                         mesh.point_data, mesh.cell_data,
    #                                                         mesh.field_data)
    #
    # however, the new dict objects have a different ordering than in older version 3.
    # we will thus take the list object and convert it to a version 3 dictionary object ourself...
    #
    points, cells, point_data, cell_data, field_data = (mesh.points, mesh.cells,
                                                        mesh.point_data, mesh.cell_data,
                                                        mesh.field_data)

    # meshio versions >= 4.0.0:
    if meshio_major_version >= 4:
        print("current meshio version >= 4")
        # creates a cells_old object as dictionaries in the same format like older version meshio 3.3.x
        if not isinstance(cells, dict):
            print("")
            print("using older meshio cell format for exporting mesh...")

            # creates new dictionary for cells and cell_data
            cells_old = {}
            nblocks = 0
            for cell_block in cells:
                # cells dictionary
                nblocks += 1
                #print("cell block",nblocks,cell_block,len(cell_block))
                if len(cell_block) == 2:
                    # must have 2 items ('name',data)
                    name = cell_block[0]
                    data = cell_block[1]
                    print("  cell block",name,len(data))
                    #for i in cell_block: print(i)
                    if name in cells_old:
                        # append data to existing key ('quad9' can occur multiple times)
                        cells_old[name] = np.concatenate((cells_old[name],data))
                    else:
                        # new entry
                        cells_old[name] = data
                else:
                    # error
                    print("Invalid cells returned",cells,cell_block)
                    sys.exit(1)
            print("")

            # replace cells object with older-format dictionary
            cells.clear()
            cells = cells_old.copy()

            #debug
            #print("cells",cells)

            # creates new dictionary for cell_data
            # format: quad       -> (gmsh:physical, 600)
            #         hexahedron -> (gmsh:physical, 1000)
            cell_data_old = {}
            for name,data in cell_data.items():
                # newer format has entries:
                #     gmsh:geometrical,
                #     gmsh:physical
                #debug
                #print("  name,data",name,data)
                # only gmsh:physical entries
                if 'gmsh:physical' in name:
                    # physical mesh infos
                    # data list should have same length as cells list entry
                    nkeys = len(cells.keys())
                    if len(data) != nkeys:
                        print("Error: number of cell_data entries and cells type differ ",len(data),nkeys)
                        sys.exit(1)

                    # create new dictionary
                    for s,array in enumerate(data):
                        # gets name from cells entry
                        data_name = list(cells.keys())[s]
                        print("  cell data block",s,len(array),"name:",data_name)
                        # dictionary has key gmsh:physical (ignoring gmsh:geometry)
                        dic = {'gmsh:physical': array}
                        # sets as new entry
                        cell_data_old[data_name] = dic
                        # debug
                        #print("  data name,dic",data_name,dic)
            print("")
            #debug
            #print("cell_data",cell_data)

            # replace cells object with older-format dictionary
            cell_data.clear()
            cell_data = cell_data_old.copy()

            # debug
            #print("new cell_data",cell_data)

    # mesh info
    print("mesh data:")
    print("number of points: ",len(points))
    print("cells     : ",len(cells),"items")
    for name,_ in cells.items():
        print("  ",name)

    print("point_data: ",len(point_data),"items")
    for name,dic in point_data.items():
        print("  ",name)
        for s,array in dic.items():
            print("    ",s,len(array))

    print("cell_data : ",len(cell_data),"items")
    for name,data in cell_data.items():
        print("  ",name)
        if isinstance(data, dict):
            # older meshio versions 3.x
            for s,array in data.items():
                print("    ",s,len(array))
        else:
            # newer meshio versions >= 4.0 (data is a list object)
            for s,array in enumerate(data):
                print("          ",s,len(array))

    print("field_data: ",len(field_data),"items")
    for name,_ in field_data.items():
        print("  ",name)

    print("")

    #debug
    #print("cells:")
    #print(cells)
    #print("cell_data:")
    #print(cell_data)
    #print("field_data:")
    #print(field_data)
    #print("")

    # checks types
    # linear    hex9  elements == 'hexahedron' (used by default)
    # quadratic hex27 elements == 'hexahedron27'
    if 'hexahedron27' in cells.keys():
        cell_type = 'hexahedron27'
    elif 'hexahedron' in cells.keys():
        cell_type = 'hexahedron'

    # linear    quad4  faces == 'quad' (used by default)
    # quadratic quad9  faces == 'quad9'
    if 'quad9' in cells.keys():
        quad_type = 'quad9'
    elif 'quad' in cells.keys():
        quad_type = 'quad'

    print("element types found: ")
    print("  hex type : ",cell_type)
    print("  quad type: ",quad_type)
    print("")
    if cell_type == '' or quad_type == '':
        print("hex and quad types not found correctly, please check mesh...")
        sys.exit(1)

    # cleans points array (removing unused points)
    points = clean_out_unused_points(points,cells)

    return points,cells,point_data,cell_data,field_data

#--------------------------------------------------------------------------------------------------


def clean_out_unused_points(points,cells):
    # for some reason, there might be extra nodes stored. we will clean them out and only use the onces needed by the hexas.
    # assumes that the surface quads won't require new nodes.
    global cell_type,quad_type

    print("cleaning points...")

    # checks if anything to do
    if not cell_type in cells.keys():
        return points

    # get array with all point entries
    cell_all = cells[cell_type].flatten()

    # remove redundant entries
    cell_points = (list(set(cell_all)))

    # point indices [0,1,2,..]
    points_index = range(len(points))

    # array with points not listed by cells (difference between the two sets)
    unused_points = list(set(points_index) - set(cell_points))

    print("  number of points                     = ",len(points))
    print("  number of points needed by hexahedra = ",len(cell_points))
    print("  unused points: ",unused_points)

    # checks if anything to do
    if len(unused_points) > 0:
        # remove unused points from lists
        points = np.delete(points, unused_points, axis=0)

        print("  new number of points                 = ",len(points))

        # not used: updates point data
        #for key in point_data:
        #    point_data[key] = numpy.delete(point_data[key],unused_points,axis=0)

    return points



#--------------------------------------------------------------------------------------------------


def get_surface_dimension(surface,points):
    # determines min/max dimensions of all points on a surface
    surf_points_x = []
    surf_points_y = []
    surf_points_z = []
    for elem in surface:
        #print(elem)
        for i in range(1,5):
            index = elem[i]
            if index < 0 or index >= len(points):
                print("invalid point index ",index,"in surface")
                sys.exit(1)
            #print(elem,index,points[index])
            surf_points_x.append(points[index][0])
            surf_points_y.append(points[index][1])
            surf_points_z.append(points[index][2])

    xmin = min(surf_points_x)
    xmax = max(surf_points_x)
    ymin = min(surf_points_y)
    ymax = max(surf_points_y)
    zmin = min(surf_points_z)
    zmax = max(surf_points_z)

    #zmax = max(p[2] for p in surf_points)

    #debug
    #print(surface)
    #print(surf_points_x)
    #print(surf_points_y)
    #print(surf_points_z)

    return xmin,xmax,ymin,ymax,zmin,zmax

#--------------------------------------------------------------------------------------------------

def detect_boundary_surfaces(points,cells,cell_data,field_data):
    """
    determines surfaces at xmin,xmax,ymin,ymax,top and bottom
    """
    global cell_type,quad_type

    print("determining boundary surfaces:")

    surf_top = []
    surf_bottom = []
    surf_xmin = []
    surf_xmax = []
    surf_ymin = []
    surf_ymax = []

    if not quad_type in cells.keys():
        print("invalid mesh: surfaces need quad cells")
        sys.exit(1)
    if not cell_type in cells.keys():
        print("invalid mesh: surfaces need hexahedron cells")
        sys.exit(1)

    quads = cells[quad_type]
    hexas = cells[cell_type]

    # checks if valid entries
    if len(quads) == 0:
        print("Error: need entries in '{}' cells".format(quad_type))
        sys.exit(1)
    if len(hexas) == 0:
        print("Error: need entries in '{}' cells".format(cell_type))
        sys.exit(1)

    # gets surface index
    index = [0,0,0,0,0,0]  # format top,bottom,xmin,xmax,ymin,ymax
    for key in field_data.keys():
        # top surface
        if key == 'top':
            # for example item: u'top': array([2, 2]) # first index is surface id
            index[0] = field_data[key][0]
        # bottom
        if key == 'bottom':
            index[1] = field_data[key][0]
        # xmin
        if key == 'xmin':
            index[2] = field_data[key][0]
        # xmax
        if key == 'xmax':
            index[3] = field_data[key][0]
        # ymin
        if key == 'ymin':
            index[4] = field_data[key][0]
        # ymax
        if key == 'ymax':
            index[5] = field_data[key][0]
    print("  surface index: ",index)
    print("")

    # gets cell data
    if not quad_type in cell_data.keys():
        print("invalid mesh: needs '{}' entry in cell_data".format(quad_type))
        sys.exit(1)
    quad_dict = cell_data[quad_type]

    # meshio creates 'gmsh:physical' and 'gmsh:geometrical' items, using physical volume associations
    if not 'gmsh:physical' in quad_dict.keys():
        print("invalid mesh: needs 'gmsh:physical' entry in cell_data for '{}'".format(quad_type))
        sys.exit(1)
    quad_mat = quad_dict['gmsh:physical']
    print("  number of quads with material = ",len(quad_mat))
    print("")

    if len(quad_mat) == 0:
        print("Error: need quads with material for boundary surfaces, please define physical surfaces using Gmsh")
        sys.exit(1)
    if len(quad_mat) != len(quads):
        print("invalid total number of quads ",len(quads),"and number of quads with cell data ",len(quad_mat))
        sys.exit(1)

    index_count = [0,0,0,0,0,0]
    for i,q in enumerate(quad_mat):
        ind = quad_mat[i]
        if ind == index[0]: index_count[0] += 1
        if ind == index[1]: index_count[1] += 1
        if ind == index[2]: index_count[2] += 1
        if ind == index[3]: index_count[3] += 1
        if ind == index[4]: index_count[4] += 1
        if ind == index[5]: index_count[5] += 1
    print("  surface index counts: ",index_count)
    print("")

    # the following assumes linear/quadratic elements, i.e., quads defined by 4 or 9 corners only
    if len(quads[0]) != 4 and len(quads[0]) != 9:
        print("Invalid element type for quads {}: need linear/quadratic elements, higher order elements not supported yet".format(len(quads[0])))
        sys.exit(1)

    # loops through quads
    for i,q in enumerate(quads):
        #print(i,q)   # i starts at 0,1,2,..; quads[id1,id2,id3,id4] indices also start at 0

        # for specfem, we will need to find corresponding hexahedral element which contains this quad
        # DK DK significantly faster way suggested by Deyu Ming
        #el = [ [j,h] for j,h in enumerate(hexas) if q[0] in h and q[1] in h and q[2] in h and q[3] in h ]
        if quad_type == 'quad9':
            el = np.where(np.count_nonzero(np.isin(hexas,q),axis=1) == 9)
        else:
            el = np.where(np.count_nonzero(np.isin(hexas,q),axis=1) == 4)

        # boundary surfaces should only associate with a single hexahedral element
        if len(el) != 1:
            print("Error: could not find hexahedral for surface quad ",i,q,"el = ",el)
            sys.exit(1)
        hex_index = el[0][0]

        # adds surface id to quad corners
        # see Gmsh ordering: http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
        if quad_type == 'quad9':
            elem = [hex_index,q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8]]
        else:
            elem = [hex_index,q[0],q[1],q[2],q[3]]

        # gets associated surface type index
        ind = quad_mat[i]

        # appends to surface
        if ind == index[0]: surf_top.append(elem)
        if ind == index[1]: surf_bottom.append(elem)
        if ind == index[2]: surf_xmin.append(elem)
        if ind == index[3]: surf_xmax.append(elem)
        if ind == index[4]: surf_ymin.append(elem)
        if ind == index[5]: surf_ymax.append(elem)

    print("  top    surface: number of quads = ",len(surf_top))
    if len(surf_top) > 0:
        xmin,xmax,ymin,ymax,zmin,zmax = get_surface_dimension(surf_top,points)
        print("    dimension:")
        print("    xmin/xmax = ",xmin,"/",xmax)
        print("    ymin/ymax = ",ymin,"/",ymax)
        print("    zmin/zmax = ",zmin,"/",zmax)

    print("  bottom surface: number of quads = ",len(surf_bottom))
    if len(surf_bottom) > 0:
        xmin,xmax,ymin,ymax,zmin,zmax = get_surface_dimension(surf_bottom,points)
        print("    dimension:")
        print("    xmin/xmax = ",xmin,"/",xmax)
        print("    ymin/ymax = ",ymin,"/",ymax)
        print("    zmin/zmax = ",zmin,"/",zmax)

    print("  xmin   surface: number of quads = ",len(surf_xmin))
    if len(surf_xmin) > 0:
        xmin,xmax,ymin,ymax,zmin,zmax = get_surface_dimension(surf_xmin,points)
        print("    dimension:")
        print("    xmin/xmax = ",xmin,"/",xmax)
        print("    ymin/ymax = ",ymin,"/",ymax)
        print("    zmin/zmax = ",zmin,"/",zmax)

    print("  xmax   surface: number of quads = ",len(surf_xmax))
    if len(surf_xmax) > 0:
        xmin,xmax,ymin,ymax,zmin,zmax = get_surface_dimension(surf_xmax,points)
        print("    dimension:")
        print("    xmin/xmax = ",xmin,"/",xmax)
        print("    ymin/ymax = ",ymin,"/",ymax)
        print("    zmin/zmax = ",zmin,"/",zmax)

    print("  ymin   surface: number of quads = ",len(surf_ymin))
    if len(surf_ymin) > 0:
        xmin,xmax,ymin,ymax,zmin,zmax = get_surface_dimension(surf_ymin,points)
        print("    dimension:")
        print("    xmin/xmax = ",xmin,"/",xmax)
        print("    ymin/ymax = ",ymin,"/",ymax)
        print("    zmin/zmax = ",zmin,"/",zmax)

    print("  ymax   surface: number of quads = ",len(surf_ymax))
    if len(surf_ymax) > 0:
        xmin,xmax,ymin,ymax,zmin,zmax = get_surface_dimension(surf_ymax,points)
        print("    dimension:")
        print("    xmin/xmax = ",xmin,"/",xmax)
        print("    ymin/ymax = ",ymin,"/",ymax)
        print("    zmin/zmax = ",zmin,"/",zmax)

    print("")

    return surf_top,surf_bottom,surf_xmin,surf_xmax,surf_ymin,surf_ymax

#--------------------------------------------------------------------------------------------------

def export_mesh(points,cells,point_data,cell_data,field_data):
    """
    creates files in MESH/ directory
    """
    global cell_type,quad_type,meshio_major_version

    print("")
    print("exporting mesh files:")

    # outputs mesh files
    os.system("mkdir -p MESH/")

    # node coordinates
    print("  number of nodes       = ",len(points))
    filename = "MESH/nodes_coords_file"
    with open(filename,"w") as f:
        f.write("%d\n" % len(points))
        for i,p in enumerate(points):
            # format: #id_node #x_coordinate #y_coordinate #z_coordinate
            f.write("%d %f %f %f\n" % (i+1,p[0],p[1],p[2]) )
    print("  exported to: ",filename)
    print("")

    # cells
    if not cell_type in cells.keys():
        print("invalid mesh: needs hexahedron cells")
        sys.exit(1)
    hexas = cells[cell_type]
    print("  number of hexahedra   = ",len(hexas))
    filename = "MESH/mesh_file"
    with open(filename,"w") as f:
        f.write("%d\n" % len(hexas))
        for i,h in enumerate(hexas):
            # format: # element_id  #id_node1 ... #id_node8
            # note: specfem assumes indices start at 1 (fortran style, not at 0 like in python or C)
            #
            # for ordering:
            # see Gmsh ordering: http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
            # see VTK  ordering: https://lorensen.github.io/VTKExamples/site/VTKBook/05Chapter5/
            #
            # note: ordering between Gmsh, VTK and SPECFEM are all different.
            #       the module meshio changes the ordering for hexahedron20 (to fit "VTK" format), but not hexahedron27 elements.
            #       thus, the ordering is still original from Gmsh for hex8 and hex27 elements.
            #
            # note: with meshio versions 4.x, ordering changed again...
            #
            if cell_type == 'hexahedron27':
                # 27-node element -> NGNOD == 27
                if len(h) != 27:
                    print("Invalid hex length {} for cell type {}, exiting...".format(len(h),cell_type))
                    sys.exit(1)
                # orders from Gmsh to VTK:
                # ordering corners indexing from (bottom) to (top): (5 4 7 6) to (1 0 3 2)
                #    mid-points (bottom edges), (top edges), (bottom-to-top edges): (16 17 19 18) (8 9 13 11) (12 10 15 14)
                #    face centers (left),(right),(front),(back),(bottom),(top): (23,22,21,24,25,20)
                #    cell center: (26)
                #    -> ordering = [5,4,7,6, 1,0,3,2, 16,17,19,18, 8,9,13,11, 12,10,15,14, 23,22,21,24,25,20, 26]
                #
                # orders from Gmsh to SPECFEM:
                # ordering corners indexing from (bottom) to (top): (5 4 7 6) to (1 0 3 2)
                #    mid-points (bottom edges),(bottom-to-top edges), (top edges), : (16 17 19 18) (12 10 15 14) (8 9 13 11)
                #    face centers (bottom),(front),(right),(back),(left),(top): 25,21,22,24,23,20
                #    cell center: (26)
                #    -> ordering = [5,4,7,6, 1,0,3,2, 16,17,19,18, 12,10,15,14, 8,9,13,11, 25,21,22,24,23,20, 26]
                #
                # newer meshio ordering, see:
                # https://github.com/nschloe/meshio/wiki/Node-ordering-in-cells
                #
                # indexing from bottom surface to top
                if meshio_major_version >= 4:
                    # orders from meshio to specfem
                    ordering = [5,4,7,6, 1,0,3,2, 12,15,14,13, 17,16,19,18, 8,11,10,9, 25,22,20,23,21,24, 26]
                else:
                    # orders from gmsh to specfem
                    ordering = [5,4,7,6, 1,0,3,2, 16,17,19,18, 12,10,15,14, 8,9,13,11, 25,21,22,24,23,20, 26]
            else:
                # 8-node element given by corner points -> NGNOD == 8
                if len(h) != 8:
                    print("Invalid hex length {} for cell type {}, exiting...".format(len(h),cell_type))
                    sys.exit(1)
                # orders from Gmsh to VTK/SPECFEM:
                # ordering corner indexing from (bottom) to (top): (5 4 7 6) to (1 0 3 2)
                # this stays the same for SPECFEM ordering.
                ordering = [5,4,7,6, 1,0,3,2]

            # creates entry string
            str = "{} ".format(i+1)  # element id
            for j in ordering: str += " {}".format(h[j]+1)  # nodes
            str += "\n"
            f.write(str)

            #f.write("%d %d %d %d %d %d %d %d %d\n" % (i+1,h[4]+1,h[5]+1,h[1]+1,h[0]+1,h[7]+1,h[6]+1,h[2]+1,h[3]+1) )

    print("  exported to: ",filename)
    print("")

    # material (or volume) associations
    if not cell_type in cell_data.keys():
        print("invalid mesh: needs hexahedron entry in cell_data")
        sys.exit(1)
    hex_dict = cell_data[cell_type]

    # meshio creates 'gmsh:physical' and 'gmsh:geometrical' items, using physical volume associations
    if not 'gmsh:physical' in hex_dict.keys():
        print("invalid mesh: needs 'gmsh:physical' entry in cell_data for hexahedron")
        sys.exit(1)
    hex_mat = hex_dict['gmsh:physical']
    print("  number of hexahedra with material = ",len(hex_mat))
    # checks if number of physical hexahedras match
    if len(hex_mat) != len(hexas):
        print("invalid numbers of hexahedras ",len(hexas),"and material associations ",len(hex_mat))
        print("please define physical volumes using Gmsh")
        sys.exit(1)
    filename = "MESH/materials_file"
    with open(filename,"w") as f:
        for i,mat in enumerate(hex_mat):
            # format: #id_element #flag
            f.write("%d %d\n" % (i+1,mat) )
    print("  exported to: ",filename)
    print("")

    # creates a default material
    # format: #(1)material_domain_id #(2)material_id  #(3)rho  #(4)vp   #(5)vs   #(6)Q_kappa   #(7)Q_mu  #(8)anisotropy_flag
    #
    # note: here, we will write a default nummaterial_velocity_file file (since SPECFEM3D requires one for external meshes).
    #       it must be overwritten with your actual values after running this script.
    #
    # sets default values
    domain_id = 2 # 1==acoustic/2==elastic
    rho = 2300.0
    vp = 2800.0
    vs = 1500.0
    q_kappa = 9999.0
    q_mu = 9999.0
    ani_flag = 0

    # checks mesh definitions
    if len(hex_mat) == 0:
        print("Error: need hexahedra with material, please define physical volumes using Gmsh")
        sys.exit(1)

    material_ids = [hex_mat[0]]
    for mat in hex_mat:
        if not mat in material_ids: material_ids.append(mat)
    print("  material ids: ",material_ids)

    filename = "MESH/nummaterial_velocity_file"
    if not os.path.isfile(filename):
        print("  creating a default nummaterial velocity file...")
        with open(filename,"w") as f:
            for mat_id in material_ids:
                # format: #(1)material_domain_id #(2)material_id  #(3)rho  #(4)vp   #(5)vs   #(6)Q_kappa   #(7)Q_mu  #(8)anisotropy_flag
                f.write("%d %d %f %f %f %f %f %d\n" % (domain_id,mat_id,rho,vp,vs,q_kappa,q_mu,ani_flag))
            # adds info
            f.write("\n")
            f.write("! note: format of nummaterial_velocity_file must be\n")
            f.write("! #(1)material_domain_id #(2)material_id  #(3)rho  #(4)vp   #(5)vs   #(6)Q_kappa   #(7)Q_mu  #(8)anisotropy_flag\n")
        print("  exported to: ",filename)
    print("")


    # surface elements
    print("surfaces:")
    if not quad_type in cells.keys():
        print("invalid mesh: surfaces need quad cells")
        sys.exit(1)
    quads = cells[quad_type]
    print("  number of quads       = ",len(quads))
    print("")

    surf_top,surf_bottom,surf_xmin,surf_xmax,surf_ymin,surf_ymax = detect_boundary_surfaces(points,cells,cell_data,field_data)

    # surfaces
    short     = ["zmax", "bottom", "xmin", "xmax", "ymin", "ymax" ]
    filenames = ["MESH/free_or_absorbing_surface_file_zmax", "MESH/absorbing_surface_file_bottom",
                 "MESH/absorbing_surface_file_xmin", "MESH/absorbing_surface_file_xmax",
                 "MESH/absorbing_surface_file_ymin", "MESH/absorbing_surface_file_ymax" ]
    surfaces  = [ surf_top, surf_bottom, surf_xmin, surf_xmax, surf_ymin, surf_ymax ]

    for name,filename,surface in zip(short,filenames,surfaces):
        # info
        print("surface {:>8}: {}".format(name,len(surface)))
        if len(surface) > 0:
            with open(filename,"w") as f:
                f.write("%d\n" % len(surface))
                for elem in surface:
                    if quad_type == 'quad9':
                        # quad 9 face -> NGNOD2D == 9 used together with hex NGNOD == 27 elements
                        # check
                        if len(elem) != 9 + 1:
                            print("Invalid quad length {} for quad type {}, exiting...".format(len(elem),quad_type))
                            sys.exit(1)
                        # creates entry string
                        str = ""
                        for j in range(len(elem)): str += " {}".format(elem[j]+1)  # nodes
                        str += "\n"
                        f.write(str)
                    else:
                        # quad 4 face corners -> NGNOD2D = 4
                        # format: #id_(element containing the face) #id_node1_face .. #id_node4_face
                        # note: specfem assumes indices start at 1 (fortran style, not at 0 like in python or C)
                        f.write("%d %d %d %d %d\n" % (elem[0]+1,elem[1]+1,elem[2]+1,elem[3]+1,elem[4]+1))
            print("  exported to: ",filename)
            print("")

    print("")
    print("  done, see files in directory: MESH/")
    print("")
    return

#--------------------------------------------------------------------------------------------------

def export2SPECFEM3D(filename):
    """
    converts Gmsh file to SPECFEM mesh files
    """
    global meshio_major_version

    # version info
    version = meshio.__version__
    meshio_major_version = int(version.split(".")[0])
    print("meshio version: ",version)
    print("")

    # reads in mesh file
    points,cells,point_data,cell_data,field_data = read_mesh_file_msh(filename)

    # exports in specfem format
    export_mesh(points,cells,point_data,cell_data,field_data)

    print("finished Gmsh2specfem successfully")
    print("")
    return

def usage():
    print("usage: ./Gmsh2specfem.py file")
    print("   where")
    print("       file - a Gmsh file, *.msh (Gmsh Ascii-format)")
    return

if __name__ == '__main__':
    # gets argument
    if len(sys.argv) != 2:
        usage()
        sys.exit(1)
    else:
        filename = sys.argv[1]

    # main routine
    export2SPECFEM3D(filename)

