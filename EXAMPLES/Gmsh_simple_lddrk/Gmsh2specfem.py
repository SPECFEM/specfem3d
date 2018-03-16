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

#--------------------------------------------------------------------------------------------------

def read_mesh_file_msh(file):
    """
    reads in Gmsh file
    """
    print("")
    print("reading file: ",file)
    print("")

    # reads *.msh file (Gmsh Ascii-format)
    points, cells, point_data, cell_data, field_data = meshio.read(file, file_format='gmsh-ascii')

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
    for name,dic in cell_data.items():
        print("  ",name)
        for s,array in dic.items():
            print("    ",s,len(array))

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

    # cleans points array (removing unused points)
    points = clean_out_unused_points(points,cells)

    return points,cells,point_data,cell_data,field_data

#--------------------------------------------------------------------------------------------------


def clean_out_unused_points(points,cells):
    # for some reason, there might be extra nodes stored. we will clean them out and only use the onces needed by the hexas.
    # assumes that the surface quads won't require new nodes.
    print("cleaning points...")

    # checks if anything to do
    if not 'hexahedron' in cells.keys():
        return points

    # get array with all point entries
    cell_all = cells['hexahedron'].flatten()

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
    print("determining boundary surfaces:")

    surf_top = []
    surf_bottom = []
    surf_xmin = []
    surf_xmax = []
    surf_ymin = []
    surf_ymax = []

    # checks
    if not 'quad' in cells.keys():
        print("invalid mesh: surfaces need 'quad' cells")
        sys.exit(1)
    if not 'hexahedron' in cells.keys():
        print("invalid mesh: surfaces need 'hexahedron' cells")
        sys.exit(1)

    quads = cells['quad']
    hexas = cells['hexahedron']

    # checks if valid entries
    if len(quads) == 0:
        print("Error: need entries in 'quad' cells")
        sys.exit(1)
    if len(hexas) == 0:
        print("Error: need entries in 'hexahedron' cells")
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
    if not 'quad' in cell_data.keys():
        print("invalid mesh: needs 'quad' entry in cell_data")
        sys.exit(1)
    quad_dict = cell_data['quad']

    # meshio creates 'gmsh:physical' and 'gmsh:geometrical' items, using physical volume associations
    if not 'gmsh:physical' in quad_dict.keys():
        print("invalid mesh: needs 'gmsh:physical' entry in cell_data for 'quads'")
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

    # the following assumes linear elements, i.e., quads defined by 4 corners only
    if len(quads[0]) != 4:
        print("Invalid element type for quads: need linear elements, higher order elements not supported yet")

    # loops through quads
    for i,q in enumerate(quads):
        #print(i,q)   # i starts at 0,1,2,..; quads[id1,id2,id3,id4] indices also start at 0

        # for specfem, we will need to find corresponding hexahedral element which contains this quad
        el = [ [j,h] for j,h in enumerate(hexas) if q[0] in h and q[1] in h and q[2] in h and q[3] in h ]

        # boundary surfaces should only associate with a single hexahedral element
        if len(el) != 1:
            print("Error: could not find hexahedral for surface quad ",i,q,"el = ",el)
            sys.exit(1)
        hex_index = el[0][0]

        # adds surface id to quad corners
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
    if not 'hexahedron' in cells.keys():
        print("invalid mesh: needs 'hexahedron' cells")
        sys.exit(1)
    hexas = cells['hexahedron']
    print("  number of hexahedra   = ",len(hexas))
    filename = "MESH/mesh_file"
    with open(filename,"w") as f:
        f.write("%d\n" % len(hexas))
        for i,h in enumerate(hexas):
            # format: # element_id  #id_node1 ... #id_node8
            # note: specfem assumes indices start at 1 (fortran style, not at 0 like in python or C)
            f.write("%d %d %d %d %d %d %d %d %d\n" % (i+1,h[0]+1,h[1]+1,h[2]+1,h[3]+1,h[4]+1,h[5]+1,h[6]+1,h[7]+1) )
    print("  exported to: ",filename)
    print("")

    # material (or volume) associations
    if not 'hexahedron' in cell_data.keys():
        print("invalid mesh: needs 'hexahedron' entry in cell_data")
        sys.exit(1)
    hex_dict = cell_data['hexahedron']

    # meshio creates 'gmsh:physical' and 'gmsh:geometrical' items, using physical volume associations
    if not 'gmsh:physical' in hex_dict.keys():
        print("invalid mesh: needs 'gmsh:physical' entry in cell_data for 'hexahedron'")
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

    # default values
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
    if not 'quad' in cells.keys():
        print("invalid mesh: surfaces need 'quad' cells")
        sys.exit(1)
    quads = cells['quad']
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
    # reads in mesh file
    points,cells,point_data,cell_data,field_data = read_mesh_file_msh(filename)

    # exports in specfem format
    export_mesh(points,cells,point_data,cell_data,field_data)


def usage():
    print("usage: ./Gmsh2specfem.py file")
    print("   where")
    print("       file - a Gmsh file, *.msh (Gmsh Ascii-format)")


if __name__ == '__main__':
    # gets argument
    if len(sys.argv) != 2:
        usage()
        sys.exit(1)
    else:
        filename = sys.argv[1]

    # main routine
    export2SPECFEM3D(filename)

