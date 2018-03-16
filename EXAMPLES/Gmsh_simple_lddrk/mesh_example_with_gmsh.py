#!/usr/bin/env python
#
# uses pygmsh:
#   https://github.com/nschloe/pygmsh
#   https://pygmsh.readthedocs.io/en/latest/
#
# install with:
#   pip install -U pygmsh
#
# this alreay includes module meshio. in case needed separately, install with:
#   pip install -U meshio
#
# in case you have already created an input file *.geo for Gmsh scripting,
# rather use:
# > from subprocess import call
# > call(["gmsh", "-2", "-order 2", "input.geo", "-o", "out.msh"])
#
# this script creates a volume with an acoustic and elastic layer.
# it uses gmsh to generate a hexahedral mesh.
#
from __future__ import print_function
import sys
import os
import math

import numpy as np

try:
    import pygmsh
except:
    print("importing module pygmsh failed, please make sure to install it via: pip install -U pygmsh")
    sys.exit(1)

import meshio

# from python file import
from Gmsh2specfem import export2SPECFEM3D


##############################################################
##
## Model Parameters
##
##############################################################

# model size
xsize = 10000.0
ysize = 8000.0
zsize = 5000.0

# model parameters
rho = 2300.0      # kg/m^3
vp = 2800.0       # m/s
vs = 1500.0       # m/s
Q_kappa = 400.0
Q_mu = 300.0
aniso_flag = 0
domain_id = 2     # 1==acoustic/ 2==elastic

# mesh size
mesh_element_size_Z = 1250.0
mesh_element_size_XY = 1000.0

##############################################################


#--------------------------------------------------------------------------------------------------
#
# helper functions
#
#--------------------------------------------------------------------------------------------------

def print_script(script=""):
    print("")
    print("-------- geo script")
    print(script)
    print("-------- end geo script")
    print("")

#--------------------------------------------------------------------------------------------------

def get_dimension(vol):
    """
    determines dimensions of volume, i.e., xmin/xmax, ymin/ymax and zmin/zmax
    """
    # initializes
    xmin = 1.e24
    xmax = -1.e24
    ymin = 1.e24
    ymax = -1.e24
    zmin = 1.e24
    zmax = -1.e24
    # surfaces
    surfaces = vol.surface_loop.surfaces
    for s in surfaces:
        # lines
        for l in s.line_loop.lines:
            # points
            for p in l.points:
                x,y,z = p.x
                if x < xmin: xmin = x
                if x > xmax: xmax = x
                if y < ymin: ymin = y
                if y > ymax: ymax = y
                if z < zmin: zmin = z
                if z > zmax: zmax = z
    return xmin,xmax,ymin,ymax,zmin,zmax

#--------------------------------------------------------------------------------------------------

def get_surface_midpoint(surface):
    """"
    determines midpoin on surface
    """
    # gets surface edges
    corners = []
    lines = surface.line_loop.lines
    #print("number of lines = ",len(lines))
    for l in lines:
        # corner points
        for p in l.points:
            #print("point ",p.x,"corners ",corners)
            x,y,z = p.x
            add_point = True
            for corner in corners:
                if corner[0] == x and corner[1] == y and corner[2] == z: add_point = False
            if add_point:
                corners.append([x,y,z])
    # checks if 4 corner points
    if len(corners) != 4:
        print("invalid number of surface corners",corners)
        sys.exit(1)
    # midpoint
    midpoint = [0.0, 0.0, 0.0]
    for corner in corners:
        midpoint[0] += corner[0]
        midpoint[1] += corner[1]
        midpoint[2] += corner[2]
    midpoint[0] = midpoint[0] * 0.25
    midpoint[1] = midpoint[1] * 0.25
    midpoint[2] = midpoint[2] * 0.25

    return midpoint

#--------------------------------------------------------------------------------------------------

def get_top_surface(vol):
    """
    determines surface at top of volume
    """
    # bounding dimensions
    xmin,xmax,ymin,ymax,zmin,zmax = get_dimension(vol)

    # checks z-coordinate of midpoints
    top_surf = None
    for s in vol.surface_loop.surfaces:
        mid = get_surface_midpoint(s)
        # check z-coord
        if mid[2] == zmax:
            # sets as top surface
            print(s.id,mid,zmax)
            top_surf = s
    # checks if found
    if top_surf == None:
        print("failed to find top surface",vol)
        sys.exit(1)

    return top_surf

#--------------------------------------------------------------------------------------------------

def stretch_to_elevation(points,zmin,zmax):
    """
    stretches z-coordinates of points uniformly from bottom zmin to top zmax
    """
    global xsize,ysize

    print("stretching vertical coordinates:")

    # layer height
    H = zmax - zmin

    # checks if anything to do
    if H <= 0.0:
        return points

    # simple sine-function as elevation variations
    print("  using sine-function")
    kx = 2.0 * 3.14159265 / xsize
    ky = 4.0 * 3.14159265 / ysize
    # amplitude
    A = 200.0

    ele_max = -1.e24
    ele_min = +1.e24

    for p in points:
        x = p[0]
        y = p[1]
        z = p[2]

        # determine elevation at x,y with respect to zmax
        elevation = A * math.sin(x*kx) * math.sin(y*ky)

        if elevation > ele_max: ele_max = elevation
        if elevation < ele_min: ele_min = elevation

        # stretch factor
        if z <= zmin:
            fac = 0.0
        else:
            fac = (z - zmin) / H  # scales between [0,1]

        # stretches z coordinate
        znew = z + fac * elevation

        p[2] = znew

    print("  elevation min/max = ",ele_min,ele_max)
    print("")

    return points


#--------------------------------------------------------------------------------------------------


def create_material_file():
    """
    overwrites material file
    """
    global rho,vp,vs,Q_kappa,Q_mu,aniso_flag,domain_id

    print("creating nummaterial velocity file...")
    print("  rho,vp,vs    = ",rho,vp,vs)
    print("  Q_mu,Q_kappa = ",Q_mu,Q_kappa)
    print("  aniso_flag   = ",aniso_flag)
    print("  domain id    = ",domain_id)
    print("")

    # default volume id gets material
    mat_id = 1

    filename = "MESH/nummaterial_velocity_file"
    with open(filename,"w") as f:
        # format: #(1)material_domain_id #(2)material_id  #(3)rho  #(4)vp   #(5)vs   #(6)Q_kappa   #(7)Q_mu  #(8)anisotropy_flag
        f.write("%d %d %f %f %f %f %f %d\n" % (domain_id,mat_id,rho,vp,vs,Q_kappa,Q_mu,aniso_flag))
        # adds info
        f.write("\n")
        f.write("! note: format of nummaterial_velocity_file must be\n")
        f.write("! #(1)material_domain_id #(2)material_id  #(3)rho  #(4)vp   #(5)vs   #(6)Q_kappa   #(7)Q_mu  #(8)anisotropy_flag\n")
    print("  written to: ",filename)
    print("")

    return

#--------------------------------------------------------------------------------------------------
#
# mesh functions
#
#--------------------------------------------------------------------------------------------------


def mesh_3D():
    """
    creates a simple model with hexahedral mesh
    """
    #
    # note: we want to create a hexahedral mesh with mostly regular shape, like a structured grid.
    #       Gmsh allows this by using an extrusion. we will thus create a surface at top and extrude it down along vertical direction.
    #
    #       however, the extrude-command will only return identifiers for the top and lateral surfaces, but not the bottom.
    #       we will thus create explicitly the bottom surface, which will get meshed and merged with the mesh from the extrusion.
    #       this will allow us to define physical surfaces for all bounding surfaces, and the mesh file will contain all surface quads.
    #
    #       for topography, using a top surface with topography and then extrude it down would lead to problems (flat bottom surface?).
    #       thus to obtain a top with some elevation, we will stretch the mesh points accordingly after the mesh was generated.
    #       it is done here in a simple way to show how one could modify the mesh.
    global xsize,ysize,zsize,mesh_element_size_XY,mesh_element_size_Z

    # output directory for mesh files
    os.system("mkdir -p MESH/")

    # dimensions
    xmin = 0.0
    xmax = xmin + xsize

    ymin = 0.0
    ymax = ymin + ysize

    z_top = 0.0
    z_bottom = -zsize

    # characteristic length of mesh elements
    lc_XY = mesh_element_size_XY

    number_of_element_layers_Z = int(zsize / mesh_element_size_Z)
    lc_Z = zsize / number_of_element_layers_Z

    print("")
    print("meshing:")
    print("characteristic length: ",lc_XY,lc_Z)

    ## geometry
    # gmsh geometry object
    geom = pygmsh.built_in.Geometry()

    ## mesh options: http://gmsh.info/doc/texinfo/gmsh.html#Mesh-options-list
    # quads
    geom.add_raw_code('Mesh.RecombineAll = 1;')
    # hex
    geom.add_raw_code('Mesh.Recombine3DAll = 1;')
    # uses DelQuad algorithm, instead of default Blossom for regular meshing
    # 1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=BAMG, 8=DelQuad
    geom.add_raw_code('Mesh.Algorithm = 8;')
    # uses frontal hex - outcommented as it will distort mesh
    # 1=Delaunay, 2=New Delaunay, 4=Frontal, 5=Frontal Delaunay, 6=Frontal Hex, 7=MMG3D, 9=R-tree
    #geom.add_raw_code('Mesh.Algorithm3D = 6;')
    # turns off mesh smoothing - outcommented as it will distort mesh
    #geom.add_raw_code('Mesh.Smoothing = 0;')
    # mesh element order 1=linear,2=quadratic,.. (example works for linear elements only)
    geom.add_raw_code('Mesh.ElementOrder = 1;')
    # add lateral surfaces to extrude output
    geom.add_raw_code('Geometry.ExtrudeReturnLateralEntities = 1;')

    ## creates geometry
    # creates top rectangle
    rec_top = geom.add_rectangle(xmin, xmax, ymin, ymax, z_top, lc_XY)
    surface_top = rec_top.surface
    print("top surface    : ",surface_top.id,"midpoint = ",get_surface_midpoint(surface_top))

    # creates rectangle for bottom
    rec_bottom = geom.add_rectangle(xmin, xmax, ymin, ymax, z_bottom, lc_XY)
    surface_bottom = rec_bottom.surface
    print("bottom surface : ",surface_bottom.id,"midpoint = ",get_surface_midpoint(surface_bottom))

    ## meshing
    ## creates a structured grid using extrude command

    # extrudes surface mesh along z-axis
    print("extrude surface: ",surface_top.id)
    # direction and length
    axis = [0,0,z_bottom]

    # uniform element-layers
    #top,vol,lat = geom.extrude(surface_top, translation_axis=axis, num_layers=number_of_element_layers_Z, recombine=True)

    # splits 2 element-layers into top (0.1 * height) layer, rest of element-layers for second (from 0.1 till bottom) layer
    layers = '{2,%d},{0.1,1}' % (number_of_element_layers_Z - 2)
    top,vol,lat = geom.extrude(surface_top, translation_axis=axis, num_layers=layers, recombine=True)

    # to make sure geometry points get merged, should be default however...
    geom.add_raw_code('Coherence;')

    # Physical Volume
    geom.add_physical_volume(vol, label="vol1")

    # Physical Surfaces
    # note: extrude only returns lateral surfaces. we thus had to create the bottom surface explicitly.
    #       geometries and mesh from extruded surface will then merge with bottom surface.
    # top and bottom
    geom.add_physical_surface(surface_top, label='top')
    geom.add_physical_surface(surface_bottom, label='bottom')
    # lateral sides
    geom.add_physical_surface(lat[0], label='ymin')
    geom.add_physical_surface(lat[1], label='xmax')
    geom.add_physical_surface(lat[2], label='ymax')
    geom.add_physical_surface(lat[3], label='xmin')

    # explicit volume meshing .. not needed, will be done when called by: gmsh -3 .. below
    #geom.add_raw_code('Mesh 3;')
    #geom.add_raw_code('Coherence Mesh;')

    # gets Gmsh script produced so far
    script = geom.get_code()
    print_script(script)

    # saves corresponding Gmsh script
    filename = "MESH/box.geo"
    with open(filename,"w") as f:
      f.write(script)
    print("script written to: ",filename)
    print("")

    # meshing
    points, cells, point_data, cell_data, field_data = pygmsh.generate_mesh(geom,verbose=True,num_quad_lloyd_steps=0)

    # mesh info
    print("")
    print("Mesh results: ")
    print("  points: ",len(points))
    print("  cells: ",len(cells))
    for name,array in cells.items():
        print("        ",name,len(array))
    print("  point_data: ",len(point_data))
    print("  cell_data: ",len(cell_data))
    for name,dic in cell_data.items():
        print("        ",name)
        for s,array in dic.items():
            print("          ",s,len(array))
    print("  field_data: ",len(field_data))
    for name,array in field_data.items():
        print("        ",name,len(array))
    print("")

    # uniform vertical stretch from bottom (at z_bottom) to a topography (with respect to z_top)
    points = stretch_to_elevation(points,z_bottom,z_top)

    # save as vtk-file
    filename = "MESH/box.vtu"
    meshio.write(filename, points, cells)
    print("VTK file written to : ",filename)

    # saves as Gmsh-file (msh mesh format)
    filename = "MESH/box.msh"
    meshio.write(filename, points, cells, point_data=point_data, cell_data=cell_data, field_data=field_data,
                 file_format='gmsh-ascii')
    print("Gmsh file written to: ",filename)
    print("")

    # export
    print("exporting Gmsh file to specfem format...")
    export2SPECFEM3D(filename)

    # overwrite material properties
    create_material_file()

    return

#--------------------------------------------------------------------------------------------------


def create_mesh():
    # version info
    print("pygmsh version: ",pygmsh.__version__)

    version = pygmsh.get_gmsh_major_version()
    print("Gmsh version  : ",version)

    # creates a simple 3D mesh
    mesh_3D()


if __name__ == '__main__':
    create_mesh()

