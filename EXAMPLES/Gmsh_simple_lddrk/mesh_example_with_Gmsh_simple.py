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
import gmsh

# from python file import
sys.path.append('../../utils/Cubit_or_Gmsh/')
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

# globals
python_major_version = 0
pygmsh_major_version = 0
meshio_major_version = 0

#--------------------------------------------------------------------------------------------------
#
# helper functions
#
#--------------------------------------------------------------------------------------------------

def print_script(geom=None):
    """
    prints out .geo script
    """
    print("")
    print("current script:")
    print("")

    if pygmsh_major_version >= 7:
        # pygmsh versions 7.x
        # note: with pygmsh versions >= 7, the *.geo script with get_code() will not be created anymore,
        #       as it uses the gmsh-python interface directly (instead of calling the gmsh-executable at the end)
        #       still, it provides a save_geometry() function which uses gmsh.write() to output mesh infos.
        #
        #       the generated .geo_unrolled script is basic and doesn't contain parameter settings or extrude calls.
        #       also, it misses physical groups if called before generate_mesh().
        #
        # saves corresponding Gmsh script (bare-bone, without mesh options)
        filename = "MESH/box.geo_unrolled"
        geom.save_geometry(filename)
        print("script written to: ",filename)
        # user output
        with open(filename,'r') as f:
            lines = f.readlines()
            script = ""
            for l in lines:
                script += "%s" % l
            print("")
            print("-------- geo unrolled script")
            print(script)
            print("-------- end geo unrolled script")
            print("")
    else:
        # pygmsh versions 6.x
        script = geom.get_code()
        print("")
        print("-------- geo script")
        print(script)
        print("-------- end geo script")
        print("")

        # saves corresponding Gmsh script
        filename = "MESH/box.geo"
        with open(filename,"w") as f:
            f.write(script)
            print("script written to: ",filename)
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
    global pygmsh_major_version

    # gets surface edges
    corners = []

    if pygmsh_major_version >= 7:
        # pygmsh versions 7.x
        lines = surface.curve_loop.curves
    else:
        # pygmsh versions 6.x
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
#
# modelling routines
#
#--------------------------------------------------------------------------------------------------

def create_gmsh_model(xmin,xmax,ymin,ymax,z_top,z_bottom,lc_XY,number_of_element_layers_Z):
    """
    creates a box-like model meshed with gmsh
    """
    global pygmsh_major_version,meshio_major_version
    global mesh_type
    global hex_type

    # gmsh geometry object
    if pygmsh_major_version >= 7:
        # pygmsh versions 7.x
        print("using pygmsh version 7...")
        print("")

        # main change in workflow with version 7: uses gmsh python interface directly
        # (instead of an explicit call to the gmsh executable with a .geo file)
        #
        # version 7.x changes usage of geometry object, and uses with-construct to initialize python object:
        #
        #     with pygmsh.geo.Geometry() as geom:
        #         geom.add_rectangle(xmin, xmax, ymin, ymax, z_top, lc_XY)
        #         ..
        #
        #     using without the with-construct, misses gmsh's initialization and an error occurs:
        #     Error   : Gmsh has not been initialized
        #
        # we'll call __enter__() and __exit__() object-functions explicitly to be able to use the
        # same lines of code below between pygmsh version 6.x and 7.x
        #
        # includes initialization:
        #    gmsh.initialize()
        geom = pygmsh.geo.Geometry()
        geom.__enter__()

    else:
        # pygmsh version 6.x
        geom = pygmsh.built_in.Geometry()

    ## mesh options: http://gmsh.info/doc/texinfo/gmsh.html#Mesh-options-list
    # quads
    if pygmsh_major_version >= 7:
        # pygmsh versions 7.x
        # no add_raw_code() anymore.. adding as gmsh options directly
        gmsh.option.setNumber("Mesh.RecombineAll", 1)
        gmsh.option.setNumber("Mesh.Recombine3DAll", 1)
        gmsh.option.setNumber("Mesh.Algorithm", 8)
        # mesh element order 1=linear,2=quadratic,.. (example works for linear elements only)
        gmsh.option.setNumber("Mesh.ElementOrder", 1)
        # add lateral surfaces to extrude output
        gmsh.option.setNumber("Geometry.ExtrudeReturnLateralEntities", 1)
    else:
        # pygmsh versions 6.x
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
    if pygmsh_major_version >= 7:
        # pygmsh versions 7.x
        surface_top_id = surface_top._id
    else:
        # pygmsh versions 6.x
        surface_top_id = surface_top.id
    print("top surface    : ",surface_top_id,"midpoint = ",get_surface_midpoint(surface_top))

    # creates rectangle for bottom
    rec_bottom = geom.add_rectangle(xmin, xmax, ymin, ymax, z_bottom, lc_XY)
    surface_bottom = rec_bottom.surface
    if pygmsh_major_version >= 7:
        # pygmsh versions 7.x
        surface_bottom_id = surface_bottom._id
    else:
        # pygmsh versions 6.x
        surface_bottom_id = surface_bottom.id
    print("bottom surface : ",surface_bottom_id,"midpoint = ",get_surface_midpoint(surface_bottom))

    ## meshing
    ## creates a structured grid using extrude command

    # extrudes surface mesh along z-axis
    print("extrude surface: ",surface_top_id)
    # direction and length
    axis = [0,0,z_bottom]

    # uniform element-layers
    #top,vol,lat = geom.extrude(surface_top, translation_axis=axis, num_layers=number_of_element_layers_Z, recombine=True)

    # splits 2 element-layers into top (0.1 * height) layer, rest of element-layers for second (from 0.1 till bottom) layer
    if pygmsh_major_version >= 7:
        # pygmsh versions 7.x
        layers = [2,number_of_element_layers_Z - 2]
        heights = [0.1,1]
        print("           layers ",layers," heights ",heights)
        top,vol,lat = geom.extrude(surface_top, translation_axis=axis, num_layers=layers, heights=heights, recombine=True)
    else:
        # pygmsh versions 6.x
        layers = '{2,%d},{0.1,1}' % (number_of_element_layers_Z - 2)
        print("           layers ",layers)
        top,vol,lat = geom.extrude(surface_top, translation_axis=axis, num_layers=layers, recombine=True)

    # to make sure geometry points get merged, should be default however...
    if pygmsh_major_version >= 7:
        # pygmsh versions 7.x
        # built-in geometry kernel Gmsh executes the Coherence command automatically, default is already set to 1
        # just to make sure...
        gmsh.option.setNumber("Geometry.AutoCoherence", 1)
    else:
        # pygmsh versions 6.x
        geom.add_raw_code('Coherence;')

    # Physical Volume
    if pygmsh_major_version < 6:
        # pygmsh versions 5.x
        geom.add_physical_volume(vol, label="vol1")
    else:
        geom.add_physical(vol, label="vol1")

    # Physical Surfaces
    # note: extrude only returns lateral surfaces. we thus had to create the bottom surface explicitly.
    #       geometries and mesh from extruded surface will then merge with bottom surface.
    # top and bottom
    if pygmsh_major_version < 6:
        # pygmsh versions 5.x
        geom.add_physical_surface(surface_top, label='top')
        geom.add_physical_surface(surface_bottom, label='bottom')
        # lateral sides
        geom.add_physical_surface(lat[0], label='ymin')
        geom.add_physical_surface(lat[1], label='xmax')
        geom.add_physical_surface(lat[2], label='ymax')
        geom.add_physical_surface(lat[3], label='xmin')
    else:
        geom.add_physical(surface_top, label='top')
        geom.add_physical(surface_bottom, label='bottom')
        # lateral sides
        geom.add_physical(lat[0], label='ymin')
        geom.add_physical(lat[1], label='xmax')
        geom.add_physical(lat[2], label='ymax')
        geom.add_physical(lat[3], label='xmin')

    # explicit volume meshing .. not needed, will be done when called by: gmsh -3 .. below
    #geom.add_raw_code('Mesh 3;')
    #geom.add_raw_code('Coherence Mesh;')

    # gets Gmsh script produced so far (also useful for debugging)
    if pygmsh_major_version < 7:
        print_script(geom)

    # meshing
    cells = {}
    cell_data = {}
    points = []
    point_data = []
    field_data = {}
    cell_sets = {}
    if pygmsh_major_version >= 7:
        # pygmsh versions 7.x
        mesh = geom.generate_mesh(verbose=True)
        # note: pygmsh versions 7.x require meshio versions >= 4.x
        #
        # note: for some reason, the mesh returned by generate_mesh() is missing cell data & point data.
        #       as a workaround, we save the geometry, then read it in again by meshio to have a complete mesh-object
        # saves original mesh
        print("")
        print("original mesh file...")
        geom.save_geometry("MESH/box_org.msh")
        #debug
        #print(""); print(mesh); print("")
        # reads in original mesh to fill meshio's mesh-object
        mesh = meshio.read("MESH/box_org.msh")
        #debug
        #print(""); print(mesh); print(""); mesh.write("MESH/box0.msh", file_format='gmsh',binary=False)
        print("")
        # gets data
        points,cells,point_data,cell_data,field_data,cell_sets = (mesh.points,mesh.cells_dict,mesh.point_data,
                                                                  mesh.cell_data_dict,mesh.field_data,mesh.cell_sets)
    else:
        # pygmsh versions 6.x
        mesh = pygmsh.generate_mesh(geom,verbose=True)
        # gets data
        # version changes:
        #   pygmsh <= 6.0.2 & meshio 3.3.1: cells and cell_data are returned as dictionary to mesh object
        #   pygmsh >= 6.0.3 & meshio >= 4.0.0: cells and cell_data are returns as lists to mesh object,
        #                                      instead use cells_dict and cell_data_dict; new provides cell_sets
        if meshio_major_version >= 4:
            points,cells,point_data,cell_data,field_data,cell_sets = (mesh.points,mesh.cells_dict,mesh.point_data,
                                                                      mesh.cell_data_dict,mesh.field_data,mesh.cell_sets)
        else:
            points,cells,point_data,cell_data,field_data = (mesh.points,mesh.cells,mesh.point_data,
                                                            mesh.cell_data,mesh.field_data)

    # gets Gmsh script .geo_unrolled file
    if pygmsh_major_version >= 7:
        print_script(geom)

    # mesh output info
    print("")
    print("Mesh results: ")
    print("  points: ",len(points))
    print("  cells: ",len(cells))
    if isinstance(cells, dict):
        # older pygmsh version <= 6.0.2
        for name,array in cells.items():
            print("        ",name,len(array))
    else:
        # newer pygmsh versions (cells is a list object)
        for cell_block in cells:
            name = cell_block[0]
            array = cell_block[1]
            print("        ",name,len(array))
    print("  point_data: ",len(point_data))
    print("  cell_data : ",len(cell_data))
    for name,data in cell_data.items():
        print("        ",name)
        if isinstance(data, dict):
            # older pygmsh version <= 6.0.2
            for s,array in data.items():
                print("          ",s,len(array))
        else:
            # newer pygmsh versions (data is a list object)
            for s,array in enumerate(data):
                print("          ",s,len(array))
    print("  field_data: ",len(field_data))
    for name,array in field_data.items():
        print("        ",name,len(array))
    print("  cell_sets : ",len(cell_sets))
    for name,array in cell_sets.items():
        print("        ",name,len(array))
    print("")
    print("")

    return mesh, points, cells, point_data, cell_data, field_data


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
    global python_major_version,pygmsh_major_version

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

    # output info
    lc_Z = zsize / number_of_element_layers_Z
    print("")
    print("meshing:")
    print("characteristic length: ",lc_XY,lc_Z)

    ## geometry
    mesh, points, cells, point_data, cell_data, field_data = create_gmsh_model(xmin,xmax,ymin,ymax,z_top,z_bottom,
                                                                               lc_XY,number_of_element_layers_Z)

    ## mesh modification
    # uniform vertical stretch from bottom (at z_bottom) to a topography (with respect to z_top)
    points = stretch_to_elevation(points,z_bottom,z_top)

    # creates new modified mesh
    # (otherwise, points have been modified by its reference pointer and thus mesh contains now modification)
    if False:
        # creates new mesh object with modified point locations
        mesh = meshio.Mesh(points, cells, point_data, cell_data, field_data)

    # save as vtk-file
    filename = "MESH/box.vtu"
    meshio.write(filename, mesh)
    print("VTK file written to : ",filename)

    # saves as Gmsh-file (msh mesh format)
    filename = "MESH/box.msh"
    if pygmsh_major_version >= 7:
        # pygmsh versions 7.x
        # note: reading back in with Gmsh2specfem.py requires msh version 2 format for now as
        #       version 4.x format will return different mesh list orderings...
        meshio.write(filename, mesh, file_format='gmsh22',binary=False)
    else:
        # pygmsh versions 6.x
        meshio.write(filename, mesh, file_format='gmsh2-ascii')
    print("Gmsh file written to: ",filename)
    print("")

    ## export
    print("exporting Gmsh file to specfem format...")
    export2SPECFEM3D(filename)

    # overwrite material properties
    create_material_file()

    return

#--------------------------------------------------------------------------------------------------


def create_mesh():
    global python_major_version,pygmsh_major_version,meshio_major_version

    # version info
    python_major_version = sys.version_info[0]
    python_minor_version = sys.version_info[1]
    print("Python version: ","{}.{}".format(python_major_version,python_minor_version))

    # meshio version
    version = meshio.__version__
    meshio_major_version = int(version.split(".")[0])
    print("meshio version: ",version)

    # pygmsh version
    version = pygmsh.__version__
    pygmsh_major_version = int(version.split(".")[0])
    print("pygmsh version: ",version)

    # Gmsh version
    if pygmsh_major_version >= 7:
        # pygmsh version >= 7.x
        version = pygmsh.__gmsh_version__
    else:
        # pygmsh version 6.x
        version = pygmsh.get_gmsh_major_version()
    print("Gmsh version  : ",version, " module version ",gmsh.__version__)

    # creates a simple 3D mesh
    mesh_3D()


if __name__ == '__main__':
    create_mesh()

