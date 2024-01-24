#!/usr/bin/env blender --background --python-use-system-env --python plot_with_blender.py --
#
# or for example on Mac:
#  /usr/bin/env /Applications/Blender.app/Contents/MacOS/Blender --background --python-use-system-env --python plot_blender_sphere.py --
#
#
# run with: > ./plot_with_blender.py --help
#
###############################################################################################

import sys
import os
import time
import datetime

# blender
try:
    import bpy
except:
    print("Failed importing bpy. Please make sure to have Blender python working properly.")
    sys.exit(1)

print("")
print("Blender: version ",bpy.app.version_string)
print("")

# blender meshing
import bmesh

# adds path
#sys.path.append('/Applications/Blender.app/Contents/Resources/3.6/python/lib/python3.10/site-packages')

try:
    import vtk
except:
    print("Failed importing vtk. Please make sure to have Blender python with vtk working properly.")
    print("system path:")
    for path in sys.path:
        print("  ",path)
    sys.exit(1)

print("VTK: version ",vtk.vtkVersion().GetVTKVersion())
print("")

try:
    import numpy as np
except:
    print("Failed importing numpy. Please make sure to have Blender python with numpy working properly.")
    sys.exit(1)

from mathutils import Vector,Matrix
from math import radians,atan2,sin,cos

###############################################################################################
## USER parameters

## renderer
# render image size
blender_img_resolution_X = 2400
blender_img_resolution_Y = 1600

# cycles rendering (for better glass effect of buildings)
use_cycles_renderer = False

## transparent sea-level plane
use_transparent_sea_level_plane = False

###############################################################################################

# Constants
PI = 3.141592653589793
DEGREE_TO_RAD = PI / 180.0

# Global variables
mesh_scale_factor = 1.0
mesh_origin = [0.0,0.0,0.0]

# baseline shift for buildings (to move above SPECFEM mesh)
shift_building_baseline = 0.0  # 0.0==no-shift default, or e.g. --shift-building-baseline=0.0005

# rendering output
suppress_renderer_output = False

# animation
use_animation_dive_in = True
use_animation_rotation = True

# camera positioning
camera_elevation_offset = 0.15
camera_y_offset = -0.7
camera_angle_depth_degree = 75.0

# center point to focus camera view on
x_center = 0.0
y_center = 0.0
z_center = 0.0

# elevation at center point
z_elevation = 0.0


# class to avoid long stdout output by renderer
# see: https://stackoverflow.com/questions/24277488/in-python-how-to-capture-the-stdout-from-a-c-shared-library-to-a-variable/29834357
class SuppressStream(object):
    def __init__(self, stream=sys.stderr,suppress=False):
        # turns on/off suppressing stdout of renderer process
        self.SUPPRESS_STDOUT = suppress

        if self.SUPPRESS_STDOUT:
            self.orig_stream_fileno = stream.fileno()

    def __enter__(self):
        if self.SUPPRESS_STDOUT:
            self.orig_stream_dup = os.dup(self.orig_stream_fileno)
            self.devnull = open(os.devnull, 'w')
            os.dup2(self.devnull.fileno(), self.orig_stream_fileno)

    def __exit__(self, type, value, traceback):
        if self.SUPPRESS_STDOUT:
            os.close(self.orig_stream_fileno)
            os.dup2(self.orig_stream_dup, self.orig_stream_fileno)
            os.close(self.orig_stream_dup)
            self.devnull.close()


def convert_vtk_to_obj(vtk_file,colormap=0,color_max=None):
    global mesh_scale_factor,mesh_origin

    # Path to your .vtu file
    print("converting vtk file: ",vtk_file)

    # check file
    if len(vtk_file) == 0:
        print("Error: no vtk file specified...")
        usage()
        sys.exit(1)

    if not os.path.exists(vtk_file):
        print("Error: vtk file specified not found...")
        sys.exit(1)

    # gets file extension
    extension = os.path.splitext(vtk_file)[1]
    # reads the vtk file
    if extension == '.vtk':
        # .vtk file
        reader = vtk.vtkDataSetReader()
        reader.SetFileName(vtk_file)
        reader.Update()
    elif extension == '.vtu':
        # .vtu
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(vtk_file)
        reader.Update()
    elif extension == '.inp':
        # AVS .inp
        #reader = vtk.vtkSimplePointsReader()
        reader = vtk.vtkAVSucdReader()
        reader.SetFileName(vtk_file)
        reader.Update()
    else:
        print("unknown vtk file extension ",extension," - exiting...")
        sys.exit(1)

    #debug
    #print(reader)

    ## scale coordinates to be in range [-1,1] for visualization
    # gets the points (vertices) from the dataset
    points = reader.GetOutput().GetPoints()

    # number of points
    num_points = points.GetNumberOfPoints()

    #print(points)
    #print(points.GetBounds())
    xmin,xmax,ymin,ymax,zmin,zmax = points.GetBounds()

    # defines the scaling factors to fit within +/- 1
    #min_coords = np.array(points.GetPoint(0))
    #max_coords = np.array(points.GetPoint(0))
    #for i in range(1, num_points):
    #    point = np.array(points.GetPoint(i))
    #    min_coords = np.minimum(min_coords, point)
    #    max_coords = np.maximum(max_coords, point)

    min_coords = np.array([xmin,ymin,zmin])
    max_coords = np.array([xmax,ymax,zmax])
    dimensions = max_coords - min_coords

    # info
    print("  minimum coordinates:", min_coords)
    print("  maximum coordinates:", max_coords)
    print("  dimensions         :", dimensions)
    print("")

    # gets data values on nodes
    data_array = None
    if reader.GetOutput().GetPointData().GetNumberOfArrays() > 0:
        data_array = reader.GetOutput().GetPointData().GetArray(0)  # Example: Accessing the first data array
        #debug
        #print(data_array)
        #info
        print("  data: ")
        print("  range = ",data_array.GetRange())
        print("")


    # determines origin of mesh to move it back to (0,0,0) and scale it between [-1,1] to better locating it in blender
    mesh_origin = min_coords + 0.5 * (max_coords - min_coords)

    # z-coordinate: leaves sea level at 0 for mesh origin
    mesh_origin[2] = 0.0

    # takes maximum size in x/y direction
    dim_max = np.maximum(dimensions[0],dimensions[1])

    if np.abs(dim_max) > 0.0:
        mesh_scale_factor = 2.0 / dim_max

    print("  mesh scaling:")
    print("  origin       :",mesh_origin)
    print("  scale factor :",mesh_scale_factor)
    print("")

    # creates an array to store scaled points
    scaled_points = vtk.vtkPoints()

    # Loop through each point, scale its coordinates, and add it to the new points array
    for i in range(num_points):
        # translation by origin
        point = np.array(points.GetPoint(i)) - mesh_origin
        # uniform scaling
        scaled_point = point * mesh_scale_factor
        # stores updated points
        scaled_points.InsertNextPoint(scaled_point)

    # creates a new polydata with the scaled points
    scaled_polydata = vtk.vtkPolyData()
    scaled_polydata.SetPoints(scaled_points)

    # vertex data
    if data_array:
        num_points = data_array.GetNumberOfTuples()
        num_components = data_array.GetNumberOfComponents()
        data_min,data_max = data_array.GetRange()

        print("  data array: ")
        print("    number of points     = ",num_points)
        print("    number of components = ",num_components)
        print("    range min/max        = ",data_min,"/",data_max)
        print("")

        # checks if fixing maximum value
        if color_max:
            # Convert VTK data array to NumPy array
            array = np.zeros((num_points, num_components))
            for i in range(num_points):
                for j in range(num_components):
                    value = data_array.GetComponent(i, j)
                    array[i,j] = value

            # limit size
            if 1 == 0:
                ## determines maximum value as a multiple of 10
                # reduce first by 10%
                total_max = abs(array).max()
                total_max = 0.9 * total_max  # sets limit at 90% of the maximum

                # get maximum value in power of 10
                if total_max != 0.0:
                    total_max = 1.0 * 10**(int(np.log10(total_max)))  # example: 2.73e-11 limits to 1.e-11
                    #total_max = 1.0 * 10**(int(np.log10(total_max))-1)  # example: 2.73e-11 limits to 1.e-11
                    #total_max = 1.0 * 10**(int(np.log10(total_max))-2)  # example: 2.73e-11 limits to 1.e-12
                else:
                    total_max = 0.0
                    print("  data: color data min/max   = ",array.min(),array.max())
                    print("  data: zero color data - nothing to show")
                    # nothing left to do
                    #sys.exit(1)

            # checks if fixing maximum value
            if color_max:
                total_max = color_max

            print("  limiting color data range:")
            print("    data min/max   = ",array.min(),array.max())
            if color_max:
                print("    data total max = ",total_max," (fixed)")
            else:
                print("    data total max = ",total_max)
            print("")

            # limits range [-total_max,total_max]
            array = np.where(array < -total_max, -total_max, array)
            array = np.where(array > total_max, total_max, array)

            # in case color-max is larger than actual range, this sets an arbitrary point to the maximum value
            # to get the correct range value when plotting the data
            if abs(array).max() < total_max:
                array[0,0] = total_max

            # for shakemaps, start range with a minimum value of 0
            if 'shaking' in vtk_file or 'shakemap' in vtk_file:
                array[1,0] = 0.0

            # sets updated range back to vtk array
            print("    new data: color data min/max   = ",array.min(),array.max())
            for i in range(num_points):
                for j in range(num_components):
                    value = array[i,j]
                    data_array.SetComponent(i, j, value)

            # Inform VTK that the array has been modified
            data_array.Modified()
            data_min,data_max = data_array.GetRange()
            print("    new data: range min/max = ",data_min,"/",data_max)
            print("")

        # sets vertex data
        scaled_polydata.GetPointData().SetScalars(data_array)

    # updates the points in the original dataset with the scaled points
    reader.GetOutput().SetPoints(scaled_points)
    reader.Update()

    # output scaled bounds
    points = reader.GetOutput().GetPoints()
    xmin,xmax,ymin,ymax,zmin,zmax = points.GetBounds()
    print("  mesh dimensions after scaling: x min/max = ",xmin,xmax)
    print("                                 y min/max = ",ymin,ymax)
    print("                                 z min/max = ",zmin,zmax)
    print("")

    # Get the unstructured grid data
    unstructured_grid = reader.GetOutput()

    # Extract colors if available
    colors_array = None
    if unstructured_grid.GetPointData().GetNumberOfArrays() > 0:
        print("  grid: point data arrays = ",unstructured_grid.GetPointData().GetNumberOfArrays())
        print("")
        colors_array = unstructured_grid.GetPointData().GetArray(0)  # Assuming colors are in the first array
        print("  colors: name = ",colors_array.GetName())
        print("          range = ",colors_array.GetRange())
        print("")

    #debug
    #print("colors_array: ",colors_array)

    # convert to .ply data file
    # Path to your generated .ply file
    obj_file = 'output.ply'

    # convert the data to polydata
    geometry_filter = vtk.vtkGeometryFilter()
    geometry_filter.SetInputData(unstructured_grid)
    geometry_filter.Update()

    # creates a new polydata object
    polydata = geometry_filter.GetOutput()

    print("  polydata: initial number of points",polydata.GetNumberOfPoints())
    print("  polydata: initial number of verts",polydata.GetNumberOfVerts())
    print("")

    # checks if we have points
    # for proc****_free_surface.vtk files, only points are stored in the .vtk file
    # and the geometry_filter won't fill the polydata object.
    # here we check that we have points & verts filled, otherwise we assume to have points only in the .vtk file
    # and we will try to get a connectivity by Delauney triangulation
    if polydata.GetNumberOfPoints() == 0:
        print("  unstructured grid: number of points",unstructured_grid.GetNumberOfPoints())
        if unstructured_grid.GetNumberOfPoints() > 0:
            points = unstructured_grid.GetPoints()
            polydata.SetPoints(points)

    if polydata.GetNumberOfVerts() == 0:
        print("  unstructured grid: number of cells",unstructured_grid.GetNumberOfCells())
        if unstructured_grid.GetNumberOfCells() > 0:
            verts = unstructured_grid.GetCells()
            polydata.SetVerts(verts)
        else:
            # Perform Delaunay triangulation to generate connectivity
            print("  getting Delaunay 2D connectivity...")
            delaunay = vtk.vtkDelaunay2D()
            delaunay.SetInputData(polydata)
            delaunay.Update()

            # Get the output triangles
            triangles = delaunay.GetOutput()
            #print(triangles)
            print("  triangles: number of verts",triangles.GetNumberOfVerts())
            print("  triangles: number of cells",triangles.GetNumberOfCells())
            polydata = triangles

        print("  polydata: number of points",polydata.GetNumberOfPoints())
        print("  polydata: number of verts",polydata.GetNumberOfVerts())
        print("  polydata: number of cells",polydata.GetNumberOfCells())
        print("  polydata: number of strips",polydata.GetNumberOfStrips())
        print("  polydata: number of data arrays",polydata.GetPointData().GetNumberOfArrays())
        print("")

    # we need to set a default lookup table for the polydata set,
    # otherwise the color float values on the points won't get stored
    #
    # set lookup table for depth values to colors
    if colors_array:
        lut = vtk.vtkLookupTable()
        lut.SetTableRange(data_min, data_max)
        lut.SetNumberOfTableValues(256) # Set the number of table values
        #lut.SetRampToLinear()
        #lut.SetRampToSQRT()

        # VTK by default maps the value range to colors from red to blue
        # determine custom type
        if colormap == 0:
            print("  color map: default VTK")
            # nothing to special to add, let's just vtk internally do it
            # colormap is going from red to white to blue
        elif colormap == 1:
            # topo
            print("  color map: topo")
            colors_rgb = [
                [0.3,  0.3,   0.3], # gray
                [0.1,  0.1,   0.4], # blue
                [0.2,  0.5,   0.2],
                [0.25, 0.625, 0.5],
                [0.0,  0.5,   0.25],
                [0.5,  0.365, 0.0],
                [0.75, 0.625, 0.25],
                [1.0,  0.75,  0.625],
                [1.0,  0.75,  0.5],
                [1,    1,     1],   # white
            ]

        elif colormap == 2:
            # Scientific Colour Map Categorical Palette
            # https://www.fabiocrameri.ch/colourmaps/
            # Crameri, F. (2018). Scientific colour maps. Zenodo. http://doi.org/10.5281/zenodo.1243862
            print("  color map: lisbon")
            # lisbon 10 Swatches
            colors_rgb255 = [
                [230, 229, 255],  #  lisbon-1 #E6E5FF
                # or start with a less white color
                #[200, 208, 237], #  lisbon-12 #C8D0ED
                [155, 175, 211],  #  lisbon-29 #9BAFD3
                [ 81, 119, 164],  #  lisbon-58 #5177A4
                [ 30,  67, 104],  #  lisbon-86 #1E4368
                [ 17,  30,  44],  #  lisbon-114 #111E2C
                [ 39,  37,  26],  #  lisbon-143 #27251A
                [ 87,  81,  52],  #  lisbon-171 #575134
                [141, 133,  86],  #  lisbon-199 #8D8556
                [201, 195, 144],  #  lisbon-228 #C9C390
                [255, 255, 217],  #  lisbon-256 #FFFFD9
            ]
            # converts the colors from 0-255 range to 0-1 range
            colors_rgb = [[comp / 255.0 for comp in color] for color in colors_rgb255]

        elif colormap == 3:
            # Scientific Colour Map Categorical Palette
            # https://www.fabiocrameri.ch/colourmaps/
            # Crameri, F. (2018). Scientific colour maps. Zenodo. http://doi.org/10.5281/zenodo.1243862
            print("  color map: lajolla")
            # lajolla 10 Swatches
            colors_rgb255 = [
                [ 25,  25,   0], #  lajolla-1 #191900
                [ 51,  34,  15], # lajolla-29 #33220F
                [ 91,  48,  35], #  lajolla-58 #5B3023
                [143,  64,  61], #  lajolla-86 #8F403D
                [199,  80,  75], #  lajolla-114 #C7504B
                [224, 114,  79], #  lajolla-143 #E0724F
                [231, 148,  82], #  lajolla-171 #E79452
                [238, 181,  85], #  lajolla-199 #EEB555
                [248, 223, 124], #  lajolla-228 #F8DF7C
                [255, 254, 203], #  lajolla-256 #FFFECB
            ]
            # converts the colors from 0-255 range to 0-1 range
            colors_rgb = [[comp / 255.0 for comp in color] for color in colors_rgb255]

        elif colormap == 4:
            # Scientific Colour Map Categorical Palette
            # https://www.fabiocrameri.ch/colourmaps/
            # Crameri, F. (2018). Scientific colour maps. Zenodo. http://doi.org/10.5281/zenodo.1243862
            print("  color map: lipari")
            # lipari 10 Swatches
            colors_rgb255 = [
                [  3,  19,  38], #  lipari-1 #031326
                [ 19,  56,  90], #  lipari-29 #13385A
                [ 71,  88, 122], #  lipari-58 #47587A
                [107,  95, 118], #  lipari-86 #6B5F76
                [142,  97, 108], #  lipari-114 #8E616C
                [188, 100,  97], #  lipari-143 #BC6461
                [229, 123,  98], #  lipari-171 #E57B62
                [231, 162, 121], #  lipari-199 #E7A279
                [233, 201, 159], #  lipari-228 #E9C99F
                [253, 245, 218], #  lipari-256 #FDF5DA
            ]
            # converts the colors from 0-255 range to 0-1 range
            colors_rgb = [[comp / 255.0 for comp in color] for color in colors_rgb255]

        elif colormap == 5:
            # Scientific Colour Map Categorical Palette
            # https://www.fabiocrameri.ch/colourmaps/
            # Crameri, F. (2018). Scientific colour maps. Zenodo. http://doi.org/10.5281/zenodo.1243862
            print("  color map: davos")
            # davos 10 Swatches
            colors_rgb255 = [
                [  0,   5,  74], #  davos-1 #00054A
                [ 17,  44, 113], #  davos-29 #112C71
                [ 41,  82, 145], #  davos-58 #295291
                [ 67, 112, 157], #  davos-86 #43709D
                [ 94, 133, 152], #  davos-114 #5E8598
                [121, 150, 141], #  davos-143 #79968D
                [153, 173, 136], #  davos-171 #99AD88
                [201, 210, 158], #  davos-199 #C9D29E
                [243, 243, 210], #  davos-228 #F3F3D2
                [254, 254, 254], #  davos-256 #FEFEFE
            ]
            # converts the colors from 0-255 range to 0-1 range
            colors_rgb = [[comp / 255.0 for comp in color] for color in colors_rgb255]

        elif colormap == 6:
            # Scientific Colour Map Categorical Palette
            # https://www.fabiocrameri.ch/colourmaps/
            # Crameri, F. (2018). Scientific colour maps. Zenodo. http://doi.org/10.5281/zenodo.1243862
            print("  color map: turku")
            # turku 10 Swatches
            colors_rgb255 = [
                [  0,   0,   0], #  turku-1 #000000
                [ 36,  36,  32], #  turku-29 #242420
                [ 66,  66,  53], #  turku-58 #424235
                [ 95,  95,  68], #  turku-86 #5F5F44
                [126, 124,  82], #  turku-114 #7E7C52
                [169, 153, 101], #  turku-143 #A99965
                [207, 166, 124], #  turku-171 #CFA67C
                [234, 173, 152], #  turku-199 #EAAD98
                [252, 199, 195], #  turku-228 #FCC7C3
                [255, 230, 230], #  turku-256 #FFE6E6
            ]
            # converts the colors from 0-255 range to 0-1 range
            colors_rgb = [[comp / 255.0 for comp in color] for color in colors_rgb255]

        elif colormap == 7:
            # Scientific Colour Map Categorical Palette
            # https://www.fabiocrameri.ch/colourmaps/
            # Crameri, F. (2018). Scientific colour maps. Zenodo. http://doi.org/10.5281/zenodo.1243862
            print("  color map: berlin")
            # berlin 10 Swatches
            colors_rgb255 = [
                [158, 176, 255], #  berlin-1 #9EB0FF
                [ 91, 164, 219], #  berlin-29 #5BA4DB
                [ 45, 117, 151], #  berlin-58 #2D7597
                [ 26,  66,  86], #  berlin-86 #1A4256
                [ 17,  25,  30], #  berlin-114 #11191E
                [ 40,  13,   1], #  berlin-143 #280D01
                [ 80,  24,   3], #  berlin-171 #501803
                [138,  63,  42], #  berlin-199 #8A3F2A
                [196, 117, 106], #  berlin-228 #C4756A
                [255, 173, 173], #  berlin-256 #FFADAD
            ]
            # converts the colors from 0-255 range to 0-1 range
            colors_rgb = [[comp / 255.0 for comp in color] for color in colors_rgb255]

        elif colormap == 8:
            # Scientific Colour Map Categorical Palette
            # https://www.fabiocrameri.ch/colourmaps/
            # Crameri, F. (2018). Scientific colour maps. Zenodo. http://doi.org/10.5281/zenodo.1243862
            print("  color map: grayC")
            colors_rgb255 = [
                [0,   0,   0],  #  grayC-1 #000000
                [35,  35,  35], #  grayC-29 #232323
                [61,  61,  61], #  grayC-58 #3D3D3D
                [86,  86,  86], #  grayC-86 #565656
                [108, 108, 108],#  grayC-114 #6C6C6C
                [130, 130, 130],#  grayC-143 #828282
                [154, 154, 154],#  grayC-171 #9A9A9A
                [182, 182, 182],#  grayC-199 #B6B6B6
                [216, 216, 216],#  grayC-228 #D8D8D8
                [255, 255, 255],#  grayC-256 #FFFFFF
            ]
            # converts the colors from 0-255 range to 0-1 range
            colors_rgb = [[comp / 255.0 for comp in color] for color in colors_rgb255]

        elif colormap == 9:
            # custom snow
            print("  color map: snow")
            colors_rgb255 = [
                [204, 204, 204], # gray
                [153, 178, 204],
                [ 71,  88, 122], #  lipari-58 #47587A
                [107,  95, 118], #  lipari-86 #6B5F76
                [142,  97, 108], #  lipari-114 #8E616C
                [188, 100,  97], #  lipari-143 #BC6461
                [229, 123,  98], #  lipari-171 #E57B62
                [231, 162, 121], #  lipari-199 #E7A279
                [255, 229, 204],
                [255, 255, 255], # white
            ]
            # converts the colors from 0-255 range to 0-1 range
            colors_rgb = [[comp / 255.0 for comp in color] for color in colors_rgb255]

        elif colormap == 10:
            # custom shakeGreen
            print("  color map: shakeGreen")
            colors_rgb = [
                [0.8,  0.8,   0.8], # gray
                [0.5,  0.5,   0.4],
                [0.5,  0.4,   0.2],
                [0.6,  0.6,   0.0], # green
                [0.72, 0.25,  0.0 ], # orange
                [0.81, 0.5,   0.0 ],
                [0.9,  0.74,  0.0 ], # yellow
                [1.0,  0.99,  0.0 ],
                [1.0,  0.99,  0.25],
                [1.0,  1.0,   1.0 ],  # white
            ]

        elif colormap == 11:
            # custom shakeRed
            print("  color map: shakeRed")
            colors_rgb = [
                [0.85, 0.85,  0.85], # gray
                [0.7,  0.7,   0.7 ],
                [0.5,  0.5,   0.5 ],
                [0.63, 0.0,   0.0 ], # red
                [0.72, 0.25,  0.0 ], # orange
                [0.81, 0.5,   0.0 ],
                [0.9,  0.74,  0.0 ], # yellow
                [1.0,  0.99,  0.0 ],
                [1.0,  0.99,  0.25],
                [1.0,  1.0,   1.0 ],  # white
            ]

        elif colormap == 12:
            # custom shakeUSGS
            # taken from a shakemap plot of the USGS
            # https://earthquake.usgs.gov/earthquakes/eventpage/us6000lqf9/shakemap/intensity
            print("  color map: shakeUSGS")
            colors_rgb = [
                [1.0,  1.0,  1.0 ], # I: white
                [0.8,  0.8,  0.8 ],
                [0.77, 0.81, 1.0 ], # II-III : light purple
                [0.5,  1.0,  0.98], # IV: turquoise
                [0.5,  1.0,  0.54], # V : green
                [1.0,  0.98, 0.0 ], # VI : yellow
                [1.0,  0.77, 0.0 ], # VII : orange
                [0.99, 0.52, 0.0 ], # VIII: darkorange
                [0.98, 0.0,  0.0 ], # IX : red
                [0.78, 0.0,  0.0 ], # X+ : darkred
            ]

        elif colormap == 13:
            # custom shakeUSGSgray
            # starts with darker gray than the default USGS
            print("  color map: shakeUSGSgray")
            colors_rgb = [
                [0.8,  0.8,  0.8 ], # gray
                [0.5,  0.5,  0.5 ],
                [0.77, 0.81, 1.0 ], # II-III : light purple
                [0.5,  1.0,  0.98], # IV: turquoise
                [0.5,  1.0,  0.54], # V : green
                [1.0,  0.98, 0.0 ], # VI : yellow
                [1.0,  0.77, 0.0 ], # VII : orange
                [0.99, 0.52, 0.0 ], # VIII: darkorange
                [0.98, 0.0,  0.0 ],  # IX : red
                [0.5, 0.0,  0.0 ],  # X+ : darkred
            ]

        else:
            print("Warning: colormap with type {} is not supported, exiting...".format(colormap))
            sys.exit(1)

        # sets lookup table entries
        if colormap != 0:
            # Create a vtkColorTransferFunction
            color_transfer_func = vtk.vtkColorTransferFunction()
            # add specific scalar values in the color transfer function
            for i, color in enumerate(colors_rgb):
                val = i / (len(colors_rgb) - 1.0)
                color_transfer_func.AddRGBPoint(val, color[0], color[1], color[2])
            # Calculate the color values for the lookup table by interpolating from the color transfer function
            for i in range(256):
                scalar = i / 255.0  # Normalized scalar value from 0 to 1
                color = color_transfer_func.GetColor(scalar)
                lut.SetTableValue(i, color[0], color[1], color[2], 1.0)
            print("")

        # build lookup table
        lut.Build()
        colors_array.SetLookupTable(lut)

    # Write the data to PLY format
    writer = vtk.vtkPLYWriter()
    writer.SetInputData(polydata)

    # Include vertex colors if available
    if colors_array:
        writer.SetArrayName(colors_array.GetName())
        writer.SetLookupTable(lut)
        # info
        print("  writer: color mode = ",writer.GetColorMode())
        print("  writer: color array name = ",writer.GetArrayName())
        print("  writer: color component  = ",writer.GetComponent())
        print("")

    os.system('rm -f output.ply')

    writer.SetFileName(obj_file)
    writer.Write()

    if not os.path.exists(obj_file):
        print("Error writing file ",obj_file)
        sys.exit(1)

    print("")
    print("  converted to: ",obj_file)
    print("")

    # work-around for .obj files
    # however, .obj file can by default only store the mesh, not the color data on the vertices
    # we thus prefer to work with the .ply file format above.
    #
    ## save mesh as .obj file
    ## Path to your generated .obj file
    #obj_file = 'output.obj'
    #
    ## Convert the data to polydata
    #geometry_filter = vtk.vtkGeometryFilter()
    #geometry_filter.SetInputConnection(reader.GetOutputPort())
    #geometry_filter.Update()
    #
    ## Write the data to .obj format
    #writer = vtk.vtkOBJWriter()
    #writer.SetInputConnection(geometry_filter.GetOutputPort())
    #
    #os.system('rm -f output.obj')
    #
    #writer.SetFileName(obj_file)  # Output .obj file path
    #writer.Write()
    #
    #if not os.path.exists(obj_file):
    #    print("Error writing file ",obj_file)
    #    sys.exit(1)
    #
    ## appends data lines
    ##..
    ##d val
    #if data_array:
    #    min_val = data_array.GetValue(0)
    #    max_val = data_array.GetValue(0)
    #
    #    with open(obj_file, 'a') as f:
    #        f.write("# Associated data values:\n")
    #        for i in range(data_array.GetNumberOfTuples()):
    #            data_value = data_array.GetValue(i)
    #            min_val = np.minimum(min_val,data_value)
    #            max_val = np.maximum(max_val,data_value)
    #            f.write(f"# vertex {i}: {data_value}\n")
    #            #f.write(f"d {data_value}\n")
    #    print("  appended data: min/max = ",min_val,max_val)
    #    data_min = min_val
    #    data_max = max_val
    #
    #print("")
    #print("  converted to: ",obj_file)

    return obj_file

def create_blender_setup(obj_file=""):
    ## Blender setup
    print("blender setup:")
    print("")

    # clears all existing objects in the scene
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete()

    # clears only default Cube object and any previous output & output_dfx object,
    # leaves Camera and Light object
    #print("  objects: ",bpy.data.objects.keys())
    #objs = bpy.data.objects
    #for obj in objs:
    #    print("  objects: ",obj.name)
    #    if obj.name == "Cube": objs.remove(obj, do_unlink=True)
    #    if obj.name.contains("output"): objs.remove(obj, do_unlink=True)

    print("  objects: ",bpy.data.objects.keys())
    objs = bpy.data.objects
    for obj in objs:
        print("  objects: ",obj.name)
    print("")

    # import mesh object into blender
    if len(obj_file) > 0:
        print("  importing mesh file: ",obj_file)
        # gets file extension
        extension = os.path.splitext(obj_file)[1]
        # reads the mesh object file
        if extension == '.ply':
            # Import .ply file
            # New experimental, but much faster
            bpy.ops.wm.ply_import(filepath=obj_file)
            # or standard/legacy ply importer
            #bpy.ops.import_mesh.ply(filepath=obj_file)
        elif extension == '.obj':
            # Import .obj file into Blender
            bpy.ops.import_scene.obj(filepath=obj_file)
        else:
            print("unknown mesh object file extension ",extension," - exiting...")
            sys.exit(1)

        print("  imported in blender: ",obj_file)
        print("")

    # blender info
    print("  scenes : ",bpy.data.scenes.keys())
    print("  objects: ",bpy.data.objects.keys())
    for obj in bpy.data.objects:
        print("    object: ",obj.name)
    print("")

    ## mesh object
    #obj = bpy.context.object
    #print(obj.name, ":", obj)
    #objs = bpy.context.selected_objects
    #print(", ".join(o.name for o in objs))
    # Select the imported object
    obj = bpy.data.objects['output']
    if obj == None:
        print("Error: no mesh object in blender available, exiting...")
        sys.exit(1)

    # object is a mesh
    mesh = obj.data

    #debug
    #print("  obj: ",obj)
    print("  obj type: ",obj.type)
    print("  obj data: ",mesh)
    print("  obj polygons: ",mesh.polygons)
    #print("  obj polygon 0: ",mesh.polygons[0])
    #print("  obj polygon 0 loop: ",mesh.polygons[0].loop_indices)
    #print("  obj data loops: ",mesh.loops)
    #print("  obj data loops: ",mesh.loops[0])
    print("  obj data vertex_colors: ",mesh.vertex_colors)
    #print("  obj data vertex_layers_float: ",mesh.vertex_layers_float)
    #print("  obj data vertex_layers_int: ",mesh.vertex_layers_int)
    #print("  obj data vertex_layers_string: ",mesh.vertex_layers_string)
    #print("  obj data vertex_colors: ",mesh.vertex_colors.get("a"))
    #for loop_index in mesh.polygons[0].loop_indices:
    #      index = mesh.loops[loop_index].vertex_index
    #      print("  loop: index",loop_index,obj.data.loops[loop_index].vertex_index)
    #      color = vertex_colors[loop_index].color
    #print("  obj data polygon_layers_float: ",mesh.polygon_layers_float)
    #print("  obj data polygon_layers_int: ",mesh.polygon_layers_int)
    #print("  obj data uv_layers: ",mesh.uv_layers)
    #print("  obj data vertices",mesh.vertices)
    #print("  obj data vertex 0",mesh.vertices[0])
    #print("  obj data vertex 0 coord",mesh.vertices[0].co)
    #print("  obj data vertex 0 keys",mesh.vertices.keys())
    #print("  obj data vertex 0 get",mesh.vertices.items())
    print("")

    #print("obj data vertex_colors 0: ",obj.data.vertex_colors[0])
    #print("obj data vertex_colors active data: ",obj.data.vertex_colors.active.data)

    if obj is not None:
        # Ensure the object has a mesh and vertex colors
        if obj.type == 'MESH' and not obj.data.vertex_colors is None:
            print("  Object 'output' has vertex colors and is a mesh.")
            # Access vertex color data
            #vertex_colors = obj.data.vertex_colors.active.data
            # Iterate through vertex color data
            #for poly in obj.data.polygons:
            #    for loop_index in poly.loop_indices:
            #        vertex_index = obj.data.loops[loop_index].vertex_index
            #        color = vertex_colors[loop_index].color
            #        # Print information about vertex colors
            #        #print(f"Vertex {vertex_index}: Color {color}")
        else:
            print("  Object 'output' does not have vertex colors or is not a mesh.")
    else:
        print("  Object 'output' not found.")
        sys.exit(1)

    print("")
    print("  mesh: setting up shader nodes...")
    print("")

    # assigns new material
    mat = bpy.data.materials.new(name="VertexColorMaterial")
    mat.use_nodes = True
    # node-graph
    nodes = mat.node_tree.nodes
    links = mat.node_tree.links

    # Clear default nodes
    #for node in nodes:
    #    nodes.remove(node)

    # Create Attribute node to fetch vertex color
    color_attribute = nodes.new(type='ShaderNodeAttribute')
    color_attribute.attribute_name = "Col"  # Use "Col" as it's the default name for vertex color
    # Set the Color Attribute node to use vertex colors
    color_attribute.attribute_type = 'GEOMETRY' # 'COLOR'

    # Create Principled BSDF shader node
    # checks default node
    bsdf = nodes["Principled BSDF"]
    if bsdf == None:
        bsdf = nodes.new(type='ShaderNodeBsdfPrincipled')

    bsdf.inputs['Metallic'].default_value = 0.4 # 0.1
    bsdf.inputs['Roughness'].default_value = 0.5 # 0.8
    bsdf.inputs['Specular'].default_value = 0.2 # 0.3

    # (default) color map from vtk file
    links.new(bsdf.inputs['Base Color'],color_attribute.outputs['Color'])

    # no emission
    bsdf.inputs['Emission Strength'].default_value = 0.0
    # w/ emission (default) for brighter colors
    #links.new(bsdf.inputs['Emission'],color_attribute.outputs['Color'])
    #links.new(bsdf.inputs['Emission Strength'],color_attribute.outputs['Fac'])

    # custom mesh coloring w/ color ramp node
    if 1 == 0:
        # creates a custom color ramp node to highlight shaking regions
        ramp = nodes.new(type='ShaderNodeValToRGB')
        ramp.color_ramp.interpolation = 'LINEAR'
        #ramp.color_ramp.interpolation = 'B_SPLINE'
        # default 2 slots
        ramp.color_ramp.elements[0].position = 0.0
        ramp.color_ramp.elements[0].color = [1,1,1,1]
        ramp.color_ramp.elements[1].position = 0.4
        ramp.color_ramp.elements[1].color = [0.5,0.5,0.5,1]  # gray
        # add color slot
        ramp.color_ramp.elements.new(0.5)
        ramp.color_ramp.elements[2].color = [0.2,0.2,0.3,1]  # dark blue
        # add color slot
        ramp.color_ramp.elements.new(0.6)
        ramp.color_ramp.elements[3].color = [0.8,0.0,0.0,1]  # red
        # add color slot
        ramp.color_ramp.elements.new(0.65)
        ramp.color_ramp.elements[4].color = [1.0,0.6,0.0,1]  # yellow
        # add color slot
        ramp.color_ramp.elements.new(0.7)
        ramp.color_ramp.elements[5].color = [1,1,1,1]       # white

        # Link Math node output to Color Ramp factor input
        links.new(ramp.inputs["Fac"],color_attribute.outputs["Fac"])
        # custom Color Ramp output to Principled BSDF node
        links.new(bsdf.inputs['Base Color'],ramp.outputs['Color'])

    # adds additional emission
    if 1 == 0:
        # Create Math node to manipulate grayscale value
        math_node = nodes.new(type='ShaderNodeMath')
        math_node.operation = 'GREATER_THAN'
        math_node.inputs[1].default_value = 0.7  # Set threshold value
        links.new(math_node.inputs['Value'],color_attribute.outputs['Fac'])

        # Create RGB to BW node to convert color to black/white float value
        #rgb_to_bw = nodes.new(type='ShaderNodeRGBToBW')
        #links.new(rgb_to_bw.inputs['Color'],color_attribute.outputs['Color'])
        #links.new(math_node.inputs['Value'],rgb_to_bw.outputs['Val'])

        # emission node
        emission = nodes.new('ShaderNodeEmission')
        emission.name = "Emission"
        links.new(emission.inputs['Color'],color_attribute.outputs['Color'])
        links.new(emission.inputs['Strength'],math_node.outputs['Value'])

        # mix light emission and main image
        mix = nodes.new('ShaderNodeMixShader')
        mix.name = "Mix Shader"
        #links.new(mix.inputs['Fac'], ramp.outputs['Alpha'])
        # takes output from main BSDF node
        links.new(mix.inputs[1], bsdf.outputs[0])
        links.new(mix.inputs[2], emission.outputs[0])

        # link mixer to final material output
        material_output = nodes["Material Output"]
        links.new(material_output.inputs["Surface"], mix.outputs["Shader"])

    # Assign the material to the object
    if obj.data.materials:
        obj.data.materials[0] = mat
    else:
        obj.data.materials.append(mat)

    print("  blender mesh done")
    print("")


def add_blender_buildings(buildings_file=""):
    """
    adds buildings given by input .ply file
    """
    global mesh_origin,mesh_scale_factor
    global use_cycles_renderer

    # checks if anything to do
    if len(buildings_file) == 0: return

    print("buildings file: ",buildings_file)
    print("")

    # check file
    if not os.path.exists(buildings_file):
        print("Error: buildings file specified not found...")
        sys.exit(1)

    # buildings need to be given as .ply (Stanford) format file
    print("  reading in .ply mesh...")

    # Enable the experimental .ply importer
    # New experimental, but much faster
    bpy.ops.wm.ply_import(filepath=buildings_file)
    # or standard/legacy ply importer
    #bpy.ops.import_mesh.ply(filepath=buildings_file)

    # Print the names of imported objects
    print("")
    for obj in bpy.context.selected_objects:
        print("  imported object:", obj.name," - type: ",obj.type)
        # select buildings mesh (name must contain dxf)
        #if obj.type == 'MESH' and 'dxf' in obj.name:
        #    bpy.context.view_layer.objects.active = obj
        #    obj_buildings = obj
        #    break
    print("")

    # gets active object
    obj = bpy.context.view_layer.objects.active

    if obj is None:
        print("  Object for buildings not found.")
        sys.exit(1)

    # note: due to different resolution of the meshes, the mesh elevation of SPECFEM mesh (e.g., AVS_shaking_map.inp)
    #       and the buildings (e.g., from SwissTopo) might be slightly off by a few meters.

    # moves and scales UTM mesh
    if 'utm' in buildings_file:
        # need to translate and scale the UTM mesh to place it within the vtk mesh
        print("  UTM mesh: moving & scaling mesh...")
        print("            mesh origin       = ",mesh_origin)
        print("            mesh scale factor = ",mesh_scale_factor)
        print("")

        # Get the mesh data
        mesh = obj.data

        # Create a BMesh object
        bm = bmesh.new()
        bm.from_mesh(mesh)

        # Translate the vertices (example: translating by (1, 0, 0))
        translation_vector = Vector(mesh_origin)
        for vert in bm.verts:
            vert.co -= translation_vector

        # Scale the vertices (example: scaling by 1.5 in all axes)
        scale_factor = mesh_scale_factor
        for vert in bm.verts:
            vert.co *= scale_factor

        # Update the mesh with the modified vertices
        bm.to_mesh(mesh)
        bm.free()

    # baseline shift
    if shift_building_baseline != 0.0:
        print("  shifting buildings up: factor = ",shift_building_baseline)
        print("                          in m  = ",shift_building_baseline / mesh_scale_factor)
        print("")
        obj.location = (0, 0, shift_building_baseline)

    # Create a new material
    mat = bpy.data.materials.new('buildingsMaterial')

    # enable node-graph edition mode
    mat.use_nodes = True

    nodes = mat.node_tree.nodes
    links = mat.node_tree.links

    # Clear default nodes
    for node in nodes:
        nodes.remove(node)

    # gets scene
    scene = bpy.context.scene

    if use_cycles_renderer:
        # Set render engine to Cycles
        scene.render.engine = 'CYCLES'

        # Create a Glass BSDF node
        mat_node = nodes.new(type='ShaderNodeBsdfGlass')
        # Set IOR (Index of Refraction) and Roughness values
        mat_node.inputs['IOR'].default_value = 1.5       # (e.g., 1.5 for typical glass)
        mat_node.inputs['Roughness'].default_value = 0.6  # Roughness value (0.0 for perfectly smooth)

    else:
        # BLENDER_EEVEE renderer
        scene.render.engine = 'BLENDER_EEVEE'

        # create a glossy node
        mat_node = nodes.new(type='ShaderNodeBsdfGlossy')
        mat_node.inputs['Roughness'].default_value = 0.5  # Roughness value (0.0 for perfectly smooth)
        mat_node.inputs['Color'].default_value = (0.8, 0.74, 0.56, 1)

        # create a glass node
        # note: transparency effect w/ eevee for buildings looks not good enough yet...
        #       glass buildings need cycles renderer.
        #       thus, i'm using glossy buildings to get some environment colors onto buildings.
        #mat_node = nodes.new(type='ShaderNodeBsdfGlass')
        # node color
        #mat_node.inputs['Color'].default_value = (0.6, 0.6, 0.6, 1)
        # enable transparency for eevee
        # material transparency
        #mat.blend_method = 'BLEND'   # 'OPAQUE', 'BLEND', ..
        #mat.shadow_method = 'HASHED'
        #mat.use_backface_culling = False
        #mat.show_transparent_back = True
        #mat.use_screen_refraction = True
        #mat.refraction_depth = 0.01

    # Create Output node
    output_node = nodes.new(type='ShaderNodeOutputMaterial')
    output_node.location = 400, 0

    # Link nodes
    links.new(mat_node.outputs['BSDF'], output_node.inputs['Surface'])

    # Assign the material to the object
    if obj.data.materials:
        obj.data.materials[0] = mat
    else:
        obj.data.materials.append(mat)

    # blender info
    print("  scenes : ",bpy.data.scenes.keys())
    print("  objects: ",bpy.data.objects.keys())
    for obj in bpy.data.objects:
        print("    object: ",obj.name)
    print("")

    print("  blender buildings mesh done")
    print("")


def get_mesh_elevation(point_of_interest):
    """
    determines elevation of object (obj) at a given point by ray intersection
    """
    # gets vtk mesh object
    obj = bpy.data.objects['output']
    if obj == None:
        print("Info: no mesh object in blender available to determine elevation...")
        return None

    # Define the origin of the ray (above the mesh)
    ray_origin = Vector((point_of_interest.x, point_of_interest.y, 10.0))  # Adjust the Z coordinate as needed

    # Get the world matrix of the object
    matrix_world = obj.matrix_world

    # Calculate the direction of the ray (pointing downwards)
    ray_direction = Vector((0, 0, -1))

    # Transform the ray direction to world space
    ray_direction.rotate(matrix_world.to_quaternion())

    # Perform the ray casting
    success, location, _, _ = obj.ray_cast(ray_origin, ray_direction)

    if success:
        # Get the elevation (Z coordinate) of the mesh at the point of interest
        return location.z
    else:
        return None


def get_mesh_elevation_at_origin():
    """
    determines elevation of mesh at origin/center point
    """
    global x_center,y_center,z_center

    print("    origin center            = {} / {} / {}".format(x_center,y_center,z_center))

    # Define the point of interest (X, Y, Z coordinates)
    point = Vector((x_center, y_center, z_center))  # origin

    # get elevation
    elevation = get_mesh_elevation(point)

    # checks if valid
    if not elevation is None:
        # got an elevation
        print("    mesh elevation at origin = ", elevation)
    else:
        # sets a default height
        elevation = 0.1
        print("    mesh elevation at origin: no intersection found")
        print("                              setting default = ",elevation)

    return elevation

def add_camera(title="",close_up_view=False):
    """
    adds camera object
    """
    global x_center,y_center,z_center
    global camera_y_offset,camera_elevation_offset,camera_angle_depth_degree
    global z_elevation

    print("  adding camera")

    # gets scene
    scene = bpy.context.scene

    # adds camera position
    bpy.ops.object.camera_add(enter_editmode=False, align='VIEW')
    # current object
    cam = bpy.context.object
    cam.name = "Camera"   #cam = bpy.data.objects["Camera"]
    scene.camera = cam

    # clip range
    # note: our objects are scaled between [-1,1].
    #       a clip end plane around 10 should be sufficient to cover the whole area.
    #       for the clip start, it is more challenging as choosing a value below 0.1 will lead to depth Z-buffer artifacts.
    #       to avoid this, for an overview scene, the default value of 0.1 is good, when zooming in into close-up views,
    #       a clip start around 0.01 is usually better, otherwise the front of the scene will be cut-off.
    scene.camera.data.clip_start = 0.1    # blender default: 0.1
    scene.camera.data.clip_end = 10.0     # blender default: 1000
    #scene.camera.data.lens = 70.0

    # determine elevation from mesh at origin point to position the camera
    z_elevation = get_mesh_elevation_at_origin()

    # elevation offset for camera positioning
    camera_elevation_offset = 0.15
    camera_y_offset = -0.7
    camera_angle_depth_degree = 75.0

    # specifics
    if title == 'Lauterbrunnen' or title == 'Zermatt':
        print("    camera setup for Lauterbrunnen/Zermatt")
        camera_elevation_offset = 0.25
        camera_angle_depth_degree = 70.0

    # center point to focus camera view on
    x_center = 0.0
    y_center = 0.0
    z_center = 0.0

    ## Rome colosseum
    if title == "Rome":
        print("    camera setup for Rome")
        camera_elevation_offset = 0.0018
        camera_y_offset = -0.013
        camera_angle_depth_degree = 82.0
        # use colosseum position as center of views
        if 1 == 1:
            # colosseum: N 41.890258 E 12.492335
            #            33 T -> utm 291957.383 4640632.710
            x_center_utm = 291957.383; y_center_utm = 4640632.710
            # need to translate and scale the UTM position to place it within the scene
            point = Vector((x_center_utm, y_center_utm, 0.0))
            translation_vector = Vector(mesh_origin)
            point -= translation_vector
            point *= mesh_scale_factor
            # set center position
            x_center = point.x
            y_center = point.y
        else:
            # directly from using cursor position in blender
            x_center = -0.01251
            y_center = -0.2263
        z_center = 0.0018

        # get elevation at x/y position
        point = Vector((x_center, y_center, 0.0))  # origin
        center_elevation = get_mesh_elevation(point)
        if not center_elevation is None:
            z_elevation = center_elevation
        else:
            z_elevation = 0.001

        #cam.location = (-0.02251, -0.24362, 0.00632)
        #cam.rotation_euler = (83.401 * DEGREE_TO_RAD, 0.000133 * DEGREE_TO_RAD, -28.0 * DEGREE_TO_RAD)

        # metallic material for buildings
        if not use_cycles_renderer:
            if 'buildingsMaterial' in bpy.data.materials.keys():
                mat = bpy.data.materials['buildingsMaterial']
                nodes = mat.node_tree.nodes
                mat_node = nodes['Glossy BSDF']
                mat_node.inputs['Roughness'].default_value = 0.25  # shinier, more metallic
                print("      using metallic buildings")

    if close_up_view:
        ## close-up view
        # adjust clip range
        scene.camera.data.clip_start = 0.01
        #scene.camera.data.clip_end = 10.0
        #scene.camera.data.lens = 70.0
        # camera location
        x_cam = x_center
        y_cam = y_center + camera_y_offset
        z_cam = z_center + z_elevation + camera_elevation_offset
        print("    close-up scene")
        print("    center point   : x/y/z = {:.6f} / {:.6f} / {:.6f}".format(x_center,y_center,z_center))
        print("    camera location: x/y/z = {:.6f} / {:.6f} / {:.6f}".format(x_cam,y_cam,z_cam))
        print("    camera location relative to mesh elevation: ",z_elevation)
        # Set camera translation
        scene.camera.location = (x_cam,y_cam,z_cam)
        # Set camera rotation in euler angles
        scene.camera.rotation_mode = 'XYZ'
        scene.camera.rotation_euler = (camera_angle_depth_degree * DEGREE_TO_RAD, 0, 0)
        # Set camera fov in degrees
        scene.camera.data.angle = float(30.0 * DEGREE_TO_RAD)
    else:
        ## overview
        # camera location
        x_cam = x_center
        y_cam = y_center - 4.0
        z_cam = z_center + 4.0
        print("    overview scene")
        print("    center point   : x/y/z = {:.6f} / {:.6f} / {:.6f}".format(x_center,y_center,z_center))
        print("    camera location: x/y/z = {:.6f} / {:.6f} / {:.6f}".format(x_cam,y_cam,z_cam))
        # Set camera translation
        scene.camera.location = (x_cam,y_cam,z_cam)
        # Set camera rotation in euler angles
        scene.camera.rotation_mode = 'XYZ'
        scene.camera.rotation_euler = (44.0 * DEGREE_TO_RAD, 0, 0)
        # Set camera fov in degrees
        scene.camera.data.angle = float(30.0 * DEGREE_TO_RAD)


def add_light():
    """
    adds a light to the scene
    """
    print("  adding light")
    # adds sun position
    bpy.ops.object.light_add(type='AREA')
    # current object
    light = bpy.context.object
    light.name = "Light"  #light = bpy.data.objects["Light"]

    light.location = (1.5, -0.5, 1.3)
    light.rotation_mode = 'XYZ'
    light.rotation_euler = (0, 40.0 * DEGREE_TO_RAD, 0)
    # sets light to rectangular (plane)
    light.data.type = 'AREA'
    light.data.shape = 'SQUARE'
    # Change the light's color
    light.data.color = (1, 1, 1)  # Set light color to white
    # intensity
    light.data.energy = 80  # W
    # Set the light's size
    light.data.size = 0.1
    light.data.use_contact_shadow = True  # Enable contact shadows


def add_plane(close_up_view=False):
    """
    adds a plane at sea-level
    """
    global use_transparent_sea_level_plane

    print("  adding sea-level plane")
    # Create a mesh plane (to capture shadows and indicate sea level)
    bpy.ops.mesh.primitive_plane_add(size=10, enter_editmode=False, location=(0, 0, 0))

    # Get the created plane object
    plane_object = bpy.context.object

    # Set the object's material to white
    mat = bpy.data.materials.new(name="White")
    # adds transparency to the plane
    if use_transparent_sea_level_plane:
        print("    using transparent sea-level")
        # Set the material to use a Principled BSDF shader
        mat.use_nodes = True
        principled_bsdf = mat.node_tree.nodes.get('Principled BSDF')
        # color
        if close_up_view:
            principled_bsdf.inputs["Base Color"].default_value = (0.8, 0.8, 0.8, 1)
        else:
            principled_bsdf.inputs["Base Color"].default_value = (0.135, 0.135, 0.135, 1) # gray for better contrast
        # Set the shader to be transparent
        principled_bsdf.inputs['Alpha'].default_value = 0.6  # Set the alpha to control transparency
        # Set the blend mode to 'Alpha Blend'
        mat.blend_method = 'BLEND'
        # Set the shadow mode to 'Alpha Hashed'
        mat.shadow_method = 'HASHED'
    else:
        if close_up_view:
            mat.diffuse_color = (0.8, 0.8, 0.8, 1)  # similar as background color to have an infinite background
        else:
            mat.diffuse_color = (0.135, 0.135, 0.135, 1)  # gray for better contrast

    plane_object.data.materials.append(mat)


def add_title(title=""):
    """
    adds a title text
    """
    if len(title) > 0:
        print("  adding text object:")
        print("    title = ",title)

        # Create a new text object
        bpy.ops.object.text_add(location=(-0.3, -1.2, 0.01))  # Adjust the location as needed
        text_object = bpy.context.object
        text_object.data.body = title  # Set the text content

        # Set text properties (font, size, etc.)
        text_object.data.size = 0.2  # Adjust the font size

        #debug
        #print("  blender fonts available: ",bpy.data.fonts.keys())

        if 'Bfont' in bpy.data.fonts:
            text_object.data.font = bpy.data.fonts['Bfont']  # Use a specific default font
        elif 'Bfont Regular' in bpy.data.fonts:
            text_object.data.font = bpy.data.fonts['Bfont Regular']  # Use a specific default font
        elif 'Arial Regular' in bpy.data.fonts:
            text_object.data.font = bpy.data.fonts['Arial Regular']

        #text_object.data.font = bpy.data.fonts.load("/path/to/your/font.ttf")  # Replace with your font path

        # Set text material
        text_material = bpy.data.materials.new(name="TextMaterial")
        text_material.use_nodes = True
        text_material.node_tree.nodes["Principled BSDF"].inputs["Base Color"].default_value = (0.5, 0.5, 0.5, 1)

        text_object.data.materials.append(text_material)

def set_scene():
    """
    sets scene options
    """
    global blender_img_resolution_X,blender_img_resolution_Y
    global use_cycles_renderer

    print("")
    print("  setting scene options")

    # gets scene
    scene = bpy.context.scene

    # renderer options
    if use_cycles_renderer:
        # Set render engine to Cycles
        scene.render.engine = 'CYCLES'

        ## adjust settings for faster rendering w/ cycles
        #scene.cycles.device = 'GPU'
        # Set tile size
        scene.cycles.tile_size = 512 # 256
        # Lower the number of samples
        scene.cycles.samples = 128 # Adjust the number of samples as needed
        # Enable denoising
        scene.cycles.use_denoising = True

        # Adjust light path bounces
        scene.cycles.max_bounces = 4  # Adjust bounce values for diffuse, glossy, transmission, etc.

        ## additional parameters to check for better performance:
        # Use Branched Path Tracing integrator
        #scene.cycles.progressive = 'BRANCHED_PATH'
        # Adjust Branched Path Tracing settings
        #scene.cycles.use_square_samples = False  # Enable square samples for the branched path tracing
        #scene.cycles.diffuse_samples = 3  # Adjust samples per type (diffuse, glossy, etc.)
        #scene.cycles.glossy_samples = 3
        #scene.cycles.transmission_samples = 3
        #scene.cycles.ao_samples = 3
        #scene.cycles.mesh_light_samples = 3
        #scene.cycles.subsurface_samples = 3

    else:
        # Set render engine to Eevee
        scene.render.engine = 'BLENDER_EEVEE'

        # ambient occlusion
        scene.eevee.use_gtao = True
        scene.eevee.gtao_distance = 0.01
        scene.eevee.gtao_factor = 1.0

        # turns on screen space reflections
        scene.eevee.use_ssr = True
        # refraction
        scene.eevee.use_ssr_refraction = True
        scene.eevee.ssr_thickness = 0.1

        # turns on bloom
        #scene.eevee.use_bloom = True

    # render resolution
    scene.render.resolution_x = blender_img_resolution_X
    scene.render.resolution_y = blender_img_resolution_Y
    # Set the render percentage (optional, for scaling down the final output)
    #scene.render.resolution_percentage = 50  # Adjust as needed

    # sets white background color
    bpy.data.worlds["World"].node_tree.nodes["Background"].inputs[0].default_value = (1, 1, 1, 1)

    print("    renderer: ",scene.render.engine)
    print("")


def render_animation():
    """
    renders a movie animation
    """
    global use_animation_dive_in,use_animation_rotation
    global suppress_renderer_output
    global camera_elevation_offset,camera_y_offset,camera_angle_depth_degree
    global z_elevation

    print("animation:")

    # gets scene
    scene = bpy.context.scene

    # initializes frames
    number_of_keyframes = 0
    total_frames = 0

    start_frame = -1
    end_frame = -1

    scene.frame_start = -1
    scene.frame_end = -1

    ## dive-in keyframes
    if use_animation_dive_in:
        # setup keyframes
        number_of_keyframes += 4
        keyframe_interval = 20

        total_frames = (number_of_keyframes-1) * keyframe_interval

        print("  dive-in:")
        print("    keyframe interval      = ",keyframe_interval)
        print("    number of keyframes    = ",number_of_keyframes)
        print("")

        # define a frame timeline
        start_frame = 0
        inter_frame_1 = keyframe_interval
        inter_frame_2 = 2*keyframe_interval
        end_frame   = total_frames

        scene.frame_start = start_frame
        scene.frame_end = end_frame

        # moves camera
        cam = bpy.data.objects["Camera"]

        # 1. start frame
        cam.keyframe_insert("location", frame=start_frame)        # (0, -4, 4)
        cam.keyframe_insert("rotation_euler", frame=start_frame)  # (44.0 * DEGREE_TO_RAD, 0, 0)

        # 2. intermediate frame
        cam.location = (0, -2.5, 1)
        cam.rotation_euler = (60.0 * DEGREE_TO_RAD, 0, 0)
        cam.keyframe_insert("location", frame=inter_frame_1)
        cam.keyframe_insert("rotation_euler", frame=inter_frame_1)
        # adjust clip range
        cam.data.clip_start = 0.1
        cam.data.keyframe_insert("clip_start", frame=inter_frame_1)
        #cam.data.clip_end = 10.0
        #cam.data.keyframe_insert("clip_end", frame=inter_frame_1)

        # 3. intermediate frame
        x_cam = 0
        y_cam = -1.2
        z_cam = z_center + z_elevation + camera_elevation_offset
        cam.location = (x_cam, y_cam, z_cam)
        cam.rotation_euler = (camera_angle_depth_degree * DEGREE_TO_RAD, 0, 0)
        cam.keyframe_insert("location", frame=inter_frame_2)
        cam.keyframe_insert("rotation_euler", frame=inter_frame_2)
        # adjust clip range
        cam.data.clip_start = 0.01
        cam.data.keyframe_insert("clip_start", frame=inter_frame_2)
        #cam.data.clip_end = 10.0
        #cam.data.keyframe_insert("clip_end", frame=inter_frame_2)

        # 4. end frame
        x_cam = x_center
        y_cam = y_center + camera_y_offset
        z_cam = z_center + z_elevation + camera_elevation_offset
        cam.location = (x_cam, y_cam, z_cam) # Ending position
        cam.rotation_euler = (camera_angle_depth_degree * DEGREE_TO_RAD, 0, 0)
        cam.keyframe_insert("location", frame=end_frame)
        cam.keyframe_insert("rotation_euler", frame=end_frame)

    ## rotation
    if use_animation_rotation:
        # appends frames for rotation
        frames = 60  # Number of frames for the animation
        keyframe_interval = 5
        number_of_keyframes += frames

        total_frames += frames * keyframe_interval

        print("  rotation:")
        print("    keyframe interval      = ",keyframe_interval)
        print("    number of keyframes    = ",frames)
        print("")

        if start_frame == -1:
            start_frame = 0
        else:
            start_frame = end_frame + 1 * keyframe_interval

        if end_frame == -1:
            end_frame = total_frames
        else:
            end_frame = total_frames + 1 * keyframe_interval

        # updates end frame count
        scene.frame_end = end_frame

        # Set the rotation parameters
        rotation_degrees = 360  # Total rotation in degrees

        # moves camera around the z-axis
        cam = bpy.data.objects['Camera']

        #if scene.frame_start == -1:
        #    # add start frame
        #    scene.frame_start = start_frame
        #    # start frame
        #    cam.keyframe_insert("location", frame=start_frame)        # (0, -4, 4)
        #    cam.keyframe_insert("rotation_euler", frame=start_frame)  # (44.0 * DEGREE_TO_RAD, 0, 0)

        # camera offset from z-axis
        x_0 = cam.location[0]
        y_0 = cam.location[1]
        z_0 = cam.location[2]

        offset_distance = ((x_0-x_center)**2 + (y_0-y_center)**2)**0.5
        angle_0 = atan2(y_0-y_center,x_0-x_center)

        # adjust clip range since rotations will be done after dive-in in more close-up views
        # to check if we zoomed in, we will use the camera's z-position which for an overview is set to +4
        if z_0 < 2.0:
            cam.data.clip_start = 0.01
            #cam.data.clip_end = 10.0

        # Set keyframes for rotation
        for frame in range(0,frames+1):
            angle = angle_0 + radians(frame * (rotation_degrees / frames))

            # Calculate camera position around the Z-axis
            x = offset_distance * cos(angle)
            y = offset_distance * sin(angle)

            # rotate around center location
            x += x_center
            y += y_center

            # new camera position
            camera_location = Vector((x, y, z_0))  # Maintain Z position

            # Set camera rotation towards the origin
            look_at = Vector((x_center, y_center, z_center + z_elevation))  # Origin point on mesh

            # gets camera rotation angle
            direction = look_at - camera_location
            rot_quat = direction.to_track_quat('-Z', 'Y')  # Point -Z axis to direction
            camera_rotation = rot_quat.to_euler()

            # there seems to be a problem with frame interpolation when the angles going past 2 PI, i.e,
            # having the angle set between [-pi,pi] from the to_euler() output and interpolating a step
            # between +pi to -pi. to avoid we always increase the angle going from [0,2pi[
            #
            # in our case here, we change & increase only the third euler angle, while the first 2 stay fixed
            # when the camera rotates around the z-axis.
            if camera_rotation[2] < 0.0: camera_rotation[2] += 2.0 * PI

            # checks difference in current camera rotation and new target rotation
            # makes sure that z-angle always increases from one frame to another,
            # otherwise the linear interpolation will between frames will lead to weird results
            if cam.rotation_euler[2] > camera_rotation[2]: camera_rotation[2] += 2.0 * PI

            # info
            if frame == 0:
                print("    initial frame: camera angles = {:.2f} / {:.2f} / {:.2f} ".format(np.rad2deg(camera_rotation.x),
                                np.rad2deg(camera_rotation.y),np.rad2deg(camera_rotation.z)))
                print("")

            frame_number = start_frame + frame * keyframe_interval
            scene.frame_set(frame_number)

            # Set camera location and keyframe for each frame
            cam.rotation_euler = camera_rotation
            cam.location = camera_location

            cam.keyframe_insert("location", frame=frame_number)
            cam.keyframe_insert("rotation_euler", frame=frame_number)

            #debug
            #print("frame: ",frame,frame_number,"angle",angle,cam.location,"rotation",cam.rotation_euler)

    print("  total number of keyframes = ",number_of_keyframes)
    print("  total number of frames    = ",total_frames)
    print("")

    # smoothing move
    #for fcurve in cam.animation_data.action.fcurves:
    #    for keyframe_point in fcurve.keyframe_points:
    #        keyframe_point.interpolation = 'LINEAR'

    # render settings
    # see: https://docs.blender.org/api/current/bpy.types.FFmpegSettings.html#bpy.types.FFmpegSettings
    # ffmpeg
    scene.render.image_settings.file_format = 'FFMPEG'   # 'AVI_JPEG', 'FFMPEG'
    # codec
    scene.render.ffmpeg.codec = 'H264'                  # compression: 'MPEG4', 'H264'
    #scene.render.ffmpeg.constant_rate_factor = 'HIGH'
    scene.render.ffmpeg.format = 'MPEG4'                # output file container format
    # frames per second (blender default is 24)
    scene.render.fps = 20

    # user output
    print("  movie fps                    = ",scene.render.fps)
    print("  movie format                 = ",scene.render.ffmpeg.format )
    print("")

    # output movie
    name = './out.anim'
    filename = name + '.mp4'
    scene.render.filepath = filename
    scene.render.use_file_extension = False

    # timing
    tic = time.perf_counter()

    # redirect stdout to null
    print("rendering animation: {} ...".format(filename))
    print("")
    # to avoid long stdout output by the renderer:
    #    Fra:1 Mem:189.14M (Peak 190.26M) | Time:00:00.68 | Syncing Sun
    #    Fra:1 Mem:189.14M (Peak 190.26M) | Time:00:00.68 | Syncing Camera
    #  ..
    suppress = suppress_renderer_output
    with SuppressStream(sys.stdout,suppress):
        # render animation
        bpy.ops.render.render(animation=True)
    print("  written to: ",filename)
    print("")

    # timing
    toc = time.perf_counter()
    if toc - tic < 100.0:
        print("elapsed time for animation render is {:0.4f} seconds\n".format(toc - tic))
    else:
        min = int((toc-tic)/60.0)
        sec = (toc - tic) - min * 60
        print("elapsed time for animation render is {} min {:0.4f} sec\n".format(min,sec))


def render_image():
    """
    renders a single image
    """
    global suppress_renderer_output

    print("image:")

    # gets scene
    scene = bpy.context.scene

    # output image
    scene.render.image_settings.file_format = 'JPEG'   # 'PNG'
    name = './out'
    filename = name + '.jpg'
    scene.render.filepath = filename
    scene.render.use_file_extension = False

    # timing
    tic = time.perf_counter()

    # redirect stdout to null
    print("  rendering image: {} ...".format(filename))
    print("")
    # to avoid long stdout output by the renderer:
    #    Fra:1 Mem:462.85M (Peak 462.85M) | Time:00:00.04 | ..
    #  ..
    suppress = suppress_renderer_output
    with SuppressStream(sys.stdout,suppress):
        # Render Scene and store the scene
        bpy.ops.render.render(write_still=True)
    print("  written to: ",filename)
    print("")

    # timing
    toc = time.perf_counter()
    if toc - tic < 100.0:
        print("elapsed time for image render is {:0.4f} seconds\n".format(toc - tic))
    else:
        min = int((toc-tic)/60.0)
        sec = (toc - tic) - min * 60
        print("elapsed time for image render is {} min {:0.4f} sec\n".format(min,sec))


def render_blender_scene(title="",animation=False,close_up_view=False):

    ## blender scene setup
    print("Setting up blender scene...")
    print("")

    ## camera
    add_camera(title,close_up_view)

    ## light
    add_light()

    ## background plane
    add_plane(close_up_view)

    # text object
    add_title(title)

    # save scene and render options
    set_scene()

    # rendering
    if animation:
      # movie rendering
      render_animation()
    else:
      # image rendering
      render_image()

    # save blend file
    dir = os.getcwd()
    name = 'out.blend'

    filename = dir + "/" + name
    bpy.ops.wm.save_as_mainfile(filepath=filename)

    print("")
    print("  saved blend file: ",filename)
    print("")


# main routine
def plot_with_blender(vtk_file="",image_title="",colormap=0,color_max=None,buildings_file="",animation=False,close_up_view=False):
    """
    renders image for (earth) sphere with textures
    """
    # set current directory, in case we need it to load files
    dir = os.getcwd()
    print("current directory: ",dir)
    print("")

    # converts .vtu to .obj file for blender to read in
    obj_file = convert_vtk_to_obj(vtk_file,colormap,color_max)

    # setup mesh node with shaders
    create_blender_setup(obj_file)

    # add buildings
    add_blender_buildings(buildings_file)

    # save blender scene
    render_blender_scene(image_title,animation,close_up_view)


def usage():
    print("usage: ./plot_with_blender.py [--vtk_file=file] [--title=my_mesh_name] [--colormap=val] [--color-max=val] [--buildings=file]")
    print("                              [--with-cycles/--no-cycles] [--closeup] [--small] [--anim] [--suppress] [--help]")
    print("  with")
    print("     --vtk_file                - input mesh file (.vtk, .vtu, .inp)")
    print("     --title                   - title text (added to image rendering)")
    print("     --colormap                - color map type: 0==VTK        / 1==topo      / 2==lisbon    / 3==lajolla       / 4==lipari")
    print("                                                 5==davos      / 6==turku     / 7==berlin    / 8==grayC         / 9==snow")
    print("                                                10==shakeGreen /11==shakeRed  /12==shakeUSGS /13==shakeUSGSgray")
    print("                                                 (default is shakeUSGSgray for shakemaps)")
    print("     --color-max               - fixes maximum value of colormap for moviedata to val, e.g., 1.e-7)")
    print("     --buildings               - mesh file (.ply) with buildings to visualize for the area")
    print("     --with-cycles/--no-cycles - turns on/off CYCLES renderer (default is off, using BLENDER_EEVEE)")
    print("     --closeup                 - sets camera view closer to center of model")
    print("     --transparent-sea-level   - turns on transparency for sea-level plane")
    print("     --small                   - turns on small images size (400x600px) for preview")
    print("     --shift-buildings         - moves buildings up, e.g., a factor 0.0005 (default is no shift==0.0)")
    print("     --anim                    - turns on movie animation (dive-in and rotation by default, use --no-rotation or --no-dive-in to turn off)")
    print("     --suppress                - suppress renderer output (default is off)")
    print("     --help                    - this help for usage...")
    sys.exit(1)


if __name__ == '__main__':
    # init
    vtk_file = ""
    image_title = ""
    color_max = None
    colormap = -1
    buildings_file = ""
    animation = False
    close_up_view = False

    # reads arguments
    #print("\nnumber of arguments: " + str(len(sys.argv)))
    i = 0
    for arg in sys.argv:
        i += 1
        #print("argument "+str(i)+": " + arg)
        # get arguments
        if "--help" in arg:
            usage()
        elif "--anim" in arg:
            animation = True
        elif "--buildings=" in arg:
            buildings_file = arg.split('=')[1]
        elif "--closeup" in arg:
            close_up_view = True
        elif "--color-max" in arg:
            color_max = float(arg.split('=')[1])
        elif "--colormap" in arg:
            colormap = int(arg.split('=')[1])
        elif "--with-cycles" in arg:
            use_cycles_renderer = True
        elif "--no-cycles" in arg:
            use_cycles_renderer = False
        elif "--no-dive-in" in arg:
            use_animation_dive_in = False
        elif "--no-rotation" in arg:
            use_animation_rotation = False
        elif "--small" in arg:
            blender_img_resolution_X = 600
            blender_img_resolution_Y = 400
        elif "--shift-buildings" in arg:
            shift_building_baseline = float(arg.split('=')[1])
        elif "--suppress" in arg:
            suppress_renderer_output = True
        elif "--title=" in arg:
            image_title = arg.split('=')[1]
        elif "--transparent-sea-level" in arg:
            use_transparent_sea_level_plane = True
        elif "--vtk_file=" in arg:
            vtk_file = arg.split('=')[1]
        elif i >= 8:
            print("argument not recognized: ",arg)

    # sets default colormap
    if colormap == -1:
        # for shakemaps
        if 'shaking' in vtk_file or 'shakemap' in vtk_file:
            colormap = 13 # shakeUSGSgray
        else:
            colormap = 0  # VTK diverging red-blue

    # logging
    cmd = " ".join(sys.argv)
    filename = './plot_with_blender.log'
    with open(filename, 'a') as f:
      print("command call --- " + str(datetime.datetime.now()),file=f)
      print(cmd,file=f)
      print("command logged to file: " + filename)

    # main routine
    plot_with_blender(vtk_file,image_title,colormap,color_max,buildings_file,animation,close_up_view)
