#!/usr/bin/env python
#
# converts DXF file to a mesh file (.ply) readable by Blender
#
# this script is mostly meant to download SwissTopo buildings for some areas in Switzerland (--mesh-area=..).
# additionally, it converts the SwissTopo LV95 coordinates to UTM coordinates that are used by SPECFEM3D.
#
############################################################################################

import sys
import os
import glob
import datetime

try:
    import numpy as np
except:
    print("Failed importing numpy. Please make sure to have python with numpy working properly.")
    sys.exit(1)

# DXF/DWG reader
try:
    import ezdxf
except:
    print("Failed importing exdxf.")
    sys.exit(1)

from ezdxf.addons.drawing import matplotlib
from ezdxf.addons import meshex
from ezdxf.render import MeshBuilder
from ezdxf.render import MeshVertexMerger
from ezdxf.render import MeshTransformer

# Coordinate projections
try:
    import pyproj
except:
    print("Failed importing pyproj.")
    sys.exit(1)

from pyproj import Transformer

# Url download requests
import json
import urllib
from zipfile import ZipFile


###############################################################################################
## USER PARAMETERS

# downloads SwissTopo buildings
use_swisstopo_buildings = True

# converts SwissTopo LV95 coordinates to UTM
use_coordinate_transform_LV95_to_UTM = True

###############################################################################################

def download_mesh_area_buildings(mesh_area=None):
    global use_swisstopo_buildings
    # downloads building tiles for specified area
    # checks if anything to
    if mesh_area is None: return []
    if not use_swisstopo_buildings: return []

    # user output
    print("downloading buildings from SwissTopo...")
    print("  mesh area: ",mesh_area)
    print("")

    # download building through swisstopo data api
    url_base = 'https://data.geo.admin.ch/api/stac/v0.9/collections/ch.swisstopo.swissbuildings3d_2'

    # mesh area min/max
    lon_min = mesh_area[0]
    lat_min = mesh_area[1]
    lon_max = mesh_area[2]
    lat_max = mesh_area[3]

    # check
    if lon_min > lon_max or lat_min > lat_max:
        print("Wrong mesh area format, please provide in a format (lon_min,lat_min,lon_max,lat_max)")
        sys.exit(1)

    print("  lat min/max: ",lat_min,lat_max)
    print("  lon min/max: ",lon_min,lon_max)
    print("")

    # bounding box with lon/lat format
    suffix = "/items?bbox={},{},{},{}".format(lon_min,lat_min,lon_max,lat_max)

    # get list of tiles
    url = url_base + suffix
    f = urllib.request.urlopen(url)

    txt = f.read().decode('utf-8')
    json_result = json.loads(txt)

    #print("json: ",json_result)

    # gets tile urls
    tile_list = []
    for item in json_result['features']:
        for k,dic in item['assets'].items():
            href = dic['href']
            # keep .dxf files
            # for example: SWISSBUILDINGS3D_2_0_CHLV95LN02_1047-32.dxf.zip
            #print("href: ",href)
            if href[-8:] == '.dxf.zip':
                # makes sure the tile names have a format like: swissbuildings3d_2_2019-04_1067-44_2056_5728.dxf.zip
                # contains 6 items when split by '_', like: ['swissbuildings3d', '2', '2019-04', '1067-44', '2056', '5728.dxf.zip']
                name_items = dic['href'].split('/')[-1].split('_')
                #print("name: ",name_items)
                if len(name_items) == 6:
                    tile_list.append(dic['href'])

    num_tiles = len(tile_list)
    print("  data tiles: number of tiles found = ",num_tiles)
    print("")

    if num_tiles == 0:
        print("No data tiles found, please check...")
        sys.exit(1)

    # download files
    datafiles = []
    for i,url in enumerate(tile_list):
        name = url.split('/')[-1]
        # zip name
        # example: ./swissbuildings3d_2_2019-04_1067-44_2056_5728.dxf.zip
        filename_zip = os.path.join("./",name)

        # check if file already downloaded
        # file name
        # example: SWISSBUILDINGS3D_2_0_CHLV95LN02_1067-44.dxf
        pref = 'SWISSBUILDINGS3D_2_0_CHLV95LN02_'
        suf = '.dxf'
        # split: swissbuildings3d _ 2 _ 2019-04 _ 1067-44 _ 2056 _ 5728.dxf.zip
        # to get 1067-44 identifier
        name_dxf = pref + name.split("_")[3] + suf
        filename = os.path.join("./",name_dxf)

        print("  tile {} out of {}".format(i+1,len(tile_list)))
        print("  zip name: ",filename_zip)
        print("  name    : ",filename)

        if not os.path.isfile(filename):
            # download file
            print("  downloading...")
            #print("  url        : ",url)
            x = urllib.request.urlopen(url)
            with open(filename_zip,'wb') as f:
                f.write(x.read())

            # de-zip file
            with ZipFile(filename_zip, 'r') as zipObj:
                # Extract all the contents of zip file in current directory
                zipObj.extractall("./")

            # remove zip file
            os.remove(filename_zip)
        else:
            print("  file exists...")
        print("")

        datafiles.append(filename)

    print("downloaded files: ")
    for file in datafiles: print("  ",file)
    print("")

    return datafiles


def read_dxf_file(filename,index):
    # reads in DXF file
    print("  reading in DXF file: ",filename)
    print("")

    # Read the DXF file using ezdxf
    doc = ezdxf.readfile(filename)

    print(f"  DXF version: {doc.dxfversion}")
    print(f"  DXF database contains {len(doc.entitydb)} entities")
    print("")

    # Recommended: audit & repair DXF document before rendering
    auditor = doc.audit()

    if auditor.has_errors:
        auditor.print_error_report(auditor.errors)
        raise exception("Error: The DXF document is damaged and can't be converted!")
        sys.exit(1)
    else:
        print("  auditor is ok")
        print("")

    # image plot (pretty slow...)
    plot_image = False
    if plot_image:
        print("  plotting as image png...")
        name = './' + "tmp_{}.png".format(index)
        matplotlib.qsave(doc.modelspace(), name)
        print("  written to: ",name)
        print("")

    # model space info
    msp = doc.modelspace()

    # query all entity types available in the model space
    entity_types = []
    for e in msp:
        if not e.dxftype() in entity_types:
            entity_types.append(e.dxftype())

    print("  model space entity types:")
    print("  ",entity_types)
    print("")

    # query lines
    # we will need polyfaces to convert 3D objects
    num_lines = len(msp.query("LINE"))
    num_polylines = len(msp.query("POLYLINE"))
    num_lwpolylines = len(msp.query("LWPOLYLINE"))
    #num_polyfaces = len(msp.query("POLYFACE"))     # swiss buildings dxf contains polylines as polyface meshes

    polyfaces = [polyline for polyline in msp.query('POLYLINE') if polyline.is_poly_face_mesh]
    num_polyfaces = len(polyfaces)

    print("  modelspace:")
    print("    lines                     : ",num_lines)
    print("    lightweight polylines     : ",num_lwpolylines)
    print("    polylines                 : ",num_polylines)
    print("    polyfaces (from polylines): ",num_polyfaces)
    print("")

    # check
    if num_polyfaces == 0:
        print("Error: no polyfaces found in file, exiting...")
        sys.exit(1)

    return polyfaces


def get_minmax_coordinates_in_LV95(mesh_area):
    # pyproj coordinate system info:
    #   SwissTopo LV95                                 == CRS / EPSG:2056
    #   WGS84                                          ==       EPSG:4326
    #
    # WGS84 (GPS) lat/lon -> LV95 (SwissTopo)
    transformer_to_lv95 = Transformer.from_crs("EPSG:4326","EPSG:2056")

    # mesh area min/max
    lon_min = mesh_area[0]
    lat_min = mesh_area[1]
    lon_max = mesh_area[2]
    lat_max = mesh_area[3]

    x_min,y_min = transformer_to_lv95.transform(lat_min,lon_min)
    x_max,y_max = transformer_to_lv95.transform(lat_max,lon_max)

    # user output
    print("  mesh area:")
    print("     lat min/max      = ",lat_min,lat_max)
    print("     lon min/max      = ",lon_min,lon_max)
    print("  -> LV95 x min/max   = ",x_min,x_max)
    print("          y min/max   = ",y_min,y_max)
    print("")

    return x_min,x_max,y_min,y_max


def convert_coordinates_LV95_to_UTM(mesh):
    # user output
    print("converting LV95 coordinates to UTM...")
    print("")
    # pyproj coordinate system info:
    #   SwissTopo LV95                                 == CRS / EPSG:2056
    #   WGS84                                          ==       EPSG:4326
    #   spherical mercator, google maps, openstreetmap ==       EPSG:3857
    #
    # we first need to determine the correct UTM zone to get the EPSG code.
    # for this, we take the mesh origin point and convert it to lat/lon (GPS) position.
    # we can then query the corresponding UTM zone for this position.

    # LV95 (SwissTopo) -> WGS84 (GPS) lat/lon
    transformer_to_latlon = Transformer.from_crs("EPSG:2056", "EPSG:4326")

    # mesh origin position
    orig_x = mesh.diagnose().bbox.extmin.x + 0.5 * mesh.diagnose().bbox.size.x
    orig_y = mesh.diagnose().bbox.extmin.y + 0.5 * mesh.diagnose().bbox.size.y

    orig_lat,orig_lon = transformer_to_latlon.transform(orig_x,orig_y)

    # user output
    print("  origin: x/y      = ",orig_x,orig_y)
    print("       -> lat/lon  = ",orig_lat,orig_lon)
    print("")

    # gets list of UTM codes
    utm_crs_list = pyproj.database.query_utm_crs_info(
                      datum_name="WGS 84",
                      area_of_interest=pyproj.aoi.AreaOfInterest(west_lon_degree=orig_lon,
                                                                 south_lat_degree=orig_lat,
                                                                 east_lon_degree=orig_lon,
                                                                 north_lat_degree=orig_lat))
    utm_code = utm_crs_list[0].code
    utm_epsg = "EPSG:{}".format(utm_code)

    print("  UTM code:", utm_epsg)

    # direct transformation from LV95 to UTM zone
    transformer_to_utm = Transformer.from_crs("EPSG:2056", utm_epsg)

    #debug
    #print(transformer_to_utm)

    # user info
    utm_x,utm_y = transformer_to_utm.transform(orig_x,orig_y)
    print("       -> utm x/y  = ",utm_x,utm_y)
    print("          backward check: orig x/y = ",transformer_to_utm.transform(utm_x,utm_y,direction='INVERSE'))
    print("")

    # converts mesh point coordinates
    vertices = mesh.vertices
    num_vertices = len(vertices)
    print("  mesh:")
    print("    number of vertices = ",num_vertices)
    print("")

    for i,vertex in enumerate(mesh.vertices):
        # gets LV95 location
        x = vertex.x
        y = vertex.y
        z = vertex.z

        # converts to UTM location (x and y)
        #point = np.array([x,y])
        x_utm,y_utm = transformer_to_utm.transform(x,y)

        # sets new coordinates
        vertex_new = ezdxf.math.Vec3(x_utm,y_utm,z)
        mesh.vertices[i] = vertex_new

    print("  new UTM mesh dimensions:")
    print("  bounding box :  ",mesh.diagnose().bbox)
    orig_x = mesh.diagnose().bbox.extmin.x + 0.5 * mesh.diagnose().bbox.size.x
    orig_y = mesh.diagnose().bbox.extmin.y + 0.5 * mesh.diagnose().bbox.size.y
    print("  origin X/Y   : ",orig_x,orig_y)
    dim_x = mesh.diagnose().bbox.size.x
    dim_y = mesh.diagnose().bbox.size.y
    dim_z = mesh.diagnose().bbox.size.z
    print("  dimensions   : ",dim_x,dim_y,dim_z)
    print("")

    return mesh


def convert_autocad_dxf_to_mesh(input_file,mesh_origin=None,mesh_scale=None,mesh_area=None):
    global use_coordinate_transform_LV95_to_UTM

    # user output
    print("")
    print("convert AutoCAD DXF file(s)")
    print("  input file  : ",input_file)
    if not mesh_origin is None:
        print("  mesh_origin : ",mesh_origin)
    if not mesh_scale is None:
        print("  mesh_scale  : ",mesh_scale)
    if not mesh_area is None:
        print("  mesh_area   : ",mesh_area)

    print("")

    # current directory
    dir = os.getcwd()
    print("  current directory: ",dir)
    print("")

    # input data
    if len(input_file) == 0:
        # download buildings
        datafiles = download_mesh_area_buildings(mesh_area)
    else:
        # check if multiple files specified with wildcard (e.g., SWISSBUILDINGS3D*.dxf)
        datafiles = []
        if '*' in input_file:
            print("  files expanded:")
            files = glob.glob(input_file)
            files.sort()
            for filename in files:
                if ".dae" in filename or ".shp" in filename or ".dbf" in filename or ".gdb" in filename:
                    # skip file name (e.g. SWISSBUILDINGS3D_2_0_CHLV95LN02_1047-32.dae)
                    continue
                else:
                    # keep file
                    print("  ",filename)
                    datafiles.append(filename)
            print("")
        else:
            datafiles.append(input_file)

    print("  number of data files: ",len(datafiles))
    print("")
    if len(datafiles) == 0:
        print("no data files found, please check input...")
        sys.exit(1)

    # checks if file(s) exists
    for filename in datafiles:
        # checks if file exists
        if not os.path.isfile(filename):
            print("Please check if file exists: ",filename)
            sys.exit(1)

    # to check if mesh buildings are in specified mesh_area
    if not mesh_area is None:
        # convert lat/lon from mesh_area to SwissTopo LV95
        area_x_min,area_x_max,area_y_min,area_y_max = get_minmax_coordinates_in_LV95(mesh_area)

    # overall document layout
    doc = ezdxf.new()
    msp = doc.modelspace()

    # create a new mesh object
    mesh = MeshBuilder()

    # reads in DXF data
    icount_file = 0
    for filename in datafiles:
        # reads in model file
        icount_file += 1
        print("")
        print("***************************************************")
        print("file {} out of {}".format(icount_file,len(datafiles)))
        print("")
        print("data file: ",filename)
        print("***************************************************")
        print("")

        polyfaces = read_dxf_file(filename,icount_file)

        # adds polyface meshes
        print("adding meshes...")
        icount_added = 0
        for i,polyface in enumerate(polyfaces):
            layer = polyface.dxf.layer
            color = polyface.dxf.color

            # user output
            # number to output in steps of 10/100/1000 depending on how large the number of polyfaces is
            nout = min(1000,int(10**np.floor(np.log10(len(polyfaces)))))
            if (i+1) % nout == 0:
                print("  polyface: {} out of {} - layer: {}".format(i+1,len(polyfaces),layer))  # " - color: ",color

            # creates object mesh
            b = MeshVertexMerger.from_polyface(polyface)
            #b.render(msp1, dxfattribs={
            #    'layer': polyface.dxf.layer,
            #    'color': polyface.dxf.color,
            #})
            b.render_polyface(msp, dxfattribs={
                'layer': layer,
                'color': color,
            })

            # checks number of vertices in mesh
            if len(b.vertices) == 0: continue

            # checks if mesh in mesh_area
            if not mesh_area is None:
                # gets bounding box min/max from this building
                bbox = b.diagnose().bbox
                x_min = bbox.extmin.x
                x_max = bbox.extmax.x
                y_min = bbox.extmin.y
                y_max = bbox.extmax.y

                if x_min < area_x_min or x_max > area_x_max or \
                   y_min < area_y_min or y_max > area_y_max:
                    # building outside of mesh area
                    #print("building: ",x_min,x_max,"/",y_min,y_max,"area: ",area_x_min,area_x_max,"/",area_y_min,area_y_max)
                    # we won't add it to the total mesh object, continue with next polyface
                    continue

            # removes faces with null normals (degenerate faces),
            # otherwise will lead to problems/crashes when saving as .stl and reading the .ply file
            faces = []
            for j,face in enumerate(b.faces):
                #print("face ",j,face)
                normal = b.get_face_normal(j)
                if normal == ezdxf.math.NULLVEC:
                    #debug
                    #print("mesh problem w/ normal from get_face_normal() ",j,normal,face)
                    # we will omit this face
                    continue
                else:
                    faces.append(face)

            # re-sets mesh faces array
            b.faces = faces

            # checks number of new faces
            if len(b.faces) == 0: continue

            # checks again face normals to avoid problems with degenerate faces when saving to file
            do_add_mesh = True
            for j,normal in enumerate(b.diagnose().face_normals):
                if normal == ezdxf.math.NULLVEC:
                    #debug
                    #print("mesh problem w/ face normal ",j,normal)
                    # we will skip this mesh
                    do_add_mesh = False
                    break

            # adds object to mesh
            if do_add_mesh:
                mesh.add_mesh(mesh=b)
                # counter
                icount_added += 1

            #debug
            #if i == 2: break
        print("  done. added {} out of {}".format(icount_added,len(polyfaces)))

    # check
    if len(mesh.vertices) == 0:
        print("")
        print("Mesh is empty, there's nothing to save. Please check your inputs...")
        print("")
        sys.exit(1)

    # user output
    print("")
    print("mesh dimensions:")
    print("  bounding box :  ",mesh.diagnose().bbox)
    orig_x = mesh.diagnose().bbox.extmin.x + 0.5 * mesh.diagnose().bbox.size.x
    orig_y = mesh.diagnose().bbox.extmin.y + 0.5 * mesh.diagnose().bbox.size.y
    print("  origin X/Y   : ",orig_x,orig_y)
    dim_x = mesh.diagnose().bbox.size.x
    dim_y = mesh.diagnose().bbox.size.y
    dim_z = mesh.diagnose().bbox.size.z
    print("  dimensions   : ",dim_x,dim_y,dim_z)
    print("")

    # transform mesh coordinate from SwissTopo to UTM
    if use_coordinate_transform_LV95_to_UTM:
        mesh = convert_coordinates_LV95_to_UTM(mesh)

        ## write out UTM mesh as binary .ply file
        name = 'output_dxf_utm.ply'
        filename = './' + name
        with open(filename, "wb") as fp:
            fp.write(meshex.ply_dumpb(mesh))
        print("  UTM mesh written: ",filename)
        print("")

    # moves mesh coordinates by origin point to position around (0,0,0) for blender
    if not mesh_origin is None:
        print("mesh: moving by origin  = ",mesh_origin)
        print("")

        translation_vector = - ezdxf.math.Vec3(mesh_origin)

        mesh2 = MeshTransformer.from_builder(mesh)
        mesh2.translate(translation_vector)

        mesh = mesh2

        print("      bounding box after : ",mesh.diagnose().bbox)
        orig_x = mesh.diagnose().bbox.extmin.x + 0.5 * mesh.diagnose().bbox.size.x
        orig_y = mesh.diagnose().bbox.extmin.y + 0.5 * mesh.diagnose().bbox.size.y
        print("      origin X/Y after   : ",orig_x,orig_y)
        dim_x = mesh.diagnose().bbox.size.x
        dim_y = mesh.diagnose().bbox.size.y
        dim_z = mesh.diagnose().bbox.size.z
        print("      dimensions after   : ",dim_x,dim_y,dim_z)
        print("")

    # scales mesh coordinates by uniform scaling factor to scale between +/- 1 for blender
    if not mesh_scale is None:
        print("mesh: scaling by factor = ",mesh_scale)
        print("")

        scale_factor = mesh_scale

        mesh2 = MeshTransformer.from_builder(mesh)
        mesh2.scale_uniform(scale_factor)

        mesh = mesh2

        print("      bounding box after : ",mesh.diagnose().bbox)
        dim_x = mesh.diagnose().bbox.size.x
        dim_y = mesh.diagnose().bbox.size.y
        dim_z = mesh.diagnose().bbox.size.z
        print("      dimensions after   : ",dim_x,dim_y,dim_z)
        print("")

    ## write out ASCII .stl file
    #name = 'output_dxf.stl'
    #filename = './' + name
    #with open(filename, "wt") as fp:
    #    fp.write(meshex.stl_dumps(mesh))

    ## write out binary .stl file
    #name = 'output_dxf.stl'
    #filename = './' + name
    #with open(filename, "wb") as fp:
    #    fp.write(meshex.stl_dumpb(mesh))

    ## write out ASCII .obj file
    #name = 'output_dxf.obj'
    #filename = './' + name
    #with open(filename, "wt") as fp:
    #    fp.write(meshex.obj_dumps(mesh))

    ## write out binary .ply file
    name = 'output_dxf.ply'
    filename = './' + name
    with open(filename, "wb") as fp:
        fp.write(meshex.ply_dumpb(mesh))
    print("written: ",filename)
    print("")



def usage():
    print("usage: ./convert_dxf_to_mesh.py [--dxf_file=file] [--mesh_origin=(x0,y0,z0)] [--mesh_scale=factor] [--mesh_area=(lon_min,lat_min,lon_max,lat_max)")
    print("  with")
    print("     --dxf_file              - input mesh file (.dxf, .dwg)")
    print("     --mesh_origin           - translate mesh by origin position")
    print("     --mesh_scale            - scale mesh by a factor")
    print("     --mesh_area             - downloads/limits buildings for area specified by (lon_min,lat_min,lon_max,lat_max)")
    sys.exit(1)


if __name__ == '__main__':
    # init
    input_file = ""
    mesh_origin = None
    mesh_scale = None
    mesh_area = None

    # reads arguments
    #print("\nnumber of arguments: " + str(len(sys.argv)))
    if len(sys.argv) <= 1: usage()

    i = 0
    for arg in sys.argv:
        i += 1
        #print("argument "+str(i)+": " + arg)
        # get arguments
        if "--help" in arg:
            usage()
        elif "--dxf_file=" in arg:
            input_file = arg.split('=')[1]
        elif "--dwg_file=" in arg:
            input_file = arg.split('=')[1]
        elif "--mesh_origin=" in arg:
            str_array = arg.split('=')[1]
            mesh_origin = np.array([float(val) for val in str_array.strip('()[]').split(',')])
        elif "--mesh_scale=" in arg:
            mesh_scale = float(arg.split('=')[1])
        elif "--mesh_area=" in arg:
            str_array = arg.split('=')[1]
            mesh_area = np.array([float(val) for val in str_array.strip('()[]').split(',')])
        elif i >= 6:
            print("argument not recognized: ",arg)

    # logging
    cmd = " ".join(sys.argv)
    filename = './convert_dxf_to_mesh.log'
    with open(filename, 'a') as f:
      print("command call --- " + str(datetime.datetime.now()),file=f)
      print(cmd,file=f)
      print("command logged to file: " + filename)

    # main routine
    convert_autocad_dxf_to_mesh(input_file,mesh_origin,mesh_scale,mesh_area)
