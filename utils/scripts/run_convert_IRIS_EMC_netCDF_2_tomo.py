#!/usr/bin/env python
#
# converts an IRIS EMC netCDF file to a tomography_model.xyz file readable by SPECFEM3D
#
# this script is loosely based on IRIS EMC tool netCDF_2_GeoCSV_3D.py,
# a Python script to read a 3D netCDF Earth model file and convert it to GeoCSV format;
# provided in: https://github.com/EarthScope/emc-tools
#
# note: in this current script version, it will only produce valid tomography model files for SPECFEM3D
#       if the EMC model uses dimensions named "latitude","longitude","depth" and has regular spaced gridding.
#       this script was tested using the EMC-FWEA18 model (https://ds.iris.edu/ds/products/emc-fwea18/) and
#       might not work properly yet with other models - something to try out in future.
#
# python modules required by this script:
# - netCDF4
# - numpy
# - pygmt
# - vtk
#
# usage example:
#
# 0. download the EMC model netCDF file:
#   for example, the model FWEA18_kmps.nc from https://ds.iris.edu/ds/products/emc-fwea18/
#   and store it say in a subfolder `IRIS_EMC/FWEA18_kmps.nc`
#
# 1. create the `tomography_model.xyz` of the IRIS EMC model file by typing:
#    > ./run_convert_IRIS_EMC_netCDF_2_tomo.py --EMC_file=IRIS_EMC/FWEA18_kmps.nc --mesh_area=136.0,36.5,138.66,38.66 --maximum_depth=80.0
#
#    and store the tomography model file `tomography_model.xyz` into the folder DATA/tomo_files/
#
# Then, to run a simulation with this tomo file, make sure the DATA/meshfem3D_files/Mesh_Par_file is setup with the target region.
# For the Noto peninsula of the Ishikawa prefecture of Japan, between lon/lat [136.0,36.5,138.66,38.66] and a depth of 80 km,
# the Mesh_Par_file for the region should have set:
#        ..
#
#        LATITUDE_MIN                    = 36.5
#        LATITUDE_MAX                    = 38.66
#        LONGITUDE_MIN                   = 136.0
#        LONGITUDE_MAX                   = 138.66
#        DEPTH_BLOCK_KM                  = 80.d0
#        UTM_PROJECTION_ZONE             = 53
#        SUPPRESS_UTM_PROJECTION         = .false.
#
#        ..
#
#        # number of materials
#        NMATERIALS                      = 1
#        # define the different materials in the model as :
#        # #material_id  #rho  #vp  #vs  #Qkappa #Qmu  #anisotropy_flag #domain_id
#        #     Q                : quality factor
#        #     anisotropy_flag  : 0=no anisotropy/ 1,2,.. check with implementation in aniso_model.f90
#        #     domain_id        : 1=acoustic / 2=elastic
#        # for tomography:
#        #   #material_id #type-keyword #domain-name #tomo-filename #tomo_id #domain_id
#        #   example:
#        #   -1 tomography elastic tomography_model.xyz 0 2
#        #
#        -1 tomography elastic tomography_model.xyz 0 2
#
#        ..
#
# as a helper script to get the topography interface setup for this region, one can use the script `run_get_simulation_topography.py`
# in folder utils/scripts/:
#    > ./run_get_simulation_topography.py 136.0 36.5 138.66 38.66 topo3
#
# Please consider improving this script by contributing to this package!
#
#########################################################################################

import sys
import os
import datetime
import math

# netCDF
try:
    import netCDF4 as nc
except:
    print("Failed importing netCDF4. Please make sure to have python with netCDF working properly.")
    sys.exit(1)

# numpy
try:
    import numpy as np
except:
    print("Failed importing numpy. Please make sure to have python with numpy working properly.")
    sys.exit(1)


# GMT python interface
try:
    import pygmt
except:
    print("Error importing module `pygmt`")
    print("install by: > pip install -U pygmt")
    sys.exit(1)

## elevation package
# not supported properly yet...
# see: https://github.com/bopen/elevation
#      http://elevation.bopen.eu/en/stable/
# install by: > sudo pip install elevation
#try:
#    import elevation
#except:
#    print("Error importing module `elevation`")
#    print("install by: > pip install -U elevation")
#    sys.exit(1)
## elevation
# check
#elevation.util.selfcheck(elevation.TOOLS)

## UTM zones
# could be used to convert to UTM X/Y instead of geo2utm() routine;
# also, could make sure that given input UTM_zone makes sense for model region
#
# see: https://github.com/Turbo87/utm
# install by: > sudo pip install utm
#try:
#    import utm
#except:
#    print("Error importing module `utm`")
#    print("install by: > pip install -U utm")
#    sys.exit(1)

## VTK visualization
try:
    import vtk
except:
    print("Failed importing vtk. Please make sure to have python with vtk working properly.")
    sys.exit(1)

############################################################################################
## Model informations
#
# SPECFEM3D tomographic model format:
#
## header infos:
#origin_x #origin_y #origin_z #end_x #end_y #end_z          - start-end dimensions
#dx #dy #dz                                                 - increments
#nx #ny #nz                                                 - number of models entries
#vp_min #vp_max #vs_min #vs_max #density_min #density_max   - min/max stats
## data records
#x1 #y1 #z1 #vp #vs #density [#Qp #Qs]                      - grid point 1
#x2 #y1 #z1 #vp #vs #density [#Qp #Qs]                      - grid point 2
#..
#
#
# IRIS EMC model:
#
# FWEA18: https://ds.iris.edu/ds/products/emc-fwea18/
#
# coordinate format  : latitude / longitude / depth (below earth surface)
# model parameters   : rho, vpv, vph, vsv, vsh, eta, qmu    - radial anisotropic
# unknown model point: missing_value = 99999.f
#
#
###########################################################################################
## PARAMETERS

# missing values in data arrays will be replaced by the average value at each depth
use_replace_missing_values_with_average = True

# surface elevation
# for example, if depth is given wrt surface elevation, use ETOPO1 elevation for each grid lat/lon point;
# however, this leads to irregular mesh increments and cannot be used for tomography_model.xyz files
use_surface_elevation = False

# create JPEG figure of surface elevation (if surface elevation is used)
create_surface_elevation_figure = True

# VTK file of model for visualization
create_vtk_file_output = True

###########################################################################################

# Global parameters
# data array sizes
ndepths = 0
nlons = 0
nlats = 0

# data array indexing
depth_index = -1
lat_index = -1
lon_index = -1

## UTM conversion
def convert_lonlat2utm(zone,lon,lat):
    """
    converts lon/lat to utm_x/utm_y for given UTM zone
    """
    # checks argument
    if zone < 1 or zone > 60:  sys.exit("error zone: zone not UTM zone")

    # converts to UTM x/y
    x,y = geo2utm(lon,lat,zone)

    return x,y

#
#----------------------------------------------------------------------------------------
#

def geo2utm(lon,lat,zone):
    """
    from utm_geo.f90

    convert geodetic longitude and latitude to UTM, and back
    use iway = ILONGLAT2UTM for long/lat to UTM, IUTM2LONGLAT for UTM to lat/long
    a list of UTM zones of the world is available at www.dmap.co.uk/utmworld.htm

    implicit none

    include "constants.h"

    -----CAMx v2.03

    UTM_GEO performs UTM to geodetic (long/lat) translation, and back.

    This is a Fortran version of the BASIC program "Transverse Mercator
    Conversion", Copyright 1986, Norman J. Berls (Stefan Musarra, 2/94)
    Based on algorithm taken from "Map Projections Used by the USGS"
    by John P. Snyder, Geological Survey Bulletin 1532, USDI.

    Input/Output arguments:

      rlon                  Longitude (deg, negative for West)
      rlat                  Latitude (deg)
      rx                    UTM easting (m)
      ry                    UTM northing (m)
      UTM_PROJECTION_ZONE  UTM zone
      iway                  Conversion type
                              ILONGLAT2UTM = geodetic to UTM
                              IUTM2LONGLAT = UTM to geodetic
    """
    PI      = math.pi
    degrad  = PI/180.0
    raddeg  = 180.0/PI
    semimaj = 6378206.4
    semimin = 6356583.8
    scfa    = 0.9996

    # some extracts about UTM:
    #
    # There are 60 longitudinal projection zones numbered 1 to 60 starting at 180 W.
    # Each of these zones is 6 degrees wide, apart from a few exceptions around Norway and Svalbard.
    # There are 20 latitudinal zones spanning the latitudes 80 S to 84 N and denoted
    # by the letters C to X, ommitting the letter O.
    # Each of these is 8 degrees south-north, apart from zone X which is 12 degrees south-north.
    #
    # To change the UTM zone and the hemisphere in which the
    # calculations are carried out, need to change the fortran code and recompile. The UTM zone is described
    # actually by the central meridian of that zone, i.e. the longitude at the midpoint of the zone, 3 degrees
    # from either zone boundary.
    # To change hemisphere need to change the "north" variable:
    #  - north=0 for northern hemisphere and
    #  - north=10000000 (10000km) for southern hemisphere. values must be in metres i.e. north=10000000.
    #
    # Note that the UTM grids are actually Mercators which
    # employ the standard UTM scale factor 0.9996 and set the
    # Easting Origin to 500,000;
    # the Northing origin in the southern
    # hemisphere is kept at 0 rather than set to 10,000,000
    # and this gives a uniform scale across the equator if the
    # normal convention of selecting the Base Latitude (origin)
    # at the equator (0 deg.) is followed.  Northings are
    # positive in the northern hemisphere and negative in the
    # southern hemisphere.
    north = 0.0
    east  = 500000.0

    IUTM2LONGLAT = 1
    ILONGLAT2UTM = 2

    # use conversion to UTM
    iway = ILONGLAT2UTM

    # zone
    UTM_PROJECTION_ZONE = zone

    # lon/lat
    rlon = lon
    rlat = lat

    rx = 0.0
    ry = 0.0

    # save original parameters
    rlon_save = rlon
    rlat_save = rlat
    rx_save = rx
    ry_save = ry

    # define parameters of reference ellipsoid
    e2 = 1.0 - (semimin/semimaj)**2
    e4 = e2 * e2
    e6 = e2 * e4
    ep2 = e2/(1.0 - e2)

    if iway == IUTM2LONGLAT:
        xx = rx
        yy = ry
    else:
        dlon = rlon
        dlat = rlat

    #----- Set Zone parameters
    zone = UTM_PROJECTION_ZONE
    cm = zone * 6.0 - 183.0       # set central meridian for this zone
    cmr = cm*degrad

    #---- Lat/Lon to UTM conversion
    if iway == ILONGLAT2UTM:

        rlon = degrad*dlon
        rlat = degrad*dlat

        delam = dlon - cm
        if delam < -180.0: delam = delam + 360.0
        if delam > 180.0: delam = delam - 360.0

        delam = delam*degrad

        f1 = (1. - e2/4.0 - 3.0*e4/64.0 - 5.0*e6/256.0)*rlat
        f2 = 3.0*e2/8.0 + 3.0*e4/32.0 + 45.0*e6/1024.0
        f2 = f2 * math.sin(2.0*rlat)
        f3 = 15.0*e4/256.0 + 45.0*e6/1024.0
        f3 = f3 * math.sin(4.0*rlat)
        f4 = 35.0*e6/3072.0
        f4 = f4 * math.sin(6.0*rlat)
        rm = semimaj*(f1 - f2 + f3 - f4)

        if dlat == 90.0 or dlat == -90.0:
            xx = 0.0
            yy = scfa*rm
        else:
            rn = semimaj/math.sqrt(1.0 - e2*math.sin(rlat)**2)
            t = math.tan(rlat)**2
            c = ep2 * math.cos(rlat)**2
            a = math.cos(rlat) * delam

            f1 = (1.0 - t + c) * a**3/6.0
            f2 = 5.0 - 18.0*t + t**2 + 72.0*c - 58.0*ep2
            f2 = f2 * a**5/120.0
            xx = scfa*rn*(a + f1 + f2)
            f1 = a**2/2.0
            f2 = 5.0 - t + 9.0*c + 4.0*c**2
            f2 = f2*a**4/24.0
            f3 = 61.0 - 58.0*t + t**2 + 600.0*c - 330.0*ep2
            f3 = f3 * a**6/720.0
            yy = scfa*(rm + rn*math.tan(rlat)*(f1 + f2 + f3))

        xx = xx + east
        yy = yy + north

    else:
        #---- UTM to Lat/Lon conversion
        xx = xx - east
        yy = yy - north
        e1 = math.sqrt(1.0 - e2)
        e1 = (1.0 - e1)/(1.0 + e1)
        rm = yy/scfa
        u = 1.0 - e2/4.0 - 3.0*e4/64.0 - 5.0*e6/256.0
        u = rm/(semimaj*u)

        f1 = 3.0*e1/2.0 - 27.0*e1**3.0/32.0
        f1 = f1*math.sin(2.0*u)
        f2 = 21.0*e1**2/16.0 - 55.0*e1**4/32.0
        f2 = f2*math.sin(4.0*u)
        f3 = 151.0*e1**3/96.0
        f3 = f3*math.sin(6.0*u)
        rlat1 = u + f1 + f2 + f3
        dlat1 = rlat1*raddeg

        if dlat1 >= 90.0 or dlat1 <= -90.0:
            dlat1 = min(dlat1,90.0)
            dlat1 = max(dlat1,-90.0)
            dlon = cm
        else:
            c1 = ep2*math.cos(rlat1)**2
            t1 = math.tan(rlat1)**2
            f1 = 1.0 - e2*math.sin(rlat1)**2
            rn1 = semimaj/math.sqrt(f1)
            r1 = semimaj*(1.0 - e2)/math.sqrt(f1**3)
            d = xx/(rn1*scfa)

            f1 = rn1*math.tan(rlat1)/r1
            f2 = d**2/2.0
            f3 = 5.0 + 3.0*t1 + 10.0*c1 - 4.0*c1**2 - 9.0*ep2
            f3 = f3*d**2*d**2/24.0
            f4 = 61.0 + 90.0*t1 + 298.0*c1 + 45.0*t1**2 - 252.0*ep2 - 3.0*c1**2
            f4 = f4*(d**2)**3/720.0
            rlat = rlat1 - f1*(f2 - f3 + f4)
            dlat = rlat*raddeg

            f1 = 1.0 + 2.0*t1 + c1
            f1 = f1*d**2*d/6.0
            f2 = 5.0 - 2.0*c1 + 28.0*t1 - 3.0*c1**2 + 8.0*ep2 + 24.0*t1**2
            f2 = f2*(d**2)**2*d/120.0
            rlon = cmr + (d - f1 + f2)/math.cos(rlat1)
            dlon = rlon*raddeg
            if dlon < -180.0: dlon = dlon + 360.0
            if dlon > 180.0: dlon = dlon - 360.0

    if iway == IUTM2LONGLAT:
        rlon = dlon
        rlat = dlat
        rx = rx_save
        ry = ry_save
        return rlon,rlat

    else:
        rx = xx
        ry = yy
        rlon = rlon_save
        rlat = rlat_save
        return rx,ry

#
#----------------------------------------------------------------------------------------
#

def display_infos(model_data):
    """
    extract and display netCDF information
    """
    # netCDF header
    print("netCDF header information:", flush=True)

    very_verbose = False

    # dimension information.
    # list of netCDF dimensions
    print("  dimensions:")
    for dim in model_data.dimensions:
        print(f"    {model_data.dimensions[dim].name} size = {model_data.dimensions[dim].size}")
    print("")

    # variable information.
    # list of nc variables
    print("  variables:")
    for var in model_data.variables:
        print(f"    {var}:")
        for attr, value in vars(model_data.variables[var]).items():
            if very_verbose:
                # print all infos
                print(f"      {attr} = {value}")
            elif attr == 'long_name':
                # only long_name info
                print(f"      {attr} = {value}")
    print("")

    # global attributes
    print("  global attributes:")
    for attr, value in vars(model_data).items():
        if isinstance(value, str):
            value = value.replace('\n', '; ')
        if very_verbose:
            # print all infos
            print(f"    {attr} = {value}")
        elif attr == 'id' or attr == 'model' or attr == 'title':
            # only id and model name
            print(f"    {attr} = {value}")
    print("")


#
#-----------------------------------------------------------------------------
#

def read_data_array(var_name,model_data):
    """
    reads in variable data as numpy array
    """
    global depth_index,lat_index,lon_index
    global use_replace_missing_values_with_average

    # gets variable data array
    var_array = np.array(model_data.variables[var_name])

    # determines indexing (which index for which dimension)
    # needed also to get index of depth for the data array averaging
    depth_index,lat_index,lon_index = get_array_indexing(var_name,var_array,model_data,depth_index,lat_index,lon_index,ndepths,nlats,nlons)

    # replaces missing values
    if use_replace_missing_values_with_average:
        # gets missing value
        # checks if there is a missing value specified for the array
        has_missing_val = True
        try:
            missing_val = model_data.variables[var_name].missing_value
        except:
            # array has no missing value specified
            has_missing_val = False

        # fill missing values
        if has_missing_val:
            var_array = replace_missing_values_with_average(var_name,var_array,model_data)
        else:
            print("  ",var_name," model has no missing values specified")

    return var_array

#
#-----------------------------------------------------------------------------
#

def replace_missing_values_with_average(var_name,var_array,model_data):
    """
    replaces missing values in array
    """
    global depth_index,lat_index,lon_index
    global ndepths,nlats,nlons

    print("  replace missing values for model parameter: ",var_name)

    # gets missing value
    try:
        missing_val = model_data.variables[var_name].missing_value
    except:
        # array has no missing value specified
        # we assume that there are no such and use a default missing value
        missing_val = 99999.0

    var_array = fill_missing_data_with_average_values(var_array,missing_val,ndepths,depth_index)

    print("  new min/max: ",var_array.min(),var_array.max())
    print("")

    return var_array

#
#-----------------------------------------------------------------------------
#

def get_array_indexing(var_name,var_array,model_data,depth_index,lat_index,lon_index,ndepths,nlats,nlons):
    """
    determines which array index is for depth/latitude/longitude
    """
    # gets dimensions
    dims = model_data.variables[var_name].get_dims()

    for i,dim in enumerate(dims):
        if dim.name == 'depth':
            # checks length
            if dim.size != ndepths:
                print("Error: invalid depth shapes ",dim,"should have length ",ndepths)
                sys.exit(1)
            # sets depth index
            # assumes that all other model parameters will have same shape
            if depth_index == -1:
                depth_index = i
            else:
                # check all model array data has the same shapes
                if depth_index != i:
                    print("Error: invalid depth index for array: ",var_name,depth_index,i)
                    sys.exit(1)

        elif dim.name == 'longitude':
            # checks length
            if dim.size != nlons:
                print("Error: invalid longitude shapes ",dim,"should have length ",nlons)
                sys.exit(1)
            # sets longitude index
            # assumes that all other model parameters will have same shape
            if lon_index == -1:
                lon_index = i
            else:
                # check all model array data has the same shapes
                if lon_index != i:
                    print("Error: invalid longitude index for array: ",var_name,lon_index,i)
                    sys.exit(1)

        elif dim.name == 'latitude':
            # checks length
            if dim.size != nlats:
                print("Error: invalid latitude shapes ",dim,"should have length ",nlats)
                sys.exit(1)
            # sets latitude index
            # assumes that all other model parameters will have same shape
            if lat_index == -1:
                lat_index = i
            else:
                # check all model array data has the same shapes
                if lat_index != i:
                    print("Error: invalid latitude index for array: ",var_name,lat_index,i)
                    sys.exit(1)

        else:
            print("Error: array dimension name not recognized ",dim.name)
            sys.exit(1)

    # checks index with array shape
    nx,ny,nz = var_array.shape

    # checks each dimension length
    # depth
    if depth_index == 0:
        if nx != ndepths:
            print("Error: invalid depth shape nx ",nx,"should be ",ndepths)
            sys.exit(1)
    elif depth_index == 1:
        if ny != ndepths:
            print("Error: invalid depth shape ny ",ny,"should be ",ndepths)
            sys.exit(1)
    else:
        if nz != ndepths:
            print("Error: invalid depth shape nz ",nz,"should be ",ndepths)
            sys.exit(1)

    # latitude
    if lat_index == 0:
        if nx != nlats:
            print("Error: invalid latitude shape nx ",nx,"should be ",nlats)
            sys.exit(1)
    elif lat_index == 1:
        if ny != nlats:
            print("Error: invalid latitude shape ny ",ny,"should be ",nlats)
            sys.exit(1)
    else:
        if nz != nlats:
            print("Error: invalid latitude shape nz ",nz,"should be ",nlats)
            sys.exit(1)

    # longitude
    if lon_index == 0:
        if nx != nlons:
            print("Error: invalid longitude shape nx ",nx,"should be ",nlons)
            sys.exit(1)
    elif lon_index == 1:
        if ny != nlons:
            print("Error: invalid longitude shape ny ",ny,"should be ",nlons)
            sys.exit(1)
    else:
        if nz != nlons:
            print("Error: invalid longitude shape nz ",nz,"should be ",nlons)
            sys.exit(1)

    # check if indexing works
    if depth_index == -1 or lat_index == -1 or lon_index == -1:
        print("Error: array indexing not working ",var_name,depth_index,lat_index,lon_index)
        sys.exit(1)

    if depth_index == lat_index or depth_index == lon_index or lat_index == lon_index:
        print("Error: array indexing not working ",var_name,depth_index,lat_index,lon_index)
        sys.exit(1)

    return depth_index,lat_index,lon_index


#
#-----------------------------------------------------------------------------
#

def fill_missing_data_with_average_values(array,missing_val,ndepths,depth_index):
    """
    fills missing data value with the average value from each depth slice
    """
    # gets average value for each depth
    for i in range(ndepths):
        # gets depth slice
        if depth_index == 0:
            depth_layer = array[i,:,:]
        elif depth_index == 1:
            depth_layer = array[:,i,:]
        else:
            depth_layer = array[:,:,i]

        # mask out entries w/ missing values
        mask = depth_layer != missing_val

        # array without missing values
        masked_layer = depth_layer[mask]

        # takes average value at this depth (ignoring missing value points)
        average_value = np.mean(masked_layer)

        #debug
        #print("debug: average in depth layer ",i,average_value)

        # fills with average values where missing values
        depth_layer[~mask] = average_value

        # replaces array
        if depth_index == 0:
            array[i,:,:] = depth_layer
        elif depth_index == 1:
            array[:,i,:] = depth_layer
        else:
            array[:,:,i] = depth_layer

    return array


#
#-----------------------------------------------------------------------------
#

def get_topo_DEM(region,res='low'):
    """
    downloads and creates tif-file with elevation data from SRTM 1-arc second or SRTM 3-arc seconds
    """
    global create_surface_elevation_figure

    print("get topo:")
    print("  region array: ",region)
    print("")

    # region format: #lon_min #lat_min #lon_max #lat_max (left bottom right top) in degrees
    lon_min = region[0]
    lat_min = region[1]
    lon_max = region[2]
    lat_max = region[3]

    # for example: region = (12.35, 41.8, 12.65, 42.0)
    print("  region:")
    print("  longitude min/max = ",lon_min,"/",lon_max,"(deg)")
    print("  latitude  min/max = ",lat_min,"/",lat_max,"(deg)")
    print("")

    # earth radius
    earth_radius = 6371000.0
    # earth circumference
    earth_d = 2.0 * math.pi * earth_radius # in m ~ 40,000 km
    # lengths
    length_deg = earth_d * 1./360 # length of 1 degree ~ 111 km
    length_arcsec = length_deg * 1.0/3600 # length of 1 arcsecond ~ 30.89 m

    # range / lengths (in km)
    range_lon = region[2]-region[0] # in degreee
    range_lat = region[3]-region[1]
    length_lon = range_lon * length_deg / 1000.0 # in km
    length_lat = range_lat * length_deg / 1000.0

    print("  longitude range: ",length_lon,"(km)")
    print("  latitude  range: ",length_lat,"(km)")
    print("")

    # maximum height of earth curvature
    # e.g., http://earthcurvature.com
    #
    #             ***
    #          *** | ***
    #       ***  h |    ***
    #    A**-------|-------**B       curvature height (h) between point A and B
    #
    # formula:  alpha = (A + B)/2            mid-distance (in degree) as angle
    #           h = R * ( 1 - cos(alpha) )   with R = earth radius 6371.0 km, assuming spherical earth
    alpha = 0.5 * max(range_lon,range_lat)
    alpha = alpha * math.pi/180.0 # in rad
    h = earth_radius * (1.0 - math.cos(alpha) )
    print("  maximum Earth curvature height = ",h,"(m)")
    print("")

    # gets elevation DEM data
    if 1 == 0:
        # not working yet with elevation module...
        # fails when too many tiles for larger region has to be downloaded.
        #
        # elevation module
        print("elevation module:")
        elevation.info(product=product)
        print("")

        # resolution
        if res == 'low':
            product = 'SRTM3'
            length = 3.0 * length_arcsec
        else:
            product = 'SRTM1'
            length = length_arcsec

        print("  resolution : ",res,product)
        print("  step length: ",length,"(m)")
        print("")

        # tiles
        if res == 'low':
            ilon,ilat = elevation.datasource.srtm3_tile_ilonlat(region[0],region[1])
            if ilon < 0 or ilat < 0:
                print("invalid tile number: ",ilon,ilat,"please check if lon/lat order in your input is correct")
                sys.exit(1)
            tiles = list(elevation.datasource.srtm3_tiles_names(region[0],region[1],region[2],region[3]))
        else:
            ilon,ilat = elevation.datasource.srtm1_tile_ilonlat(region[0],region[1])
            if ilon < 0 or ilat < 0:
                print("invalid tile number: ",ilon,ilat,"please check if lon/lat order in your input is correct")
                sys.exit(1)
            tiles = list(elevation.datasource.srtm1_tiles_names(region[0],region[1],region[2],region[3]))
        print("tiles:",len(tiles))
        for name in tiles:
            print("  ",name)
        print("")

        # note: in case elevation fails to download files, check in folder:
        #       /opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/elevation/datasource.py
        #       SRTM 3 files have new address: http://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/tiff
        #       test with:
        #       > curl -s -o srtm_16_09.zip.html http://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/tiff/srtm_16_09.zip

        # tif file
        filename = "./ptopo-DEM.tif"

        # get data, bounds: left bottom right top
        print("  creating topography DEM tif file:")
        print("  region : ",region)
        print("  path   : ",filename)
        print("  product: ",product)
        elevation.clip(bounds=region,output=filename,product=product)

        # cleans cache ( ~/Library/Caches/elevation/)
        elevation.clean(product=product)

        # check
        if not os.path.isfile(filename):
            print("error getting topo DEM file: ",filename)
            sys.exit(1)

        # not working yet...
        # use pygmt instead for now
        #
        # using raster io
        import rasterio
        # Load the elevation data using rasterio
        src = rasterio.open(filename)
        topo_grid = src.read(1)  # Read the elevation band (assuming single-band data)
        # debug
        #lat = 40.0
        #lon = 100.0
        ## Convert geographic coordinates to pixel coordinates in the elevation dataset
        #row, col = src.index(target_longitude, target_latitude)
        ## Access elevation at the specified pixel coordinates
        #elevation = topo_grid[row, col]
        #print("debug: lat/lon/elevation = ",lat,lon,elevation)
        return topo_grid

    # GMT python interface
    # https://www.pygmt.org/dev/api/generated/pygmt.datasets.load_earth_relief.html
    #
    # GMT region format: e.g. -R123.0/132.0/31.0/40.0
    gmt_region = [lon_min,lon_max,lat_min,lat_max]

    # ETOPO1 (1-arc minute)
    # (topo_grid is an xarray data object)
    topo_grid = pygmt.datasets.load_earth_relief(resolution="01m",
                                                 region=gmt_region,
                                                 registration="gridline")

    # checks units are in m
    if not topo_grid.units == 'meters':
        print("Error: topo grid has invalid unit ",topo_grid.units)
        sys.exit(1)

    # debug
    #lat = 40.0
    #lon = 100.0
    #elevation = get_topo_elevation(lat,lon,topo_grid)
    #print("debug: lat/lon/elevation = ",lat,lon,elevation)

    # show/write figure plot
    if create_surface_elevation_figure:
        print("  creating surface elevation figure...")
        fig = pygmt.Figure()
        fig.grdimage(grid=topo_grid,
                     cmap="haxby",
                     projection="M10c",
                     frame=True,
                    )
        fig.colorbar(frame=["x+lelevation", "y+lm"])
        #fig.show()

        # save figure as jpeg image
        name = "./tomography_model.region-plot.jpg"
        fig.savefig(name, crop=True, dpi=720)
        print("  figure plotted to: ",name)
        print("")

    return topo_grid

#
#-----------------------------------------------------------------------------
#

def get_topo_elevation(lat,lon,topo_grid):
    # gets elevation interpolated by nearest neighbor method
    elevation = topo_grid.interp(lon=lon,lat=lat,method="nearest")

    # extract simple float value from returned xarray object
    elevation_val = elevation.data

    return elevation_val

#
#-----------------------------------------------------------------------------
#

def get_Voigt_average_model(vpv,vph):
    """
    returns isotropic velocities (vp or vs) based on Voigt average for input radial anisotropic velocities (vpv,vph) or (vsv,vsh)
    """
    # Voigt average
    vp = np.sqrt( (2.0 * vpv*vpv + vph*vph)/3.0 )

    # formula would be the same for vs
    #vs = np.sqrt( (2.0 * vsv*vsv + vsh*vsh)/3.0 )

    return vp

#
#-----------------------------------------------------------------------------
#

def scale_rho_from_vp_Brocher(vp_in):
    """
    returns density scaled from vp according to Brocher(2005)'s relationship
    """
    # Brocher 2005, Empirical Relations between Elastic Wavespeeds and Density in the Earth's Crust, BSSA
    # factors from eq. (1)
    fac1 = 1.6612
    fac2 = -0.4721
    fac3 = 0.0671
    fac4 = -0.0043
    fac5 = 0.000106

    # scaling requires vp in km/s
    # (assumes input vp is given in m/s)

    # Vp (in km/s)
    vp = vp_in * 1.0/1000.0

    vp_p2 = vp * vp
    vp_p3 = vp * vp_p2
    vp_p4 = vp * vp_p3
    vp_p5 = vp * vp_p4

    # scaling relation: eq.(1)
    # (rho given in g/cm^3)
    rho = fac1 * vp + fac2 * vp_p2 + fac3 * vp_p3 + fac4 * vp_p4 + fac5 * vp_p5

    # file output needs rho in kg/m^3
    # density scaling for rho in kg/m^3: converts g/cm^3 -> kg/m^3
    # rho [kg/m^3] = rho * 1000 [g/cm^3]
    rho *= 1000.0

    return rho

#
#-----------------------------------------------------------------------------
#

def scale_vs_from_vp_Brocher(vp_in):
    """
    returns vs scaled from vp according to Brocher(2005)'s relationship
    """
    # Brocher 2005, Empirical Relations between Elastic Wavespeeds and Density in the Earth's Crust, BSSA
    # factors from eq. (1)
    fac1 = 0.7858
    fac2 = -1.2344
    fac3 = 0.7949
    fac4 = -0.1238
    fac5 = 0.0064

    # scaling requires vp in km/s
    # (assumes input vp is given in m/s)

    # Vp (in km/s)
    vp = vp_in * 1.0/1000.0

    vp_p2 = vp * vp
    vp_p3 = vp * vp_p2
    vp_p4 = vp * vp_p3

    # scaling relation: eq.(1)
    # (Vs given in km/s)
    vs = fac1 + fac2 * vp + fac3 * vp_p2 + fac4 * vp_p3 + fac5 * vp_p4

    # file output needs vs in m/s
    # Vs scaling for vs in m/s: converts km/s -> m/s
    vs *= 1000.0

    return vs


#
#-----------------------------------------------------------------------------
#

def scale_vp_from_vs_Brocher(vs_in):
    """
    returns vp scaled from vs according to Brocher(2005)'s relationship
    """
    # Brocher 2005, Empirical Relations between Elastic Wavespeeds and Density in the Earth's Crust, BSSA
    # factors from eq. (1)
    fac1 = 0.9409
    fac2 = 2.0947
    fac3 = -0.8206
    fac4 = 0.2683
    fac5 = -0.0251

    # scaling requires vs in km/s
    # (assumes input vs is given in m/s)

    # Vs (in km/s)
    vs = vs_in * 1.0/1000.0

    vs_p2 = vs * vs
    vs_p3 = vs * vs_p2
    vs_p4 = vs * vs_p3

    # scaling relation: eq.(1)
    # (Vp given in km/s)
    vp = fac1 + fac2 * vs + fac3 * vs_p2 + fac4 * vs_p3 + fac5 * vs_p4

    # file output needs vp in m/s
    # Vp scaling for vp in m/s: converts km/s -> m/s
    vp *= 1000.0

    return vp

#
#-----------------------------------------------------------------------------
#

def scale_qp_from_qs(qs_in,vs_in,vp_in):
    """
    returns Qp scaled from Qs
    """
    # (empirical) scaling factor to scale Qp from Qs
    SCALING_FACTOR_QP_FROM_QS = 1.5

    # enforces ratio Qs/Qp >= L factor from Anderson & Hart (1978)
    # IMPORTANT: this flag applies only if USE_OLSEN_ATTENUATION is true
    use_Anderson_criteria = True

    # scales Qp from Qs
    if use_Anderson_criteria:
        # Anderson & Hart (1978), Q of the Earth, JGR, 83, No. B12
        # conversion between (Qp,Qs) and (Qkappa,Qmu)
        # factor L
        L_val = 4.0/3.0 * (vs_in/vp_in)**2

        # Anderson & Hart criteria: Qs/Qp >= L (eq. 4) since Qmu and Qkappa must be positive
        # enforces (eq. 4) from Anderson & Hart
        qp = qs_in / L_val

    else:
        # scales with an empirical scaling factor
        qp = SCALING_FACTOR_QP_FROM_QS * qs_in

    return qp

#
#-----------------------------------------------------------------------------
#

def netCDF_2_tomo(input_file,UTM_zone=None,mesh_area=None,maximum_depth=None):
    """
    create tomographic_model file from a netCDF EMC model file
    """
    global ndepths,nlats,nlons
    global use_replace_missing_values_with_average
    global use_surface_elevation
    global create_vtk_file_output

    # user output
    print("")
    print("convert netCDF file")
    print("  input file   : ",input_file)
    if not UTM_zone is None:
        print("  UTM zone     : ",UTM_zone)
    if not mesh_area is None:
        print("  mesh_area    : ",mesh_area)
    if not maximum_depth is None:
        print("  maximum depth: ",maximum_depth)
    print("")

    # target mesh area min/max
    if not mesh_area is None:
        mesh_lon_min = mesh_area[0]
        mesh_lat_min = mesh_area[1]
        mesh_lon_max = mesh_area[2]
        mesh_lat_max = mesh_area[3]

        print("  target mesh area: latitude  min/max = {} / {}".format(mesh_lat_min,mesh_lat_max))
        print("                    longitude min/max = {} / {}".format(mesh_lon_min,mesh_lon_max))
        print("")

        # check
        if mesh_lon_min > mesh_lon_max or mesh_lat_min > mesh_lat_max:
            print("Wrong mesh area format, please provide in a format (lon_min,lat_min,lon_max,lat_max)")
            sys.exit(1)

    # current directory
    dir = os.getcwd()
    print("  current directory: ",dir)
    print("")

    # checks if file exists
    if not os.path.isfile(input_file):
        print("Please check if file exists: ",input_file)
        sys.exit(1)

    # reads in netCDF file
    model_data = nc.Dataset(input_file)

    # model info
    display_infos(model_data)

    # checks if lon/lat/depth are coordinates in the netCDF file
    if model_data.variables['longitude'] == None or \
       model_data.variables['latitude'] == None or \
       model_data.variables['depth'] == None:
        print("Invalid netCDF file: need longitude/latitude/depth as coordinates")
        sys.exit(1)

    global_vars = vars(model_data)
    if 'geospatial_lat_min' in global_vars:
        model_lat_min = model_data.geospatial_lat_min
        model_lat_max = model_data.geospatial_lat_max
        if model_data.geospatial_lat_resolution:
            model_dlat = model_data.geospatial_lat_resolution

        model_lon_min = model_data.geospatial_lon_min
        model_lon_max = model_data.geospatial_lon_max
        if model_data.geospatial_lon_resolution:
            model_dlon = model_data.geospatial_lon_resolution

        model_depth_min = model_data.geospatial_vertical_min
        model_depth_max = model_data.geospatial_vertical_max

        print("  model range:")
        print("    lon   min/max = {} / {}".format(model_lon_min,model_lon_max))
        print("    lat   min/max = {} / {}".format(model_lat_min,model_lat_max))
        print("    depth min/max = {} / {}".format(model_depth_min,model_depth_max))
        print("")

        # check w/ maximum depth
        if not maximum_depth is None:
            if maximum_depth > model_depth_max:
                print("  Warning: maximum depth of model {} is smaller than target maximum depth {}".format(model_depth_max,maximum_depth))
                print("           -> model output will be cutoff at the maximum depth of the model")
                print("")

    # coordinates lon / lat / depth
    x = model_data.variables['longitude'][:]
    y = model_data.variables['latitude'][:]
    z = model_data.variables['depth'][:]

    # determines scaling factor for depth if given in km
    if model_data.variables['depth'].units == "km":
        # given in km -> m
        depth_unit_scale = 1000
    elif model_data.variables['depth'].units == "m":
        # given in m
        depth_unit_scale = 1
    else:
        print("Error: depth array has invalid unit ",model_data.variables['depth'].units)
        sys.exit(1)
    # scales depth to m
    if depth_unit_scale != 1:
        z *= depth_unit_scale

    nlons = len(x)
    nlats = len(y)
    ndepths = len(z)

    print("  number of entries:")
    print("    lon/lat/depth (nx/ny/nz) = ",nlons,nlats,ndepths)
    print("")

    if nlons == 0 or nlats == 0 or ndepths == 0:
        print("Error: Invalid model dimensions ",nlons,nlats,ndepths)
        sys.exit(1)

    # tomo file needs regular spacing
    print("  checking regular spacing:")
    # numerical tolerance in spacing increments
    TOL_SPACING = 1.e-4
    # checks longitudes
    dx0 = x[1] - x[0]
    for i in range(1,nlons):
        dx = x[i] - x[i-1]
        if np.abs(dx - dx0) > TOL_SPACING:
            print("Error: irregular spacing in longitudes: ",i,dx,dx0,"step: ",x[i],x[i-1])
            print("       using a tolerance ",TOL_SPACING)
            print("       Please select an EMC model with regular spacing, exiting...")
            sys.exit(1)
    if model_dlon:
        dlon = model_dlon
        if np.abs(dlon - dx0) > TOL_SPACING:
            print("Error: irregular spacing in longitudes: ",dx0," instead of model specification: ",dlon)
            print("       using a tolerance ",TOL_SPACING)
            print("       Please select an EMC model with regular spacing, exiting...")
            sys.exit(1)
    else:
        dlon = dx0
    print("    has regular longitude spacing ",dlon,"(deg)")

    # checks latitudes
    dy0 = y[1] - y[0]
    for i in range(1,nlats):
        dy = y[i] - y[i-1]
        if np.abs(dy - dy0) > TOL_SPACING:
            print("Error: irregular spacing in latitudes: ",i,dy,dy0,"step: ",y[i],y[i-1])
            print("       using a tolerance ",TOL_SPACING)
            print("       Please select an EMC model with regular spacing, exiting...")
            sys.exit(1)
    if model_dlat:
        dlat = model_dlat
        if np.abs(dlat - dy0) > TOL_SPACING:
            print("Error: irregular spacing in latitudes: ",dy0," instead of model specification: ",dlat)
            print("       using a tolerance ",TOL_SPACING)
            print("       Please select an EMC model with regular spacing, exiting...")
            sys.exit(1)
    else:
        dlat = dy0
    print("    has regular latitude  spacing ",dlat,"(deg)")

    # checks depths
    dz0 = z[1] - z[0]
    dz_increments = np.array([dz0])
    has_irregular_depths = False
    for i in range(1,ndepths):
        dz = z[i] - z[i-1]
        if np.abs(dz - dz0) > TOL_SPACING:
            has_irregular_depths = True
            # adds increment to list of increments
            if dz not in dz_increments:
                dz_increments = np.append(dz_increments,dz)

    if has_irregular_depths:
        print("    has irregular depth   spacing (m): ",dz_increments)
        # chooses minimum spacing
        # Convert float numbers to integers (multiply by a common factor)
        common_factor = 10.0
        integer_numbers = np.round(dz_increments * common_factor).astype(int)
        # greatest common divisor: numpy's gcd function calculates the GCD for all numbers (integers)
        gcd = np.gcd.reduce(integer_numbers)

        # regular depth increment (float)
        ddepth = float(gcd / common_factor)

        print("")
        print("    tomography model will use a regular spacing: ",ddepth,"(m)")
        print("    converting model depths to regular depths...")
        # convert irregular depths to a regular gridding
        # number of regular depth increments
        ndepths_regular = np.round(z.max() / ddepth).astype(int) + 1
        z_regular = np.linspace(z.min(),z.max(),ndepths_regular)

    else:
        # has regular spacing
        ddepth = dz0
        print("  has regular depth     spacing ",ddepth,"(km)")

    print("")
    print("    dlon/dlat/ddepth = {} (deg) / {} (deg) / {} (km)".format(dlon,dlat,ddepth/1000.0))
    print("    spacing ok")
    print("")

    # checks parameterization
    nc_vars = [var for var in model_data.variables]
    #print("debug: ",nc_vars)

    if ("vpv" in nc_vars and "vph" in nc_vars) or ("vsv" in nc_vars and "vsh" in nc_vars):
        has_radial_anisotropy = True
    else:
        has_radial_anisotropy = False

    if ("vp" in nc_vars) or ("vs" in nc_vars):
        has_isotropy = True
    else:
        has_isotropy = False

    if "rho" in nc_vars:
        has_density = True
    else:
        has_density = False

    if ("qmu" in nc_vars) or ("qs" in nc_vars) :
        has_shear_attenuation = True
    else:
        has_shear_attenuation = False

    print("model: ",model_data.id)
    print("  has radial anisotropy   : ",has_radial_anisotropy)
    print("  has isotropic velocities: ",has_isotropy)
    print("  has shear attenuation   : ",has_shear_attenuation)
    print("  has density             : ",has_density)
    print("")

    # checks if model can be used
    if (not has_radial_anisotropy) and (not has_isotropy):
        print("Error: given model is incomplete for converting to a tomographic model with (vp,vs,rho) parameterization")
        sys.exit(1)


    # checks if model range covers target region
    if not mesh_area is None:
        # latitude range
        if mesh_lat_min > model_lat_min and mesh_lat_min < model_lat_max and \
            mesh_lat_max > model_lat_min and mesh_lat_max < model_lat_max:
            print("  target latitude range {} / {} within model range {} / {}".format(mesh_lat_min,mesh_lat_max,model_lat_min,model_lat_max))
        else:
            print("")
            print("ERROR: target latitude range {} / {} outside model range {} / {}".format(mesh_lat_min,mesh_lat_max,model_lat_min,model_lat_max))
            print("       nothing to extract...")
            print("")
            print("       Please make sure that the EMC model covers your target region, exiting...")
            print("")
            sys.exit(1)

        # longitude range
        if mesh_lon_min > model_lon_min and mesh_lon_min < model_lon_max and \
            mesh_lon_max > model_lon_min and mesh_lon_max < model_lon_max:
            print("  target longitude range {} / {} within model range {} / {}".format(mesh_lon_min,mesh_lon_max,model_lon_min,model_lon_max))
        else:
            print("")
            print("ERROR: target longitude range {} / {} outside model range {} / {}".format(mesh_lon_min,mesh_lon_max,model_lon_min,model_lon_max))
            print("       nothing to extract...")
            print("")
            print("       Please make sure that the EMC model covers your target mesh area, exiting...")
            print("")
            sys.exit(1)
        print("")

    # checks depth
    if not maximum_depth is None:
        if maximum_depth > model_depth_max:
            print("WARNING: target maximum depth {} below model maximum depth {}".format(maximum_depth,model_depth_max))
            print("         tomographic model output will be cut-off at model depth {}".format(model_depth_max))
            print("")

    # fill missing values with average value from each depth
    if use_replace_missing_values_with_average:
        print("  replacing missing values with depth average")
        print("")

    # get model values
    model = dict()

    # reads available velocity values
    if has_radial_anisotropy:
        # Vp
        if "vpv" in nc_vars and "vph" in nc_vars:
            vpv = read_data_array('vpv',model_data)
            vph = read_data_array('vph',model_data)

            # file output requires velocities in m/s
            # determines scaling factor
            # (based on vpv array, assuming all other velocity arrays have the same units)
            if model_data.variables['vpv'].units in ["km.s-1","km/s"]:
                # given in km/s -> m/s
                unit_scale = 1000
            elif model_data.variables['vpv'].units in ["m.s-1","m/s"]:
                # given in m/s
                unit_scale = 1
            else:
                print("Error: vpv array has invalid unit ",model_data.variables['vpv'].units)
                sys.exit(1)
            # scales velocities to m/s
            if unit_scale != 1:
                vpv *= unit_scale
                vph *= unit_scale

            # convert to isotropic model using Voigt average
            vp = get_Voigt_average_model(vpv,vph)
            model['vp'] = vp

        # Vs
        if "vsv" in nc_vars and "vsh" in nc_vars:
            vsv = read_data_array('vsv',model_data)
            vsh = read_data_array('vsh',model_data)

            # file output requires velocities in m/s
            # determines scaling factor
            # (based on vpv array, assuming all other velocity arrays have the same units)
            if model_data.variables['vsv'].units in ["km.s-1","km/s"]:
                # given in km/s -> m/s
                unit_scale = 1000
            elif model_data.variables['vsv'].units in ["m.s-1","m/s"]:
                # given in m/s
                unit_scale = 1
            else:
                print("Error: vsv array has invalid unit ",model_data.variables['vsv'].units)
                sys.exit(1)
            # scales velocities to m/s
            if unit_scale != 1:
                vsv *= unit_scale
                vsh *= unit_scale

            # convert to isotropic model using Voigt average
            vs = get_Voigt_average_model(vsv,vsh)
            model['vs'] = vs


        # not used any further...
        #if "eta" in nc_vars:
        #    eta = read_data_array('eta',model_data)

    elif has_isotropy:
        if "vp" in nc_vars:
            # vp given
            vp = read_data_array('vp',model_data)

            # file output requires velocities in m/s
            # determines scaling factor
            if model_data.variables['vp'].units in ["km.s-1","km/s"]:
                # given in km/s -> m/s
                unit_scale = 1000
            elif model_data.variables['vp'].units in ["m.s-1","m/s"]:
                # given in m/s
                unit_scale = 1
            else:
                print("Error: vp array has invalid unit ",model_data.variables['vp'].units)
                sys.exit(1)
            # scales velocities to m/s
            if unit_scale != 1:
                vp *= unit_scale

            model['vp'] = vp

        if "vs" in nc_vars:
            # vs given
            vs = read_data_array('vs',model_data)

            # file output requires velocities in m/s
            # determines scaling factor
            if model_data.variables['vs'].units in ["km.s-1","km/s"]:
                # given in km/s -> m/s
                unit_scale = 1000
            elif model_data.variables['vs'].units in ["m.s-1","m/s"]:
                # given in m/s
                unit_scale = 1
            else:
                print("Error: vs array has invalid unit ",model_data.variables['vs'].units)
                sys.exit(1)
            # scales velocities to m/s
            if unit_scale != 1:
                vs *= unit_scale

            model['vs'] = vs
    print("")


    # parameter scaling
    # checks missing velocity parameters
    if ('vp' in model) and (not 'vs' in model):
        vs = scale_vs_from_vp_Brocher(vp)
        model['vs'] = vs

    if ('vs' in model) and (not 'vp' in model):
        vp = scale_vp_from_vs_Brocher(vs)
        model['vp'] = vp

    # adds density
    if has_density:
        rho = read_data_array('rho',model_data)
        model['rho'] = rho

        # file output requires density in kg/m^3
        # determines scaling factor
        if model_data.variables['rho'].units in ["g.cm-3","g/cm3"]:
            # given in g/cm^3 -> kg/m^3
            # rho [kg/m^3] = rho * 1000 [g/cm^3]
            unit_scale = 1000
        elif model_data.variables['rho'].units in ["kg.cm-3","kg/cm3"]:
            # given in kg/cm^3 -> kg/m^3
            # rho [kg/m^3] = rho * 1000000 [kg/cm^3]
            unit_scale = 1000000
        elif model_data.variables['rho'].units in ["kg.m-3","kg/m3"]:
            # given in kg/m^3
            unit_scale = 1
        else:
            print("Error: rho array has invalid unit ",model_data.variables['rho'].units)
            sys.exit(1)
        # converts density to default kg/m^3
        if unit_scale != 1:
            model['rho'] *= unit_scale
    else:
        # scale density from vp using Brocher relationship
        rho = scale_rho_from_vp_Brocher(model['vp'])
        model['rho'] = rho

    # adds attenuation
    if has_shear_attenuation:
        if "qmu" in nc_vars:
            # Qs == qmu are equal
            qs = read_data_array('qmu',model_data)
            model['qs'] = qs
        if "qs" in nc_vars:
            qs = read_data_array('qs',model_data)
            model['qs'] = qs

        # checks missing Q arrays
        if not 'qs' in model:
            print("Error: misses Qs in model, no shear attenuation available")
            sys.exit(1)

        if "qp" in nc_vars:
            qp = read_data_array('qp',model_data)
            # Qs == qmu are equal
            model['qp'] = qp
        else:
            # scales Qp from Qs
            qp = scale_qp_from_qs(qs,model['vs'],model['vp'])
            model['qp'] = qp


    # checks if depth given with respect to earth surface (or sea level)
    depth_descr = model_data.variables['depth'].long_name
    # user info
    print("depth:")
    if "below earth surface" in depth_descr:
        print("  depth given with respect to earth surface")
    elif "below sea level" in depth_descr:
        print("  depth given with respect to sea level")
    elif "below mean Earth radius of 6371 km" in depth_descr:
        print("  depth given with respect to mean Earth radius of 6371 km")
    else:
        print("  description given: ",depth_descr)
    print("")

    # elevation
    if use_surface_elevation:
        print("using surface elevation data...")
        print("")

        # region format: #lon_min #lat_min #lon_max #lat_max (left bottom right top) in degrees
        # for example: region = (12.35, 41.8, 12.65, 42.0)
        region = (model_lon_min, model_lat_min, model_lon_max, model_lat_max)

        # surface elevation data
        topo_grid = get_topo_DEM(region,res='low')

        # loops over latitudes
        print("  calculating elevations...")

        # too slow...
        #model_elevation = np.empty((nlats,nlons),dtype=float)
        #for i, lat in enumerate(y):
        #    # loops over longitudes
        #    for j, lon in enumerate(x):
        #        elevation = get_topo_elevation(lat,lon,topo_grid)
        #        model_elevation[i,j] = elevation
        #        print("debug: elevation ",i,j,lat,lon,elevation)
        #        #debug
        #        if j >= 10: break
        #    #debug
        #    if i >= 2: break
        #
        # faster using numpy arrays
        lon = np.array(x)
        lat = np.array(y)

        model_elevation = get_topo_elevation(lat,lon,topo_grid)

        print("  elevation min/max = {} / {}".format(model_elevation.min(),model_elevation.max()))
        print("")

    # vtk file output
    if create_vtk_file_output:
        print("creating vtk file for visualization...")
        print("")

        # Create points
        points = vtk.vtkPoints()

        # Create model data
        vtk_vp = vtk.vtkFloatArray()
        vtk_vp.SetNumberOfComponents(1)  # 1 components (vp)
        vtk_vp.SetName("Vp")

        # Create model data
        vtk_vs = vtk.vtkFloatArray()
        vtk_vs.SetNumberOfComponents(1)  # 1 components (vs)
        vtk_vs.SetName("Vs")

        # Create model data
        vtk_rho = vtk.vtkFloatArray()
        vtk_rho.SetNumberOfComponents(1)  # 1 components (rho)
        vtk_rho.SetName("Density")

    data_header = list()
    output_data = list()

    # initializes header variables
    header_origin_x = sys.float_info.max
    header_origin_y = sys.float_info.max
    header_origin_z = sys.float_info.max

    header_end_x = -sys.float_info.max
    header_end_y = -sys.float_info.max
    header_end_z = -sys.float_info.max

    header_vp_min = sys.float_info.max
    header_vp_max = -sys.float_info.max

    header_vs_min = sys.float_info.max
    header_vs_max = -sys.float_info.max

    header_rho_min = sys.float_info.max
    header_rho_max = -sys.float_info.max

    print("creating tomography model...")
    print("")

    # counters for actual dimension of target model
    ndim_depths = 0
    ndim_lats = 0
    ndim_lons = 0

    # switches depth array to loop over regular depths
    if has_irregular_depths:
        # original, irregular depth z-array
        z_irreg = z.copy()
        # switches z to regular depth z-array
        z = z_regular

    # loops over depth
    for k, depth in enumerate(z):
        # determines target area
        if not maximum_depth is None:
            # skips depths larger than target maximum
            if depth > maximum_depth * depth_unit_scale + ddepth:
                continue

        # counter
        ndim_depths += 1

        # user output
        if depth is not None:
            # Show the progress.
            if k == 0:
                zero_depth = depth
            elif k == 1:
                print(f'[PROCESSING] Depth range: {z[0]/1000.0} to {z[-1]/1000.0} (km)', flush=True)
                print(f"{zero_depth/1000.0}, {depth/1000.0}", end=' ', flush=True)
            else:
                print(f", {depth/1000.0}", end=' ', flush=True)

        # Go through each latitude
        for i, lat in enumerate(y):
            # determines target area
            if not mesh_area is None:
                # mesh area min/max
                # skips latitudes outside of target min/max
                if (lat < mesh_lat_min - dlat) or (lat > mesh_lat_max + dlat):
                    continue

            # counter
            if ndim_depths == 1: ndim_lats += 1

            # Go through each longitude.
            for j, lon in enumerate(x):
                # determines target area
                if not mesh_area is None:
                    # mesh area min/max
                    # skips longitudes outside of target min/max
                    if (lon < mesh_lon_min - dlon) or (lon > mesh_lon_max + dlon):
                        continue

                # counter
                if ndim_depths == 1 and ndim_lats == 1: ndim_lons += 1

                # determines x/y coordinates
                if not UTM_zone is None:
                    # UTM x/y
                    # convert lon/lat/depth to UTM
                    utm_x,utm_y = convert_lonlat2utm(UTM_zone,lon,lat)

                    #debug
                    #print("debug: lon/lat/depth = ",lon,lat,depth," UTM x/y = ",utm_x,utm_y)

                    x_val = utm_x
                    y_val = utm_y
                else:
                    # lon/lat
                    x_val = lon       # x-direction for longitudes
                    y_val = lat       # y-direction for latitudes

                # depth -> Z coordinate (Z coordinates: positive up)
                if use_surface_elevation:
                    elevation = model_elevation[i,j]
                    z_val = elevation - depth
                else:
                    # no surface elevation (assumes top being at 0.0)
                    z_val = depth    # z-direction (positive down)

                # tuple to access data arrays in correct order
                ijk = [0,0,0]
                ijk[lat_index] = i
                ijk[lon_index] = j

                if has_irregular_depths:
                    # find index in irregular depth array
                    # Find the matching index of the depth value in the irregular z_irreg array
                    k_irreg = np.max(np.where(z_irreg <= depth))
                    ijk[depth_index] = k_irreg
                else:
                    # regular depth intervals, simply use index k from loop
                    ijk[depth_index] = k

                # gets tuple for array indexing
                index = tuple(ijk)

                # model parameters
                vp_val = model['vp'][index]
                vs_val = model['vs'][index]
                rho_val = model['rho'][index]

                # data line format:
                if has_shear_attenuation:
                    qp_val = model['qp'][index]
                    qs_val = model['qs'][index]
                    #x #y #z #vp #vs #density #Q_p #Q_s
                    output_data.append("{} {} {} {} {} {} {} {}\n".format(x_val,y_val,z_val,vp_val,vs_val,rho_val,qp_val,qs_val))
                else:
                    #x #y #z #vp #vs #density
                    output_data.append("{} {} {} {} {} {}\n".format(x_val,y_val,z_val,vp_val,vs_val,rho_val))

                # header stats
                header_origin_x = min(header_origin_x,x_val)
                header_origin_y = min(header_origin_y,y_val)
                header_origin_z = min(header_origin_z,z_val)

                header_end_x = max(header_end_x,x_val)
                header_end_y = max(header_end_y,y_val)
                header_end_z = max(header_end_z,z_val)

                header_vp_min = min(header_vp_min,vp_val)
                header_vp_max = max(header_vp_max,vp_val)

                header_vs_min = min(header_vs_min,vs_val)
                header_vs_max = max(header_vs_max,vs_val)

                header_rho_min = min(header_rho_min,rho_val)
                header_rho_max = max(header_rho_max,rho_val)

                # vtk file output
                if create_vtk_file_output:
                    # adds grid point
                    if not UTM_zone is None:
                        # uses utm coordinates
                        points.InsertNextPoint(x_val, y_val, z_val)
                    else:
                        # uses grid point indices (as using lon/lat/depth gives very skewed meshes)
                        points.InsertNextPoint(j, i, k)
                    # adds model values
                    vtk_vp.InsertNextValue(vp_val)
                    vtk_vs.InsertNextValue(vs_val)
                    vtk_rho.InsertNextValue(rho_val)

    print("")
    print("")
    print("tomographic model statistics:")
    print("  vp  min/max : {:.2f} / {:.2f} (m/s)".format(header_vp_min,header_vp_max))
    print("  vs  min/max : {:.2f} / {:.2f} (m/s)".format(header_vs_min,header_vs_max))
    print("  rho min/max : {:.2f} / {:.2f} (kg/m3)".format(header_rho_min,header_rho_max))
    print("")

    if header_vp_min <= 0.0:
        print("WARNING: Vp has invalid entries with a minimum of zero!")
        print("         The provided output model is likely invalid.")
        print("         Please check with your inputs...")
        print("")

    if header_vs_min <= 0.0:
        print("WARNING: Vs has entries with a minimum of zero!")
        print("         The provided output model is likely invalid.")
        print("         Please check with your inputs...")
        print("")

    if header_rho_min <= 0.0:
        print("WARNING: Density has entries with a minimum of zero!")
        print("         The provided output model is likely invalid.")
        print("         Please check with your inputs...")
        print("")

    # header infos
    header_dx = dlon        # increment x-direction for longitudes
    header_dy = dlat        # increment y-direction for latitudes
    header_dz = ddepth      # increment z-direction for depths (positive down)

    header_nx = ndim_lons   # x-direction
    header_ny = ndim_lats   # y-direction
    header_nz = ndim_depths # z-direction

    ## SPECFEM tomographic model format
    ## file header
    data_header.append("# tomography model - converted using script run_convert_IRIS_EMC_netCDF_2_tomo.py\n")
    data_header.append("#\n")

    # providence
    data_header.append("# providence\n")
    data_header.append("# original netCDF file: {}\n".format(input_file))
    data_header.append("# created             : {}\n".format(str(datetime.datetime.now())))
    data_header.append("# command             : {}\n".format(" ".join(sys.argv)))
    data_header.append("#\n")

    # tomographic model format
    data_header.append("# model format\n")
    data_header.append("# model type          : IRIS EMC model\n")

    # note: the following comment line is important as a header comment for this file type!
    #       it will be recognized by SPECFEM3D when reading in the tomography_model.xyz file,
    #       and in case the "coordinate format" together with the "lon / lat / depth" string is found,
    #       SPECFEM3D will set a flag to read the coordinate positions at lon/lat/depth instead of x/y/z.
    #
    if not UTM_zone is None:
        if use_surface_elevation:
            data_header.append("# coordinate format   : UTM_X / UTM_Y / Z     # z-direction (positive up) using surface elevation\n")
        else:
            data_header.append("# coordinate format   : UTM_X / UTM_Y / depth # z-direction (positive down)\n")
    else:
        if use_surface_elevation:
            data_header.append("# coordinate format   : lon / lat / Z      # z-direction (positive up) using surface elevation\n")
        else:
            data_header.append("# coordinate format   : lon / lat / depth  # z-direction (positive down)\n")
    data_header.append("#\n")

    # tomography model header infos
    #origin_x #origin_y #origin_z #end_x #end_y #end_z          - start-end dimensions
    data_header.append("#origin_x #origin_y #origin_z #end_x #end_y #end_z\n")
    data_header.append("{} {} {} {} {} {}\n".format(header_origin_x,header_origin_y,header_origin_z,header_end_x,header_end_y,header_end_z))
    #dx #dy #dz                                                 - increments
    data_header.append("#dx #dy #dz\n")
    data_header.append("{} {} {}\n".format(header_dx,header_dy,header_dz))
    #nx #ny #nz                                                 - number of models entries
    data_header.append("#nx #ny #nz\n")
    data_header.append("{} {} {}\n".format(header_nx,header_ny,header_nz))
    #vp_min #vp_max #vs_min #vs_max #density_min #density_max   - min/max stats
    data_header.append("#vp_min #vp_max #vs_min #vs_max #density_min #density_max\n")
    data_header.append("{} {} {} {} {} {}\n".format(header_vp_min,header_vp_max,header_vs_min,header_vs_max,header_rho_min,header_rho_max))
    data_header.append("#\n")

    # data record info
    data_header.append("# data records - format:\n")
    # full format: x #y #z #vp #vs #density #Qp #Qs
    # coordinate format
    if not UTM_zone is None:
        if use_surface_elevation:
            data_header.append("#UTM_X #UTM_Y #Z ")
        else:
            data_header.append("#UTM_X #UTM_Y #depth ")
    else:
        if use_surface_elevation:
            data_header.append("#lon #lat #Z ")
        else:
            data_header.append("#lon #lat #depth ")
    # parameter format
    if has_shear_attenuation:
        #x #y #z #vp #vs #density #Qp #Qs
        data_header.append("#vp #vs #density #Qp #Qs\n")
    else:
        #x #y #z #vp #vs #density
        data_header.append("#vp #vs #density\n")

    ## writes output file
    filename = "./tomography_model.xyz"
    with open(filename, "w") as fp:
        fp.write(''.join(data_header))
        fp.write(''.join(output_data))
        fp.close()

    print("tomography model written to: ",filename)
    print("")


    # vtk file output
    if create_vtk_file_output:
        # finish vtk file
        # structured grid file
        # Define grid dimensions
        dimensions = (ndim_lons, ndim_lats, ndim_depths)

        # Create a structured grid
        grid = vtk.vtkStructuredGrid()
        grid.SetDimensions(dimensions)

        grid.SetPoints(points)

        # Add model data to the grid
        grid.GetPointData().AddArray(vtk_vp)
        grid.GetPointData().AddArray(vtk_vs)
        grid.GetPointData().AddArray(vtk_rho)

        # Write grid to a VTK file
        filename = "./tomography_model.vts"
        writer = vtk.vtkXMLStructuredGridWriter()
        writer.SetFileName(filename)
        writer.SetInputData(grid)
        writer.Write()

        print("vtk file written to: ",filename)
        print("")


#
#-----------------------------------------------------------------------------
#

def usage():
    print("usage: ./run_convert_IRIS_EMC_netCDF_2_tomo.py --EMC_file=file [--mesh_area=lon_min,lat_min,lon_max,lat_max] [--UTM_zone=zone] [--maximum_depth=depth]")
    print("   where")
    print("       file                            - IRIS EMC model file (netCDF format .nc)")
    print("       lon_min,..,lat_max              - (optional) target output area defined by lon/lat minimum and lon/lat maximum coordinates (in degrees);")
    print("                                          if not provided, output will cover whole region defined in model file")
    print("       zone                            - (optional) UTM zone (1-60) to convert longitude/latitude to specific UTM x/y coordinates")
    print("       depth                           - (optional) target model maximum depth (to cutoff model output; in same unit as netCDF model, e.g., in km)")
    sys.exit(1)


if __name__ == '__main__':
    # init
    input_file = ""
    UTM_zone = None
    mesh_area = None
    maximum_depth = None

    #debug
    #print("\nnumber of arguments: " + str(len(sys.argv)))

    # reads arguments
    if len(sys.argv) <= 1: usage()
    i = 0
    for arg in sys.argv:
        i += 1
        #print("argument "+str(i)+": " + arg)
        # get arguments
        if "--help" in arg:
            usage()
        elif "--EMC_file=" in arg:
            input_file = arg.split('=')[1]
        elif "--UTM_zone=" in arg:
            str_array = arg.split('=')[1]
            UTM_zone = int(arg.split('=')[1])
        elif "--mesh_area=" in arg:
            str_array = arg.split('=')[1]
            mesh_area = np.array([float(val) for val in str_array.strip('()[]').split(',')])
        elif "--maximum_depth=" in arg:
            str_array = arg.split('=')[1]
            maximum_depth = float(arg.split('=')[1])
        elif i >= 6:
            print("argument not recognized: ",arg)

    # logging
    cmd = " ".join(sys.argv)
    filename = './run_convert_IRIS_EMC_netCDF_2_tomo.log'
    with open(filename, 'a') as f:
      print("command call --- " + str(datetime.datetime.now()),file=f)
      print(cmd,file=f)
      print("command logged to file: " + filename)

    # converts netCDF model to tomography model for SPECFEM3D
    netCDF_2_tomo(input_file,UTM_zone,mesh_area,maximum_depth)

