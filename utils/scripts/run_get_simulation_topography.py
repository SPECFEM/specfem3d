#!/usr/bin/env python
#
# script to setup topography for a given region
#
# required python modules:
#   - utm       - version 0.4.2
#   - elevation - version 1.0.5
#
# required external packages:
#   - GMT       - version 5.4.2
#   - gdal      - version 2.2.3
#
from __future__ import print_function

import os
import sys
import subprocess
import math
import datetime

## elevation package
# see: https://github.com/bopen/elevation
#      http://elevation.bopen.eu/en/stable/
# install by: > sudo pip install elevation
try:
    import elevation
except:
    print("Error importing module `elevation`")
    print("install by: > pip install -U elevation")
    sys.exit(1)

## elevation
# check
elevation.util.selfcheck(elevation.TOOLS)

## UTM zones
# see: https://github.com/Turbo87/utm
# install by: > sudo pip install utm
try:
    import utm
except:
    print("Error importing module `utm`")
    print("install by: > pip install -U utm")
    sys.exit(1)

## DEM analysis (slope, topographic wetness index, ..)
# http://conference.scipy.org/proceedings/scipy2015/pdfs/mattheus_ueckermann.pdf
# install by: > sudo pip install pyDEM


## GMT
## http://gmt.soest.hawaii.edu
## install by: > sudo apt install gmt
## setup gmt functions by:
## > source /usr/share/gmt/tools/gmt_functions.sh
try:
    cmd = 'gmt --version'
    print("> ",cmd)
    version = subprocess.check_output(cmd, shell=True)
except:
    print("Error using `gmt`")
    print("install by: > sudo apt install gmt")
    sys.exit(1)
# avoid bytes string issues with strings like b'Hello', converts to text string
if isinstance(version, (bytes, bytearray)): version = version.decode("utf-8")
version = version.strip()
print("GMT version: %s" % (version))
print("")
# get version numbers for later (grdconvert command format changes between version 5.3 and 5.4)
elem = version.split(".")
gmt_major = int(elem[0])
gmt_minor = int(elem[1])

# GMT python interface
# todo: not used yet, calling gmt commands directly as shell commands...
#
#try:
#    import pygmt
#except:
#    print("Error importing module `pygmt`")
#    print("install by: > pip install -U pygmt")
#    sys.exit(1)


## GDAL/OGR
## https://www.gdal.org
try:
    cmd = 'gdalinfo --version'
    print("> ",cmd)
    version = subprocess.check_output(cmd, shell=True)
except:
    print("Error using `gdalinfo`")
    print("install by: > sudo apt install gdal-bin")
    sys.exit(1)
# avoid bytes string issues with strings like b'Hello', converts to text string
if isinstance(version, (bytes, bytearray)): version = version.decode("utf-8")
version = version.strip()
print("GDAL version: %s" % (version.strip()))
print("")

# GDAL python interface
# todo: not used yet, calling gdal commands directly as shell commands...
#
#try:
#    from osgeo import gdal
#except:
#    print("Error importing module `gdal`")
#    print("install by: > pip install -U gdal")
#    sys.exit(1)


#########################################################################
## USER PARAMETERS

## SRTM data
# type: 'low' == SRTM 90m / 'high' == SRTM 30m / else e.g. 'etopo' == topo30 (30-arc seconds)
SRTM_type = 'low'

## GMT grid sampling (in degrees)
# to convert increment to km: incr_dx * math.pi/180.0 * 6371.0 -> 1 degree = 111.1949 km
# sampling coarse (~5.5km)
#incr_dx = 0.05
# sampling fine (~1.1km)
#incr_dx = 0.01
# sampling fine (~500m)
incr_dx = 0.0045
# sampling fine (~110m)
#incr_dx = 0.001

# topography shift for 2nd interface (shifts original topography downwards)
toposhift = 8000.0
#toposhift = 1000.0 (small, local meshes)

# scaling factor for topography variations
toposcale = 0.1

## local data directory
datadir = 'topo_data'

#########################################################################

# globals
utm_zone = 0
gmt_region = ""


def get_topo_DEM(region,filename_path,res='low'):
    """
    downloads and creates tif-file with elevation data from SRTM 1-arc second or SRTM 3-arc seconds
    """
    # check
    if len(filename_path) == 0:
        print("error invalid filename",filename_path)
        sys.exit(1)

    # region format: #lon_min #lat_min #lon_max #lat_max (left bottom right top) in degrees
    # for example: region = (12.35, 41.8, 12.65, 42.0)
    print("region:")
    print("  longitude min/max = ",region[0],"/",region[2],"(deg)")
    print("  latitude  min/max = ",region[1],"/",region[3],"(deg)")
    print("")

    # earth radius
    earth_radius = 6371000.0
    # earth circumference
    earth_d = 2.0 * 3.14159 * earth_radius # in m ~ 40,000 km
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

    # resolution
    if res == 'low':
        product = 'SRTM3'
        length = 3.0 * length_arcsec
    else:
        product = 'SRTM1'
        length = length_arcsec

    print("resolution: ",res,product)
    print("  step length: ",length,"(m)")
    print("")
    print("elevation module:")
    elevation.info(product=product)
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

    # get data, bounds: left bottom right top
    print("getting topography DEM file:")
    print("  region : ",region)
    print("  path   : ",filename_path)
    print("  product: ",product)
    elevation.clip(bounds=region,output=filename_path,product=product)

    # cleans cache ( ~/Library/Caches/elevation/)
    elevation.clean(product=product)

    # check
    if not os.path.isfile(filename_path):
        print("error getting topo DEM file: ",filename_path)
        sys.exit(1)

    print("  done")
    print("")


#
#-----------------------------------------------------------------------------
#

def get_topo(lon_min,lat_min,lon_max,lat_max):
    """
    gets topography data for given region and stores it in file ptopo.xyz
    """
    global incr_dx
    global SRTM_type
    global gmt_region

    # region format: #lon_min #lat_min #lon_max #lat_max (left bottom right top) in degrees
    # for example: region = (12.35, 41.8, 12.65, 42.0)
    region = (lon_min, lat_min, lon_max, lat_max)

    print("*******************************")
    print("get topo:")
    print("*******************************")
    print("  region: ",region)
    print("")

    ## gmt
    # region format: e.g. -R123.0/132.0/31.0/40.0
    gmt_region = '-R' + str(lon_min) + '/' + str(lon_max) + '/' + str(lat_min) + '/' + str(lat_max)
    # sampling fine (~1.1km)
    gmt_interval = '-I' + str(incr_dx) + '/' + str(incr_dx)

    # current directory
    dir = os.getcwd()

    name = 'ptopo-DEM.tif'
    filename = dir + '/' + name

    print("  current directory:",dir)
    print("  topo file name   :",filename)
    print("")

    ## get topo file
    if SRTM_type == 'low' or SRTM_type == 'high':
        # enlarge data region slightly to fit desired locations
        #eps = 0.0001
        #region_extended = ( lon_min - eps, lat_min - eps, lon_max + eps, lat_max + eps)
        # SRTM data
        get_topo_DEM(region,filename,res=SRTM_type)
    elif SRTM_type == 'etopo2' \
      or SRTM_type == 'topo30' \
      or SRTM_type == 'srtm30s' \
      or SRTM_type == 'srtm_1km' \
      or SRTM_type == 'topo15' \
      or SRTM_type == 'srtm15s' \
      or SRTM_type == 'srtm_500m' \
      or SRTM_type == 'topo3' \
      or SRTM_type == 'srtm3s' \
      or SRTM_type == 'srtm_100m' \
      or SRTM_type == 'topo1' \
      or SRTM_type == 'srtm1s' \
      or SRTM_type == 'srtm_30m' :
        # gmt grid
        gridfile = 'ptopo-DEM.grd'
        if gmt_major >= 6:
            # new version uses grdcut and earth relief grids from server
            # http://gmt.soest.hawaii.edu/doc/latest/datasets.html
            if SRTM_type == 'etopo2':
                # ETOPO2 (2-arc minutes)
                cmd = 'gmt grdcut @earth_relief_02m ' + gmt_region + ' -G' + gridfile
                # interpolation?
                incr_dx = 0.05 # coarse (~5.5km)
                gmt_interval = '-I' + str(incr_dx) + '/' + str(incr_dx)
                #cmd += '; gmt blockmean ' + gridfile + ' -I2m -G' + gridfile
            elif SRTM_type == 'etopo1':
                # ETOPO1 (1-arc minute)
                cmd = 'gmt grdcut @earth_relief_01m ' + gmt_region + ' -G' + gridfile
                # interpolation?
                incr_dx = 0.01 # fine (~1.1km)
                gmt_interval = '-I' + str(incr_dx) + '/' + str(incr_dx)
                #cmd += '; gmt blockmean ' + gridfile + ' -I0.00833 -G' + gridfile
            elif SRTM_type == 'topo30' or SRTM_type == 'srtm30s' or SRTM_type == 'srtm_1km':
                # srtm 30s (30-arc seconds) ~ 1km resolution
                cmd = 'gmt grdcut @earth_relief_30s ' + gmt_region + ' -G' + gridfile
                # interpolation?
                incr_dx = 0.009 # fine (~1km)
                gmt_interval = '-I' + str(incr_dx) + '/' + str(incr_dx)
                #cmd += '; gmt blockmean ' + gridfile + ' -I0.00833 -G' + gridfile
            elif SRTM_type == 'topo15' or SRTM_type == 'srtm15s' or SRTM_type == 'srtm_500m':
                # srtm 15s (15-arc seconds) ~ 0.5km resolution
                cmd = 'gmt grdcut @earth_relief_15s ' + gmt_region + ' -G' + gridfile
                # interpolation?
                incr_dx = 0.0045 # fine (~500m)
                gmt_interval = '-I' + str(incr_dx) + '/' + str(incr_dx)
                #cmd += '; gmt blockmean ' + gridfile + ' -I0.004166 -G' + gridfile
            elif SRTM_type == 'topo3' or SRTM_type == 'srtm3s' or SRTM_type == 'srtm_100m':
                # srtm 3s (3-arc seconds) ~ 100m resolution
                cmd = 'gmt grdcut @earth_relief_03s ' + gmt_region + ' -G' + gridfile
                # interpolation?
                incr_dx = 0.00083 # fine (~100m)
                gmt_interval = '-I' + str(incr_dx) + '/' + str(incr_dx)
                #cmd += '; gmt blockmean ' + gridfile + ' -I0.0008 -G' + gridfile
            elif SRTM_type == 'topo1' or SRTM_type == 'srtm1s' or SRTM_type == 'srtm_30m':
                # srtm 1s (1-arc seconds) ~ 30m resolution
                cmd = 'gmt grdcut @earth_relief_01s ' + gmt_region + ' -G' + gridfile
                # interpolation?
                incr_dx = 0.0003 # fine (~30m)
                gmt_interval = '-I' + str(incr_dx) + '/' + str(incr_dx)
                #cmd += '; gmt blockmean ' + gridfile + ' -I0.0003 -G' + gridfile
            else:
                print("Error invalid SRTM_type " + SRTM_type)
                sys.exit(1)
        else:
            # older version < 6, like 5.4 ...
            # or use grdraster for lower resolution
            if SRTM_type == 'etopo2':
                # ETOPO2 (2-arc minutes)
                cmd = 'gmt grdraster ETOPO2 ' + gmt_region + ' -I2m -G' + gridfile
            elif SRTM_type == 'topo30':
                # topo30 (30-arc seconds) ~ 1km resolution
                cmd = 'gmt grdraster topo30 ' + gmt_region + ' -I0.00833 -G' + gridfile
            elif SRTM_type == 'srtm_500m':
                # srtm 500m resolution
                cmd = 'gmt grdraster srtm_500m ' + gmt_region + ' -I0.004166 -G' + gridfile
            else:
                print("Error invalid SRTM_type " + SRTM_type)
                sys.exit(1)
        print("  > ",cmd)
        status = subprocess.call(cmd, shell=True)
        check_status(status)

        print("")
        # convert gmt grid-file to gdal GTiff for shading
        if gmt_major >= 5 and gmt_minor >= 4:
            # version > 5.4
            cmd = 'gmt grdconvert ' + gridfile + ' -G' + filename + '=gd:Gtiff'
        else:
            # older version < 5.3
            cmd = 'gmt grdconvert ' + gridfile + ' ' + filename + '=gd:Gtiff'
        print("  > ",cmd)
        status = subprocess.call(cmd, shell=True)
        check_status(status)
        print("")

    else:
        print("Error invalid SRTM_type " + SRTM_type)
        sys.exit(1)

    print("  GMT:")
    print("  region  : ",gmt_region)
    print("  interval: ",gmt_interval)
    print("")

    # topography info
    cmd = 'gmt grdinfo ' + filename
    print("  topography file info:")
    print("  > ",cmd)
    status = subprocess.call(cmd, shell=True)
    check_status(status)
    print("")

    # hillshade
    gif_file = filename + '.hillshaded.gif'
    cmd = 'gdaldem hillshade ' + filename + ' ' + gif_file + ' -of GTiff'
    print("  hillshade figure:")
    print("  > ",cmd)
    status = subprocess.call(cmd, shell=True)
    check_status(status)
    print("")

    # defines model region
    ## resampling
    # even sampling
    # note: the region might be slightly off, then this resampling can lead to wrong values.
    #       check output of gmt command and see what region it is using
    #
    # e.g. grdsample etopo.grd $region -I$dx/$dx -Getopo.sampled.grd
    #cmd = 'grdsample ' + filename + ' ' + gmt_region + ' ' + gmt_interval + ' -Getopo.sampled.grd'
    # uses region specified from grid file
    gridfile = 'ptopo.sampled.grd'
    cmd = 'gmt grdsample ' + filename + ' ' + gmt_interval + ' -G' + gridfile
    print("  resampling topo data")
    print("  > ",cmd)
    status = subprocess.call(cmd, shell=True)
    check_status(status)
    print("")

    cmd = 'gmt grdinfo ' + gridfile
    print("  > ",cmd)
    status = subprocess.call(cmd, shell=True)
    check_status(status)
    print("")

    # converts to xyz format
    xyz_file = 'ptopo.xyz'
    cmd = 'gmt grd2xyz ' + gridfile + ' > ' + xyz_file
    print("  converting to xyz")
    print("  > ",cmd)
    status = subprocess.call(cmd, shell=True)
    check_status(status)
    print("")

    # gif
    cmd = 'gdaldem hillshade ' + xyz_file + ' ' + xyz_file + '.hillshaded1.gif -of GTiff'
    print("  creating gif-image ...")
    print("  > ",cmd)
    status = subprocess.call(cmd, shell=True)
    check_status(status)
    print("")

    # mesh interpolation (only needed when different gridding interval)
    mean_file = 'ptopo.mean.xyz'
    # note: can lead to problems, if mean interval is similar to grid interval
    #cmd = 'gmt blockmean ' + gmt_interval + ' ptopo.xyz ' + gmt_region + ' > ptopo.mean.xyz'
    gmt_interval2 = '-I' + str(incr_dx/10.0) + '/' + str(incr_dx/10.0)

    cmd = 'gmt blockmean ' + gmt_interval2 + ' ' + xyz_file + ' ' + gmt_region + ' > ' + mean_file
    print("  mesh interpolation")
    print("  > ",cmd)
    status = subprocess.call(cmd, shell=True)
    check_status(status)
    cmd = 'mv -v ptopo.xyz ptopo.xyz.org; mv -v ptopo.mean.xyz ptopo.xyz'
    print("  > ",cmd)
    status = subprocess.call(cmd, shell=True)
    check_status(status)
    print("")

    # gif
    cmd = 'gdaldem hillshade ' + xyz_file + ' ' + xyz_file + '.hillshaded2.gif -of GTiff'
    print("  creating gif-image ...")
    print("  > ",cmd)
    status = subprocess.call(cmd, shell=True)
    check_status(status)
    print("")

    # map
    plot_map(gmt_region)

    return xyz_file

#
#-----------------------------------------------------------------------------
#

def plot_map(gmt_region,gridfile="ptopo.sampled.grd"):
    global datadir

    print("*******************************")
    print("plotting map ...")
    print("*******************************")

    # current directory
    dir = os.getcwd()
    print("  current directory:",dir)
    print("")

    #cmd = 'cd ' + datadir + '/' + ';'
    #status = subprocess.call(cmd, shell=True)
    #check_status(status)

    ps_file = "map.ps"
    pdf_file = "map.pdf"

    # gmt plotting
    cmd = 'gmt pscoast ' + gmt_region + ' -JM6i -Dh -G220 -P -K > ' + ps_file + ';'
    # topography shading
    #makecpt -Cgray -T0/1/0.01 > topo.cpt
    #cmd += 'makecpt -Cglobe -T-2500/2500/100 > topo.cpt' + ';'
    #cmd += 'makecpt -Cterra -T-2500/2500/100 > ptopo.cpt' + ';'
    #cmd += 'makecpt -Ctopo -T-2500/2500/100 > ptopo.cpt' + ';'
    cmd += 'gmt makecpt -Crelief -T-2500/2500/100 > ptopo.cpt' + ';'
    cmd += 'gmt grdgradient ' + gridfile + ' -Nt1 -A45 -Gptopogradient.grd -V' + ';'
    cmd += 'gmt grdimage ' + gridfile + ' -Iptopogradient.grd -J -R -Cptopo.cpt -V -O -K >> ' + ps_file + ';'
    cmd += 'gmt pscoast -R -J -Di -N1/1.5p,gray40 -A1000 -W1 -O -K >> ' + ps_file + ';'
    cmd += 'gmt psbasemap -O -R -J -Ba1g1:"Map": -P -V  >> ' + ps_file + ';'
    status = subprocess.call(cmd, shell=True)
    check_status(status)
    print("  map plotted in file: ",ps_file)

    # imagemagick converts ps to pdf
    cmd = 'convert ' + ps_file + ' ' + pdf_file
    print("  > ",cmd)
    status = subprocess.call(cmd, shell=True)
    check_status(status)
    print("  map plotted in file: ",pdf_file)
    print("")
    return

#
#-----------------------------------------------------------------------------
#

def create_AVS_file():
    global datadir
    global utm_zone
    global gmt_region

    print("*******************************")
    print("creating AVS border file ...")
    print("*******************************")

    # current directory
    dir = os.getcwd()
    #print("current directory:",dir)
    #print("")

    # GMT segment file
    name = "map_segment.dat"
    cmd = 'gmt pscoast ' + gmt_region + ' -Dh -W -M > ' + name + ';'
    status = subprocess.call(cmd, shell=True)
    check_status(status)
    print("  GMT segment file plotted in file: ",name)

    # note: GMT segment file has format
    # > Shore Bin # .., Level ..
    #   lon1 lat1
    #   lon2 lat2
    # > Shore Bin # .., Level ..
    #   ..

    print("  Getting boundaries from file %s ..." % name)

    # reads gmt boundary file
    with open(name,'r') as f:
        content = f.readlines()

    if len(content) == 0:
        print("")
        print("  INFO: no boundaries in file")
        print("")
        return

    # counts segments and points
    numsegments = 0
    numpoints = 0
    for line in content:
        if ">" in line:
            # new segment
            numsegments += 1
        else:
            # point
            numpoints += 1

    print("  There are %i contours" % numsegments)
    print("  There are %i data points" % numpoints)

    # read the GMT file to get the number of individual line segments
    currentelem = 0
    previous_was_comment = 1
    for line in content:
        line = line.strip()
        #debug
        #print("currentelem %i %i %s" % (currentelem,previous_was_comment,line))
        # skip comment lines
        if line[0:1] == "#": continue
        # get line marker (comment in file)
        if ">" in line:
            previous_was_comment = 1
        else:
            if previous_was_comment == 0: currentelem += 1
            previous_was_comment = 0

    num_individual_lines = currentelem
    print("  There are %i individual line segments" % num_individual_lines)

    print("")
    print("converting to .inp format:")

    avsfile = "AVS_boundaries_utm.inp"
    with open(avsfile,'w') as f:
        # write header for AVS (with point data)
        f.write("%i %i 1 0 0\n" % (numpoints,num_individual_lines))

        # read the GMT file to get the points
        currentpoint = 0
        for line in content:
            line = line.strip()
            # skip comment lines
            if line[0:1] == "#": continue

            #   get point only if line is not a comment
            if ">" not in line:
                currentpoint += 1
                elem = line.split()

                ## global lon/lat coordinates
                # longitude is the number before the white space
                lon = float(elem[0])
                # latitude is the number after the white space
                lat = float(elem[1])

                # perl example
                # convert geographic latitude to geocentric colatitude and convert to radians
                # $pi = 3.14159265;
                # $theta = $pi/2. - atan2(0.99329534 * tan($latitude * $pi / 180.),1) ;
                # $phi = $longitude * $pi / 180. ;
                # compute the Cartesian position of the receiver (ignore ellipticity for AVS)
                # assume a sphere of radius one
                # $r_target = 1. ;
                ## DK DK make the radius a little bit bigger to make sure it is
                ## DK DK correctly superimposed to the mesh in final AVS figure
                # $r_target = 1.015 ;
                # $x_target = $r_target*sin($theta)*cos($phi) ;
                # $y_target = $r_target*sin($theta)*sin($phi) ;
                # $z_target = $r_target*cos($theta) ;

                ## UTM
                # utm_x is the number before the white space
                #utm_x = float(elem[0])
                # utm_y is the number after the white space
                #utm_y = float(elem[1])
                #x_target = utm_x
                #y_target = utm_y

                # converts to UTM x/y
                x,y = geo2utm(lon,lat,utm_zone)

                # location
                x_target = x
                y_target = y
                z_target = 0.0     # assuming models use depth in negative z-direction

                f.write("%i %f %f %f\n" % (currentpoint,x_target,y_target,z_target))

        # read the GMT file to get the lines
        currentline = 0
        currentelem = 0
        currentpoint = 0
        previous_was_comment = 1
        for line in content:
            line = line.strip()
            # skip comment lines
            if line[0:1] == "#": continue
            #   get line marker (comment in file)
            if ">" in line:
                # check if previous was line was also a segment
                # for example: there can be empty segment lines
                #  > Shore Bin # 4748, Level 1
                #  > Shore Bin # 4748, Level 1
                #  > Shore Bin # 4748, Level 1
                #  136.117036698   36.2541237507
                #  136.121248188   36.2533302815
                #  ..
                if currentline > 0 and previous_was_comment :
                    continue
                else:
                    currentline += 1
                    currentpoint += 1
                    previous_was_comment = 1
                #print("processing contour %i named %s" % (currentline,line))
            else:
                if previous_was_comment == 0:
                    previouspoint = currentpoint
                    currentelem  +=1
                    currentpoint  +=1
                    # new line
                    f.write("%i %i line %i %i\n" % (currentelem,currentline,previouspoint,currentpoint))
                previous_was_comment = 0

        # dummy variable names
        f.write(" 1 1\n")
        f.write(" Zcoord, meters\n")
        # create data values for the points
        for currentpoint in range(1,numpoints+1):
            f.write("%i 255.\n" % (currentpoint))

    # check
    if numpoints != currentpoint:
        print("  WARNING:")
        print("    possible format corruption: total number of points ",numpoints," should match last line point id ",currentpoint)
        print("")

    print("  see file: %s" % avsfile)
    print("")
    return

#
#-----------------------------------------------------------------------------
#

def topo_extract(filename):
    global toposhift
    global toposcale

    # ./topo_extract.sh ptopo.mean.xyz
    #cmd = './topo_extract.sh ptopo.mean.xyz'
    print("*******************************")
    print("extracting interface data for xmeshfem3D ...")
    print("*******************************")

    import numpy

    ## shift/downscale topography
    print("  topo shift   = ",toposhift,"(m)")
    print("  scale factor = ",toposcale)
    print("")

    file1 = filename + '.1.dat'
    file2 = filename + '.2.dat'

    # statistics
    cmd = 'gmt gmtinfo ' + filename
    print("  > ",cmd)
    status = subprocess.call(cmd, shell=True)
    check_status(status)
    print("")

    # cleanup
    cmd = 'rm -f ' + file1 + ';'
    cmd += 'rm -f ' + file2 + ';'
    print("  > ",cmd)
    status = subprocess.call(cmd, shell=True)
    check_status(status)
    print("")

    # reads in lon/lat/elevation
    print("  reading file " + filename + " ...")
    data = numpy.loadtxt(filename)
    #debug
    #print(data)
    elevation = data[:,2]

    # extracts only elevation data
    with open(file1,'w') as f:
        for i in range(0,len(elevation)):
            f.write("%f\n" % (elevation[i]) )

    # shifts topography surface down,
    with open(file2,'w') as f:
        for i in range(0,len(elevation)):
            f.write("%f\n" % (elevation[i] * toposcale - toposhift) )

    print("")
    print("  check: ",file1,file2)
    print("")

    cmd = 'gmt gmtinfo ' + file1 + ';'
    cmd += 'gmt gmtinfo ' + file2 + ';'
    print("  > ",cmd)
    status = subprocess.call(cmd, shell=True)
    check_status(status)
    print("")

    print("  number of points along x (NXI) and y (NETA):")
    x0 = data[0,0]
    y0 = data[0,1]
    i0 = 1
    nx = 0; ny = 0
    xmin = x0; xmax = x0
    ymin = y0; ymax = y0

    for i in range(1,len(data)):
      x = data[i,0]
      y = data[i,1]
      dx = x - x0
      x0 = x
      if x < xmin: xmin = x
      if x > xmax: xmax = x
      if y < ymin: ymin = y
      if y > ymax: ymax = y
      if dx < 0.0:
          ii = i + 1
          if nx > 0 and ii - i0 != nx:
              print("  non-regular nx: ",nx,ii-i0,"on line ",i+1)
          nx = ii - i0
          ny += 1
          deltay = y - y0
          y0 = y
          i0 = ii
      else:
          deltax = dx

    ii = len(data) + 1
    if nx > 0 and ii - i0 != nx:
        print("  non-regular nx: ",nx,ii-i0,"on line ",ii)
    nx = ii - i0
    ny += 1
    print("  --------------------------------------------")
    print("  NXI  = ",nx)
    print("  NETA = ",ny)
    print("  xmin/xmax = ",xmin,xmax)
    print("  ymin/ymax = ",ymin,ymax)
    print("  deltax = ",deltax,"average = ",(xmax-xmin)/(nx-1))
    print("  deltay = ",deltay,"average = ",(ymax-ymin)/(ny-1))
    print("  --------------------------------------------")
    print("")
    return nx,ny,deltax,deltay



#
#-----------------------------------------------------------------------------
#

def check_status(status):
    if status != 0:
        print("error: status returned ",status)
        sys.exit(status)
    return

#
#-----------------------------------------------------------------------------
#

def update_Mesh_Par_file(dir,lon_min,lat_min,lon_max,lat_max,nx,ny,dx,dy,xyz_file):
    global datadir
    global utm_zone

    # change working directory back to DATA/
    path = dir + '/' + 'DATA/meshfem3D_files/'
    os.chdir(path)

    print("*******************************")
    print("updating Mesh_Par_file ...")
    print("*******************************")
    print("  working directory: ",os.getcwd())
    print("  min lon/lat      : ",lon_min,"/",lat_min)
    print("  max lon/lat      : ",lon_max,"/",lat_max)
    print("")

    cmd = 'echo "";'
    cmd += 'echo "setting lat min/max = ' + str(lat_min) + ' ' + str(lat_max) + '";'
    cmd += 'echo "        lon min/max = ' + str(lon_min) + ' ' + str(lon_max) + '";'
    cmd += 'echo "        utm zone = ' + str(utm_zone) + '";'
    cmd += 'echo "";'
    cmd += 'sed -i "s:^LATITUDE_MIN .*:LATITUDE_MIN                    = ' + str(lat_min) + ':" Mesh_Par_file' + ';'
    cmd += 'sed -i "s:^LATITUDE_MAX .*:LATITUDE_MAX                    = ' + str(lat_max) + ':" Mesh_Par_file' + ';'
    cmd += 'sed -i "s:^LONGITUDE_MIN .*:LONGITUDE_MIN                   = ' + str(lon_min) + ':" Mesh_Par_file' + ';'
    cmd += 'sed -i "s:^LONGITUDE_MAX .*:LONGITUDE_MAX                   = ' + str(lon_max) + ':" Mesh_Par_file' + ';'
    cmd += 'sed -i "s:^UTM_PROJECTION_ZONE .*:UTM_PROJECTION_ZONE             = ' + str(utm_zone) + ':" Mesh_Par_file' + ';'
    print("  > ",cmd)
    status = subprocess.call(cmd, shell=True)
    check_status(status)
    print("")

    # interfaces.dat file
    nxi = nx
    neta = ny

    dxi = dx
    deta = dy

    if dx < 0.0:
        # negative increment, starts from maximum location
        lon = lon_max
    else:
        # positive increment, starts from minimum location
        lon = lon_min

    if dy < 0.0:
        # negative increment, starts from maximum location
        lat = lat_max
    else:
        # positive increment, starts from minimum location
        lat = lat_min

    # format:
    # #SUPPRESS_UTM_PROJECTION #NXI #NETA #LONG_MIN #LAT_MIN #SPACING_XI #SPACING_ETA
    # .false. 901 901   123.d0 40.0d0 0.01 -0.01
    cmd = 'echo "";'
    cmd += 'echo "interfaces nxi = ' + str(nxi) + ' neta = ' + str(neta) + '";'
    cmd += 'echo "           long = ' + str(lon) + ' lat = ' + str(lat) + '";'
    cmd += 'echo "           spacing_xi = ' + str(dxi) + ' spacing_eta = ' + str(deta) + '";'
    cmd += 'echo "";'
    line = '.false. ' + str(nxi) + ' ' + str(neta) + ' ' + str(lon) + ' ' + str(lat) + ' ' + str(dxi) + ' ' + str(deta)
    cmd += 'sed -i "s:^.false. .*:' + line + ':g" interfaces.dat' + ';'
    print("  > ",cmd)
    status = subprocess.call(cmd, shell=True)
    check_status(status)
    print("")

    # link to topography files
    topodir = '../../' + datadir
    xyz_file1 = xyz_file + '.1.dat'
    xyz_file2 = xyz_file + '.2.dat'
    cmd = 'rm -f ' + xyz_file1 + ' ' + xyz_file2 + ';'
    cmd += 'ln -s ' + topodir + '/' + xyz_file1 + ';'
    cmd += 'ln -s ' + topodir + '/' + xyz_file2 + ';'
    print("  > ",cmd)
    status = subprocess.call(cmd, shell=True)
    check_status(status)
    print("")

    print("  file updated: %s/Mesh_Par_file" % path)
    print("")

    return

#
#-----------------------------------------------------------------------------
#


def update_Par_file(dir):
    global datadir
    global utm_zone

    # change working directory back to DATA/
    path = dir + '/' + 'DATA/'
    os.chdir(path)

    print("*******************************")
    print("updating Par_file ...")
    print("*******************************")
    print("  working directory: ",os.getcwd())
    print("  utm_zone = ",utm_zone)
    print("")

    cmd = 'sed -i "s:^UTM_PROJECTION_ZONE .*:UTM_PROJECTION_ZONE             = ' + str(utm_zone) + ':" Par_file' + ';'
    print("  > ",cmd)
    status = subprocess.call(cmd, shell=True)
    check_status(status)
    print("")

    print("  file updated: %s/Par_file" % path)
    print("")

    return


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

    #---------------------------------------------------------------
    # use conversion to UTM
    iway = ILONGLAT2UTM

    # zone
    UTM_PROJECTION_ZONE = zone

    # lon/lat
    rlon = lon
    rlat = lat

    rx = 0.0
    ry = 0.0

    #---------------------------------------------------------------


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

def convert_lonlat2utm(file_in,zone,file_out):
    """
    converts file with lon/lat/elevation to output file utm_x/utm_y/elevation
    """
    print("converting lon/lat to UTM: " + file_in + " ...")
    print("  zone: %i " % zone)

    # checks argument
    if zone < 1 or zone > 60:  sys.exit("error zone: zone not UTM zone")

    # grab all the locations in file
    with open(file_in,'r') as f:
        content = f.readlines()

    nlines = len(content)
    print("  number of lines: %i" % nlines)
    print("")
    print("  output file: " + file_out)
    print("  format: #UTM_x #UTM_y #elevation")

    # write out converted file
    with open(file_out,'w') as f:
        for line in content:
            # reads lon/lat coordinates
            items = line.split()
            lon = float(items[0])
            lat = float(items[1])
            ele = float(items[2])

            # converts to UTM x/y
            x,y = geo2utm(lon,lat,zone)
            #print("%18.8f\t%18.8f\t%18.8f" % (x,y,ele))
            f.write("%18.8f\t%18.8f\t%18.8f\n" % (x,y,ele))

#
#-----------------------------------------------------------------------------
#


def setup_simulation(lon_min,lat_min,lon_max,lat_max):
    """
    sets up directory for a SPECFEM3D simulation
    """
    global datadir
    global utm_zone
    global incr_dx,toposhift,toposcale

    print("")
    print("*******************************")
    print("setup simulation topography")
    print("*******************************")
    print("")
    print("  topo                  : ",SRTM_type)
    print("  grid sampling interval: ",incr_dx,"(deg) ",incr_dx * math.pi/180.0 * 6371.0, "(km)")
    print("  topo down shift       : ",toposhift)
    print("  topo down scaling     : ",toposcale)
    print("")

    # current directory
    dir = os.getcwd()

    # creates data directory
    cmd = 'mkdir -p ' + datadir
    print("  > ",cmd)
    status = subprocess.call(cmd, shell=True)
    check_status(status)
    print("")

    # change working directory to ./topo_data/
    path = dir + '/' + datadir
    os.chdir(path)
    print("  working directory     : ",os.getcwd())
    print("")

    # check min/max is given in correct order
    if lat_min > lat_max:
        tmp = lat_min
        lat_min = lat_max
        lat_max = tmp

    if lon_min > lon_max:
        tmp = lon_min
        lon_min = lon_max
        lon_max = tmp

    # get topography data
    xyz_file = get_topo(lon_min,lat_min,lon_max,lat_max)

    ## UTM zone
    print("*******************************")
    print("determining UTM coordinates...")
    print("*******************************")

    midpoint_lat = (lat_min + lat_max)/2.0
    midpoint_lon = (lon_min + lon_max)/2.0

    utm_zone = utm.latlon_to_zone_number(midpoint_lat,midpoint_lon)
    print("  region midpoint lat/lon: %f / %f " %(midpoint_lat,midpoint_lon))
    print("  UTM zone: %d" % utm_zone)
    print("")

    # converting to utm
    utm_file = 'ptopo.utm'
    convert_lonlat2utm(xyz_file,utm_zone,utm_file)
    #script version
    #cmd = '../topo/convert_lonlat2utm.py ' + xyz_file + ' ' + str(utm_zone) + ' > ' + utm_file
    #print("converting lon/lat to UTM: ptopo.utm ...")
    #print("> ",cmd)
    #status = subprocess.call(cmd, shell=True)
    #check_status(status)
    print("")

    # AVS UCD file with region borders
    create_AVS_file()

    # extracts interface data for xmeshfem3D
    # uses file with degree increments to determine dx,dy in degrees, as needed for interfaces.dat
    nx,ny,dx,dy = topo_extract(xyz_file)

    # creates parameter files
    # main Par_file
    update_Par_file(dir)

    # mesher Mesh_Par_file
    update_Mesh_Par_file(dir,lon_min,lat_min,lon_max,lat_max,nx,ny,dx,dy,xyz_file)

    print("")
    print("topo output in directory: ",datadir)
    print("")
    print("all done")
    print("")
    return

#
#-----------------------------------------------------------------------------
#

def usage():
    global incr_dx,toposhift,toposcale

    # default increment in km
    incr_dx_km = incr_dx * math.pi/180.0 * 6371.0

    print("usage: ./run_get_simulation_topography.py lon_min lat_min lon_max lat_max [--SRTM=SRTM] [--dx=incr_dx] [--toposhift=toposhift] [--toposcale=toposcale]")
    print("   where")
    print("       lon_min lat_min lon_max lat_max - region given by points: left bottom right top")
    print("                                         for example: 12.35 42.0 12.65 41.8 (Rome)")
    print("       SRTM                            - (optional) name options are:")
    print("                                          'low'  == SRTM 90m ")
    print("                                          'high' == SRTM 30m")
    print("                                          'etopo2' == ETOPO2 (2-arc minutes)")
    print("                                          'topo30' / 'srtm30s' / 'srtm_1km' == SRTM topo 30s (30-arc seconds)")
    print("                                          'topo15' / 'srtm15s' / 'srtm_500m' == SRTM topo 15s (15-arc seconds)")
    print("                                          'topo3' / 'srtm3s' / 'srtm_100m' == SRTM topo 3s (3-arc seconds)")
    print("                                          'topo1' / 'srtm1s' / 'srtm_30m' == SRTM topo 1s (1-arc seconds)")
    print("       incr_dx                         - (optional) GMT grid sampling (in degrees)")
    print("                                         e.g., 0.01 == (~1.1km) [default %f degrees ~ %f km]" % (incr_dx,incr_dx_km))
    print("       toposhift                       - (optional) topography shift for 2nd interface (in m)")
    print("                                         to shift original topography downwards [default %f m]" % toposhift)
    print("       toposcale                       - (optional) scalefactor to topography for shifted 2nd interface (e.g., 0.1) [default %f]" %toposcale)
    sys.exit(1)

if __name__ == '__main__':
    # gets arguments
    if len(sys.argv) < 5:
        usage()

    else:
        lon_min = float(sys.argv[1])
        lat_min = float(sys.argv[2])
        lon_max = float(sys.argv[3])
        lat_max = float(sys.argv[4])

        i = 0
        for arg in sys.argv:
            #print("argument "+str(i)+": " + arg)
            # get arguments
            if "--help" in arg:
                usage()
            elif "--SRTM=" in arg:
                # type: 'low' == SRTM 90m / 'high' == SRTM 30m / else e.g. 'etopo' == topo30 (30-arc seconds)
                SRTM_type = arg.split('=')[1]
            elif "--dx=" in arg:
                # GMT grid sampling interval
                incr_dx = float(arg.split('=')[1])
            elif "--toposhift=" in arg:
                # topography shift for 2nd interface
                toposhift = float(arg.split('=')[1])
            elif "--toposcale=" in arg:
                # topography scalefactor
                toposcale = float(arg.split('=')[1])
            elif i >= 5:
                print("argument not recognized: ",arg)
                usage()
            i += 1

    # logging
    cmd = " ".join(sys.argv)
    filename = './run_get_simulation_topography.log'
    with open(filename, 'a') as f:
      print("command call --- " + str(datetime.datetime.now()),file=f)
      print(cmd,file=f)
      print("command logged to file: " + filename)

    # initializes
    setup_simulation(lon_min,lat_min,lon_max,lat_max)

