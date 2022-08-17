#!/usr/bin/env python
#
# converts between UTM and lon/lat coordinates
#
# usage example:
# convert from lon/lat to UTM: (zone 31)
#   > ./convert_utm2geo.py 2.6 51.5 31 2
# convert from UTM x/y to lon/lat:
#   > ./convert_utm2geo.py 472234.9543355125 5705505.016786354 31 1
#
from __future__ import print_function

import sys
import math

IUTM2LONGLAT = 1
ILONGLAT2UTM = 2

#
#----------------------------------------------------------------------------------------
#

def utm_geo(lon,lat,zone,iway,ellipsoid=23):
    """
    from utm_geo.f90

    convert geodetic longitude and latitude to UTM, and back
    use iway = ILONGLAT2UTM for long/lat to UTM, IUTM2LONGLAT for UTM to lat/long
    a list of UTM zones of the world is available at www.dmap.co.uk/utmworld.htm

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
    global IUTM2LONGLAT,ILONGLAT2UTM

    PI      = math.pi
    degrad  = PI/180.0
    raddeg  = 180.0/PI

    # default ellipsoid WGS-84
    semimaj = 6378137.0
    semimin = 6356752.314245

    # Clarke 1866
    if ellipsoid == 5:
        semimaj = 6378206.4
        semimin = 6356583.8
    # WGS-84 (World Geodetic System 1984)
    if ellipsoid == 23:
        semimaj = 6378137.0
        semimin = 6356752.314245

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

    # define parameters of reference ellipsoid
    e2 = 1.0 - (semimin/semimaj)**2
    e4 = e2 * e2
    e6 = e2 * e4
    ep2 = e2/(1.0 - e2)

    #---------------------------------------------------------------

    # lon/lat
    if iway == IUTM2LONGLAT:
        xx = lon
        yy = lat
    else:
        dlon = lon
        dlat = lat


    #----- Set Zone parameters
    # zone
    UTM_PROJECTION_ZONE = zone

    lsouth = False
    if UTM_PROJECTION_ZONE < 0: lsouth = True
    zone = abs(UTM_PROJECTION_ZONE)

    cm = zone * 6.0 - 183.0       # set central meridian for this zone
    cmr = cm * degrad

    #---- Lat/Lon to UTM conversion
    if iway == ILONGLAT2UTM:

        rlon = degrad*dlon
        rlat = degrad*dlat

        delam = dlon - cm
        if delam < -180.0: delam = delam + 360.0
        if delam > 180.0: delam = delam - 360.0

        delam = delam*degrad

        # page 61, eq. 3-21 for M
        f1 = (1.0 - e2 / 4.0 - 3.0 * e4 / 64.0 - 5.0 * e6 / 256.0)*rlat
        f2 = 3.0 *e2 / 8.0 + 3.0 * e4 / 32.0 + 45.0 * e6 / 1024.0
        f2 = f2 * math.sin(2.0*rlat)
        # corrected: using .. + 45 e6 / 1024 instead of .. * 45 e6 / 1024
        ##wrong: f3 = 15.0 * e4 / 256.0 * 45.0 * e6 /1024.0
        f3 = 15.0 * e4 / 256.0 + 45.0 * e6 / 1024.0
        f3 = f3 * math.sin(4.0*rlat)
        f4 = 35.0 * e6 / 3072.0
        f4 = f4 * math.sin(6.0*rlat)
        rm = semimaj * (f1 - f2 + f3 - f4)

        if dlat == 90.0 or dlat == -90.0:
            xx = 0.0
            yy = scfa*rm
        else:
            # page 61, eq. 4-20
            rn = semimaj/math.sqrt(1.0 - e2 * math.sin(rlat)**2)
            # page 61, eq. 8-13
            t = math.tan(rlat)**2
            # page 61, eq. 8-14
            c = ep2 * math.cos(rlat)**2
            a = math.cos(rlat) * delam

            # page 61, eq. 8-9 for x
            f1 = (1.0 - t + c) * a**3 / 6.0
            f2 = 5.0 - 18.0*t + t**2 + 72.0*c - 58.0*ep2
            f2 = f2 * a**5 / 120.0
            xx = scfa * rn * (a + f1 + f2)
            # page 61, eq. 8-10 for y
            f1 = a**2 / 2.0
            f2 = 5.0 - t + 9.0 * c + 4.0 * c**2
            f2 = f2 * a**4 / 24.0
            f3 = 61.0 - 58.0*t + t**2 + 600.0*c - 330.0*ep2
            f3 = f3 * a**6 / 720.0
            yy = scfa * (rm + rn * math.tan(rlat) * (f1 + f2 + f3))

        xx = xx + east
        yy = yy + north

    else:
        #---- UTM to Lat/Lon conversion
        xx = xx - east
        yy = yy - north

        ## inverse formulas
        # page 63, eq. 3-24 for e_1
        e1 = math.sqrt(1.0 - e2)
        e1 = (1.0 - e1) / (1.0 + e1)
        rm = yy / scfa
        # page 63, eq. 7-19 for mu
        u = 1.0 - e2 / 4.0 - 3.0 * e4 / 64.0 - 5.0 * e6 / 256.0
        u = rm / (semimaj*u)

        #  page 63, eq. 3-26 for phi_1
        f1 = 3.0 * e1 / 2.0 - 27.0 * e1**3 / 32.0
        f1 = f1 * math.sin(2.0*u)
        f2 = 21.0 * e1**2 / 16.0 - 55.0 * e1**4 / 32.0
        f2 = f2 * math.sin(4.0*u)
        f3 = 151.0 * e1**3 / 96.0
        f3 = f3 * math.sin(6.0*u)
        rlat1 = u + f1 + f2 + f3
        dlat1 = rlat1 * raddeg

        if dlat1 >= 90.0 or dlat1 <= -90.0:
            dlat1 = min(dlat1,90.0)
            dlat1 = max(dlat1,-90.0)
            dlon = cm
        else:
            c1 = ep2 * math.cos(rlat1)**2
            t1 = math.tan(rlat1)**2
            f1 = 1.0 - e2 * math.sin(rlat1)**2
            rn1 = semimaj / math.sqrt(f1)
            r1 = semimaj * (1.0 - e2)/math.sqrt(f1**3)
            d = xx / (rn1*scfa)

            # page 63, eq. 8-17 for phi
            f1 = rn1 * math.tan(rlat1) / r1
            f2 = d**2 / 2.0
            # corrected: factor 5 + 3 T1 + .. instead of 5 * 3 T1 ..
            ##wrong: f3 = 5.0 * 3.0 * t1 + 10.0 * c1 - 4.0 *c1**2 - 9.0 * ep2
            f3 = 5.0 + 3.0 * t1 + 10.0 * c1 - 4.0 *c1**2 - 9.0 * ep2
            f3 = f3 * d**4 / 24.0
            f4 = 61.0 + 90.0 * t1 + 298.0 * c1 + 45.0 * t1**2 - 252.0 * ep2 - 3.0 * c1**2
            f4 = f4 * d**6 / 720.0
            rlat = rlat1 - f1 * (f2 - f3 + f4)
            dlat = rlat * raddeg

            # page 63, eq. 8-18 for lambda
            f1 = 1.0 + 2.0 * t1 + c1
            f1 = f1 * d**3 / 6.0
            f2 = 5.0 - 2.0 * c1 + 28.0 * t1 - 3.0 * c1**2 + 8.0 * ep2 + 24.0 * t1**2
            f2 = f2 * d**5 / 120.0
            rlon = cmr + (d - f1 + f2) / math.cos(rlat1)
            dlon = rlon * raddeg

            if dlon < -180.0: dlon = dlon + 360.0
            if dlon > 180.0: dlon = dlon - 360.0

    if iway == IUTM2LONGLAT:
        rlon = dlon
        rlat = dlat
        return rlon,rlat

    else:
        rx = xx
        if lsouth: yy = yy + 1.e7
        ry = yy
        return rx,ry


#
#----------------------------------------------------------------------------------------
#

# copy from: CUBIT_GEOCUBIT/geocubitlib/LatLongUTMconversion.py
#
# or to avoid copy import directly from geocubit:
#   sys.path.append('../../CUBIT_GEOCUBIT/geocubitlib')
#   from LatLongUTMconversion import LLtoUTM,UTMtoLL

from math import pi, sin, cos, tan, sqrt

# LatLong- UTM conversion..h
# definitions for lat/long to UTM and UTM to lat/lng conversions
# include <string.h>

_deg2rad = pi / 180.0
_rad2deg = 180.0 / pi

_EquatorialRadius = 2
_eccentricitySquared = 3

_ellipsoid = [
    #  id, Ellipsoid name, Equatorial Radius, square of eccentricity
    # first once is a placeholder only, To allow array indices to match id
    # numbers
    [-1, "Placeholder", 0, 0],
    [1, "Airy", 6377563, 0.00667054],
    [2, "Australian National", 6378160, 0.006694542],
    [3, "Bessel 1841", 6377397, 0.006674372],
    [4, "Bessel 1841 (Nambia] ", 6377484, 0.006674372],
    [5, "Clarke 1866", 6378206, 0.006768658],
    [6, "Clarke 1880", 6378249, 0.006803511],
    [7, "Everest", 6377276, 0.006637847],
    [8, "Fischer 1960 (Mercury] ", 6378166, 0.006693422],
    [9, "Fischer 1968", 6378150, 0.006693422],
    [10, "GRS 1967", 6378160, 0.006694605],
    [11, "GRS 1980", 6378137, 0.00669438],
    [12, "Helmert 1906", 6378200, 0.006693422],
    [13, "Hough", 6378270, 0.00672267],
    [14, "International", 6378388, 0.00672267],
    [15, "Krassovsky", 6378245, 0.006693422],
    [16, "Modified Airy", 6377340, 0.00667054],
    [17, "Modified Everest", 6377304, 0.006637847],
    [18, "Modified Fischer 1960", 6378155, 0.006693422],
    [19, "South American 1969", 6378160, 0.006694542],
    [20, "WGS 60", 6378165, 0.006693422],
    [21, "WGS 66", 6378145, 0.006694542],
    [22, "WGS-72", 6378135, 0.006694318],
    [23, "WGS-84", 6378137, 0.00669438]
]

# Reference ellipsoids derived from Peter H. Dana's website-
# http://www.utexas.edu/depts/grg/gcraft/notes/datum/elist.html
# Department of Geography, University of Texas at Austin
# Internet: pdana@mail.utexas.edu
# 3/22/95

# Source
# Defense Mapping Agency. 1987b. DMA Technical Report: Supplement to
# Department of Defense World Geodetic System
# 1984 Technical Report. Part I and II. Washington, DC: Defense Mapping Agency


def LLtoUTM(ReferenceEllipsoid, Lat, Long):
    # converts lat/long to UTM coords.  Equations from USGS Bulletin 1532
    # East Longitudes are positive, West longitudes are negative.
    # North latitudes are positive, South latitudes are negative
    # Lat and Long are in decimal degrees
    # Written by Chuck Gantz- chuck.gantz@globalstar.com

    a = _ellipsoid[ReferenceEllipsoid][_EquatorialRadius]
    eccSquared = _ellipsoid[ReferenceEllipsoid][_eccentricitySquared]
    k0 = 0.9996

# Make sure the longitude is between -180.00 .. 179.9
    LongTemp = (Long + 180) - int((Long + 180) / 360) * \
        360 - 180  # -180.00 .. 179.9

    LatRad = Lat * _deg2rad
    LongRad = LongTemp * _deg2rad

    ZoneNumber = int((LongTemp + 180) / 6) + 1

    if Lat >= 56.0 and Lat < 64.0 and LongTemp >= 3.0 and LongTemp < 12.0:
        ZoneNumber = 32

    # Special zones for Svalbard
    if Lat >= 72.0 and Lat < 84.0:
        if LongTemp >= 0.0 and LongTemp < 9.0:
            ZoneNumber = 31
        elif LongTemp >= 9.0 and LongTemp < 21.0:
            ZoneNumber = 33
        elif LongTemp >= 21.0 and LongTemp < 33.0:
            ZoneNumber = 35
        elif LongTemp >= 33.0 and LongTemp < 42.0:
            ZoneNumber = 37

    # +3 puts origin in middle of zone
    LongOrigin = (ZoneNumber - 1) * 6 - 180 + 3
    LongOriginRad = LongOrigin * _deg2rad

    # compute the UTM Zone from the latitude and longitude
    UTMZone = "%d%c" % (ZoneNumber, _UTMLetterDesignator(Lat))

    eccPrimeSquared = (eccSquared) / (1 - eccSquared)
    N = a / sqrt(1 - eccSquared * sin(LatRad) * sin(LatRad))
    T = tan(LatRad) * tan(LatRad)
    C = eccPrimeSquared * cos(LatRad) * cos(LatRad)
    A = cos(LatRad) * (LongRad - LongOriginRad)

    M = a * ((1 - eccSquared / 4 - 3 * eccSquared * eccSquared / 64 - 5 *
              eccSquared * eccSquared * eccSquared / 256) *
             LatRad - (3 * eccSquared / 8 + 3 * eccSquared * eccSquared / 32 +
                       45 * eccSquared * eccSquared * eccSquared / 1024) *
             sin(2 * LatRad) + (15 * eccSquared * eccSquared / 256 + 45 *
                                eccSquared * eccSquared * eccSquared / 1024) *
             sin(4 * LatRad) -
             (35 * eccSquared * eccSquared * eccSquared / 3072) *
             sin(6 * LatRad))

    UTMEasting = (k0 * N *
                  (A + (1 - T + C) * A * A * A / 6 +
                   (5 - 18 * T + T * T + 72 * C - 58 * eccPrimeSquared) *
                   A * A * A * A * A / 120) +
                  500000.0)

    UTMNorthing = (k0 *
                   (M + N * tan(LatRad) *
                    (A * A / 2 + (5 - T + 9 * C + 4 * C * C) *
                     A * A * A * A / 24 +
                     (61 - 58 * T + T * T + 600 * C - 330 * eccPrimeSquared) *
                     A * A * A * A * A * A / 720)))

    if Lat < 0:
        # 10000000 meter offset for southern hemisphere
        UTMNorthing = UTMNorthing + 10000000.0
    return (UTMZone, UTMEasting, UTMNorthing)


def _UTMLetterDesignator(Lat):
    # This routine determines the correct UTM letter designator
    # for the given latitude
    # returns 'Z' if latitude is outside the UTM limits of 84N to 80S
    # Written by Chuck Gantz- chuck.gantz@globalstar.com

    if 84 >= Lat >= 72:
        return 'X'
    elif 72 > Lat >= 64:
        return 'W'
    elif 64 > Lat >= 56:
        return 'V'
    elif 56 > Lat >= 48:
        return 'U'
    elif 48 > Lat >= 40:
        return 'T'
    elif 40 > Lat >= 32:
        return 'S'
    elif 32 > Lat >= 24:
        return 'R'
    elif 24 > Lat >= 16:
        return 'Q'
    elif 16 > Lat >= 8:
        return 'P'
    elif 8 > Lat >= 0:
        return 'N'
    elif 0 > Lat >= -8:
        return 'M'
    elif -8 > Lat >= -16:
        return 'L'
    elif -16 > Lat >= -24:
        return 'K'
    elif -24 > Lat >= -32:
        return 'J'
    elif -32 > Lat >= -40:
        return 'H'
    elif -40 > Lat >= -48:
        return 'G'
    elif -48 > Lat >= -56:
        return 'F'
    elif -56 > Lat >= -64:
        return 'E'
    elif -64 > Lat >= -72:
        return 'D'
    elif -72 > Lat >= -80:
        return 'C'
    else:
        return 'Z'  # if the Latitude is outside the UTM limits


def UTMtoLL(ReferenceEllipsoid, northing, easting, zone):

    # converts UTM coords to lat/long.  Equations from USGS Bulletin 1532
    # East Longitudes are positive, West longitudes are negative.
    # North latitudes are positive, South latitudes are negative
    # Lat and Long are in decimal degrees.
    # Written by Chuck Gantz- chuck.gantz@globalstar.com
    # Converted to Python by Russ Nelson <nelson@crynwr.com>

    k0 = 0.9996
    a = _ellipsoid[ReferenceEllipsoid][_EquatorialRadius]
    eccSquared = _ellipsoid[ReferenceEllipsoid][_eccentricitySquared]
    e1 = (1 - sqrt(1 - eccSquared)) / (1 + sqrt(1 - eccSquared))
    # NorthernHemisphere; //1 for northern hemispher, 0 for southern

    x = easting - 500000.0  # remove 500,000 meter offset for longitude
    y = northing

    ZoneLetter = zone[-1]
    ZoneNumber = int(zone[:-1])
    if ZoneLetter >= 'N':
        NorthernHemisphere = 1  # point is in northern hemisphere
    else:
        NorthernHemisphere = 0  # point is in southern hemisphere
        y -= 10000000.0         # remove 10,000,000 meter offset used for S

    # +3 puts origin in middle of zone
    LongOrigin = (ZoneNumber - 1) * 6 - 180 + 3

    eccPrimeSquared = (eccSquared) / (1 - eccSquared)

    M = y / k0
    mu = M / (a * (1 - eccSquared / 4 - 3 * eccSquared * eccSquared /
                   64 - 5 * eccSquared * eccSquared * eccSquared / 256))

    phi1Rad = (mu + (3 * e1 / 2 - 27 * e1 * e1 * e1 / 32) * sin(2 * mu) +
               (21 * e1 * e1 / 16 - 55 * e1 * e1 * e1 * e1 / 32) *
               sin(4 * mu) + (151 * e1 * e1 * e1 / 96) * sin(6 * mu))
    phi1 = phi1Rad * _rad2deg

    N1 = a / sqrt(1 - eccSquared * sin(phi1Rad) * sin(phi1Rad))
    T1 = tan(phi1Rad) * tan(phi1Rad)
    C1 = eccPrimeSquared * cos(phi1Rad) * cos(phi1Rad)
    R1 = a * (1 - eccSquared) / pow(1 - eccSquared *
                                    sin(phi1Rad) * sin(phi1Rad), 1.5)
    D = x / (N1 * k0)

    Lat = phi1Rad - (N1 * tan(phi1Rad) / R1) *\
        (D * D / 2 -
         (5 + 3 * T1 + 10 * C1 - 4 * C1 * C1 - 9 * eccPrimeSquared) *
         D * D * D * D / 24 + (61 + 90 * T1 + 298 * C1 + 45 * T1 * T1 - 252 *
                               eccPrimeSquared - 3 * C1 * C1) *
         D * D * D * D * D * D / 720)
    Lat = Lat * _rad2deg

    Long = (D - (1 + 2 * T1 + C1) * D * D * D / 6 +
            (5 - 2 * C1 + 28 * T1 - 3 * C1 * C1 + 8 *
             eccPrimeSquared + 24 * T1 * T1) *
            D * D * D * D * D / 120) / cos(phi1Rad)
    Long = LongOrigin + Long * _rad2deg
    return (Lat, Long)


#
#----------------------------------------------------------------------------------------
#


def usage():
    print("Usage: convert_geo2utm.py [lon/utm_x] [lat/utm_y] [zone] [iway] ")
    print("")
    print("  with: zone - UTM zone")
    print("        iway - 1 == UTM2LONGLAT / 2 == LONGLAT2UTM")
    sys.exit(1)


if __name__ == '__main__':
    # gets input arguments
    if len(sys.argv) != 5: usage()

    lon_utmx = float(sys.argv[1])
    lat_utmy = float(sys.argv[2])
    zone = int(sys.argv[3])
    iway = int(sys.argv[4])

    print("input:")
    if iway == ILONGLAT2UTM:
        lon = lon_utmx
        lat = lat_utmy
        print("  lon   : ",lon)
        print("  lat   : ",lat)
        print("  zone  : ",zone)
        print("  iway  : ",iway," (LONGLAT2UTM)")
        # converts lon/lat to UTM x/y
        x,y = utm_geo(lon,lat,zone,iway)
    else:
        utmx = lon_utmx
        utmy = lat_utmy
        print("  utm_x : ",utmx)
        print("  utm_y : ",utmy)
        print("  zone  : ",zone)
        print("  iway  : ",iway," (UTM2LONGLAT)")
        # converts UTM x/y to lon/lat
        x,y = utm_geo(utmx,utmy,zone,iway)

    print("")
    print("result:")
    #print("  %18.8f\t%18.8f" % (x,y))
    if iway == ILONGLAT2UTM:
        print("  utm_x / utm_y = ",x,y)
    else:
        print("  lon   / lat   = ",x,y)
    print("")

    # from GEOCUBIT
    print("result geocubitlib:")
    if iway == ILONGLAT2UTM:
        #x,y = geo2utm(lon,lat,unit='geo')
        #print("  geocubitlib: geo2utm              utm_x / utm_y = ",x,y)
        ellipsoid = 23
        (zone, x, y) = LLtoUTM(ellipsoid, lat, lon)
        print("  LLtoUTM utm_x / utm_y = ",x,y," zone = ",zone)
    else:
        ellipsoid = 23
        z = str(zone) + "N"
        (lat, lon) = UTMtoLL(ellipsoid, utmy, utmx, z)  # UTMtoLL(ellipsoid, n, e, z)
        print("  UTMtoLL lon / lat = ",lon,lat)
    print("")
