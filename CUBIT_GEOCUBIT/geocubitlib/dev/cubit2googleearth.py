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
#
try:
    cubit.cmd('comment')
except:
    import cubit

KMLstart='<kml xmlns="http://earth.google.com/kml/2.2">\n<Folder><name>prova</name>\n'
placemarkstart='<Placemark>\n      <name></name>\n      <description></description>\n<TimeStamp><begin>%s</begin><end>%s</end></TimeStamp>\n        <MultiGeometry>\n            <altitudeMode>absolute</altitudeMode>\n'
placemarkend='        </MultiGeometry>\n</Placemark>\n'
templatepoly='                <Polygon>\n                    <altitudeMode>absolute</altitudeMode>\n                    <outerBoundaryIs>\n                        <LinearRing>\n                            <coordinates>\n                            %s,%s,%s\n                            %s,%s,%s\n                            %s,%s,%s\n                            %s,%s,%s\n                            %s,%s,%s \n                            </coordinates>\n                        </LinearRing>\n                    </outerBoundaryIs>\n                </Polygon>\n'

KMLend='    </Folder>\n</kml>\n'


import datetime
now=datetime.datetime.now().year

import random

# Lat Long - UTM, UTM - Lat Long conversions

from math import pi, sin, cos, tan, sqrt

#LatLong- UTM conversion..h
#definitions for lat/long to UTM and UTM to lat/lng conversions
#include <string.h>

_deg2rad = pi / 180.0
_rad2deg = 180.0 / pi

_EquatorialRadius = 2
_eccentricitySquared = 3

_ellipsoid = [[ -1, "Placeholder", 0, 0],[ 1, "Airy", 6377563, 0.00667054],[ 2, "Australian National", 6378160, 0.006694542],[ 3, "Bessel 1841", 6377397, 0.006674372],[ 4, "Bessel 1841 (Nambia] ", 6377484, 0.006674372],[ 5, "Clarke 1866", 6378206, 0.006768658],[ 6, "Clarke 1880", 6378249, 0.006803511],[ 7, "Everest", 6377276, 0.006637847],[ 8, "Fischer 1960 (Mercury] ", 6378166, 0.006693422],[ 9, "Fischer 1968", 6378150, 0.006693422],[ 10, "GRS 1967", 6378160, 0.006694605],[ 11, "GRS 1980", 6378137, 0.00669438],[ 12, "Helmert 1906", 6378200, 0.006693422],[ 13, "Hough", 6378270, 0.00672267],[ 14, "International", 6378388, 0.00672267],[ 15, "Krassovsky", 6378245, 0.006693422],[ 16, "Modified Airy", 6377340, 0.00667054],[ 17, "Modified Everest", 6377304, 0.006637847],[ 18, "Modified Fischer 1960", 6378155, 0.006693422],[ 19, "South American 1969", 6378160, 0.006694542],[ 20, "WGS 60", 6378165, 0.006693422],[ 21, "WGS 66", 6378145, 0.006694542],[ 22, "WGS-72", 6378135, 0.006694318],[ 23, "WGS-84", 6378137, 0.00669438]]

#Reference ellipsoids derived from Peter H. Dana's website- 
#http://www.utexas.edu/depts/grg/gcraft/notes/datum/elist.htmlV
#Department of Geography, University of Texas at Austin
#Internet: pdana@mail.utexas.edu
#3/22/95

#Source
#Defense Mapping Agency. 1987b. DMA Technical Report: Supplement to Department of Defense World Geodetic System
#1984 Technical Report. Part I and II. Washington, DC: Defense Mapping Agency

def LLtoUTM(ReferenceEllipsoid, Lat, Long):
    #converts lat/long to UTM coords.  Equations from USGS Bulletin 1532 
    #East Longitudes are positive, West longitudes are negative. 
    #North latitudes are positive, South latitudes are negative
    #Lat and Long are in decimal degrees
    #Written by Chuck Gantz- chuck.gantz@globalstar.com
    a = _ellipsoid[ReferenceEllipsoid][_EquatorialRadius]
    eccSquared = _ellipsoid[ReferenceEllipsoid][_eccentricitySquared]
    k0 = 0.9996#Make sure the longitude is between -180.00 .. 179.9
    LongTemp = (Long+180)-int((Long+180)/360)*360-180 # -180.00 .. 179.9
    LatRad = Lat*_deg2rad
    LongRad = LongTemp*_deg2rad
    ZoneNumber = int((LongTemp + 180)/6) + 1
    if Lat >= 56.0 and Lat < 64.0 and LongTemp >= 3.0 and LongTemp < 12.0:
        ZoneNumber = 32
    # Special zones for Svalbard
    if Lat >= 72.0 and Lat < 84.0:
        if  LongTemp >= 0.0  and LongTemp <  9.0:ZoneNumber = 31
        elif LongTemp >= 9.0  and LongTemp < 21.0: ZoneNumber = 33
        elif LongTemp >= 21.0 and LongTemp < 33.0: ZoneNumber = 35
        elif LongTemp >= 33.0 and LongTemp < 42.0: ZoneNumber = 37
    LongOrigin = (ZoneNumber - 1)*6 - 180 + 3 #+3 puts origin in middle of zone
    LongOriginRad = LongOrigin * _deg2rad
    #compute the UTM Zone from the latitude and longitude
    UTMZone = "%d%c" % (ZoneNumber, _UTMLetterDesignator(Lat))
    eccPrimeSquared = (eccSquared)/(1-eccSquared)
    N = a/sqrt(1-eccSquared*sin(LatRad)*sin(LatRad))
    T = tan(LatRad)*tan(LatRad)
    C = eccPrimeSquared*cos(LatRad)*cos(LatRad)
    A = cos(LatRad)*(LongRad-LongOriginRad)
    M = a*((1
            - eccSquared/4
            - 3*eccSquared*eccSquared/64
            - 5*eccSquared*eccSquared*eccSquared/256)*LatRad 
           - (3*eccSquared/8
              + 3*eccSquared*eccSquared/32
              + 45*eccSquared*eccSquared*eccSquared/1024)*sin(2*LatRad)
           + (15*eccSquared*eccSquared/256 + 45*eccSquared*eccSquared*eccSquared/1024)*sin(4*LatRad) 
           - (35*eccSquared*eccSquared*eccSquared/3072)*sin(6*LatRad))
    UTMEasting = (k0*N*(A+(1-T+C)*A*A*A/6
                        + (5-18*T+T*T+72*C-58*eccPrimeSquared)*A*A*A*A*A/120)
                  + 500000.0)
    UTMNorthing = (k0*(M+N*tan(LatRad)*(A*A/2+(5-T+9*C+4*C*C)*A*A*A*A/24
                                        + (61
                                           -58*T
                                           +T*T
                                           +600*C
                                           -330*eccPrimeSquared)*A*A*A*A*A*A/720)))
    if Lat < 0:
        UTMNorthing = UTMNorthing + 10000000.0; #10000000 meter offset for southern hemisphere
    return (UTMZone, UTMEasting, UTMNorthing)


def _UTMLetterDesignator(Lat):
#This routine determines the correct UTM letter designator for the given latitude
#returns 'Z' if latitude is outside the UTM limits of 84N to 80S
#Written by Chuck Gantz- chuck.gantz@globalstar.com
    if 84 >= Lat >= 72: return 'X'
    elif 72 > Lat >= 64: return 'W'
    elif 64 > Lat >= 56: return 'V'
    elif 56 > Lat >= 48: return 'U'
    elif 48 > Lat >= 40: return 'T'
    elif 40 > Lat >= 32: return 'S'
    elif 32 > Lat >= 24: return 'R'
    elif 24 > Lat >= 16: return 'Q'
    elif 16 > Lat >= 8: return 'P'
    elif  8 > Lat >= 0: return 'N'
    elif  0 > Lat >= -8: return 'M'
    elif -8> Lat >= -16: return 'L'
    elif -16 > Lat >= -24: return 'K'
    elif -24 > Lat >= -32: return 'J'
    elif -32 > Lat >= -40: return 'H'
    elif -40 > Lat >= -48: return 'G'
    elif -48 > Lat >= -56: return 'F'
    elif -56 > Lat >= -64: return 'E'
    elif -64 > Lat >= -72: return 'D'
    elif -72 > Lat >= -80: return 'C'
    else: return 'Z'	# if the Latitude is outside the UTM limits

#void UTMtoLL(int ReferenceEllipsoid, const double UTMNorthing, const double UTMEasting, const char* UTMZone,
#			  double& Lat,  double& Long )

def UTMtoLL(ReferenceEllipsoid, northing, easting, zone):
    #converts UTM coords to lat/long.  Equations from USGS Bulletin 1532 
    #East Longitudes are positive, West longitudes are negative. 
    #North latitudes are positive, South latitudes are negative
    #Lat and Long are in decimal degrees. 
    #Written by Chuck Gantz- chuck.gantz@globalstar.com
    #Converted to Python by Russ Nelson <nelson@crynwr.com>
    k0 = 0.9996
    a = _ellipsoid[ReferenceEllipsoid][_EquatorialRadius]
    eccSquared = _ellipsoid[ReferenceEllipsoid][_eccentricitySquared]
    e1 = (1-sqrt(1-eccSquared))/(1+sqrt(1-eccSquared))
    #NorthernHemisphere; //1 for northern hemispher, 0 for southern
    x = easting - 500000.0 #remove 500,000 meter offset for longitude
    y = northing
    ZoneLetter = zone[-1]
    ZoneNumber = int(zone[:-1])
    if ZoneLetter >= 'N':
        NorthernHemisphere = 1  # point is in northern hemisphere
    else:
        NorthernHemisphere = 0  # point is in southern hemisphere
        y -= 10000000.0         # remove 10,000,000 meter offset used for southern hemisphere
    LongOrigin = (ZoneNumber - 1)*6 - 180 + 3  # +3 puts origin in middle of zone
    eccPrimeSquared = (eccSquared)/(1-eccSquared)
    M = y / k0
    mu = M/(a*(1-eccSquared/4-3*eccSquared*eccSquared/64-5*eccSquared*eccSquared*eccSquared/256))
    phi1Rad = (mu + (3*e1/2-27*e1*e1*e1/32)*sin(2*mu) 
               + (21*e1*e1/16-55*e1*e1*e1*e1/32)*sin(4*mu)
               +(151*e1*e1*e1/96)*sin(6*mu))
    phi1 = phi1Rad*_rad2deg;
    N1 = a/sqrt(1-eccSquared*sin(phi1Rad)*sin(phi1Rad))
    T1 = tan(phi1Rad)*tan(phi1Rad)
    C1 = eccPrimeSquared*cos(phi1Rad)*cos(phi1Rad)
    R1 = a*(1-eccSquared)/pow(1-eccSquared*sin(phi1Rad)*sin(phi1Rad), 1.5)
    D = x/(N1*k0)
    Lat = phi1Rad - (N1*tan(phi1Rad)/R1)*(D*D/2-(5+3*T1+10*C1-4*C1*C1-9*eccPrimeSquared)*D*D*D*D/24
                                          +(61+90*T1+298*C1+45*T1*T1-252*eccPrimeSquared-3*C1*C1)*D*D*D*D*D*D/720)
    Lat = Lat * _rad2deg
    Long = (D-(1+2*T1+C1)*D*D*D/6+(5-2*C1+28*T1-3*C1*C1+8*eccPrimeSquared+24*T1*T1)
            *D*D*D*D*D/120)/cos(phi1Rad)
    Long = LongOrigin + Long * _rad2deg
    return (Lat, Long)

class cubit2GE(object):
    def __init__(self,unit='utm',x_origin=None,y_origin=None,z_origin=None):
        super(cubit2GE, self).__init__()
        self.unit=unit
        self.xo=x_origin
        self.yo=y_origin
        self.zo=z_origin
        cubit.cmd('compress all')
    def __repr__(self):
        pass
    def create_blocks(self,list_entity=None):
        txt=' face in surface '                                 
        if list_entity:                                         
            if not isinstance(list_entity,list):                
                list_entity=[list_entity]
        iblock=cubit.get_next_block_id()                       
        for entity in list_entity:                                                  
            command = "block "+str(iblock)+ txt +str(entity)    
            cubit.cmd(command)
        command = "block "+str(iblock)+ "name 'toGE'"
        cubit.cmd(command)
        self.iblock=iblock 
    def geo2utm(self,lon,lat,unit,ellipsoid=23):
        if unit == 'geo' :          
           (zone, x, y) = LLtoUTM(ellipsoid, lat, lon)
        elif unit == 'utm' : 
           x=lon
           y=lat
           zone=0
        return zone,x,y
    def utm2geo(self,lon,lat,zone,unit='utm'):
        if unit == 'utm' and zone:          
           (y, x) = UTMtoLL(23, lat, lon,zone)
        else: 
           x=lon
           y=lat
        return x,y                            
    def write(self,GEname=None,t1=0,t2=1):
        cubit.cmd('set info off')
        cubit.cmd('set echo off')
        cubit.cmd('set journal off')
        from sets import Set
        if not GEname: GEname='cubit2GE.kml'
        if self.xo and self.yo:
            zone_utm,xo_utm,yo_utm=self.geo2utm(self.xo,self.yo,'geo')
        else:
            xo_utm=0.
            yo_utm=0.
            zone_utm=None
        if self.zo:
            zo=self.zo
        else:
            zo=0
        #
        GE=open(GEname,'w')
        print 'Writing '+GEname+'.....'
        #
        #
        GE.write(KMLstart)
        for it in range(int(t1),int(t2),1):
            tbegin=now+seconds*it
            tafter=now+(seconds*(it+1))
            print tbegin.isoformat()
            GE.write(placemarkstart%(str(tbegin.isoformat()),str(tafter.isoformat())))
            #
            block=self.iblock
            nameb=cubit.get_exodus_entity_name('block',block)
            print nameb,block
            quads_all=cubit.get_block_faces(block)
            print len(quads_all)
            for f in quads_all:
                nodes=cubit.get_connectivity('Face',f)
                coord=[]
                for n in [nodes[0],nodes[1],nodes[2],nodes[3],nodes[0]]:
                    x,y,z=cubit.get_nodal_coordinates(n)
                    if self.unit == 'geo':
                        zone,x,y=self.geo2utm(x,y,'geo')
                    new_x=x+xo_utm+sin(it)*3
                    new_y=y+yo_utm+cos(it)*3
                    new_z=z+zo+sin(it)*3
                    lon,lat=self.utm2geo(new_x,new_y,zone_utm,unit='utm')
                    coord.append(str(lon))
                    coord.append(str(lat))
                    coord.append(str(new_z))
                GE.write(templatepoly % tuple(coord))
            GE.write(placemarkend)
        GE.write(KMLend)
        GE.close()
        #
        print 'Ok'
        cubit.cmd('set info on')
        cubit.cmd('set echo on')

def toGE(name='c2GE.kml',list_entity=None,unit='utm',x_origin=None,y_origin=None,z_origin=None,t1=0,t2=5):
    GE=cubit2GE(unit=unit,x_origin=x_origin,y_origin=y_origin,z_origin=z_origin)
    GE.create_blocks(list_entity=list_entity)
    GE.write(GEname=name,t1=t1,t2=t2)

name='/Users/emanuele/Desktop/c2GE.kml'
toGE(name,list_entity='1 to 493',x_origin=8.767022,y_origin= 44.460383, z_origin=220)



