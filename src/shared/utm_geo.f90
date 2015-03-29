!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

!=====================================================================
!
!  UTM (Universal Transverse Mercator) projection from the USGS
!
!=====================================================================

! taken from the open-source package CAMx at http://www.camx.com/download
! converted by Dimitri Komatitsch to Fortran90 and slightly modified to add one parameter to the subroutine call
! and change the order of the arguments for compatibility with SPECFEM calls.
! Also converted to double precision for all calculations and all results.
! Also defined the UTM easting and northing in meters instead of kilometers for compatibility with SPECFEM calls.

! convert geodetic longitude and latitude to UTM, and back
! use iway = ILONGLAT2UTM for long/lat to UTM, IUTM2LONGLAT for UTM to lat/long
! a list of UTM zones of the world is available at www.dmap.co.uk/utmworld.htm

  subroutine utm_geo(rlon4,rlat4,rx4,ry4,UTM_PROJECTION_ZONE,iway,SUPPRESS_UTM_PROJECTION)

!
!---- CAMx v6.10 2014/04/02
!
!     UTM_GEO performs UTM to geodetic (long/lat) translation, and back.
!
!     This is a Fortran version of the BASIC program "Transverse Mercator
!     Conversion", Copyright 1986, Norman J. Berls (Stefan Musarra, 2/94)
!     Based on algorithm taken from "Map Projections Used by the USGS"
!     by John P. Snyder, Geological Survey Bulletin 1532, USDI.
!
!     Portions Copyright 1996 - 2014
!     ENVIRON International Corporation
!
!     Modifications:
!      2012/12/02   Added logic for southern hemisphere UTM zones
!                   Use zones +1 to +60 for the Northern hemisphere, -1 to -60 for the Southern hemisphere
!                   Equator is defined as 0 km North for the Northern hemisphere, 10,000 km North for the Southern hemisphere
!
!     Input/Output arguments:
!
!        rlon4                 Longitude (degrees, negative for West)
!        rlat4                 Latitude (degrees)
!        rx4                   UTM easting (meters)
!        ry4                   UTM northing (meters)
!        UTM_PROJECTION_ZONE   UTM zone
!                              The Northern hemisphere corresponds to zones +1 to +60
!                              The Southern hemisphere corresponds to zones -1 to -60
!        iway                  Conversion type
!                              ILONGLAT2UTM = geodetic to UTM
!                              IUTM2LONGLAT = UTM to geodetic

! Some general information about UTM:
! (for more details see e.g. http://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system )
!
! There are 60 longitudinal projection zones numbered 1 to 60 starting at 180 degrees W.
! Each of these zones is 6 degrees wide, apart from a few exceptions around Norway and Svalbard.
! There are 20 latitudinal zones spanning the latitudes 80 degrees S to 84 degrees N and denoted
! by the letters C to X, ommitting the letter O.
! Each of these is 8 degrees south-north, apart from zone X which is 12 degrees south-north.
!
! The UTM zone is described by the central meridian of that zone, i.e. the longitude at the midpoint of the zone,
! 3 degrees away from both zone boundary.
!
  use constants, only: PI,ILONGLAT2UTM,IUTM2LONGLAT

  implicit none

! input/output parameters
  double precision :: rx4,ry4,rlon4,rlat4
  integer, intent(in) :: UTM_PROJECTION_ZONE,iway
  logical, intent(in) :: SUPPRESS_UTM_PROJECTION

! local parameters

! From http://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system :
! The Universal Transverse Mercator coordinate system was developed by the United States Army Corps of Engineers in the 1940s.
! The system was based on an ellipsoidal model of Earth. For areas within the contiguous United States
! the Clarke Ellipsoid of 1866 was used. For the remaining areas of Earth, including Hawaii, the International Ellipsoid was used.
! The WGS84 ellipsoid is now generally used to model the Earth in the UTM coordinate system,
! which means that current UTM northing at a given point can be 200+ meters different from the old one.
! For different geographic regions, other datum systems (e.g.: ED50, NAD83) can be used.

! Clarke 1866
! double precision, parameter :: SEMI_MAJOR_AXIS = 6378206.4d0, SEMI_MINOR_AXIS = 6356583.8d0

! WGS84 (World Geodetic System 1984)
  double precision, parameter :: SEMI_MAJOR_AXIS = 6378137.0d0, SEMI_MINOR_AXIS = 6356752.314245d0

! Note that the UTM grids are actually Mercators which
! employ the standard UTM scale factor 0.9996 and set the Easting Origin to 500,000.
  double precision, parameter :: scfa=0.9996d0
  double precision, parameter :: north=0.d0, east=500000.d0

  double precision, parameter :: DEGREES_TO_RADIANS=PI/180.d0, RADIANS_TO_DEGREES=180.d0/PI

! local variables
  integer :: zone
  double precision :: rlon,rlat
  double precision :: e2,e4,e6,ep2,xx,yy,dlat,dlon,cm,cmr,delam
  double precision :: f1,f2,f3,f4,rm,rn,t,c,a,e1,u,rlat1,dlat1,c1,t1,rn1,r1,d
  double precision :: rx_save,ry_save,rlon_save,rlat_save
  logical :: lsouth

  ! checks if conversion to UTM has to be done
  if (SUPPRESS_UTM_PROJECTION) then
    if (iway == ILONGLAT2UTM) then
      rx4 = rlon4
      ry4 = rlat4
    else
      rlon4 = rx4
      rlat4 = ry4
    endif
    return
  endif

! save original parameters
  rlon_save = rlon4
  rlat_save = rlat4
  rx_save = rx4
  ry_save = ry4

  e2=1.d0-(SEMI_MINOR_AXIS/SEMI_MAJOR_AXIS)**2
  e4=e2*e2
  e6=e2*e4
  ep2=e2/(1.d0-e2)

!
!---- Set Zone parameters
!

  lsouth = .false.
  if (UTM_PROJECTION_ZONE < 0) lsouth = .true.
  zone = abs(UTM_PROJECTION_ZONE)
  cm = zone*6.0d0 - 183.d0
  cmr = cm*DEGREES_TO_RADIANS

  if (iway == IUTM2LONGLAT) then
    xx = rx4
    yy = ry4
    if (lsouth) yy = yy - 1.d7
  else
    dlat = rlat4
    dlon = rlon4
  endif

!
!---- Lat/Lon to UTM conversion
!
  if (iway == ILONGLAT2UTM) then

    rlon = DEGREES_TO_RADIANS*dlon
    rlat = DEGREES_TO_RADIANS*dlat

    delam = dlon - cm
    if (delam < -180.d0) delam = delam + 360.d0
    if (delam > 180.d0) delam = delam - 360.d0
    delam = delam*DEGREES_TO_RADIANS

    f1 = (1.d0 - e2/4.d0 - 3.d0*e4/64.d0 - 5.d0*e6/256d0)*rlat
    f2 = 3.d0*e2/8.d0 + 3.d0*e4/32.d0 + 45.d0*e6/1024.d0
    f2 = f2*sin(2.d0*rlat)
    f3 = 15.d0*e4/256.d0*45.d0*e6/1024.d0
    f3 = f3*sin(4.d0*rlat)
    f4 = 35.d0*e6/3072.d0
    f4 = f4*sin(6.d0*rlat)
    rm = SEMI_MAJOR_AXIS*(f1 - f2 + f3 - f4)
    if (dlat == 90.d0 .or. dlat == -90.d0) then
      xx = 0.d0
      yy = scfa*rm
    else
      rn = SEMI_MAJOR_AXIS/sqrt(1.d0 - e2*sin(rlat)**2)
      t = tan(rlat)**2
      c = ep2*cos(rlat)**2
      a = cos(rlat)*delam

      f1 = (1.d0 - t + c)*a**3/6.d0
      f2 = 5.d0 - 18.d0*t + t**2 + 72.d0*c - 58.d0*ep2
      f2 = f2*a**5/120.d0
      xx = scfa*rn*(a + f1 + f2)
      f1 = a**2/2.d0
      f2 = 5.d0 - t + 9.d0*c + 4.d0*c**2
      f2 = f2*a**4/24.d0
      f3 = 61.d0 - 58.d0*t + t**2 + 600.d0*c - 330.d0*ep2
      f3 = f3*a**6/720.d0
      yy = scfa*(rm + rn*tan(rlat)*(f1 + f2 + f3))
    endif
    xx = xx + east
    yy = yy + north

!
!---- UTM to Lat/Lon conversion
!
  else

    xx = xx - east
    yy = yy - north
    e1 = sqrt(1.d0 - e2)
    e1 = (1.d0 - e1)/(1.d0 + e1)
    rm = yy/scfa
    u = 1.d0 - e2/4.d0 - 3.d0*e4/64.d0 - 5.d0*e6/256.d0
    u = rm/(SEMI_MAJOR_AXIS*u)

    f1 = 3.d0*e1/2.d0 - 27.d0*e1**3.d0/32.d0
    f1 = f1*sin(2.d0*u)
    f2 = 21.d0*e1**2/16.d0 - 55.d0*e1**4/32.d0
    f2 = f2*sin(4.d0*u)
    f3 = 151.d0*e1**3.d0/96.d0
    f3 = f3*sin(6.d0*u)
    rlat1 = u + f1 + f2 + f3
    dlat1 = rlat1*RADIANS_TO_DEGREES
    if (dlat1 >= 90.d0 .or. dlat1 <= -90.d0) then
      dlat1 = dmin1(dlat1,90.d0)
      dlat1 = dmax1(dlat1,-90.d0)
      dlon = cm
    else
      c1 = ep2*cos(rlat1)**2
      t1 = tan(rlat1)**2
      f1 = 1.d0 - e2*sin(rlat1)**2
      rn1 = SEMI_MAJOR_AXIS/sqrt(f1)
      r1 = SEMI_MAJOR_AXIS*(1.d0 - e2)/sqrt(f1**3)
      d = xx/(rn1*scfa)

      f1 = rn1*tan(rlat1)/r1
      f2 = d**2/2.d0
      f3 = 5.d0*3.d0*t1 + 10.d0*c1 - 4.d0*c1**2 - 9.d0*ep2
      f3 = f3*d**2*d**2/24.d0
      f4 = 61.d0 + 90.d0*t1 + 298.d0*c1 + 45.d0*t1**2 - 252.d0*ep2 - 3.d0*c1**2
      f4 = f4*(d**2)**3.d0/720.d0
      rlat = rlat1 - f1*(f2 - f3 + f4)
      dlat = rlat*RADIANS_TO_DEGREES

      f1 = 1.d0 + 2.d0*t1 + c1
      f1 = f1*d**2*d/6.d0
      f2 = 5.d0 - 2.d0*c1 + 28.d0*t1 - 3.d0*c1**2 + 8.d0*ep2 + 24.d0*t1**2
      f2 = f2*(d**2)**2*d/120.d0
      rlon = cmr + (d - f1 + f2)/cos(rlat1)
      dlon = rlon*RADIANS_TO_DEGREES
      if (dlon < -180.d0) dlon = dlon + 360.d0
      if (dlon > 180.d0) dlon = dlon - 360.d0
    endif
  endif

!
!----- output
!
  if (iway == IUTM2LONGLAT) then
    rlon4 = dlon
    rlat4 = dlat
    rx4 = rx_save
    ry4 = ry_save
  else
    rx4 = xx
    if (lsouth) yy = yy + 1.d7
    ry4 = yy
    rlon4 = rlon_save
    rlat4 = rlat_save
  endif

  end subroutine utm_geo

