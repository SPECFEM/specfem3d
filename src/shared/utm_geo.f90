!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
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

  subroutine utm_geo(rlon,rlat,rx,ry,UTM_PROJECTION_ZONE,iway,SUPPRESS_UTM_PROJECTION)

! convert geodetic longitude and latitude to UTM, and back
! use iway = ILONGLAT2UTM for long/lat to UTM, IUTM2LONGLAT for UTM to lat/long
! a list of UTM zones of the world is available at www.dmap.co.uk/utmworld.htm

  implicit none

  include "constants.h"

!
!-----CAMx v2.03
!
!     UTM_GEO performs UTM to geodetic (long/lat) translation, and back.
!
!     This is a Fortran version of the BASIC program "Transverse Mercator
!     Conversion", Copyright 1986, Norman J. Berls (Stefan Musarra, 2/94)
!     Based on algorithm taken from "Map Projections Used by the USGS"
!     by John P. Snyder, Geological Survey Bulletin 1532, USDI.
!
!     Input/Output arguments:
!
!        rlon                  Longitude (deg, negative for West)
!        rlat                  Latitude (deg)
!        rx                    UTM easting (m)
!        ry                    UTM northing (m)
!        UTM_PROJECTION_ZONE  UTM zone
!        iway                  Conversion type
!                              ILONGLAT2UTM = geodetic to UTM
!                              IUTM2LONGLAT = UTM to geodetic
!

  integer UTM_PROJECTION_ZONE,iway
  double precision rx,ry,rlon,rlat
  logical SUPPRESS_UTM_PROJECTION

  double precision, parameter :: degrad=PI/180.d0, raddeg=180.d0/PI
  double precision, parameter :: semimaj=6378206.4d0, semimin=6356583.8d0
  double precision, parameter :: scfa=0.9996d0

! some extracts about UTM:
!
! There are 60 longitudinal projection zones numbered 1 to 60 starting at 180°W.
! Each of these zones is 6 degrees wide, apart from a few exceptions around Norway and Svalbard.
! There are 20 latitudinal zones spanning the latitudes 80°S to 84°N and denoted
! by the letters C to X, ommitting the letter O.
! Each of these is 8 degrees south-north, apart from zone X which is 12 degrees south-north.
!
! To change the UTM zone and the hemisphere in which the
! calculations are carried out, need to change the fortran code and recompile. The UTM zone is described
! actually by the central meridian of that zone, i.e. the longitude at the midpoint of the zone, 3 degrees
! from either zone boundary.
! To change hemisphere need to change the "north" variable:
!  - north=0 for northern hemisphere and
!  - north=10000000 (10000km) for southern hemisphere. values must be in metres i.e. north=10000000.
!
! Note that the UTM grids are actually Mercators which
! employ the standard UTM scale factor 0.9996 and set the
! Easting Origin to 500,000;
! the Northing origin in the southern
! hemisphere is kept at 0 rather than set to 10,000,000
! and this gives a uniform scale across the equator if the
! normal convention of selecting the Base Latitude (origin)
! at the equator (0 deg.) is followed.  Northings are
! positive in the northern hemisphere and negative in the
! southern hemisphere.
  double precision, parameter :: north=0.d0
  double precision, parameter :: east=500000.d0

  double precision e2,e4,e6,ep2,xx,yy,dlat,dlon,zone,cm,cmr,delam
  double precision f1,f2,f3,f4,rm,rn,t,c,a,e1,u,rlat1,dlat1,c1,t1,rn1,r1,d
  double precision rx_save,ry_save,rlon_save,rlat_save

  ! checks if conversion to utm has to be done
  if(SUPPRESS_UTM_PROJECTION) then
    if (iway == ILONGLAT2UTM) then
      rx = rlon
      ry = rlat
    else
      rlon = rx
      rlat = ry
    endif
    return
  endif

! save original parameters
  rlon_save = rlon
  rlat_save = rlat
  rx_save = rx
  ry_save = ry

  xx = 0.d0
  yy = 0.d0
  dlat = 0.d0
  dlon = 0.d0

! define parameters of reference ellipsoid
  e2=1.0-(semimin/semimaj)**2.0
  e4=e2*e2
  e6=e2*e4
  ep2=e2/(1.-e2)

  if (iway == IUTM2LONGLAT) then
    xx = rx
    yy = ry
  else
    dlon = rlon
    dlat = rlat
  endif
!
!----- Set Zone parameters
!
  zone = dble(UTM_PROJECTION_ZONE)
  ! sets central meridian for this zone
  cm = zone*6.0 - 183.0
  cmr = cm*degrad
!
!---- Lat/Lon to UTM conversion
!
  if (iway == ILONGLAT2UTM) then

  rlon = degrad*dlon
  rlat = degrad*dlat

  delam = dlon - cm
  if (delam < -180.) delam = delam + 360.
  if (delam > 180.) delam = delam - 360.
  delam = delam*degrad

  f1 = (1. - e2/4. - 3.*e4/64. - 5.*e6/256)*rlat
  f2 = 3.*e2/8. + 3.*e4/32. + 45.*e6/1024.
  f2 = f2*sin(2.*rlat)
  f3 = 15.*e4/256.*45.*e6/1024.
  f3 = f3*sin(4.*rlat)
  f4 = 35.*e6/3072.
  f4 = f4*sin(6.*rlat)
  rm = semimaj*(f1 - f2 + f3 - f4)
  if (dlat == 90. .or. dlat == -90.) then
    xx = 0.
    yy = scfa*rm
  else
    rn = semimaj/sqrt(1. - e2*sin(rlat)**2)
    t = tan(rlat)**2
    c = ep2*cos(rlat)**2
    a = cos(rlat)*delam

    f1 = (1. - t + c)*a**3/6.
    f2 = 5. - 18.*t + t**2 + 72.*c - 58.*ep2
    f2 = f2*a**5/120.
    xx = scfa*rn*(a + f1 + f2)
    f1 = a**2/2.
    f2 = 5. - t + 9.*c + 4.*c**2
    f2 = f2*a**4/24.
    f3 = 61. - 58.*t + t**2 + 600.*c - 330.*ep2
    f3 = f3*a**6/720.
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
  e1 = sqrt(1. - e2)
  e1 = (1. - e1)/(1. + e1)
  rm = yy/scfa
  u = 1. - e2/4. - 3.*e4/64. - 5.*e6/256.
  u = rm/(semimaj*u)

  f1 = 3.*e1/2. - 27.*e1**3./32.
  f1 = f1*sin(2.*u)
  f2 = 21.*e1**2/16. - 55.*e1**4/32.
  f2 = f2*sin(4.*u)
  f3 = 151.*e1**3./96.
  f3 = f3*sin(6.*u)
  rlat1 = u + f1 + f2 + f3
  dlat1 = rlat1*raddeg
  if (dlat1 >= 90. .or. dlat1 <= -90.) then
    dlat1 = dmin1(dlat1,dble(90.) )
    dlat1 = dmax1(dlat1,dble(-90.) )
    dlon = cm
  else
    c1 = ep2*cos(rlat1)**2.
    t1 = tan(rlat1)**2.
    f1 = 1. - e2*sin(rlat1)**2.
    rn1 = semimaj/sqrt(f1)
    r1 = semimaj*(1. - e2)/sqrt(f1**3)
    d = xx/(rn1*scfa)

    f1 = rn1*tan(rlat1)/r1
    f2 = d**2/2.
    f3 = 5.*3.*t1 + 10.*c1 - 4.*c1**2 - 9.*ep2
    f3 = f3*d**2*d**2/24.
    f4 = 61. + 90.*t1 + 298.*c1 + 45.*t1**2. - 252.*ep2 - 3.*c1**2
    f4 = f4*(d**2)**3./720.
    rlat = rlat1 - f1*(f2 - f3 + f4)
    dlat = rlat*raddeg

    f1 = 1. + 2.*t1 + c1
    f1 = f1*d**2*d/6.
    f2 = 5. - 2.*c1 + 28.*t1 - 3.*c1**2 + 8.*ep2 + 24.*t1**2.
    f2 = f2*(d**2)**2*d/120.
    rlon = cmr + (d - f1 + f2)/cos(rlat1)
    dlon = rlon*raddeg
    if (dlon < -180.) dlon = dlon + 360.
    if (dlon > 180.) dlon = dlon - 360.
  endif
  endif

  if (iway == IUTM2LONGLAT) then
    rlon = dlon
    rlat = dlat
    rx = rx_save
    ry = ry_save
  else
    rx = xx
    ry = yy
    rlon = rlon_save
    rlat = rlat_save
  endif

  end subroutine utm_geo

