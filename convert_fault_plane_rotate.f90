
! convert ASCII data based on tsurf to set of CMTSOLUTION files

  program convert_tsurf_to_CMTSOLUTION

  implicit none

  include "constants.h"

! UTM zone for L.A.
  integer, parameter :: UTM_PROJECTION_ZONE = 11

 double precision horiz_dist_fault,time_shift,time_shift_min,time_shift_max

 integer ipoin,ispec,isource,NSOURCES,isourceshiftmin,isource_current
 integer iglob1_store,iglob2_store,iglob3_store,iglob_dummy

 logical exclude_source

 double precision horizdistval,TOLERANCE
 double precision area_current,area_min,area_max,area_sum
 double precision nx,ny,nz,norm
 double precision x_center,y_center
 double precision Mxx,Myy,Mzz,Mxy,Mxz,Myz,long,lat

!! DK DK coordinates of fault top center from Wald et al. 1996
 double precision, parameter :: long_center = -118.55,lat_center = 34.344,depth_center = -5000.

! fault edges
 double precision x_begin,x_end,y_begin,y_end
 double precision xmin,xmax,ymin,ymax

! slip vector Cartesian components
 double precision ex,ey,ez,real_slip_length

 double precision, external :: area_triangle

! for parameter file
  integer NEX_ETA,NEX_XI,NPROC_ETA,NPROC_XI,NSTEP
  double precision DEPTH_BLOCK_KM,LAT_MIN,LAT_MAX,LONG_MIN,LONG_MAX,DT

! definition of fault plane
  integer, parameter :: npoin = 4
  integer, parameter :: nspec = 1

  double precision x1,y1,z1
  double precision x2,y2,z2
  double precision x3,y3,z3
  double precision x4,y4,z4

  double precision xold1,yold1,zold1
  double precision xold2,yold2,zold2
  double precision xold3,yold3,zold3
  double precision xold4,yold4,zold4

  double precision xveryold1,yveryold1,zveryold1
  double precision xveryold2,yveryold2,zveryold2
  double precision xveryold3,yveryold3,zveryold3
  double precision xveryold4,yveryold4,zveryold4

  double precision thetadip,dip
  double precision thetastrike,strike

!! DK DK convert coordinates of top center of fault
  call utm_geo(long_center,lat_center,x_center,y_center,UTM_PROJECTION_ZONE,ILONGLAT2UTM)

!! DK DK
  xveryold1 = 0.
  xveryold2 = 0.
  xveryold3 = 0.
  xveryold4 = 0.

  yveryold1 = -9000.
  yveryold2 = +9000.
  yveryold3 = +9000.
  yveryold4 = -9000.

  zveryold1 = 0.
  zveryold2 = 0.
  zveryold3 = -24000.
  zveryold4 = -24000.

!! DK DK rotation to implement dip of fault
  dip = 40.
  thetadip = - (90. - dip) * PI / 180.
  xold1 =  cos(thetadip)*xveryold1 + sin(thetadip)*zveryold1
  zold1 = -sin(thetadip)*xveryold1 + cos(thetadip)*zveryold1
  yold1 = yveryold1

  xold2 =  cos(thetadip)*xveryold2 + sin(thetadip)*zveryold2
  zold2 = -sin(thetadip)*xveryold2 + cos(thetadip)*zveryold2
  yold2 = yveryold2

  xold3 =  cos(thetadip)*xveryold3 + sin(thetadip)*zveryold3
  zold3 = -sin(thetadip)*xveryold3 + cos(thetadip)*zveryold3
  yold3 = yveryold3

  xold4 =  cos(thetadip)*xveryold4 + sin(thetadip)*zveryold4
  zold4 = -sin(thetadip)*xveryold4 + cos(thetadip)*zveryold4
  yold4 = yveryold4

!! DK DK then rotation to implement strike of fault
  strike = 122.
  thetastrike = + strike * PI / 180.
  x1 =  cos(thetastrike)*xold1 + sin(thetastrike)*yold1   + x_center
  y1 = -sin(thetastrike)*xold1 + cos(thetastrike)*yold1   + y_center
  z1 = zold1     + depth_center

  x2 =  cos(thetastrike)*xold2 + sin(thetastrike)*yold2   + x_center
  y2 = -sin(thetastrike)*xold2 + cos(thetastrike)*yold2   + y_center
  z2 = zold2     + depth_center

  x3 =  cos(thetastrike)*xold3 + sin(thetastrike)*yold3   + x_center
  y3 = -sin(thetastrike)*xold3 + cos(thetastrike)*yold3   + y_center
  z3 = zold3     + depth_center

  x4 =  cos(thetastrike)*xold4 + sin(thetastrike)*yold4   + x_center
  y4 = -sin(thetastrike)*xold4 + cos(thetastrike)*yold4   + y_center
  z4 = zold4     + depth_center

!! DK DK UGLY
  print *,'distance between p1 and p2 = ',sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
  print *,'distance between p3 and p4 = ',sqrt((x3-x4)**2 + (y3-y4)**2 + (z3-z4)**2)

!====================================

! write DX file with normals, to check orientation

! write result to DX file
  open(unit=11,file='DX_plane_4.dx',status='unknown')

! write points
   write(11,*) 'object 1 class array type float rank 1 shape 3 items ',npoin,' data follows'

       write(11,*) sngl(x1),sngl(y1),sngl(z1)
       write(11,*) sngl(x2),sngl(y2),sngl(z2)
       write(11,*) sngl(x3),sngl(y3),sngl(z3)
       write(11,*) sngl(x4),sngl(y4),sngl(z4)

! write elements
   write(11,*) 'object 2 class array type int rank 1 shape 4 items ',nspec,' data follows'

       write(11,*) '0 3 1 2'

   write(11,*) 'attribute "element type" string "quads"'
   write(11,*) 'attribute "ref" string "positions"'
   write(11,*) 'object 3 class array type float rank 0 items ',nspec,' data follows'

       write(11,*) '200'

   write(11,*) 'attribute "dep" string "connections"'
   write(11,*) 'object "irregular positions irregular connections" class field'
   write(11,*) 'component "positions" value 1'
   write(11,*) 'component "connections" value 2'
   write(11,*) 'component "data" value 3'
   write(11,*) 'end'

  close(11)

  end program convert_tsurf_to_CMTSOLUTION

!=====================================================================
!
!  UTM (Universal Transverse Mercator) projection from the USGS
!
!=====================================================================

  subroutine utm_geo(rlon,rlat,rx,ry,UTM_PROJECTION_ZONE,iway)

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

  double precision, parameter :: degrad=PI/180., raddeg=180./PI
  double precision, parameter :: semimaj=6378206.4d0, semimin=6356583.8d0
  double precision, parameter :: scfa=.9996d0
  double precision, parameter :: north=0.d0, east=500000.d0

  double precision e2,e4,e6,ep2,xx,yy,dlat,dlon,zone,cm,cmr,delam
  double precision f1,f2,f3,f4,rm,rn,t,c,a,e1,u,rlat1,dlat1,c1,t1,rn1,r1,d
  double precision rx_save,ry_save,rlon_save,rlat_save

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

