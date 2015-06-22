!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2008
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

! get backazimuth baz from event and station coordinates the, phe, ths and phs
  subroutine get_backazimuth(the,phe,ths,phs,baz)

  implicit none

  double precision the, phe
  double precision ths, phs
  double precision az,baz,xdeg

  double precision a, a1, b, b1, c, c1
  double precision d, d1, e, e1
  double precision ec2, eps, f, f1, g, g1, h, h1, onemec2, pherad
  double precision phsrad, sc, sd, ss
  double precision temp, therad, thg, thsrad

  double precision, parameter :: rad = 6378.160
  double precision, parameter :: fl = 0.00335293
  double precision, parameter :: twopideg = 360.
  double precision, parameter :: c00 = 1.
  double precision, parameter :: c01 = 0.25
  double precision, parameter :: c02 = -4.6875e-02
  double precision, parameter :: c03 = 1.953125e-02
  double precision, parameter :: c21 = -0.125
  double precision, parameter :: c22 = 3.125e-02
  double precision, parameter :: c23 = -1.46484375e-02
  double precision, parameter :: c42 = -3.90625e-03
  double precision, parameter :: c43 = 2.9296875e-03
  double precision, parameter :: degtokm = 111.3199
  double precision, parameter :: pi = 3.141592654
  double precision, parameter :: TORAD = pi/180.
  double precision, parameter :: TODEG = 1./TORAD


  !=====================================================================
  ! PURPOSE:  To compute the distance and azimuth between locations.
  !=====================================================================
  ! INPUT ARGUMENTS:
  !    THE:     Event latitude in decimal degrees, North positive. [r]
  !    PHE:     Event longitude, East positive. [r]
  !    THS:     Array of station latitudes. [r]
  !    PHS:     Array of station longitudes. [r]
  !    NS:      Length of THS and PHS. [i]
  !=====================================================================
  ! OUTPUT ARGUMENTS:
  !    DIST:    Array of epicentral distances in km. [r]
  !    AZ:      Array of azimuths in degrees. [r]
  !    BAZ:     Array of back azimuths. [r]
  !    XDEG:    Array of great circle arc lengths. [r]
  !    NERR:    Error flag:
  !             =    0   No error.
  !             = 0904   Calculation failed internal consistency checks.
  !=====================================================================
  ! MODULE/LEVEL:  DFM/4
  !=====================================================================
  ! GLOBAL INPUT:
  !    MACH:
  !=====================================================================
  ! SUBROUTINES CALLED:
  !    SACLIB:  SETMSG, APCMSG
  !=====================================================================
  ! LOCAL VARIABLES:
  !=====================================================================
  ! KNOWN ERRORS:
  ! - Problem with equation for distance. See discussion below.
  !=====================================================================
  ! PROCEDURE:
  ! - Calculations are based upon the reference spheroid of 1968 and
  !   are defined by the major radius (RAD) and the flattening (FL).
  ! - Initialize.
  !nerr = 0

  ec2 = 2.*fl - fl*fl
  onemec2 = 1. - ec2
  eps = 1. + ec2/onemec2

  ! - Convert event location to radians.
  !   (Equations are unstable for latidudes of exactly 0 degrees.)

  temp = the
  if( temp == 0. ) temp = 1.0e-08
  therad = TORAD*temp
  pherad = TORAD*phe

  ! - Must convert from geographic to geocentric coordinates in order
  !   to use the spherical trig equations.  This requires a latitude
  !   correction given by: 1-EC2=1-2*FL+FL*FL

  if ( the == 90 .or. the == -90 ) then         ! special attention at the poles
              thg = the*TORAD                   ! ... to avoid division by zero.
  else
              thg = atan( onemec2*tan( therad ) )
  endif

  d = sin( pherad )
  e = -cos( pherad )
  f = -cos( thg )
  c = sin( thg )
  a = f*e
  b = -f*d
  g = -c*e
  h = c*d


  ! -- Convert to radians.
  temp = Ths
  if( temp == 0. ) temp = 1.0e-08
  thsrad = TORAD*temp
  phsrad = TORAD*Phs

  ! -- Calculate some trig constants.
  if ( Ths == 90 .or. Ths == -90 ) then
        thg = Ths * TORAD
  else
        thg = atan( onemec2*tan( thsrad ) )
  endif

  d1 = sin( phsrad )
  e1 = -cos( phsrad )
  f1 = -cos( thg )
  c1 = sin( thg )
  a1 = f1*e1
  b1 = -f1*d1
  g1 = -c1*e1
  h1 = c1*d1
  sc = a*a1 + b*b1 + c*c1

  ! - Spherical trig relationships used to compute angles.

  sd = 0.5*sqrt( ((a - a1)**2 + (b - b1)**2 + (c - &
   c1)**2)*((a + a1)**2 + (b + b1)**2 + (c + c1)**2) )
  Xdeg = atan2( sd, sc )*TODEG
  if( Xdeg < 0. ) &
      Xdeg = Xdeg + twopideg

  ss = (a1 - d)**2 + (b1 - e)**2 + (c1)**2 - 2.
  sc = (a1 - g)**2 + (b1 - h)**2 + (c1 - f)**2 - 2.
  Az = atan2( ss, sc )*TODEG
  if( Az < 0. ) &
      Az = Az + twopideg

  ss = (a - d1)**2 + (b - e1)**2 + (c)**2 - 2.
  sc = (a - g1)**2 + (b - h1)**2 + (c - f1)**2 - 2.
  Baz = atan2( ss, sc )*TODEG
  if( Baz < 0. ) &
      Baz = Baz + twopideg

end subroutine get_backazimuth
