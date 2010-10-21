!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  1 . 4
!               ---------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
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

! recompute 3D jacobian at a given point for a 8-node element

  subroutine recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
                   xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

  implicit none

  include "constants.h"

  double precision x,y,z,xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
  double precision xi,eta,gamma,jacobian

! coordinates of the control points
  double precision xelm(NGNOD),yelm(NGNOD),zelm(NGNOD)

! 3D shape functions and their derivatives at receiver
  double precision shape3D(NGNOD)
  double precision dershape3D(NDIM,NGNOD)

  double precision xxi,yxi,zxi
  double precision xeta,yeta,zeta
  double precision xgamma,ygamma,zgamma
  double precision ra1,ra2,rb1,rb2,rc1,rc2

  integer ia

! for 8-node element
  double precision, parameter :: ONE_EIGHTH = 0.125d0

! recompute jacobian for any (xi,eta,gamma) point, not necessarily a GLL point

! check that the parameter file is correct
  if(NGNOD /= 8) stop 'elements should have 8 control nodes'

! ***
! *** create the 3D shape functions and the Jacobian for an 8-node element
! ***

!--- case of an 8-node 3D element (Dhatt-Touzot p. 115)

  ra1 = one + xi
  ra2 = one - xi

  rb1 = one + eta
  rb2 = one - eta

  rc1 = one + gamma
  rc2 = one - gamma

  shape3D(1) = ONE_EIGHTH*ra2*rb2*rc2
  shape3D(2) = ONE_EIGHTH*ra1*rb2*rc2
  shape3D(3) = ONE_EIGHTH*ra1*rb1*rc2
  shape3D(4) = ONE_EIGHTH*ra2*rb1*rc2
  shape3D(5) = ONE_EIGHTH*ra2*rb2*rc1
  shape3D(6) = ONE_EIGHTH*ra1*rb2*rc1
  shape3D(7) = ONE_EIGHTH*ra1*rb1*rc1
  shape3D(8) = ONE_EIGHTH*ra2*rb1*rc1

  dershape3D(1,1) = - ONE_EIGHTH*rb2*rc2
  dershape3D(1,2) = ONE_EIGHTH*rb2*rc2
  dershape3D(1,3) = ONE_EIGHTH*rb1*rc2
  dershape3D(1,4) = - ONE_EIGHTH*rb1*rc2
  dershape3D(1,5) = - ONE_EIGHTH*rb2*rc1
  dershape3D(1,6) = ONE_EIGHTH*rb2*rc1
  dershape3D(1,7) = ONE_EIGHTH*rb1*rc1
  dershape3D(1,8) = - ONE_EIGHTH*rb1*rc1

  dershape3D(2,1) = - ONE_EIGHTH*ra2*rc2
  dershape3D(2,2) = - ONE_EIGHTH*ra1*rc2
  dershape3D(2,3) = ONE_EIGHTH*ra1*rc2
  dershape3D(2,4) = ONE_EIGHTH*ra2*rc2
  dershape3D(2,5) = - ONE_EIGHTH*ra2*rc1
  dershape3D(2,6) = - ONE_EIGHTH*ra1*rc1
  dershape3D(2,7) = ONE_EIGHTH*ra1*rc1
  dershape3D(2,8) = ONE_EIGHTH*ra2*rc1

  dershape3D(3,1) = - ONE_EIGHTH*ra2*rb2
  dershape3D(3,2) = - ONE_EIGHTH*ra1*rb2
  dershape3D(3,3) = - ONE_EIGHTH*ra1*rb1
  dershape3D(3,4) = - ONE_EIGHTH*ra2*rb1
  dershape3D(3,5) = ONE_EIGHTH*ra2*rb2
  dershape3D(3,6) = ONE_EIGHTH*ra1*rb2
  dershape3D(3,7) = ONE_EIGHTH*ra1*rb1
  dershape3D(3,8) = ONE_EIGHTH*ra2*rb1

! compute coordinates and jacobian matrix
  x=ZERO
  y=ZERO
  z=ZERO
  xxi=ZERO
  xeta=ZERO
  xgamma=ZERO
  yxi=ZERO
  yeta=ZERO
  ygamma=ZERO
  zxi=ZERO
  zeta=ZERO
  zgamma=ZERO

  do ia=1,NGNOD
    x=x+shape3D(ia)*xelm(ia)
    y=y+shape3D(ia)*yelm(ia)
    z=z+shape3D(ia)*zelm(ia)

    xxi=xxi+dershape3D(1,ia)*xelm(ia)
    xeta=xeta+dershape3D(2,ia)*xelm(ia)
    xgamma=xgamma+dershape3D(3,ia)*xelm(ia)
    yxi=yxi+dershape3D(1,ia)*yelm(ia)
    yeta=yeta+dershape3D(2,ia)*yelm(ia)
    ygamma=ygamma+dershape3D(3,ia)*yelm(ia)
    zxi=zxi+dershape3D(1,ia)*zelm(ia)
    zeta=zeta+dershape3D(2,ia)*zelm(ia)
    zgamma=zgamma+dershape3D(3,ia)*zelm(ia)
  enddo

  jacobian = xxi*(yeta*zgamma-ygamma*zeta) - xeta*(yxi*zgamma-ygamma*zxi) + xgamma*(yxi*zeta-yeta*zxi)

  if(jacobian <= ZERO) stop '3D Jacobian undefined'

! invert the relation (Fletcher p. 50 vol. 2)
  xix=(yeta*zgamma-ygamma*zeta)/jacobian
  xiy=(xgamma*zeta-xeta*zgamma)/jacobian
  xiz=(xeta*ygamma-xgamma*yeta)/jacobian
  etax=(ygamma*zxi-yxi*zgamma)/jacobian
  etay=(xxi*zgamma-xgamma*zxi)/jacobian
  etaz=(xgamma*yxi-xxi*ygamma)/jacobian
  gammax=(yxi*zeta-yeta*zxi)/jacobian
  gammay=(xeta*zxi-xxi*zeta)/jacobian
  gammaz=(xxi*yeta-xeta*yxi)/jacobian

  end subroutine recompute_jacobian

