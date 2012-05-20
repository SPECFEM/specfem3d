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

! recompute 3D jacobian at a given point for 27-node elements

  subroutine recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
                   xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

  implicit none

  include "constants.h"

  double precision x,y,z,xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
  double precision xi,eta,gamma,jacobian

! coordinates of the control points of the surface element
  double precision xelm(NGNOD),yelm(NGNOD),zelm(NGNOD)

! 3D shape functions and their derivatives at receiver
  double precision shape3D(NGNOD)
  double precision dershape3D(NDIM,NGNOD)

  double precision l1xi,l2xi,l3xi
  double precision l1eta,l2eta,l3eta
  double precision l1gamma,l2gamma,l3gamma
  double precision l1pxi,l2pxi,l3pxi
  double precision l1peta,l2peta,l3peta
  double precision l1pgamma,l2pgamma,l3pgamma

  double precision xxi,yxi,zxi
  double precision xeta,yeta,zeta
  double precision xgamma,ygamma,zgamma

  integer ia

! recompute jacobian for any given (xi,eta,gamma) point
! not necessarily a GLL point

! check that the parameter file is correct
  if(NGNOD /= 27) stop 'elements should have 27 control nodes'

  l1xi=HALF*xi*(xi-ONE)
  l2xi=ONE-xi**2
  l3xi=HALF*xi*(xi+ONE)

  l1pxi=xi-HALF
  l2pxi=-TWO*xi
  l3pxi=xi+HALF

  l1eta=HALF*eta*(eta-ONE)
  l2eta=ONE-eta**2
  l3eta=HALF*eta*(eta+ONE)

  l1peta=eta-HALF
  l2peta=-TWO*eta
  l3peta=eta+HALF

  l1gamma=HALF*gamma*(gamma-ONE)
  l2gamma=ONE-gamma**2
  l3gamma=HALF*gamma*(gamma+ONE)

  l1pgamma=gamma-HALF
  l2pgamma=-TWO*gamma
  l3pgamma=gamma+HALF

! corner nodes

  shape3D(1)=l1xi*l1eta*l1gamma
  shape3D(2)=l3xi*l1eta*l1gamma
  shape3D(3)=l3xi*l3eta*l1gamma
  shape3D(4)=l1xi*l3eta*l1gamma
  shape3D(5)=l1xi*l1eta*l3gamma
  shape3D(6)=l3xi*l1eta*l3gamma
  shape3D(7)=l3xi*l3eta*l3gamma
  shape3D(8)=l1xi*l3eta*l3gamma

  dershape3D(1,1)=l1pxi*l1eta*l1gamma
  dershape3D(1,2)=l3pxi*l1eta*l1gamma
  dershape3D(1,3)=l3pxi*l3eta*l1gamma
  dershape3D(1,4)=l1pxi*l3eta*l1gamma
  dershape3D(1,5)=l1pxi*l1eta*l3gamma
  dershape3D(1,6)=l3pxi*l1eta*l3gamma
  dershape3D(1,7)=l3pxi*l3eta*l3gamma
  dershape3D(1,8)=l1pxi*l3eta*l3gamma

  dershape3D(2,1)=l1xi*l1peta*l1gamma
  dershape3D(2,2)=l3xi*l1peta*l1gamma
  dershape3D(2,3)=l3xi*l3peta*l1gamma
  dershape3D(2,4)=l1xi*l3peta*l1gamma
  dershape3D(2,5)=l1xi*l1peta*l3gamma
  dershape3D(2,6)=l3xi*l1peta*l3gamma
  dershape3D(2,7)=l3xi*l3peta*l3gamma
  dershape3D(2,8)=l1xi*l3peta*l3gamma

  dershape3D(3,1)=l1xi*l1eta*l1pgamma
  dershape3D(3,2)=l3xi*l1eta*l1pgamma
  dershape3D(3,3)=l3xi*l3eta*l1pgamma
  dershape3D(3,4)=l1xi*l3eta*l1pgamma
  dershape3D(3,5)=l1xi*l1eta*l3pgamma
  dershape3D(3,6)=l3xi*l1eta*l3pgamma
  dershape3D(3,7)=l3xi*l3eta*l3pgamma
  dershape3D(3,8)=l1xi*l3eta*l3pgamma

! midside nodes

  shape3D(9)=l2xi*l1eta*l1gamma
  shape3D(10)=l3xi*l2eta*l1gamma
  shape3D(11)=l2xi*l3eta*l1gamma
  shape3D(12)=l1xi*l2eta*l1gamma
  shape3D(13)=l1xi*l1eta*l2gamma
  shape3D(14)=l3xi*l1eta*l2gamma
  shape3D(15)=l3xi*l3eta*l2gamma
  shape3D(16)=l1xi*l3eta*l2gamma
  shape3D(17)=l2xi*l1eta*l3gamma
  shape3D(18)=l3xi*l2eta*l3gamma
  shape3D(19)=l2xi*l3eta*l3gamma
  shape3D(20)=l1xi*l2eta*l3gamma

  dershape3D(1,9)=l2pxi*l1eta*l1gamma
  dershape3D(1,10)=l3pxi*l2eta*l1gamma
  dershape3D(1,11)=l2pxi*l3eta*l1gamma
  dershape3D(1,12)=l1pxi*l2eta*l1gamma
  dershape3D(1,13)=l1pxi*l1eta*l2gamma
  dershape3D(1,14)=l3pxi*l1eta*l2gamma
  dershape3D(1,15)=l3pxi*l3eta*l2gamma
  dershape3D(1,16)=l1pxi*l3eta*l2gamma
  dershape3D(1,17)=l2pxi*l1eta*l3gamma
  dershape3D(1,18)=l3pxi*l2eta*l3gamma
  dershape3D(1,19)=l2pxi*l3eta*l3gamma
  dershape3D(1,20)=l1pxi*l2eta*l3gamma

  dershape3D(2,9)=l2xi*l1peta*l1gamma
  dershape3D(2,10)=l3xi*l2peta*l1gamma
  dershape3D(2,11)=l2xi*l3peta*l1gamma
  dershape3D(2,12)=l1xi*l2peta*l1gamma
  dershape3D(2,13)=l1xi*l1peta*l2gamma
  dershape3D(2,14)=l3xi*l1peta*l2gamma
  dershape3D(2,15)=l3xi*l3peta*l2gamma
  dershape3D(2,16)=l1xi*l3peta*l2gamma
  dershape3D(2,17)=l2xi*l1peta*l3gamma
  dershape3D(2,18)=l3xi*l2peta*l3gamma
  dershape3D(2,19)=l2xi*l3peta*l3gamma
  dershape3D(2,20)=l1xi*l2peta*l3gamma

  dershape3D(3,9)=l2xi*l1eta*l1pgamma
  dershape3D(3,10)=l3xi*l2eta*l1pgamma
  dershape3D(3,11)=l2xi*l3eta*l1pgamma
  dershape3D(3,12)=l1xi*l2eta*l1pgamma
  dershape3D(3,13)=l1xi*l1eta*l2pgamma
  dershape3D(3,14)=l3xi*l1eta*l2pgamma
  dershape3D(3,15)=l3xi*l3eta*l2pgamma
  dershape3D(3,16)=l1xi*l3eta*l2pgamma
  dershape3D(3,17)=l2xi*l1eta*l3pgamma
  dershape3D(3,18)=l3xi*l2eta*l3pgamma
  dershape3D(3,19)=l2xi*l3eta*l3pgamma
  dershape3D(3,20)=l1xi*l2eta*l3pgamma

! side center nodes

  shape3D(21)=l2xi*l2eta*l1gamma
  shape3D(22)=l2xi*l1eta*l2gamma
  shape3D(23)=l3xi*l2eta*l2gamma
  shape3D(24)=l2xi*l3eta*l2gamma
  shape3D(25)=l1xi*l2eta*l2gamma
  shape3D(26)=l2xi*l2eta*l3gamma

  dershape3D(1,21)=l2pxi*l2eta*l1gamma
  dershape3D(1,22)=l2pxi*l1eta*l2gamma
  dershape3D(1,23)=l3pxi*l2eta*l2gamma
  dershape3D(1,24)=l2pxi*l3eta*l2gamma
  dershape3D(1,25)=l1pxi*l2eta*l2gamma
  dershape3D(1,26)=l2pxi*l2eta*l3gamma

  dershape3D(2,21)=l2xi*l2peta*l1gamma
  dershape3D(2,22)=l2xi*l1peta*l2gamma
  dershape3D(2,23)=l3xi*l2peta*l2gamma
  dershape3D(2,24)=l2xi*l3peta*l2gamma
  dershape3D(2,25)=l1xi*l2peta*l2gamma
  dershape3D(2,26)=l2xi*l2peta*l3gamma

  dershape3D(3,21)=l2xi*l2eta*l1pgamma
  dershape3D(3,22)=l2xi*l1eta*l2pgamma
  dershape3D(3,23)=l3xi*l2eta*l2pgamma
  dershape3D(3,24)=l2xi*l3eta*l2pgamma
  dershape3D(3,25)=l1xi*l2eta*l2pgamma
  dershape3D(3,26)=l2xi*l2eta*l3pgamma

! center node

  shape3D(27)=l2xi*l2eta*l2gamma

  dershape3D(1,27)=l2pxi*l2eta*l2gamma
  dershape3D(2,27)=l2xi*l2peta*l2gamma
  dershape3D(3,27)=l2xi*l2eta*l2pgamma

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

  jacobian = xxi*(yeta*zgamma-ygamma*zeta) - xeta*(yxi*zgamma-ygamma*zxi) + &
             xgamma*(yxi*zeta-yeta*zxi)

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

