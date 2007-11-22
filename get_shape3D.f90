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

! 3D shape functions for 8-node element

  subroutine get_shape3D(myrank,shape3D,dershape3D,xigll,yigll,zigll)

  implicit none

  include "constants.h"

  integer myrank

! Gauss-Lobatto-Legendre points of integration
  double precision xigll(NGLLX)
  double precision yigll(NGLLY)
  double precision zigll(NGLLZ)

! 3D shape functions and their derivatives
  double precision shape3D(NGNOD,NGLLX,NGLLY,NGLLZ)
  double precision dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ)

  integer i,j,k,ia

! location of the nodes of the 3D quadrilateral elements
  double precision xi,eta,gamma
  double precision ra1,ra2,rb1,rb2,rc1,rc2

! for checking the 3D shape functions
  double precision sumshape,sumdershapexi,sumdershapeeta,sumdershapegamma

  double precision, parameter :: ONE_EIGHTH = 0.125d0

! check that the parameter file is correct
  if(NGNOD /= 8) call exit_MPI(myrank,'elements should have 8 control nodes')

! ***
! *** create 3D shape functions and jacobian
! ***

!--- case of a 3D 8-node element (Dhatt-Touzot p. 115)

  do i=1,NGLLX
  do j=1,NGLLY
  do k=1,NGLLZ

  xi = xigll(i)
  eta = yigll(j)
  gamma = zigll(k)

  ra1 = one + xi
  ra2 = one - xi

  rb1 = one + eta
  rb2 = one - eta

  rc1 = one + gamma
  rc2 = one - gamma

  shape3D(1,i,j,k) = ONE_EIGHTH*ra2*rb2*rc2
  shape3D(2,i,j,k) = ONE_EIGHTH*ra1*rb2*rc2
  shape3D(3,i,j,k) = ONE_EIGHTH*ra1*rb1*rc2
  shape3D(4,i,j,k) = ONE_EIGHTH*ra2*rb1*rc2
  shape3D(5,i,j,k) = ONE_EIGHTH*ra2*rb2*rc1
  shape3D(6,i,j,k) = ONE_EIGHTH*ra1*rb2*rc1
  shape3D(7,i,j,k) = ONE_EIGHTH*ra1*rb1*rc1
  shape3D(8,i,j,k) = ONE_EIGHTH*ra2*rb1*rc1

  dershape3D(1,1,i,j,k) = - ONE_EIGHTH*rb2*rc2
  dershape3D(1,2,i,j,k) = ONE_EIGHTH*rb2*rc2
  dershape3D(1,3,i,j,k) = ONE_EIGHTH*rb1*rc2
  dershape3D(1,4,i,j,k) = - ONE_EIGHTH*rb1*rc2
  dershape3D(1,5,i,j,k) = - ONE_EIGHTH*rb2*rc1
  dershape3D(1,6,i,j,k) = ONE_EIGHTH*rb2*rc1
  dershape3D(1,7,i,j,k) = ONE_EIGHTH*rb1*rc1
  dershape3D(1,8,i,j,k) = - ONE_EIGHTH*rb1*rc1

  dershape3D(2,1,i,j,k) = - ONE_EIGHTH*ra2*rc2
  dershape3D(2,2,i,j,k) = - ONE_EIGHTH*ra1*rc2
  dershape3D(2,3,i,j,k) = ONE_EIGHTH*ra1*rc2
  dershape3D(2,4,i,j,k) = ONE_EIGHTH*ra2*rc2
  dershape3D(2,5,i,j,k) = - ONE_EIGHTH*ra2*rc1
  dershape3D(2,6,i,j,k) = - ONE_EIGHTH*ra1*rc1
  dershape3D(2,7,i,j,k) = ONE_EIGHTH*ra1*rc1
  dershape3D(2,8,i,j,k) = ONE_EIGHTH*ra2*rc1

  dershape3D(3,1,i,j,k) = - ONE_EIGHTH*ra2*rb2
  dershape3D(3,2,i,j,k) = - ONE_EIGHTH*ra1*rb2
  dershape3D(3,3,i,j,k) = - ONE_EIGHTH*ra1*rb1
  dershape3D(3,4,i,j,k) = - ONE_EIGHTH*ra2*rb1
  dershape3D(3,5,i,j,k) = ONE_EIGHTH*ra2*rb2
  dershape3D(3,6,i,j,k) = ONE_EIGHTH*ra1*rb2
  dershape3D(3,7,i,j,k) = ONE_EIGHTH*ra1*rb1
  dershape3D(3,8,i,j,k) = ONE_EIGHTH*ra2*rb1

  enddo
  enddo
  enddo

!--- check the shape functions and their derivatives

  do i=1,NGLLX
  do j=1,NGLLY
  do k=1,NGLLZ

  sumshape = ZERO
  sumdershapexi = ZERO
  sumdershapeeta = ZERO
  sumdershapegamma = ZERO

  do ia=1,NGNOD
    sumshape = sumshape + shape3D(ia,i,j,k)
    sumdershapexi = sumdershapexi + dershape3D(1,ia,i,j,k)
    sumdershapeeta = sumdershapeeta + dershape3D(2,ia,i,j,k)
    sumdershapegamma = sumdershapegamma + dershape3D(3,ia,i,j,k)
  enddo

! sum of shape functions should be one
! sum of derivative of shape functions should be zero
  if(abs(sumshape-one) >  TINYVAL) call exit_MPI(myrank,'error shape functions')
  if(abs(sumdershapexi) >  TINYVAL) call exit_MPI(myrank,'error derivative xi shape functions')
  if(abs(sumdershapeeta) >  TINYVAL) call exit_MPI(myrank,'error derivative eta shape functions')
  if(abs(sumdershapegamma) >  TINYVAL) call exit_MPI(myrank,'error derivative gamma shape functions')

  enddo
  enddo
  enddo

  end subroutine get_shape3D

