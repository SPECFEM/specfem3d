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

  subroutine get_shape2D(myrank,shape2D,dershape2D,xigll,yigll,NGLLA,NGLLB)

  implicit none

  include "constants.h"

! generic routine that accepts any polynomial degree in each direction

  integer NGLLA,NGLLB,myrank

  double precision xigll(NGLLA)
  double precision yigll(NGLLB)

! 2D shape functions and their derivatives
  double precision shape2D(NGNOD2D,NGLLA,NGLLB)
  double precision dershape2D(NDIM2D,NGNOD2D,NGLLA,NGLLB)

  integer i,j,ia

! location of the nodes of the 2D quadrilateral elements
  double precision xi,eta
  double precision xi_map,eta_map

! for checking the 2D shape functions
  double precision sumshape,sumdershapexi,sumdershapeeta

! check that the parameter file is correct
  if(NGNOD /= 8) call exit_MPI(myrank,'elements should have 8 control nodes')
  if(NGNOD2D /= 4) call exit_MPI(myrank,'surface elements should have 4 control nodes')

! generate the 2D shape functions and their derivatives (4 nodes)
  do i=1,NGLLA

  xi=xigll(i)

  do j=1,NGLLB

    eta=yigll(j)

! map coordinates to [0,1]
    xi_map = (xi + 1.) / 2.
    eta_map = (eta + 1.) / 2.

! corner nodes
    shape2D(1,i,j) = (1 - xi_map)*(1 - eta_map)
    shape2D(2,i,j) = xi_map*(1 - eta_map)
    shape2D(3,i,j) = xi_map*eta_map
    shape2D(4,i,j) = (1 - xi_map)*eta_map

    dershape2D(1,1,i,j) = (eta - 1.) / 4.
    dershape2D(2,1,i,j) = (xi - 1.) / 4.

    dershape2D(1,2,i,j) = (1. - eta) / 4.
    dershape2D(2,2,i,j) = (-1. - xi) / 4.

    dershape2D(1,3,i,j) = (1. + eta) / 4.
    dershape2D(2,3,i,j) = (1. + xi) / 4.

    dershape2D(1,4,i,j) = (- 1. - eta) / 4.
    dershape2D(2,4,i,j) = (1. - xi) / 4.

    enddo
  enddo

! check the 2D shape functions
  do i=1,NGLLA
    do j=1,NGLLB

    sumshape=ZERO

    sumdershapexi=ZERO
    sumdershapeeta=ZERO

    do ia=1,NGNOD2D
      sumshape=sumshape+shape2D(ia,i,j)

      sumdershapexi=sumdershapexi+dershape2D(1,ia,i,j)
      sumdershapeeta=sumdershapeeta+dershape2D(2,ia,i,j)
    enddo

!   the sum of the shape functions should be 1
    if(abs(sumshape-ONE)>TINYVAL) call exit_MPI(myrank,'error in 2D shape functions')

!   the sum of the derivatives of the shape functions should be 0
    if(abs(sumdershapexi)>TINYVAL) &
      call exit_MPI(myrank,'error in xi derivatives of 2D shape function')

    if(abs(sumdershapeeta)>TINYVAL) &
      call exit_MPI(myrank,'error in eta derivatives of 2D shape function')

    enddo
  enddo

  end subroutine get_shape2D

