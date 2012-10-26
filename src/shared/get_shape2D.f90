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

  subroutine get_shape2D(myrank,shape2D,dershape2D,xigll,yigll,NGLLA,NGLLB,NGNOD,NGNOD2D)

  implicit none

  include "constants.h"

! generic routine that accepts any polynomial degree in each direction

  integer NGLLA,NGLLB,NGNOD,NGNOD2D,myrank

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
  if(NGNOD /= 8 .and. NGNOD /= 27) &
       call exit_MPI(myrank,'volume elements should have 8 or 27 control nodes')
  if(NGNOD2D /= 4 .and. NGNOD2D /= 9) &
       call exit_MPI(myrank,'surface elements should have 4 or 9 control nodes')

  if(NGNOD2D == 4) then

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

  else

    ! note: put further initialization for ngnod2d == 9 into subroutine
    !       to avoid compilation errors in case ngnod2d == 4
    call get_shape2D_9(NGNOD2D,shape2D,dershape2D,xigll,yigll,NGLLA,NGLLB)

  endif

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

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_shape2D_9(NGNOD2D,shape2D,dershape2D,xigll,yigll,NGLLA,NGLLB)

  implicit none

  include "constants.h"

! generic routine that accepts any polynomial degree in each direction

  integer :: NGNOD2D
  integer :: NGLLA,NGLLB

  double precision xigll(NGLLA)
  double precision yigll(NGLLB)

! 2D shape functions and their derivatives
  double precision shape2D(NGNOD2D,NGLLA,NGLLB)
  double precision dershape2D(NDIM2D,NGNOD2D,NGLLA,NGLLB)

  integer i,j

! location of the nodes of the 2D quadrilateral elements
  double precision xi,eta
  double precision l1xi,l2xi,l3xi,l1eta,l2eta,l3eta
  double precision l1pxi,l2pxi,l3pxi,l1peta,l2peta,l3peta

  ! check that the parameter file is correct
  if( NGNOD2D /= 9 ) stop 'surface elements should have 9 control nodes'

  ! generate the 2D shape functions and their derivatives (9 nodes)
  do i=1,NGLLA

    xi=xigll(i)

    l1xi=HALF*xi*(xi-ONE)
    l2xi=ONE-xi**2
    l3xi=HALF*xi*(xi+ONE)

    l1pxi=xi-HALF
    l2pxi=-TWO*xi
    l3pxi=xi+HALF

    do j=1,NGLLB

      eta=yigll(j)

      l1eta=HALF*eta*(eta-ONE)
      l2eta=ONE-eta**2
      l3eta=HALF*eta*(eta+ONE)

      l1peta=eta-HALF
      l2peta=-TWO*eta
      l3peta=eta+HALF

!   corner nodes

      shape2D(1,i,j)=l1xi*l1eta
      shape2D(2,i,j)=l3xi*l1eta
      shape2D(3,i,j)=l3xi*l3eta
      shape2D(4,i,j)=l1xi*l3eta

      dershape2D(1,1,i,j)=l1pxi*l1eta
      dershape2D(1,2,i,j)=l3pxi*l1eta
      dershape2D(1,3,i,j)=l3pxi*l3eta
      dershape2D(1,4,i,j)=l1pxi*l3eta

      dershape2D(2,1,i,j)=l1xi*l1peta
      dershape2D(2,2,i,j)=l3xi*l1peta
      dershape2D(2,3,i,j)=l3xi*l3peta
      dershape2D(2,4,i,j)=l1xi*l3peta

!   midside nodes

      shape2D(5,i,j)=l2xi*l1eta
      shape2D(6,i,j)=l3xi*l2eta
      shape2D(7,i,j)=l2xi*l3eta
      shape2D(8,i,j)=l1xi*l2eta

      dershape2D(1,5,i,j)=l2pxi*l1eta
      dershape2D(1,6,i,j)=l3pxi*l2eta
      dershape2D(1,7,i,j)=l2pxi*l3eta
      dershape2D(1,8,i,j)=l1pxi*l2eta

      dershape2D(2,5,i,j)=l2xi*l1peta
      dershape2D(2,6,i,j)=l3xi*l2peta
      dershape2D(2,7,i,j)=l2xi*l3peta
      dershape2D(2,8,i,j)=l1xi*l2peta

!   center node

      shape2D(9,i,j)=l2xi*l2eta

      dershape2D(1,9,i,j)=l2pxi*l2eta
      dershape2D(2,9,i,j)=l2xi*l2peta

    enddo
  enddo

  end subroutine get_shape2D_9

