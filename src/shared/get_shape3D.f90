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

! 3D shape functions for 8-node element

  subroutine get_shape3D(myrank,shape3D,dershape3D,xigll,yigll,zigll,NGNOD)

  implicit none

  include "constants.h"

  integer NGNOD,myrank

! Gauss-Lobatto-Legendre points of integration
  double precision xigll(NGLLX)
  double precision yigll(NGLLY)
  double precision zigll(NGLLZ)

! 3D shape functions and their derivatives
  double precision shape3D(NGNOD,NGLLX,NGLLY,NGLLZ)
  double precision dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ)

  integer i,j,k,ia

! location of the nodes of the 3D hexahedra elements
  double precision xi,eta,gamma
  double precision ra1,ra2,rb1,rb2,rc1,rc2

! for checking the 3D shape functions
  double precision sumshape,sumdershapexi,sumdershapeeta,sumdershapegamma

  double precision, parameter :: ONE_EIGHTH = 0.125d0

! check that the parameter file is correct
  if(NGNOD /= 8 .and. NGNOD /= 27) &
       call exit_MPI(myrank,'volume elements should have 8 or 27 control nodes')

! ***
! *** create 3D shape functions and jacobian
! ***

  do i=1,NGLLX
  do j=1,NGLLY
  do k=1,NGLLZ

  xi = xigll(i)
  eta = yigll(j)
  gamma = zigll(k)

!--- case of a 3D 8-node element (Dhatt-Touzot p. 115)
  if(NGNOD == 8) then

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

    else

    ! note: put further initialization for NGNOD == 27 into subroutine
    !       to avoid compilation errors in case NGNOD == 8
      call get_shape3D_27(NGNOD,shape3D,dershape3D,xi,eta,gamma,i,j,k)

    endif

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
  if(abs(sumshape-one) >  TINYVAL) call exit_MPI(myrank,'error in 3D shape functions')
  if(abs(sumdershapexi) >  TINYVAL) call exit_MPI(myrank,'error in xi derivative of 3D shape functions')
  if(abs(sumdershapeeta) >  TINYVAL) call exit_MPI(myrank,'error in eta derivative of 3D shape functions')
  if(abs(sumdershapegamma) >  TINYVAL) call exit_MPI(myrank,'error in gamma derivative of 3D shape functions')

  enddo
  enddo
  enddo

  end subroutine get_shape3D

!
!-------------------------------------------------------------------------------------------------
!

! 3D shape functions for given, single xi/eta/gamma location

  subroutine eval_shape3D_single(myrank,shape3D,xi,eta,gamma,NGNOD)

  implicit none

  include "constants.h"

  integer :: myrank,NGNOD

  ! 3D shape functions
  double precision :: shape3D(NGNOD)

  ! location
  double precision :: xi,eta,gamma

  ! local parameters
  double precision, parameter :: ONE_EIGHTH = 0.125d0
  double precision :: ra1,ra2,rb1,rb2,rc1,rc2
  double precision :: l1xi,l2xi,l3xi,l1eta,l2eta,l3eta,l1gamma,l2gamma,l3gamma
  double precision :: sumshape
  integer :: ia

! check that the parameter file is correct
  if(NGNOD /= 8 .and. NGNOD /= 27) &
       call exit_MPI(myrank,'volume elements should have 8 or 27 control nodes')

  ! shape functions

!--- case of a 3D 8-node element (Dhatt-Touzot p. 115)
  if(NGNOD == 8) then

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

  else

  l1xi=HALF*xi*(xi-ONE)
  l2xi=ONE-xi**2
  l3xi=HALF*xi*(xi+ONE)

    l1eta=HALF*eta*(eta-ONE)
    l2eta=ONE-eta**2
    l3eta=HALF*eta*(eta+ONE)

      l1gamma=HALF*gamma*(gamma-ONE)
      l2gamma=ONE-gamma**2
      l3gamma=HALF*gamma*(gamma+ONE)

!     corner nodes
      shape3D(1)=l1xi*l1eta*l1gamma
      shape3D(2)=l3xi*l1eta*l1gamma
      shape3D(3)=l3xi*l3eta*l1gamma
      shape3D(4)=l1xi*l3eta*l1gamma
      shape3D(5)=l1xi*l1eta*l3gamma
      shape3D(6)=l3xi*l1eta*l3gamma
      shape3D(7)=l3xi*l3eta*l3gamma
      shape3D(8)=l1xi*l3eta*l3gamma

!     midside nodes
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

!     side center nodes
      shape3D(21)=l2xi*l2eta*l1gamma
      shape3D(22)=l2xi*l1eta*l2gamma
      shape3D(23)=l3xi*l2eta*l2gamma
      shape3D(24)=l2xi*l3eta*l2gamma
      shape3D(25)=l1xi*l2eta*l2gamma
      shape3D(26)=l2xi*l2eta*l3gamma

!     center node
      shape3D(27)=l2xi*l2eta*l2gamma

  endif

  ! check the shape functions
  sumshape = ZERO
  do ia=1,NGNOD
    sumshape = sumshape + shape3D(ia)
  enddo

  ! sum of shape functions should be one
  ! sum of derivative of shape functions should be zero
  if(abs(sumshape-one) >  TINYVAL) call exit_MPI(myrank,'error single shape functions')

  end subroutine eval_shape3D_single

!
!-------------------------------------------------------------------------------------------------
!

  subroutine eval_shape3D_element_corners(xelm,yelm,zelm,ispec,&
       ibool,xstore,ystore,zstore,NSPEC_AB,NGLOB_AB)

  implicit none

  include "constants.h"

  integer :: ispec
  integer :: NSPEC_AB,NGLOB_AB

  real(kind=CUSTOM_REAL),dimension(NGNOD_EIGHT_CORNERS),intent(out) :: xelm,yelm,zelm

  ! mesh coordinates
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB) :: xstore,ystore,zstore
  integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! 8 node corners
  xelm(1)=xstore(ibool(1,1,1,ispec))
  yelm(1)=ystore(ibool(1,1,1,ispec))
  zelm(1)=zstore(ibool(1,1,1,ispec))

  xelm(2)=xstore(ibool(NGLLX,1,1,ispec))
  yelm(2)=ystore(ibool(NGLLX,1,1,ispec))
  zelm(2)=zstore(ibool(NGLLX,1,1,ispec))

  xelm(3)=xstore(ibool(NGLLX,NGLLY,1,ispec))
  yelm(3)=ystore(ibool(NGLLX,NGLLY,1,ispec))
  zelm(3)=zstore(ibool(NGLLX,NGLLY,1,ispec))

  xelm(4)=xstore(ibool(1,NGLLY,1,ispec))
  yelm(4)=ystore(ibool(1,NGLLY,1,ispec))
  zelm(4)=zstore(ibool(1,NGLLY,1,ispec))

  xelm(5)=xstore(ibool(1,1,NGLLZ,ispec))
  yelm(5)=ystore(ibool(1,1,NGLLZ,ispec))
  zelm(5)=zstore(ibool(1,1,NGLLZ,ispec))

  xelm(6)=xstore(ibool(NGLLX,1,NGLLZ,ispec))
  yelm(6)=ystore(ibool(NGLLX,1,NGLLZ,ispec))
  zelm(6)=zstore(ibool(NGLLX,1,NGLLZ,ispec))

  xelm(7)=xstore(ibool(NGLLX,NGLLY,NGLLZ,ispec))
  yelm(7)=ystore(ibool(NGLLX,NGLLY,NGLLZ,ispec))
  zelm(7)=zstore(ibool(NGLLX,NGLLY,NGLLZ,ispec))

  xelm(8)=xstore(ibool(1,NGLLY,NGLLZ,ispec))
  yelm(8)=ystore(ibool(1,NGLLY,NGLLZ,ispec))
  zelm(8)=zstore(ibool(1,NGLLY,NGLLZ,ispec))

  end subroutine eval_shape3D_element_corners

!
!-------------------------------------------------------------------------------------------------
!

!--- case of a 3D 27-node element

  subroutine get_shape3D_27(NGNOD,shape3D,dershape3D,xi,eta,gamma,i,j,k)

  implicit none

  include "constants.h"

  integer :: NGNOD,i,j,k

! 3D shape functions and their derivatives
  double precision shape3D(NGNOD,NGLLX,NGLLY,NGLLZ)
  double precision dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ)

! location of the nodes of the 3D hexahedra elements
  double precision xi,eta,gamma
  double precision l1xi,l2xi,l3xi,l1eta,l2eta,l3eta,l1gamma,l2gamma,l3gamma
  double precision l1pxi,l2pxi,l3pxi,l1peta,l2peta,l3peta,l1pgamma,l2pgamma,l3pgamma

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

!     corner nodes
      shape3D(1,i,j,k)=l1xi*l1eta*l1gamma
      shape3D(2,i,j,k)=l3xi*l1eta*l1gamma
      shape3D(3,i,j,k)=l3xi*l3eta*l1gamma
      shape3D(4,i,j,k)=l1xi*l3eta*l1gamma
      shape3D(5,i,j,k)=l1xi*l1eta*l3gamma
      shape3D(6,i,j,k)=l3xi*l1eta*l3gamma
      shape3D(7,i,j,k)=l3xi*l3eta*l3gamma
      shape3D(8,i,j,k)=l1xi*l3eta*l3gamma

      dershape3D(1,1,i,j,k)=l1pxi*l1eta*l1gamma
      dershape3D(1,2,i,j,k)=l3pxi*l1eta*l1gamma
      dershape3D(1,3,i,j,k)=l3pxi*l3eta*l1gamma
      dershape3D(1,4,i,j,k)=l1pxi*l3eta*l1gamma
      dershape3D(1,5,i,j,k)=l1pxi*l1eta*l3gamma
      dershape3D(1,6,i,j,k)=l3pxi*l1eta*l3gamma
      dershape3D(1,7,i,j,k)=l3pxi*l3eta*l3gamma
      dershape3D(1,8,i,j,k)=l1pxi*l3eta*l3gamma

      dershape3D(2,1,i,j,k)=l1xi*l1peta*l1gamma
      dershape3D(2,2,i,j,k)=l3xi*l1peta*l1gamma
      dershape3D(2,3,i,j,k)=l3xi*l3peta*l1gamma
      dershape3D(2,4,i,j,k)=l1xi*l3peta*l1gamma
      dershape3D(2,5,i,j,k)=l1xi*l1peta*l3gamma
      dershape3D(2,6,i,j,k)=l3xi*l1peta*l3gamma
      dershape3D(2,7,i,j,k)=l3xi*l3peta*l3gamma
      dershape3D(2,8,i,j,k)=l1xi*l3peta*l3gamma

      dershape3D(3,1,i,j,k)=l1xi*l1eta*l1pgamma
      dershape3D(3,2,i,j,k)=l3xi*l1eta*l1pgamma
      dershape3D(3,3,i,j,k)=l3xi*l3eta*l1pgamma
      dershape3D(3,4,i,j,k)=l1xi*l3eta*l1pgamma
      dershape3D(3,5,i,j,k)=l1xi*l1eta*l3pgamma
      dershape3D(3,6,i,j,k)=l3xi*l1eta*l3pgamma
      dershape3D(3,7,i,j,k)=l3xi*l3eta*l3pgamma
      dershape3D(3,8,i,j,k)=l1xi*l3eta*l3pgamma

!     midside nodes
      shape3D(9,i,j,k)=l2xi*l1eta*l1gamma
      shape3D(10,i,j,k)=l3xi*l2eta*l1gamma
      shape3D(11,i,j,k)=l2xi*l3eta*l1gamma
      shape3D(12,i,j,k)=l1xi*l2eta*l1gamma
      shape3D(13,i,j,k)=l1xi*l1eta*l2gamma
      shape3D(14,i,j,k)=l3xi*l1eta*l2gamma
      shape3D(15,i,j,k)=l3xi*l3eta*l2gamma
      shape3D(16,i,j,k)=l1xi*l3eta*l2gamma
      shape3D(17,i,j,k)=l2xi*l1eta*l3gamma
      shape3D(18,i,j,k)=l3xi*l2eta*l3gamma
      shape3D(19,i,j,k)=l2xi*l3eta*l3gamma
      shape3D(20,i,j,k)=l1xi*l2eta*l3gamma

      dershape3D(1,9,i,j,k)=l2pxi*l1eta*l1gamma
      dershape3D(1,10,i,j,k)=l3pxi*l2eta*l1gamma
      dershape3D(1,11,i,j,k)=l2pxi*l3eta*l1gamma
      dershape3D(1,12,i,j,k)=l1pxi*l2eta*l1gamma
      dershape3D(1,13,i,j,k)=l1pxi*l1eta*l2gamma
      dershape3D(1,14,i,j,k)=l3pxi*l1eta*l2gamma
      dershape3D(1,15,i,j,k)=l3pxi*l3eta*l2gamma
      dershape3D(1,16,i,j,k)=l1pxi*l3eta*l2gamma
      dershape3D(1,17,i,j,k)=l2pxi*l1eta*l3gamma
      dershape3D(1,18,i,j,k)=l3pxi*l2eta*l3gamma
      dershape3D(1,19,i,j,k)=l2pxi*l3eta*l3gamma
      dershape3D(1,20,i,j,k)=l1pxi*l2eta*l3gamma

      dershape3D(2,9,i,j,k)=l2xi*l1peta*l1gamma
      dershape3D(2,10,i,j,k)=l3xi*l2peta*l1gamma
      dershape3D(2,11,i,j,k)=l2xi*l3peta*l1gamma
      dershape3D(2,12,i,j,k)=l1xi*l2peta*l1gamma
      dershape3D(2,13,i,j,k)=l1xi*l1peta*l2gamma
      dershape3D(2,14,i,j,k)=l3xi*l1peta*l2gamma
      dershape3D(2,15,i,j,k)=l3xi*l3peta*l2gamma
      dershape3D(2,16,i,j,k)=l1xi*l3peta*l2gamma
      dershape3D(2,17,i,j,k)=l2xi*l1peta*l3gamma
      dershape3D(2,18,i,j,k)=l3xi*l2peta*l3gamma
      dershape3D(2,19,i,j,k)=l2xi*l3peta*l3gamma
      dershape3D(2,20,i,j,k)=l1xi*l2peta*l3gamma

      dershape3D(3,9,i,j,k)=l2xi*l1eta*l1pgamma
      dershape3D(3,10,i,j,k)=l3xi*l2eta*l1pgamma
      dershape3D(3,11,i,j,k)=l2xi*l3eta*l1pgamma
      dershape3D(3,12,i,j,k)=l1xi*l2eta*l1pgamma
      dershape3D(3,13,i,j,k)=l1xi*l1eta*l2pgamma
      dershape3D(3,14,i,j,k)=l3xi*l1eta*l2pgamma
      dershape3D(3,15,i,j,k)=l3xi*l3eta*l2pgamma
      dershape3D(3,16,i,j,k)=l1xi*l3eta*l2pgamma
      dershape3D(3,17,i,j,k)=l2xi*l1eta*l3pgamma
      dershape3D(3,18,i,j,k)=l3xi*l2eta*l3pgamma
      dershape3D(3,19,i,j,k)=l2xi*l3eta*l3pgamma
      dershape3D(3,20,i,j,k)=l1xi*l2eta*l3pgamma

!     side center nodes
      shape3D(21,i,j,k)=l2xi*l2eta*l1gamma
      shape3D(22,i,j,k)=l2xi*l1eta*l2gamma
      shape3D(23,i,j,k)=l3xi*l2eta*l2gamma
      shape3D(24,i,j,k)=l2xi*l3eta*l2gamma
      shape3D(25,i,j,k)=l1xi*l2eta*l2gamma
      shape3D(26,i,j,k)=l2xi*l2eta*l3gamma

      dershape3D(1,21,i,j,k)=l2pxi*l2eta*l1gamma
      dershape3D(1,22,i,j,k)=l2pxi*l1eta*l2gamma
      dershape3D(1,23,i,j,k)=l3pxi*l2eta*l2gamma
      dershape3D(1,24,i,j,k)=l2pxi*l3eta*l2gamma
      dershape3D(1,25,i,j,k)=l1pxi*l2eta*l2gamma
      dershape3D(1,26,i,j,k)=l2pxi*l2eta*l3gamma

      dershape3D(2,21,i,j,k)=l2xi*l2peta*l1gamma
      dershape3D(2,22,i,j,k)=l2xi*l1peta*l2gamma
      dershape3D(2,23,i,j,k)=l3xi*l2peta*l2gamma
      dershape3D(2,24,i,j,k)=l2xi*l3peta*l2gamma
      dershape3D(2,25,i,j,k)=l1xi*l2peta*l2gamma
      dershape3D(2,26,i,j,k)=l2xi*l2peta*l3gamma

      dershape3D(3,21,i,j,k)=l2xi*l2eta*l1pgamma
      dershape3D(3,22,i,j,k)=l2xi*l1eta*l2pgamma
      dershape3D(3,23,i,j,k)=l3xi*l2eta*l2pgamma
      dershape3D(3,24,i,j,k)=l2xi*l3eta*l2pgamma
      dershape3D(3,25,i,j,k)=l1xi*l2eta*l2pgamma
      dershape3D(3,26,i,j,k)=l2xi*l2eta*l3pgamma

!     center node
      shape3D(27,i,j,k)=l2xi*l2eta*l2gamma

      dershape3D(1,27,i,j,k)=l2pxi*l2eta*l2gamma
      dershape3D(2,27,i,j,k)=l2xi*l2peta*l2gamma
      dershape3D(3,27,i,j,k)=l2xi*l2eta*l2pgamma

  end subroutine get_shape3D_27

