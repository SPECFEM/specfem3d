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

!
!-------------------------------------------------------------------------------------------------
!

! 3D shape functions for given, single xi/eta/gamma location

  subroutine get_shape3D_single(myrank,shape3D,xi,eta,gamma)

  implicit none

  include "constants.h"

  integer :: myrank

  ! 3D shape functions 
  double precision :: shape3D(NGNOD)

  ! location 
  double precision :: xi,eta,gamma
  
  ! local parameters
  double precision :: ra1,ra2,rb1,rb2,rc1,rc2
  double precision, parameter :: ONE_EIGHTH = 0.125d0
  double precision :: sumshape
  integer :: ia

! check that the parameter file is correct
  if(NGNOD /= 8) call exit_MPI(myrank,'elements should have 8 control nodes')

!--- case of a 3D 8-node element (Dhatt-Touzot p. 115)
  ra1 = one + xi
  ra2 = one - xi

  rb1 = one + eta
  rb2 = one - eta

  rc1 = one + gamma
  rc2 = one - gamma

  ! shape functions
  shape3D(1) = ONE_EIGHTH*ra2*rb2*rc2
  shape3D(2) = ONE_EIGHTH*ra1*rb2*rc2
  shape3D(3) = ONE_EIGHTH*ra1*rb1*rc2
  shape3D(4) = ONE_EIGHTH*ra2*rb1*rc2
  shape3D(5) = ONE_EIGHTH*ra2*rb2*rc1
  shape3D(6) = ONE_EIGHTH*ra1*rb2*rc1
  shape3D(7) = ONE_EIGHTH*ra1*rb1*rc1
  shape3D(8) = ONE_EIGHTH*ra2*rb1*rc1

  ! check the shape functions
  sumshape = ZERO
  do ia=1,NGNOD
    sumshape = sumshape + shape3D(ia)
  enddo

  ! sum of shape functions should be one
  ! sum of derivative of shape functions should be zero
  if(abs(sumshape-one) >  TINYVAL) call exit_MPI(myrank,'error single shape functions')

  end subroutine get_shape3D_single

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_shape3D_element_corners(xelm,yelm,zelm,ispec,&
                        ibool,xstore,ystore,zstore,NSPEC_AB,NGLOB_AB)

  implicit none

  include "constants.h"

  integer :: ispec
  integer :: NSPEC_AB,NGLOB_AB

  real(kind=CUSTOM_REAL),dimension(NGNOD),intent(out) :: xelm,yelm,zelm
  
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

  end subroutine get_shape3D_element_corners
  
  