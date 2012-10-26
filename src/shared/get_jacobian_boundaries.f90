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

  subroutine get_jacobian_boundary_face(myrank,nspec, &
                        xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob,&
                        dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,&
                        ispec,iface,jacobian2Dw_face,normal_face,NGLLA,NGLLB,NGNOD2D)

! returns jacobian2Dw_face and normal_face (pointing outwards of element)

  implicit none

  include "constants.h"

  integer :: nspec,myrank,nglob

! arrays with the mesh
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore_dummy,ystore_dummy,zstore_dummy

! face information
  integer :: iface,ispec,NGLLA,NGLLB
  real(kind=CUSTOM_REAL), dimension(NGLLA,NGLLB) :: jacobian2Dw_face
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLA,NGLLB) :: normal_face

  integer :: NGNOD2D
  double precision dershape2D_x(NDIM2D,NGNOD2D,NGLLY,NGLLZ)
  double precision dershape2D_y(NDIM2D,NGNOD2D,NGLLX,NGLLZ)
  double precision dershape2D_bottom(NDIM2D,NGNOD2D,NGLLX,NGLLY)
  double precision dershape2D_top(NDIM2D,NGNOD2D,NGLLX,NGLLY)

  double precision, dimension(NGLLX,NGLLY) :: wgllwgll_xy
  double precision, dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  double precision, dimension(NGLLY,NGLLZ) :: wgllwgll_yz

! local parameters
! face corners
  double precision xelm(NGNOD2D),yelm(NGNOD2D),zelm(NGNOD2D)

! check that the parameter file is correct
  if(NGNOD2D /= 4 .and. NGNOD2D /= 9) &
       call exit_MPI(myrank,'surface elements should have 4 or 9 control nodes')

  select case ( iface )
  ! on reference face: xmin
  case(1)
    xelm(1)=xstore_dummy( ibool(1,1,1,ispec) )
    yelm(1)=ystore_dummy( ibool(1,1,1,ispec) )
    zelm(1)=zstore_dummy( ibool(1,1,1,ispec) )
    xelm(2)=xstore_dummy( ibool(1,NGLLY,1,ispec) )
    yelm(2)=ystore_dummy( ibool(1,NGLLY,1,ispec) )
    zelm(2)=zstore_dummy( ibool(1,NGLLY,1,ispec) )
    xelm(3)=xstore_dummy( ibool(1,NGLLY,NGLLZ,ispec) )
    yelm(3)=ystore_dummy( ibool(1,NGLLY,NGLLZ,ispec) )
    zelm(3)=zstore_dummy( ibool(1,NGLLY,NGLLZ,ispec) )
    xelm(4)=xstore_dummy( ibool(1,1,NGLLZ,ispec) )
    yelm(4)=ystore_dummy( ibool(1,1,NGLLZ,ispec) )
    zelm(4)=zstore_dummy( ibool(1,1,NGLLZ,ispec) )

    if(NGNOD2D == 9) then
          xelm(5)=xstore_dummy( ibool(1,(NGLLY+1)/2,1,ispec) )
          yelm(5)=ystore_dummy( ibool(1,(NGLLY+1)/2,1,ispec) )
          zelm(5)=zstore_dummy( ibool(1,(NGLLY+1)/2,1,ispec) )
          xelm(6)=xstore_dummy( ibool(1,NGLLY,(NGLLZ+1)/2,ispec) )
          yelm(6)=ystore_dummy( ibool(1,NGLLY,(NGLLZ+1)/2,ispec) )
          zelm(6)=zstore_dummy( ibool(1,NGLLY,(NGLLZ+1)/2,ispec) )
          xelm(7)=xstore_dummy( ibool(1,(NGLLY+1)/2,NGLLZ,ispec) )
          yelm(7)=ystore_dummy( ibool(1,(NGLLY+1)/2,NGLLZ,ispec) )
          zelm(7)=zstore_dummy( ibool(1,(NGLLY+1)/2,NGLLZ,ispec) )
          xelm(8)=xstore_dummy( ibool(1,1,(NGLLZ+1)/2,ispec) )
          yelm(8)=ystore_dummy( ibool(1,1,(NGLLZ+1)/2,ispec) )
          zelm(8)=zstore_dummy( ibool(1,1,(NGLLZ+1)/2,ispec) )
          xelm(9)=xstore_dummy( ibool(1,(NGLLY+1)/2,(NGLLZ+1)/2,ispec) )
          yelm(9)=ystore_dummy( ibool(1,(NGLLY+1)/2,(NGLLZ+1)/2,ispec) )
          zelm(9)=zstore_dummy( ibool(1,(NGLLY+1)/2,(NGLLZ+1)/2,ispec) )
    endif

    call compute_jacobian_2D_face(myrank,xelm,yelm,zelm, &
                  dershape2D_x,wgllwgll_yz, &
                  jacobian2Dw_face,normal_face,NGLLY,NGLLZ,NGNOD2D)

! on boundary: xmax
  case(2)
    xelm(1)=xstore_dummy( ibool(NGLLX,1,1,ispec) )
    yelm(1)=ystore_dummy( ibool(NGLLX,1,1,ispec) )
    zelm(1)=zstore_dummy( ibool(NGLLX,1,1,ispec) )
    xelm(2)=xstore_dummy( ibool(NGLLX,NGLLY,1,ispec) )
    yelm(2)=ystore_dummy( ibool(NGLLX,NGLLY,1,ispec) )
    zelm(2)=zstore_dummy( ibool(NGLLX,NGLLY,1,ispec) )
    xelm(3)=xstore_dummy( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    yelm(3)=ystore_dummy( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    zelm(3)=zstore_dummy( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    xelm(4)=xstore_dummy( ibool(NGLLX,1,NGLLZ,ispec) )
    yelm(4)=ystore_dummy( ibool(NGLLX,1,NGLLZ,ispec) )
    zelm(4)=zstore_dummy( ibool(NGLLX,1,NGLLZ,ispec) )

    if(NGNOD2D == 9) then
          xelm(5)=xstore_dummy( ibool(NGLLX,(NGLLY+1)/2,1,ispec) )
          yelm(5)=ystore_dummy( ibool(NGLLX,(NGLLY+1)/2,1,ispec) )
          zelm(5)=zstore_dummy( ibool(NGLLX,(NGLLY+1)/2,1,ispec) )
          xelm(6)=xstore_dummy( ibool(NGLLX,NGLLY,(NGLLZ+1)/2,ispec) )
          yelm(6)=ystore_dummy( ibool(NGLLX,NGLLY,(NGLLZ+1)/2,ispec) )
          zelm(6)=zstore_dummy( ibool(NGLLX,NGLLY,(NGLLZ+1)/2,ispec) )
          xelm(7)=xstore_dummy( ibool(NGLLX,(NGLLY+1)/2,NGLLZ,ispec) )
          yelm(7)=ystore_dummy( ibool(NGLLX,(NGLLY+1)/2,NGLLZ,ispec) )
          zelm(7)=zstore_dummy( ibool(NGLLX,(NGLLY+1)/2,NGLLZ,ispec) )
          xelm(8)=xstore_dummy( ibool(NGLLX,1,(NGLLZ+1)/2,ispec) )
          yelm(8)=ystore_dummy( ibool(NGLLX,1,(NGLLZ+1)/2,ispec) )
          zelm(8)=zstore_dummy( ibool(NGLLX,1,(NGLLZ+1)/2,ispec) )
          xelm(9)=xstore_dummy( ibool(NGLLX,(NGLLY+1)/2,(NGLLZ+1)/2,ispec) )
          yelm(9)=ystore_dummy( ibool(NGLLX,(NGLLY+1)/2,(NGLLZ+1)/2,ispec) )
          zelm(9)=zstore_dummy( ibool(NGLLX,(NGLLY+1)/2,(NGLLZ+1)/2,ispec) )
    endif

    call compute_jacobian_2D_face(myrank,xelm,yelm,zelm, &
                  dershape2D_x,wgllwgll_yz, &
                  jacobian2Dw_face,normal_face,NGLLY,NGLLZ,NGNOD2D)

! on boundary: ymin
  case(3)
    xelm(1)=xstore_dummy( ibool(1,1,1,ispec) )
    yelm(1)=ystore_dummy( ibool(1,1,1,ispec) )
    zelm(1)=zstore_dummy( ibool(1,1,1,ispec) )
    xelm(2)=xstore_dummy( ibool(NGLLX,1,1,ispec) )
    yelm(2)=ystore_dummy( ibool(NGLLX,1,1,ispec) )
    zelm(2)=zstore_dummy( ibool(NGLLX,1,1,ispec) )
    xelm(3)=xstore_dummy( ibool(NGLLX,1,NGLLZ,ispec) )
    yelm(3)=ystore_dummy( ibool(NGLLX,1,NGLLZ,ispec) )
    zelm(3)=zstore_dummy( ibool(NGLLX,1,NGLLZ,ispec) )
    xelm(4)=xstore_dummy( ibool(1,1,NGLLZ,ispec) )
    yelm(4)=ystore_dummy( ibool(1,1,NGLLZ,ispec) )
    zelm(4)=zstore_dummy( ibool(1,1,NGLLZ,ispec) )

    if(NGNOD2D == 9) then
          xelm(5)=xstore_dummy( ibool((NGLLX+1)/2,1,1,ispec) )
          yelm(5)=ystore_dummy( ibool((NGLLX+1)/2,1,1,ispec) )
          zelm(5)=zstore_dummy( ibool((NGLLX+1)/2,1,1,ispec) )
          xelm(6)=xstore_dummy( ibool(NGLLX,1,(NGLLZ+1)/2,ispec) )
          yelm(6)=ystore_dummy( ibool(NGLLX,1,(NGLLZ+1)/2,ispec) )
          zelm(6)=zstore_dummy( ibool(NGLLX,1,(NGLLZ+1)/2,ispec) )
          xelm(7)=xstore_dummy( ibool((NGLLX+1)/2,1,NGLLZ,ispec) )
          yelm(7)=ystore_dummy( ibool((NGLLX+1)/2,1,NGLLZ,ispec) )
          zelm(7)=zstore_dummy( ibool((NGLLX+1)/2,1,NGLLZ,ispec) )
          xelm(8)=xstore_dummy( ibool(1,1,(NGLLZ+1)/2,ispec) )
          yelm(8)=ystore_dummy( ibool(1,1,(NGLLZ+1)/2,ispec) )
          zelm(8)=zstore_dummy( ibool(1,1,(NGLLZ+1)/2,ispec) )
          xelm(9)=xstore_dummy( ibool((NGLLX+1)/2,1,(NGLLZ+1)/2,ispec) )
          yelm(9)=ystore_dummy( ibool((NGLLX+1)/2,1,(NGLLZ+1)/2,ispec) )
          zelm(9)=zstore_dummy( ibool((NGLLX+1)/2,1,(NGLLZ+1)/2,ispec) )
    endif

    call compute_jacobian_2D_face(myrank,xelm,yelm,zelm, &
                  dershape2D_y,wgllwgll_xz, &
                  jacobian2Dw_face,normal_face,NGLLX,NGLLZ,NGNOD2D)

! on boundary: ymax
  case(4)
    xelm(1)=xstore_dummy( ibool(1,NGLLY,1,ispec) )
    yelm(1)=ystore_dummy( ibool(1,NGLLY,1,ispec) )
    zelm(1)=zstore_dummy( ibool(1,NGLLY,1,ispec) )
    xelm(2)=xstore_dummy( ibool(NGLLX,NGLLY,1,ispec) )
    yelm(2)=ystore_dummy( ibool(NGLLX,NGLLY,1,ispec) )
    zelm(2)=zstore_dummy( ibool(NGLLX,NGLLY,1,ispec) )
    xelm(3)=xstore_dummy( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    yelm(3)=ystore_dummy( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    zelm(3)=zstore_dummy( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    xelm(4)=xstore_dummy( ibool(1,NGLLY,NGLLZ,ispec) )
    yelm(4)=ystore_dummy( ibool(1,NGLLY,NGLLZ,ispec) )
    zelm(4)=zstore_dummy( ibool(1,NGLLY,NGLLZ,ispec) )

    if(NGNOD2D == 9) then
          xelm(5)=xstore_dummy( ibool((NGLLX+1)/2,NGLLY,1,ispec) )
          yelm(5)=ystore_dummy( ibool((NGLLX+1)/2,NGLLY,1,ispec) )
          zelm(5)=zstore_dummy( ibool((NGLLX+1)/2,NGLLY,1,ispec) )
          xelm(6)=xstore_dummy( ibool(NGLLX,NGLLY,(NGLLZ+1)/2,ispec) )
          yelm(6)=ystore_dummy( ibool(NGLLX,NGLLY,(NGLLZ+1)/2,ispec) )
          zelm(6)=zstore_dummy( ibool(NGLLX,NGLLY,(NGLLZ+1)/2,ispec) )
          xelm(7)=xstore_dummy( ibool((NGLLX+1)/2,NGLLY,NGLLZ,ispec) )
          yelm(7)=ystore_dummy( ibool((NGLLX+1)/2,NGLLY,NGLLZ,ispec) )
          zelm(7)=zstore_dummy( ibool((NGLLX+1)/2,NGLLY,NGLLZ,ispec) )
          xelm(8)=xstore_dummy( ibool(1,NGLLY,(NGLLZ+1)/2,ispec) )
          yelm(8)=ystore_dummy( ibool(1,NGLLY,(NGLLZ+1)/2,ispec) )
          zelm(8)=zstore_dummy( ibool(1,NGLLY,(NGLLZ+1)/2,ispec) )
          xelm(9)=xstore_dummy( ibool((NGLLX+1)/2,NGLLY,(NGLLZ+1)/2,ispec) )
          yelm(9)=ystore_dummy( ibool((NGLLX+1)/2,NGLLY,(NGLLZ+1)/2,ispec) )
          zelm(9)=zstore_dummy( ibool((NGLLX+1)/2,NGLLY,(NGLLZ+1)/2,ispec) )
    endif

    call compute_jacobian_2D_face(myrank,xelm,yelm,zelm, &
                  dershape2D_y, wgllwgll_xz, &
                  jacobian2Dw_face,normal_face,NGLLX,NGLLZ,NGNOD2D)


! on boundary: bottom
  case(5)
    xelm(1)=xstore_dummy( ibool(1,1,1,ispec) )
    yelm(1)=ystore_dummy( ibool(1,1,1,ispec) )
    zelm(1)=zstore_dummy( ibool(1,1,1,ispec) )
    xelm(2)=xstore_dummy( ibool(NGLLX,1,1,ispec) )
    yelm(2)=ystore_dummy( ibool(NGLLX,1,1,ispec) )
    zelm(2)=zstore_dummy( ibool(NGLLX,1,1,ispec) )
    xelm(3)=xstore_dummy( ibool(NGLLX,NGLLY,1,ispec) )
    yelm(3)=ystore_dummy( ibool(NGLLX,NGLLY,1,ispec) )
    zelm(3)=zstore_dummy( ibool(NGLLX,NGLLY,1,ispec) )
    xelm(4)=xstore_dummy( ibool(1,NGLLY,1,ispec) )
    yelm(4)=ystore_dummy( ibool(1,NGLLY,1,ispec) )
    zelm(4)=zstore_dummy( ibool(1,NGLLY,1,ispec) )

    if(NGNOD2D == 9) then
          xelm(5)=xstore_dummy( ibool((NGLLX+1)/2,1,1,ispec) )
          yelm(5)=ystore_dummy( ibool((NGLLX+1)/2,1,1,ispec) )
          zelm(5)=zstore_dummy( ibool((NGLLX+1)/2,1,1,ispec) )
          xelm(6)=xstore_dummy( ibool(NGLLX,(NGLLY+1)/2,1,ispec) )
          yelm(6)=ystore_dummy( ibool(NGLLX,(NGLLY+1)/2,1,ispec) )
          zelm(6)=zstore_dummy( ibool(NGLLX,(NGLLY+1)/2,1,ispec) )
          xelm(7)=xstore_dummy( ibool((NGLLX+1)/2,NGLLY,1,ispec) )
          yelm(7)=ystore_dummy( ibool((NGLLX+1)/2,NGLLY,1,ispec) )
          zelm(7)=zstore_dummy( ibool((NGLLX+1)/2,NGLLY,1,ispec) )
          xelm(8)=xstore_dummy( ibool(1,(NGLLY+1)/2,1,ispec) )
          yelm(8)=ystore_dummy( ibool(1,(NGLLY+1)/2,1,ispec) )
          zelm(8)=zstore_dummy( ibool(1,(NGLLY+1)/2,1,ispec) )
          xelm(9)=xstore_dummy( ibool((NGLLX+1)/2,(NGLLY+1)/2,1,ispec) )
          yelm(9)=ystore_dummy( ibool((NGLLX+1)/2,(NGLLY+1)/2,1,ispec) )
          zelm(9)=zstore_dummy( ibool((NGLLX+1)/2,(NGLLY+1)/2,1,ispec) )
    endif

    call compute_jacobian_2D_face(myrank,xelm,yelm,zelm,&
                  dershape2D_bottom,wgllwgll_xy, &
                  jacobian2Dw_face,normal_face,NGLLX,NGLLY,NGNOD2D)

! on boundary: top
  case(6)
    xelm(1)=xstore_dummy( ibool(1,1,NGLLZ,ispec) )
    yelm(1)=ystore_dummy( ibool(1,1,NGLLZ,ispec) )
    zelm(1)=zstore_dummy( ibool(1,1,NGLLZ,ispec) )
    xelm(2)=xstore_dummy( ibool(NGLLX,1,NGLLZ,ispec) )
    yelm(2)=ystore_dummy( ibool(NGLLX,1,NGLLZ,ispec) )
    zelm(2)=zstore_dummy( ibool(NGLLX,1,NGLLZ,ispec) )
    xelm(3)=xstore_dummy( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    yelm(3)=ystore_dummy( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    zelm(3)=zstore_dummy( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    xelm(4)=xstore_dummy( ibool(1,NGLLY,NGLLZ,ispec) )
    yelm(4)=ystore_dummy( ibool(1,NGLLY,NGLLZ,ispec) )
    zelm(4)=zstore_dummy( ibool(1,NGLLY,NGLLZ,ispec) )

    if(NGNOD2D == 9) then
          xelm(5)=xstore_dummy( ibool((NGLLX+1)/2,1,NGLLZ,ispec) )
          yelm(5)=ystore_dummy( ibool((NGLLX+1)/2,1,NGLLZ,ispec) )
          zelm(5)=zstore_dummy( ibool((NGLLX+1)/2,1,NGLLZ,ispec) )
          xelm(6)=xstore_dummy( ibool(NGLLX,(NGLLY+1)/2,NGLLZ,ispec) )
          yelm(6)=ystore_dummy( ibool(NGLLX,(NGLLY+1)/2,NGLLZ,ispec) )
          zelm(6)=zstore_dummy( ibool(NGLLX,(NGLLY+1)/2,NGLLZ,ispec) )
          xelm(7)=xstore_dummy( ibool((NGLLX+1)/2,NGLLY,NGLLZ,ispec) )
          yelm(7)=ystore_dummy( ibool((NGLLX+1)/2,NGLLY,NGLLZ,ispec) )
          zelm(7)=zstore_dummy( ibool((NGLLX+1)/2,NGLLY,NGLLZ,ispec) )
          xelm(8)=xstore_dummy( ibool(1,(NGLLY+1)/2,NGLLZ,ispec) )
          yelm(8)=ystore_dummy( ibool(1,(NGLLY+1)/2,NGLLZ,ispec) )
          zelm(8)=zstore_dummy( ibool(1,(NGLLY+1)/2,NGLLZ,ispec) )
          xelm(9)=xstore_dummy( ibool((NGLLX+1)/2,(NGLLY+1)/2,NGLLZ,ispec) )
          yelm(9)=ystore_dummy( ibool((NGLLX+1)/2,(NGLLY+1)/2,NGLLZ,ispec) )
          zelm(9)=zstore_dummy( ibool((NGLLX+1)/2,(NGLLY+1)/2,NGLLZ,ispec) )
    endif

    call compute_jacobian_2D_face(myrank,xelm,yelm,zelm,&
                  dershape2D_top, wgllwgll_xy, &
                  jacobian2Dw_face,normal_face,NGLLX,NGLLY,NGNOD2D)

  case default
    stop 'error 2D jacobian'
  end select

  end subroutine get_jacobian_boundary_face

!
! ------------------------------------------------------------------------------------------------
!

  subroutine compute_jacobian_2D_face(myrank,xelm,yelm,zelm, &
                                dershape2D,wgllwgll, &
                                jacobian2Dw_face,normal_face,NGLLA,NGLLB,NGNOD2D)

  implicit none

  include "constants.h"

! generic routine that accepts any polynomial degree in each direction
! returns 2D jacobian and normal for this face only

  integer NGLLA,NGLLB,NGNOD2D,myrank

  double precision xelm(NGNOD2D),yelm(NGNOD2D),zelm(NGNOD2D)
  double precision dershape2D(NDIM2D,NGNOD2D,NGLLA,NGLLB)
  double precision, dimension(NGLLA,NGLLB) :: wgllwgll

  real(kind=CUSTOM_REAL), dimension(NGLLA,NGLLB) :: jacobian2Dw_face
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLA,NGLLB) :: normal_face

  integer i,j,ia
  double precision xxi,xeta,yxi,yeta,zxi,zeta
  double precision unx,uny,unz,jacobian

  do j=1,NGLLB
    do i=1,NGLLA

      xxi=ZERO
      xeta=ZERO
      yxi=ZERO
      yeta=ZERO
      zxi=ZERO
      zeta=ZERO
      do ia=1,NGNOD2D
        xxi=xxi+dershape2D(1,ia,i,j)*xelm(ia)
        xeta=xeta+dershape2D(2,ia,i,j)*xelm(ia)
        yxi=yxi+dershape2D(1,ia,i,j)*yelm(ia)
        yeta=yeta+dershape2D(2,ia,i,j)*yelm(ia)
        zxi=zxi+dershape2D(1,ia,i,j)*zelm(ia)
        zeta=zeta+dershape2D(2,ia,i,j)*zelm(ia)
      enddo

      !   calculate the unnormalized normal to the boundary
      unx=yxi*zeta-yeta*zxi
      uny=zxi*xeta-zeta*xxi
      unz=xxi*yeta-xeta*yxi
      jacobian=dsqrt(unx*unx+uny*uny+unz*unz)
      if(jacobian <= ZERO) call exit_MPI(myrank,'2D Jacobian undefined')

      !   normalize normal vector and store weighted surface jacobian

      ! distinguish if single or double precision for reals
      if(CUSTOM_REAL == SIZE_REAL) then
        jacobian2Dw_face(i,j) = sngl(jacobian * wgllwgll(i,j))
        normal_face(1,i,j)=sngl(unx/jacobian)
        normal_face(2,i,j)=sngl(uny/jacobian)
        normal_face(3,i,j)=sngl(unz/jacobian)
      else
        jacobian2Dw_face(i,j) = jacobian * wgllwgll(i,j)
        normal_face(1,i,j)=unx/jacobian
        normal_face(2,i,j)=uny/jacobian
        normal_face(3,i,j)=unz/jacobian
      endif

    enddo
  enddo

  end subroutine compute_jacobian_2D_face

