!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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

  subroutine get_jacobian_boundary_face(nspec, &
                                        xstore,ystore,zstore,ibool,nglob, &
                                        dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
                                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                                        ispec,iface,jacobian2Dw_face,normal_face,NGLLA,NGLLB,NGNOD2D)

! returns jacobian2Dw_face and normal_face (pointing outwards of element)

  use constants, only: CUSTOM_REAL,NDIM,NGLLX,NGLLY,NGLLZ,MIDX,MIDY,MIDZ,NDIM2D,myrank

  implicit none

  integer,intent(in) :: nspec,nglob

  ! arrays with the mesh
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob),intent(in) :: xstore,ystore,zstore

  ! face information
  integer,intent(in) :: iface,ispec,NGLLA,NGLLB
  real(kind=CUSTOM_REAL), dimension(NGLLA,NGLLB),intent(inout) :: jacobian2Dw_face
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLA,NGLLB),intent(inout) :: normal_face

  integer,intent(in) :: NGNOD2D
  double precision,intent(in) :: dershape2D_x(NDIM2D,NGNOD2D,NGLLY,NGLLZ)
  double precision,intent(in) :: dershape2D_y(NDIM2D,NGNOD2D,NGLLX,NGLLZ)
  double precision,intent(in) :: dershape2D_bottom(NDIM2D,NGNOD2D,NGLLX,NGLLY)
  double precision,intent(in) :: dershape2D_top(NDIM2D,NGNOD2D,NGLLX,NGLLY)

  double precision, dimension(NGLLX,NGLLY),intent(in) :: wgllwgll_xy
  double precision, dimension(NGLLX,NGLLZ),intent(in) :: wgllwgll_xz
  double precision, dimension(NGLLY,NGLLZ),intent(in) :: wgllwgll_yz

! local parameters
  ! face corners
  double precision xelm(NGNOD2D),yelm(NGNOD2D),zelm(NGNOD2D)

  ! check that the parameter file is correct
  if (NGNOD2D /= 4 .and. NGNOD2D /= 9) &
    call exit_MPI(myrank,'surface elements should have 4 or 9 control nodes')

  select case (iface)
  ! on reference face: xmin
  case (1)
    xelm(1)=xstore( ibool(1,1,1,ispec) )
    yelm(1)=ystore( ibool(1,1,1,ispec) )
    zelm(1)=zstore( ibool(1,1,1,ispec) )
    xelm(2)=xstore( ibool(1,NGLLY,1,ispec) )
    yelm(2)=ystore( ibool(1,NGLLY,1,ispec) )
    zelm(2)=zstore( ibool(1,NGLLY,1,ispec) )
    xelm(3)=xstore( ibool(1,NGLLY,NGLLZ,ispec) )
    yelm(3)=ystore( ibool(1,NGLLY,NGLLZ,ispec) )
    zelm(3)=zstore( ibool(1,NGLLY,NGLLZ,ispec) )
    xelm(4)=xstore( ibool(1,1,NGLLZ,ispec) )
    yelm(4)=ystore( ibool(1,1,NGLLZ,ispec) )
    zelm(4)=zstore( ibool(1,1,NGLLZ,ispec) )

    if (NGNOD2D == 9) then
      xelm(5)=xstore( ibool(1,MIDY,1,ispec) )
      yelm(5)=ystore( ibool(1,MIDY,1,ispec) )
      zelm(5)=zstore( ibool(1,MIDY,1,ispec) )
      xelm(6)=xstore( ibool(1,NGLLY,MIDZ,ispec) )
      yelm(6)=ystore( ibool(1,NGLLY,MIDZ,ispec) )
      zelm(6)=zstore( ibool(1,NGLLY,MIDZ,ispec) )
      xelm(7)=xstore( ibool(1,MIDY,NGLLZ,ispec) )
      yelm(7)=ystore( ibool(1,MIDY,NGLLZ,ispec) )
      zelm(7)=zstore( ibool(1,MIDY,NGLLZ,ispec) )
      xelm(8)=xstore( ibool(1,1,MIDZ,ispec) )
      yelm(8)=ystore( ibool(1,1,MIDZ,ispec) )
      zelm(8)=zstore( ibool(1,1,MIDZ,ispec) )
      xelm(9)=xstore( ibool(1,MIDY,MIDZ,ispec) )
      yelm(9)=ystore( ibool(1,MIDY,MIDZ,ispec) )
      zelm(9)=zstore( ibool(1,MIDY,MIDZ,ispec) )
    endif

    call compute_jacobian_2D_face(xelm,yelm,zelm, &
                                  dershape2D_x,wgllwgll_yz, &
                                  jacobian2Dw_face,normal_face,NGLLY,NGLLZ,NGNOD2D)

  ! on boundary: xmax
  case (2)
    xelm(1)=xstore( ibool(NGLLX,1,1,ispec) )
    yelm(1)=ystore( ibool(NGLLX,1,1,ispec) )
    zelm(1)=zstore( ibool(NGLLX,1,1,ispec) )
    xelm(2)=xstore( ibool(NGLLX,NGLLY,1,ispec) )
    yelm(2)=ystore( ibool(NGLLX,NGLLY,1,ispec) )
    zelm(2)=zstore( ibool(NGLLX,NGLLY,1,ispec) )
    xelm(3)=xstore( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    yelm(3)=ystore( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    zelm(3)=zstore( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    xelm(4)=xstore( ibool(NGLLX,1,NGLLZ,ispec) )
    yelm(4)=ystore( ibool(NGLLX,1,NGLLZ,ispec) )
    zelm(4)=zstore( ibool(NGLLX,1,NGLLZ,ispec) )

    if (NGNOD2D == 9) then
      xelm(5)=xstore( ibool(NGLLX,MIDY,1,ispec) )
      yelm(5)=ystore( ibool(NGLLX,MIDY,1,ispec) )
      zelm(5)=zstore( ibool(NGLLX,MIDY,1,ispec) )
      xelm(6)=xstore( ibool(NGLLX,NGLLY,MIDZ,ispec) )
      yelm(6)=ystore( ibool(NGLLX,NGLLY,MIDZ,ispec) )
      zelm(6)=zstore( ibool(NGLLX,NGLLY,MIDZ,ispec) )
      xelm(7)=xstore( ibool(NGLLX,MIDY,NGLLZ,ispec) )
      yelm(7)=ystore( ibool(NGLLX,MIDY,NGLLZ,ispec) )
      zelm(7)=zstore( ibool(NGLLX,MIDY,NGLLZ,ispec) )
      xelm(8)=xstore( ibool(NGLLX,1,MIDZ,ispec) )
      yelm(8)=ystore( ibool(NGLLX,1,MIDZ,ispec) )
      zelm(8)=zstore( ibool(NGLLX,1,MIDZ,ispec) )
      xelm(9)=xstore( ibool(NGLLX,MIDY,MIDZ,ispec) )
      yelm(9)=ystore( ibool(NGLLX,MIDY,MIDZ,ispec) )
      zelm(9)=zstore( ibool(NGLLX,MIDY,MIDZ,ispec) )
    endif

    call compute_jacobian_2D_face(xelm,yelm,zelm, &
                                  dershape2D_x,wgllwgll_yz, &
                                  jacobian2Dw_face,normal_face,NGLLY,NGLLZ,NGNOD2D)

  ! on boundary: ymin
  case (3)
    xelm(1)=xstore( ibool(1,1,1,ispec) )
    yelm(1)=ystore( ibool(1,1,1,ispec) )
    zelm(1)=zstore( ibool(1,1,1,ispec) )
    xelm(2)=xstore( ibool(NGLLX,1,1,ispec) )
    yelm(2)=ystore( ibool(NGLLX,1,1,ispec) )
    zelm(2)=zstore( ibool(NGLLX,1,1,ispec) )
    xelm(3)=xstore( ibool(NGLLX,1,NGLLZ,ispec) )
    yelm(3)=ystore( ibool(NGLLX,1,NGLLZ,ispec) )
    zelm(3)=zstore( ibool(NGLLX,1,NGLLZ,ispec) )
    xelm(4)=xstore( ibool(1,1,NGLLZ,ispec) )
    yelm(4)=ystore( ibool(1,1,NGLLZ,ispec) )
    zelm(4)=zstore( ibool(1,1,NGLLZ,ispec) )

    if (NGNOD2D == 9) then
      xelm(5)=xstore( ibool(MIDX,1,1,ispec) )
      yelm(5)=ystore( ibool(MIDX,1,1,ispec) )
      zelm(5)=zstore( ibool(MIDX,1,1,ispec) )
      xelm(6)=xstore( ibool(NGLLX,1,MIDZ,ispec) )
      yelm(6)=ystore( ibool(NGLLX,1,MIDZ,ispec) )
      zelm(6)=zstore( ibool(NGLLX,1,MIDZ,ispec) )
      xelm(7)=xstore( ibool(MIDX,1,NGLLZ,ispec) )
      yelm(7)=ystore( ibool(MIDX,1,NGLLZ,ispec) )
      zelm(7)=zstore( ibool(MIDX,1,NGLLZ,ispec) )
      xelm(8)=xstore( ibool(1,1,MIDZ,ispec) )
      yelm(8)=ystore( ibool(1,1,MIDZ,ispec) )
      zelm(8)=zstore( ibool(1,1,MIDZ,ispec) )
      xelm(9)=xstore( ibool(MIDX,1,MIDZ,ispec) )
      yelm(9)=ystore( ibool(MIDX,1,MIDZ,ispec) )
      zelm(9)=zstore( ibool(MIDX,1,MIDZ,ispec) )
    endif

    call compute_jacobian_2D_face(xelm,yelm,zelm, &
                                  dershape2D_y,wgllwgll_xz, &
                                  jacobian2Dw_face,normal_face,NGLLX,NGLLZ,NGNOD2D)

  ! on boundary: ymax
  case (4)
    xelm(1)=xstore( ibool(1,NGLLY,1,ispec) )
    yelm(1)=ystore( ibool(1,NGLLY,1,ispec) )
    zelm(1)=zstore( ibool(1,NGLLY,1,ispec) )
    xelm(2)=xstore( ibool(NGLLX,NGLLY,1,ispec) )
    yelm(2)=ystore( ibool(NGLLX,NGLLY,1,ispec) )
    zelm(2)=zstore( ibool(NGLLX,NGLLY,1,ispec) )
    xelm(3)=xstore( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    yelm(3)=ystore( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    zelm(3)=zstore( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    xelm(4)=xstore( ibool(1,NGLLY,NGLLZ,ispec) )
    yelm(4)=ystore( ibool(1,NGLLY,NGLLZ,ispec) )
    zelm(4)=zstore( ibool(1,NGLLY,NGLLZ,ispec) )

    if (NGNOD2D == 9) then
      xelm(5)=xstore( ibool(MIDX,NGLLY,1,ispec) )
      yelm(5)=ystore( ibool(MIDX,NGLLY,1,ispec) )
      zelm(5)=zstore( ibool(MIDX,NGLLY,1,ispec) )
      xelm(6)=xstore( ibool(NGLLX,NGLLY,MIDZ,ispec) )
      yelm(6)=ystore( ibool(NGLLX,NGLLY,MIDZ,ispec) )
      zelm(6)=zstore( ibool(NGLLX,NGLLY,MIDZ,ispec) )
      xelm(7)=xstore( ibool(MIDX,NGLLY,NGLLZ,ispec) )
      yelm(7)=ystore( ibool(MIDX,NGLLY,NGLLZ,ispec) )
      zelm(7)=zstore( ibool(MIDX,NGLLY,NGLLZ,ispec) )
      xelm(8)=xstore( ibool(1,NGLLY,MIDZ,ispec) )
      yelm(8)=ystore( ibool(1,NGLLY,MIDZ,ispec) )
      zelm(8)=zstore( ibool(1,NGLLY,MIDZ,ispec) )
      xelm(9)=xstore( ibool(MIDX,NGLLY,MIDZ,ispec) )
      yelm(9)=ystore( ibool(MIDX,NGLLY,MIDZ,ispec) )
      zelm(9)=zstore( ibool(MIDX,NGLLY,MIDZ,ispec) )
    endif

    call compute_jacobian_2D_face(xelm,yelm,zelm, &
                                  dershape2D_y, wgllwgll_xz, &
                                  jacobian2Dw_face,normal_face,NGLLX,NGLLZ,NGNOD2D)


  ! on boundary: bottom
  case (5)
    xelm(1)=xstore( ibool(1,1,1,ispec) )
    yelm(1)=ystore( ibool(1,1,1,ispec) )
    zelm(1)=zstore( ibool(1,1,1,ispec) )
    xelm(2)=xstore( ibool(NGLLX,1,1,ispec) )
    yelm(2)=ystore( ibool(NGLLX,1,1,ispec) )
    zelm(2)=zstore( ibool(NGLLX,1,1,ispec) )
    xelm(3)=xstore( ibool(NGLLX,NGLLY,1,ispec) )
    yelm(3)=ystore( ibool(NGLLX,NGLLY,1,ispec) )
    zelm(3)=zstore( ibool(NGLLX,NGLLY,1,ispec) )
    xelm(4)=xstore( ibool(1,NGLLY,1,ispec) )
    yelm(4)=ystore( ibool(1,NGLLY,1,ispec) )
    zelm(4)=zstore( ibool(1,NGLLY,1,ispec) )

    if (NGNOD2D == 9) then
      xelm(5)=xstore( ibool(MIDX,1,1,ispec) )
      yelm(5)=ystore( ibool(MIDX,1,1,ispec) )
      zelm(5)=zstore( ibool(MIDX,1,1,ispec) )
      xelm(6)=xstore( ibool(NGLLX,MIDY,1,ispec) )
      yelm(6)=ystore( ibool(NGLLX,MIDY,1,ispec) )
      zelm(6)=zstore( ibool(NGLLX,MIDY,1,ispec) )
      xelm(7)=xstore( ibool(MIDX,NGLLY,1,ispec) )
      yelm(7)=ystore( ibool(MIDX,NGLLY,1,ispec) )
      zelm(7)=zstore( ibool(MIDX,NGLLY,1,ispec) )
      xelm(8)=xstore( ibool(1,MIDY,1,ispec) )
      yelm(8)=ystore( ibool(1,MIDY,1,ispec) )
      zelm(8)=zstore( ibool(1,MIDY,1,ispec) )
      xelm(9)=xstore( ibool(MIDX,MIDY,1,ispec) )
      yelm(9)=ystore( ibool(MIDX,MIDY,1,ispec) )
      zelm(9)=zstore( ibool(MIDX,MIDY,1,ispec) )
    endif

    call compute_jacobian_2D_face(xelm,yelm,zelm, &
                                  dershape2D_bottom,wgllwgll_xy, &
                                  jacobian2Dw_face,normal_face,NGLLX,NGLLY,NGNOD2D)

  ! on boundary: top
  case (6)
    xelm(1)=xstore( ibool(1,1,NGLLZ,ispec) )
    yelm(1)=ystore( ibool(1,1,NGLLZ,ispec) )
    zelm(1)=zstore( ibool(1,1,NGLLZ,ispec) )
    xelm(2)=xstore( ibool(NGLLX,1,NGLLZ,ispec) )
    yelm(2)=ystore( ibool(NGLLX,1,NGLLZ,ispec) )
    zelm(2)=zstore( ibool(NGLLX,1,NGLLZ,ispec) )
    xelm(3)=xstore( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    yelm(3)=ystore( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    zelm(3)=zstore( ibool(NGLLX,NGLLY,NGLLZ,ispec) )
    xelm(4)=xstore( ibool(1,NGLLY,NGLLZ,ispec) )
    yelm(4)=ystore( ibool(1,NGLLY,NGLLZ,ispec) )
    zelm(4)=zstore( ibool(1,NGLLY,NGLLZ,ispec) )

    if (NGNOD2D == 9) then
      xelm(5)=xstore( ibool(MIDX,1,NGLLZ,ispec) )
      yelm(5)=ystore( ibool(MIDX,1,NGLLZ,ispec) )
      zelm(5)=zstore( ibool(MIDX,1,NGLLZ,ispec) )
      xelm(6)=xstore( ibool(NGLLX,MIDY,NGLLZ,ispec) )
      yelm(6)=ystore( ibool(NGLLX,MIDY,NGLLZ,ispec) )
      zelm(6)=zstore( ibool(NGLLX,MIDY,NGLLZ,ispec) )
      xelm(7)=xstore( ibool(MIDX,NGLLY,NGLLZ,ispec) )
      yelm(7)=ystore( ibool(MIDX,NGLLY,NGLLZ,ispec) )
      zelm(7)=zstore( ibool(MIDX,NGLLY,NGLLZ,ispec) )
      xelm(8)=xstore( ibool(1,MIDY,NGLLZ,ispec) )
      yelm(8)=ystore( ibool(1,MIDY,NGLLZ,ispec) )
      zelm(8)=zstore( ibool(1,MIDY,NGLLZ,ispec) )
      xelm(9)=xstore( ibool(MIDX,MIDY,NGLLZ,ispec) )
      yelm(9)=ystore( ibool(MIDX,MIDY,NGLLZ,ispec) )
      zelm(9)=zstore( ibool(MIDX,MIDY,NGLLZ,ispec) )
    endif

    call compute_jacobian_2D_face(xelm,yelm,zelm, &
                                  dershape2D_top, wgllwgll_xy, &
                                  jacobian2Dw_face,normal_face,NGLLX,NGLLY,NGNOD2D)

  case default
    stop 'error 2D jacobian'
  end select

  end subroutine get_jacobian_boundary_face

!
! ------------------------------------------------------------------------------------------------
!

  subroutine compute_jacobian_2D_face(xelm,yelm,zelm, &
                                      dershape2D,wgllwgll, &
                                      jacobian2Dw_face,normal_face,NGLLA,NGLLB,NGNOD2D)

  use constants

  implicit none

! generic routine that accepts any polynomial degree in each direction
! returns 2D jacobian and normal for this face only

  integer :: NGLLA,NGLLB,NGNOD2D

  double precision :: xelm(NGNOD2D),yelm(NGNOD2D),zelm(NGNOD2D)
  double precision :: dershape2D(NDIM2D,NGNOD2D,NGLLA,NGLLB)
  double precision, dimension(NGLLA,NGLLB) :: wgllwgll

  real(kind=CUSTOM_REAL), dimension(NGLLA,NGLLB) :: jacobian2Dw_face
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLA,NGLLB) :: normal_face

  integer :: i,j,ia
  double precision :: xxi,xeta,yxi,yeta,zxi,zeta
  double precision :: unx,uny,unz,jacobian

  do j = 1,NGLLB
    do i = 1,NGLLA

      xxi = ZERO
      xeta = ZERO
      yxi = ZERO
      yeta = ZERO
      zxi = ZERO
      zeta = ZERO
      do ia = 1,NGNOD2D
        xxi = xxi+dershape2D(1,ia,i,j)*xelm(ia)
        xeta = xeta+dershape2D(2,ia,i,j)*xelm(ia)

        yxi = yxi+dershape2D(1,ia,i,j)*yelm(ia)
        yeta = yeta+dershape2D(2,ia,i,j)*yelm(ia)

        zxi = zxi+dershape2D(1,ia,i,j)*zelm(ia)
        zeta = zeta+dershape2D(2,ia,i,j)*zelm(ia)
      enddo

      !   calculate the unnormalized normal to the boundary
      unx = yxi*zeta-yeta*zxi
      uny = zxi*xeta-zeta*xxi
      unz = xxi*yeta-xeta*yxi
      jacobian = dsqrt(unx*unx+uny*uny+unz*unz)
      if (jacobian <= ZERO) call exit_MPI(myrank,'2D Jacobian undefined')

      !   normalize normal vector and store weighted surface jacobian
      ! distinguish if single or double precision for reals
      jacobian2Dw_face(i,j) = real(jacobian * wgllwgll(i,j),kind=CUSTOM_REAL)

      normal_face(1,i,j) = real(unx/jacobian,kind=CUSTOM_REAL)
      normal_face(2,i,j) = real(uny/jacobian,kind=CUSTOM_REAL)
      normal_face(3,i,j) = real(unz/jacobian,kind=CUSTOM_REAL)

    enddo
  enddo

  end subroutine compute_jacobian_2D_face

