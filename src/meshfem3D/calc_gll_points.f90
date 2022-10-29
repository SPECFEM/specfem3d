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


  subroutine calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)

  implicit none

  integer,intent(in) :: NGNOD,NGLLX,NGLLY,NGLLZ
  double precision,intent(in) :: shape3D(NGNOD,NGLLX,NGLLY,NGLLZ)
  double precision,intent(inout) :: xstore(NGLLX,NGLLY,NGLLZ),ystore(NGLLX,NGLLY,NGLLZ),zstore(NGLLX,NGLLY,NGLLZ)
  double precision,intent(in) :: xelm(NGNOD),yelm(NGNOD),zelm(NGNOD)

  ! local parameters
  double precision :: xmesh,ymesh,zmesh
  double precision, parameter :: ZERO = 0.d0
  integer :: ia,i,j,k

  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX
        ! compute mesh coordinates
        xmesh = ZERO
        ymesh = ZERO
        zmesh = ZERO
        do ia = 1,NGNOD
          xmesh = xmesh + shape3D(ia,i,j,k)*xelm(ia)
          ymesh = ymesh + shape3D(ia,i,j,k)*yelm(ia)
          zmesh = zmesh + shape3D(ia,i,j,k)*zelm(ia)
        enddo
        xstore(i,j,k) = xmesh
        ystore(i,j,k) = ymesh
        zstore(i,j,k) = zmesh
      enddo
    enddo
  enddo

  end subroutine calc_gll_points
