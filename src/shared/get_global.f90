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

  subroutine get_global(npointot,xp,yp,zp,iglob,locval,ifseg,nglob,UTM_X_MIN,UTM_X_MAX)

! this routine MUST be in double precision to avoid sensitivity
! to roundoff errors in the coordinates of the points

! non-structured global numbering software provided by Paul F. Fischer

! leave sorting subroutines in same source file to allow for inlining

  use constants, only: SMALLVAL_TOL

  implicit none

  integer :: npointot
  integer :: nglob
  integer :: iglob(npointot),locval(npointot)
  logical :: ifseg(npointot)
  double precision :: xp(npointot),yp(npointot),zp(npointot)
  double precision :: UTM_X_MIN,UTM_X_MAX

  integer :: ier

  integer, dimension(:), allocatable :: ninseg,idummy

! geometry tolerance parameter to calculate number of independent grid points
! small value for double precision and to avoid sensitivity to roundoff
  double precision :: SMALLVALTOL

  ! note: the tolerance value is important to adapt according to the mesh dimensions
  !       having the tolerance too small, e.g., for double precision runs, will lead to numerical artifacts
  !       when different partitions have slightly different node positions due to numerical round-offs.
  !       we adapt the tolerance to the dimensions of the mesh (using the x-direction as a typical dimension value),
  !       to avoid round-off problems for other kind of meshes
  !       (e.g., UTM meshes have often dimensions in ~ 10 km, local meshes for engineering ~ 1 m, ultrasonic ~ 1mm range).
  !
  !       one way to avoid this mesh size dependence would be to normalize the coordinates as done in the global version.
  !       again, that would need a typical dimension of the mesh geometry used. to check for the future...

  ! define geometrical tolerance based upon typical size of the model
  SMALLVALTOL = SMALLVAL_TOL * dabs(UTM_X_MAX - UTM_X_MIN)

  ! dynamically allocate arrays
  allocate(ninseg(npointot), &
           idummy(npointot),stat=ier)
  if (ier /= 0) stop 'Error allocating ninseg arrays'

  call sort_array_coordinates(npointot,xp,yp,zp,idummy,iglob,locval,ifseg, &
                              nglob,ninseg,SMALLVALTOL)

  ! deallocate arrays
  deallocate(ninseg)
  deallocate(idummy)

  end subroutine get_global

! ------------------------------------------------------------------


  subroutine get_global_indirect_addressing(nspec,nglob,ibool)

!
!- we can create a new indirect addressing to reduce cache misses
! (put into this subroutine but compiler keeps on complaining that it cannot vectorize loops...)
!! DK DK
!! DK DK answer from Dimitri, April 2014: that is normal because the nested loops have a dependency
!! DK DK (they can write to the same memory location through the mask_ibool() array) and thus
!! DK DK they cannot be vectorized. Thus the compiler is right.
!! DK DK

  use constants, only: NGLLX,NGLLY,NGLLZ

  implicit none

  integer :: nspec,nglob
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

! mask to sort ibool
  integer, dimension(:), allocatable :: mask_ibool
  integer, dimension(:,:,:,:), allocatable :: copy_ibool_ori
  integer :: inumber
  integer:: i,j,k,ispec,ier

! copies original array
  allocate(copy_ibool_ori(NGLLX,NGLLY,NGLLZ,nspec), &
           mask_ibool(nglob),stat=ier)
  if (ier /= 0) stop 'Error in allocating copy_ibool_ori arrays'

  mask_ibool(:) = -1
  copy_ibool_ori(:,:,:,:) = ibool(:,:,:,:)

! reduces misses
  inumber = 0
  do ispec=1,nspec
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          if (mask_ibool(copy_ibool_ori(i,j,k,ispec)) == -1) then
! create a new point
            inumber = inumber + 1
            ibool(i,j,k,ispec) = inumber
            mask_ibool(copy_ibool_ori(i,j,k,ispec)) = inumber
          else
! use an existing point created previously
            ibool(i,j,k,ispec) = mask_ibool(copy_ibool_ori(i,j,k,ispec))
          endif
        enddo
      enddo
    enddo
  enddo

! cleanup
  deallocate(copy_ibool_ori,mask_ibool)

  end subroutine get_global_indirect_addressing

