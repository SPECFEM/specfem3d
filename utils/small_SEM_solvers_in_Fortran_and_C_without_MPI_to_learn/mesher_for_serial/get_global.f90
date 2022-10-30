!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2008
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

  subroutine get_global(npointot,xp,yp,zp,iglob,loc,ifseg,nglob)

! this routine MUST be in double precision to avoid sensitivity
! to roundoff errors in the coordinates of the points

! non-structured global numbering software provided by Paul F. Fischer

! leave sorting subroutines in same source file to allow for inlining

  implicit none

  include "constants.h"

! parameters
  integer, intent(in) :: npointot
  double precision, intent(in) :: xp(npointot),yp(npointot),zp(npointot)

  integer, intent(out) :: iglob(npointot),loc(npointot)
  logical, intent(out) :: ifseg(npointot)
  integer, intent(out) :: nglob

! variables
  integer ier

  integer, dimension(:), allocatable :: ninseg,idummy

! dynamically allocate arrays
  allocate(ninseg(npointot),stat=ier)
  if (ier /= 0) stop 'error allocating array ninseg'
  allocate(idummy(npointot),stat=ier)
  if (ier /= 0) stop 'error allocating array idummy'

  call sort_array_coordinates(npointot,xp,yp,zp,idummy,iglob,loc,ifseg, &
                              nglob,ninseg)

! deallocate arrays
  deallocate(ninseg)
  deallocate(idummy)

  end subroutine get_global
