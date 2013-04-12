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

  subroutine store_coords(xstore,ystore,zstore,xelm,yelm,zelm,ispec,nspec)

  implicit none

  include "constants.h"
  include "constants_meshfem3D.h"

  integer ispec,nspec

  double precision, dimension(NGNOD_EIGHT_CORNERS) :: xelm,yelm,zelm

  double precision, dimension(NGLLX_M,NGLLY_M,NGLLZ_M,nspec) :: xstore,ystore,zstore

  xstore(1,1,1,ispec) = xelm(1)
  ystore(1,1,1,ispec) = yelm(1)
  zstore(1,1,1,ispec) = zelm(1)

  xstore(2,1,1,ispec) = xelm(2)
  ystore(2,1,1,ispec) = yelm(2)
  zstore(2,1,1,ispec) = zelm(2)

  xstore(2,2,1,ispec) = xelm(3)
  ystore(2,2,1,ispec) = yelm(3)
  zstore(2,2,1,ispec) = zelm(3)

  xstore(1,2,1,ispec) = xelm(4)
  ystore(1,2,1,ispec) = yelm(4)
  zstore(1,2,1,ispec) = zelm(4)

  xstore(1,1,2,ispec) = xelm(5)
  ystore(1,1,2,ispec) = yelm(5)
  zstore(1,1,2,ispec) = zelm(5)

  xstore(2,1,2,ispec) = xelm(6)
  ystore(2,1,2,ispec) = yelm(6)
  zstore(2,1,2,ispec) = zelm(6)

  xstore(2,2,2,ispec) = xelm(7)
  ystore(2,2,2,ispec) = yelm(7)
  zstore(2,2,2,ispec) = zelm(7)

  xstore(1,2,2,ispec) = xelm(8)
  ystore(1,2,2,ispec) = yelm(8)
  zstore(1,2,2,ispec) = zelm(8)

  end subroutine store_coords

