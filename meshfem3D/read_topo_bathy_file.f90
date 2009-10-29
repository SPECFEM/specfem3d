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

  subroutine read_topo_bathy_file(itopo_bathy,NX_TOPO,NY_TOPO,topo_file)
!
!---- read topography and bathymetry file once and for all
!
  implicit none

  include "constants.h"

  integer NX_TOPO,NY_TOPO

! use integer array to store topography values
  integer itopo_bathy(NX_TOPO,NY_TOPO)

  character(len=100) topo_file

  integer ix,iy

  itopo_bathy(:,:) = 0

  open(unit=13,file=topo_file,status='old',action='read')
  do iy=1,NY_TOPO
    do ix=1,NX_TOPO
      read(13,*) itopo_bathy(ix,iy)
    enddo
  enddo
  close(13)

  end subroutine read_topo_bathy_file

