!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and University of Pau / CNRS / INRIA
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            November 2010
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

  subroutine read_moho_map(imoho_depth)
!
!---- read Lupei Zhu's Moho map of Southern California once and for all
!
  implicit none

  include "constants.h"

! use integer array to store Moho depth
  integer imoho_depth(NX_MOHO,NY_MOHO)

  integer ix,iy

  double precision long,lat,depth_km

  character(len=256) MOHO_MAP_FILE

  imoho_depth(:,:) = 0

  call get_value_string(MOHO_MAP_FILE, &
                        'model.MOHO_MAP_FILE', &
                        'DATA/moho_map/moho_lupei_zhu.dat')
  open(unit=13,file=MOHO_MAP_FILE,status='old',action='read')
! file starts from North-West corner
  do iy=NY_MOHO,1,-1
    do ix=1,NX_MOHO
      read(13,*) long,lat,depth_km
! convert depth to meters
      imoho_depth(ix,iy) = nint(depth_km * 1000.d0)
    enddo
  enddo
  close(13)

  end subroutine read_moho_map

