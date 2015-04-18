!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2008
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

  subroutine get_topo_bathy(xlat,xlon,value,ibathy_topo)

!
!---- get elevation or ocean depth in meters at a given latitude and longitude
!

  implicit none

  include "constants.h"

! use integer array to store values
  integer, dimension(NX_BATHY,NY_BATHY) :: ibathy_topo

  double precision xlat,xlon,value

  integer iadd1,iel1
  double precision samples_per_degree_topo
  double precision xlo

  xlo = xlon
  if(xlon < 0.d0) xlo = xlo + 360.d0

! compute number of samples per degree
  samples_per_degree_topo = dble(RESOLUTION_TOPO_FILE) / 60.d0

! compute offset in data file and avoid edge effects
  iadd1 = 1 + int((90.d0-xlat)/samples_per_degree_topo)
  if(iadd1 < 1) iadd1 = 1
  if(iadd1 > NY_BATHY) iadd1 = NY_BATHY

  iel1 = int(xlo/samples_per_degree_topo)
  if(iel1 <= 0 .or. iel1 > NX_BATHY) iel1 = NX_BATHY

! convert integer value to double precision
  value = dble(ibathy_topo(iel1,iadd1))

  end subroutine get_topo_bathy

! -------------------------------------------

  subroutine read_topo_bathy_file(ibathy_topo)
!
!---- read topography and bathymetry file once and for all
!
  implicit none

  include "constants.h"

  character(len=150) topo_bathy_file

! use integer array to store values
  integer, dimension(NX_BATHY,NY_BATHY) :: ibathy_topo

  integer itopo_x,itopo_y

  call get_value_string(topo_bathy_file, 'model.topoBathy.PATHNAME_TOPO_FILE', PATHNAME_TOPO_FILE)

  open(unit=13,file=topo_bathy_file,status='old',action='read')

  do itopo_y=1,NY_BATHY
    do itopo_x=1,NX_BATHY

      read(13,*) ibathy_topo(itopo_x,itopo_y)

! impose maximum height of mountains, to suppress oscillations in Himalaya etc.
  if(USE_MAXIMUM_HEIGHT_TOPO .and. ibathy_topo(itopo_x,itopo_y) > MAXIMUM_HEIGHT_TOPO) &
    ibathy_topo(itopo_x,itopo_y) = MAXIMUM_HEIGHT_TOPO

! impose maximum depth of oceans, to suppress oscillations near deep trenches
  if(USE_MAXIMUM_DEPTH_OCEANS .and. ibathy_topo(itopo_x,itopo_y) < MAXIMUM_DEPTH_OCEANS) &
    ibathy_topo(itopo_x,itopo_y) = MAXIMUM_DEPTH_OCEANS

    enddo
  enddo

  close(13)

  end subroutine read_topo_bathy_file

