!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and University of Pau / CNRS / INRIA
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            April 2011
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
!
! United States and French Government Sponsorship Acknowledged.


! obsolete...

!  subroutine read_topography_bathymetry()
!
!  use specfem_par
!  implicit none
!  integer :: ier
!
!! read topography and bathymetry file
!
!  if( OCEANS .and. TOPOGRAPHY ) then
!
!    NX_TOPO = NX_TOPO_FILE
!    NY_TOPO = NY_TOPO_FILE
!    allocate(itopo_bathy(NX_TOPO,NY_TOPO),stat=ier)
!    if( ier /= 0 ) stop 'error allocating array itopo_bathy'
!
!    call read_topo_bathy_file(itopo_bathy,NX_TOPO,NY_TOPO)
!
!    if(myrank == 0) then
!      write(IMAIN,*)
!      write(IMAIN,*) 'regional topography file read ranges in m from ', &
!        minval(itopo_bathy),' to ',maxval(itopo_bathy)
!      write(IMAIN,*)
!    endif
!
!  else
!    NX_TOPO = 1
!    NY_TOPO = 1
!    allocate(itopo_bathy(NX_TOPO,NY_TOPO),stat=ier)
!    if( ier /= 0 ) stop 'error allocating dummy array itopo_bathy'
!
!  endif
!
!  end subroutine read_topography_bathymetry
