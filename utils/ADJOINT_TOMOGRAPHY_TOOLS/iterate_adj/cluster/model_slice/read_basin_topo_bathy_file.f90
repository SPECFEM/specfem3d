!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine read_basin_topo_bathy_file(itopo_bathy_basin,NX_TOPO,NY_TOPO,topo_file)
!
!---- read basin topography and bathymetry file once and for all
!
  implicit none

  include "constants.h"

  integer NX_TOPO,NY_TOPO

! use integer array to store topography values
  integer itopo_bathy_basin(NX_TOPO,NY_TOPO)

  character(len=100) topo_file

  integer ix,iy

  itopo_bathy_basin(:,:) = 0

  open(unit=13,file=topo_file,status='old',action='read')
  do iy=1,NY_TOPO
    do ix=1,NX_TOPO
      read(13,*) itopo_bathy_basin(ix,iy)
    enddo
  enddo
  close(13)

  end subroutine read_basin_topo_bathy_file

