!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 2
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology July 2004
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine read_basin_topo_bathy_file(itopo_bathy_basin)
!
!---- read basin topography and bathymetry file once and for all
!
  implicit none

  include "constants.h"

! use integer array to store topography values
  integer itopo_bathy_basin(NX_TOPO,NY_TOPO)

  integer ix,iy

  itopo_bathy_basin(:,:) = 0

!! DK DK UGLY LACQ  open(unit=13,file='DATA/la_topography/topo_bathy_final.dat',status='old')
  open(unit=13,file='DATA/lacq_thomas/mnt_Lacq_Lambert_final_dimitri.dat',status='old')
  do iy=1,NY_TOPO
    do ix=1,NX_TOPO
      read(13,*) itopo_bathy_basin(ix,iy)
    enddo
  enddo
  close(13)

  end subroutine read_basin_topo_bathy_file

