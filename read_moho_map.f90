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

  imoho_depth(:,:) = 0

  open(unit=13,file='DATA/moho_map/moho_lupei_zhu.dat',status='old')
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

