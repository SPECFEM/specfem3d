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

  subroutine socal_model(idoubling,zmesh,rho,vp,vs)

  implicit none

  include "constants.h"

  integer idoubling
  double precision zmesh
  double precision rho,vp,vs

  if(idoubling == IFLAG_HALFSPACE_MOHO) then
        vp=7.8d0
        vs=4.5d0
        rho=3.0d0

  else if(idoubling == IFLAG_MOHO_16km) then
        vp=6.7d0
        vs=3.87d0
        rho=2.8d0

  else if(idoubling == IFLAG_ONE_LAYER_TOPOGRAPHY .or. idoubling == IFLAG_BASEMENT_TOPO) then
        vp=5.5d0
        vs=3.18d0
        rho=2.4d0

  else
        vp=6.3d0
        vs=3.64d0
        rho=2.67d0
  endif

! scale to standard units
  vp = vp * 1000.d0
  vs = vs * 1000.d0
  rho = rho * 1000.d0

  end subroutine socal_model

