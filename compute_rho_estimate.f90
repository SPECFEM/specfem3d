!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 4
!          --------------------------------------------------
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

  subroutine compute_rho_estimate(rho,vp)

! compute rho estimate in Gocad block and in Hauksson's model
! based upon Vp

  implicit none

  include "constants.h"

  double precision rho,vp

! scale density - use empirical rule from Christiane
  rho = 0.33d0 * vp + 1280.d0

! make sure density estimate is reasonable
  if(rho > DENSITY_MAX) rho = DENSITY_MAX
  if(rho < DENSITY_MIN) rho = DENSITY_MIN

  end subroutine compute_rho_estimate

