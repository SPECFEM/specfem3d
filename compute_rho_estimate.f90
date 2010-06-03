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

  subroutine compute_rho_estimate(rho,vp)

! compute rho estimate in Gocad block and in Hauksson model from Vp
! -- this may be called in create_regions_mesh.f90 to apply the
!    Vp-rho scaling throughout the entire volume

  implicit none

  include "constants.h"

  double precision rho,vp
  double precision :: P0, P1, P2, P3, P4, P5

!!$! Vp-rho empirical rule from Stidham et al. (2001)
!!$! -- used in Komatitsch et al. (2004) simulations
!!$  P1 = 0.33d0 
!!$  P0 = 1280.d0
!!$  rho = P1*vp + P0
!!$!  rho = 0.33d0 * vp + 1280.d0

! Vp-rho empirical rule from Ludwig-Nafe-Drake (1970),
! which is listed in Brocher (2005a)
! -- used in Tape et al. (2009) simulations
   P5 =  1.0600d-16
   P4 = -4.3000d-12
   P3 =  6.7100d-08
   P2 = -4.7210d-04
   P1 =  1.6612d0
   P0 =  0.0d0
   rho = P5*vp**5 + P4*vp**4 + P3*vp**3 + P2*vp**2 + P1*vp + P0

! make sure density estimate is reasonable
  if(rho > DENSITY_MAX) rho = DENSITY_MAX
  if(rho < DENSITY_MIN) rho = DENSITY_MIN

  end subroutine compute_rho_estimate

