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

! compute rho estimate in Gocad block and in Hauksson's model
! based upon Vp

  implicit none

!  include "constants.h"
  include "constants_gocad.h"

  double precision rho,vp

! scale density - use empirical rule from Christiane
  rho = 0.33d0 * vp + 1280.d0

! make sure density estimate is reasonable
  if(rho > DENSITY_MAX) rho = DENSITY_MAX
  if(rho < DENSITY_MIN) rho = DENSITY_MIN

  end subroutine compute_rho_estimate

