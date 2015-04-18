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

  subroutine xyz_2_rthetaphi(x,y,z,r,theta,phi)

! convert x y z to r theta phi, single precision call

  implicit none

  include "constants.h"

  real(kind=CUSTOM_REAL) x,y,z,r,theta,phi
  double precision xmesh,ymesh,zmesh

! distinguish between single and double precision for reals
  if(CUSTOM_REAL == SIZE_REAL) then

    xmesh = dble(x)
    ymesh = dble(y)
    zmesh = dble(z)

    if(zmesh > -SMALL_VAL_ANGLE .and. zmesh <= ZERO) zmesh = -SMALL_VAL_ANGLE
    if(zmesh < SMALL_VAL_ANGLE .and. zmesh >= ZERO) zmesh = SMALL_VAL_ANGLE
    theta = sngl(datan2(dsqrt(xmesh*xmesh+ymesh*ymesh),zmesh))
    if(xmesh > -SMALL_VAL_ANGLE .and. xmesh <= ZERO) xmesh = -SMALL_VAL_ANGLE
    if(xmesh < SMALL_VAL_ANGLE .and. xmesh >= ZERO) xmesh = SMALL_VAL_ANGLE
    phi = sngl(datan2(ymesh,xmesh))

    r = sngl(dsqrt(xmesh**2 + ymesh**2 + zmesh**2))

  else

    xmesh = x
    ymesh = y
    zmesh = z

    if(zmesh > -SMALL_VAL_ANGLE .and. zmesh <= ZERO) zmesh = -SMALL_VAL_ANGLE
    if(zmesh < SMALL_VAL_ANGLE .and. zmesh >= ZERO) zmesh = SMALL_VAL_ANGLE
    theta = datan2(dsqrt(xmesh*xmesh+ymesh*ymesh),zmesh)
    if(xmesh > -SMALL_VAL_ANGLE .and. xmesh <= ZERO) xmesh = -SMALL_VAL_ANGLE
    if(xmesh < SMALL_VAL_ANGLE .and. xmesh >= ZERO) xmesh = SMALL_VAL_ANGLE
    phi = datan2(ymesh,xmesh)

    r = dsqrt(xmesh**2 + ymesh**2 + zmesh**2)

  endif

  end subroutine xyz_2_rthetaphi

!-------------------------------------------------------------

  subroutine xyz_2_rthetaphi_dble(x,y,z,r,theta,phi)

! convert x y z to r theta phi, double precision call

  implicit none

  include "constants.h"

  double precision x,y,z,r,theta,phi
  double precision xmesh,ymesh,zmesh

  xmesh = x
  ymesh = y
  zmesh = z

  if(zmesh > -SMALL_VAL_ANGLE .and. zmesh <= ZERO) zmesh = -SMALL_VAL_ANGLE
  if(zmesh < SMALL_VAL_ANGLE .and. zmesh >= ZERO) zmesh = SMALL_VAL_ANGLE
  theta = datan2(dsqrt(xmesh*xmesh+ymesh*ymesh),zmesh)
  if(xmesh > -SMALL_VAL_ANGLE .and. xmesh <= ZERO) xmesh = -SMALL_VAL_ANGLE
  if(xmesh < SMALL_VAL_ANGLE .and. xmesh >= ZERO) xmesh = SMALL_VAL_ANGLE
  phi = datan2(ymesh,xmesh)

  r = dsqrt(xmesh**2 + ymesh**2 + zmesh**2)

  end subroutine xyz_2_rthetaphi_dble

!-------------------------------------------------------------

  subroutine rthetaphi_2_xyz(x,y,z,r,theta,phi)

! convert r theta phi to x y z

  implicit none

  include "constants.h"

  real(kind=CUSTOM_REAL) x,y,z,r,theta,phi

  x = r * sin(theta) * cos(phi)
  y = r * sin(theta) * sin(phi)
  z = r * cos(theta)

  end subroutine rthetaphi_2_xyz

