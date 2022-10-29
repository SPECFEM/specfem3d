!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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

!--------------------------------------------------------------------------------------------------
!
! 1D Southern California model
!
! model is the standard model used in southern California:
!   Kanamori and Hadley (1975), Dreger and Helmberger (1990), Wald-Hutton,Given (1995)
!
!--------------------------------------------------------------------------------------------------

  subroutine model_1D_socal(xmesh,ymesh,zmesh,rho,vp,vs,qmu_atten,qkappa_atten)

! given a GLL point, returns super-imposed velocity model values

  use create_regions_mesh_ext_par
  implicit none

  ! GLL point location
  double precision, intent(in) :: xmesh,ymesh,zmesh

  ! density, Vp and Vs
  real(kind=CUSTOM_REAL),intent(inout) :: vp,vs,rho,qmu_atten,qkappa_atten

  ! local parameters
  real(kind=CUSTOM_REAL) :: depth,x,y,z

  ! mesh point location
  x = xmesh
  y = ymesh
  z = zmesh

  ! depth in m
  depth = -zmesh

  ! assigns model parameters
  if (depth >= 32000.0) then
    ! moho
    vp=7.8_CUSTOM_REAL
    vs=4.5_CUSTOM_REAL
    rho=3.0_CUSTOM_REAL
  else if (depth > 16000.0) then
    ! moho - 16km
    vp=6.7_CUSTOM_REAL
    vs=3.87_CUSTOM_REAL
    rho=2.8_CUSTOM_REAL
  else if (depth > 5500.0) then
    ! basement
    vp=6.3_CUSTOM_REAL
    vs=3.64_CUSTOM_REAL
    rho=2.67_CUSTOM_REAL
  else
    ! up to topo surface
    vp=5.5_CUSTOM_REAL
    vs=3.18_CUSTOM_REAL
    rho=2.4_CUSTOM_REAL
  endif

  ! scale to standard units
  vp = vp * 1000._CUSTOM_REAL
  vs = vs * 1000._CUSTOM_REAL
  rho = rho * 1000._CUSTOM_REAL

  ! no attenuation information
  ! using PREM crustal values instead
  qmu_atten = 600.0
  qkappa_atten = 57827.0

  end subroutine model_1D_socal
