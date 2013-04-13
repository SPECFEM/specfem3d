!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
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

!--------------------------------------------------------------------------------------------------
!
! 1D model profile for Cascadia region
!
! by Gian Matharu
!--------------------------------------------------------------------------------------------------

  subroutine model_1D_cascadia(xmesh,ymesh,zmesh,rho,vp,vs,qmu_atten)

! given a GLL point, returns super-imposed velocity model values

  use generate_databases_par,only: nspec => NSPEC_AB,ibool
  use create_regions_mesh_ext_par
  implicit none

  ! GLL point location
  double precision, intent(in) :: xmesh,ymesh,zmesh

  ! density, Vp and Vs
  real(kind=CUSTOM_REAL),intent(inout) :: vp,vs,rho,qmu_atten

  ! local parameters
  real(kind=CUSTOM_REAL) :: x,y,z
  real(kind=CUSTOM_REAL) :: depth
  real(kind=CUSTOM_REAL) :: elevation,distmin

  ! converts GLL point location to real
  x = xmesh
  y = ymesh
  z = zmesh

  ! get approximate topography elevation at target coordinates
  distmin = HUGEVAL
  elevation = 0.0
  call get_topo_elevation_free_closest(x,y,elevation,distmin, &
                    nspec,nglob_dummy,ibool,xstore_dummy,ystore_dummy,zstore_dummy, &
                    num_free_surface_faces,free_surface_ispec,free_surface_ijk)

  ! depth in Z-direction
  if( distmin < HUGEVAL ) then
    depth = elevation - z
  else
    depth = - z
  endif

  ! depth in km
  depth = depth / 1000.0

  ! 1D profile Cascadia

  ! super-imposes values
  if( depth < 1.0 ) then
    ! vp in m/s
    vp = 5000.0
    ! vs in m/s
    vs = 2890.0
    ! density in kg/m**3
    rho = 2800.0
  else if( depth < 6.0 ) then
    vp = 6000.0
    vs = 3460.0
    rho = 2800.0
  else if( depth < 30.0 ) then
    vp = 6700.0
    vs = 3870.0
    rho = 3200.0
  else if( depth < 45.0 ) then
    vp = 7100.0
    vs = 4100.0
    rho = 3200.0
  else if( depth < 65.0 ) then
    vp = 7750.0
    vs = 4470.0
    rho = 3200.0
  else
    vp = 8100.0
    vs = 4670.0
    rho = 3200.0
  endif

  ! attenuation: PREM crust value
  qmu_atten=600.0d0

  end subroutine model_1D_cascadia
