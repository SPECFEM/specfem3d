!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
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
! generic model file
!
! note: the idea is to super-impose velocity model values on the GLL points,
!          additional to the ones assigned on the CUBIT mesh
!
! most of the routines here are place-holders, please add/implement your own routines
!
!--------------------------------------------------------------------------------------------------


  module external_model

!---
!
! ADD YOUR MODEL HERE
!
!---
  integer, parameter :: dummy_size = 1

  ! only here to illustrate an example
  type model_external_variables
      sequence
      double precision :: dvs(dummy_size)
    end type model_external_variables
  type (model_external_variables) MEXT_V

  end module external_model

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_external_broadcast()

! standard routine to setup model

  use external_model

  use constants

  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE,MESH_A_CHUNK_OF_THE_EARTH

  implicit none

  ! local parameters
  integer :: idummy

  ! dummy to ignore compiler warnings
  idummy = myrank

  ! safety check
  if (COUPLE_WITH_INJECTION_TECHNIQUE .or. MESH_A_CHUNK_OF_THE_EARTH) then
    print *,'Coupling with injection technique or mesh a chunk of the earth requires in Par_file: MODEL = coupled'
    stop 'Error model external'
  endif


!---
!
! ADD YOUR MODEL HERE
!
!---

  ! the variables read are declared and stored in structure MEXT_V
  if (myrank == 0) call read_external_model()

  ! broadcast the information read on the main to the nodes
  call bcast_all_dp(MEXT_V%dvs, size(MEXT_V%dvs))

  end subroutine model_external_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_external_model()

  use external_model

  use constants

  implicit none

!---
!
! ADD YOUR MODEL HERE
!
!---
  MEXT_V%dvs(:) = 0.d0

  end subroutine read_external_model


!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_external_values(xmesh,ymesh,zmesh,ispec,rho,vp,vs,qkappa_atten,qmu_atten,iflag_aniso,idomain_id )

! given a GLL point, returns super-imposed velocity model values

  !use constants, only: MIDX,MIDY,MIDZ

  use generate_databases_par, only: nspec => NSPEC_AB,ibool,HUGEVAL,TINYVAL,IDOMAIN_ELASTIC,CUSTOM_REAL

  use create_regions_mesh_ext_par, only: xstore_unique,ystore_unique,zstore_unique, &
    num_free_surface_faces,free_surface_ijk,free_surface_ispec,nglob_unique

  implicit none

  ! GLL point
  double precision, intent(in) :: xmesh,ymesh,zmesh
  integer, intent(in) :: ispec

  ! density, Vp and Vs
  real(kind=CUSTOM_REAL) :: vp,vs,rho

  ! attenuation flag
  real(kind=CUSTOM_REAL) :: qkappa_atten,qmu_atten

  ! anisotropy flag
  integer :: iflag_aniso

  ! acoustic/elastic/.. domain flag ( 1 = acoustic / 2 = elastic / ... )
  integer :: idomain_id

  ! local parameters
  double precision :: x,y,z
  real(kind=CUSTOM_REAL) :: xmin,xmax,ymin,ymax,zmin,zmax,x_target,y_target
  real(kind=CUSTOM_REAL) :: depth
  real(kind=CUSTOM_REAL) :: elevation,distmin
  integer :: idummy
  !integer :: iglob
  !real(kind=CUSTOM_REAL) :: mid_x,mid_y,mid_z

!---
!
! ADD YOUR MODEL HERE
!
!---

  ! GLL point location converted to real
  x = xmesh
  y = ymesh
  z = zmesh


  ! gets element midpoint location (if needed)
  idummy = ispec
  !iglob = ibool(MIDX,MIDY,MIDZ,ispec)
  !mid_x = xstore_unique(iglob)
  !mid_y = ystore_unique(iglob)
  !mid_z = zstore_unique(iglob)

  ! note: z coordinate will be negative below surface
  !          convention is z-axis points up

  ! model dimensions
  xmin = 0._CUSTOM_REAL ! minval(xstore_unique)
  xmax = 134000._CUSTOM_REAL ! maxval(xstore_unique)
  ymin = 0._CUSTOM_REAL  !minval(ystore_unique)
  ymax = 134000._CUSTOM_REAL ! maxval(ystore_unique)
  zmin = 0._CUSTOM_REAL ! minval(zstore_unique)
  zmax = 60000._CUSTOM_REAL ! maxval(zstore_unique)
  x_target = x
  y_target = y

  ! get approximate topography elevation at target coordinates from free surface
  call get_topo_elevation_free_closest(x_target,y_target,elevation,distmin, &
                                       nspec,nglob_unique,ibool,xstore_unique,ystore_unique,zstore_unique, &
                                       num_free_surface_faces,free_surface_ispec,free_surface_ijk)

  ! depth in Z-direction
  if (distmin < HUGEVAL) then
    depth = elevation - z
  else
    depth = zmin - z
  endif

  ! normalizes depth between 0 and 1
  if (abs( zmax - zmin ) > TINYVAL) depth = depth / (zmax - zmin)

  ! initial values (in m/s and kg/m^3)
  rho = 2691.0_CUSTOM_REAL
  vp = 4187.5_CUSTOM_REAL
  vs = 2151.9_CUSTOM_REAL

  ! adds a velocity depth gradient
  ! (e.g. from PREM mantle gradients:
  !     vp : 3.9382*6371/5.5
  !     vs : 2.3481*6371/5.5
  !     rho : 0.6924*6371/5.5 )
  rho = rho + 802._CUSTOM_REAL * depth
  vp = vp + 4562._CUSTOM_REAL * depth
  vs = vs + 2720._CUSTOM_REAL * depth

  ! attenuation: PREM crust value
  qmu_atten=600._CUSTOM_REAL

  ! no Q_Kappa in this model, use a dummy very high value of 9999. as a flag for no QKappa attenuation
  qkappa_atten = 9999._CUSTOM_REAL

  ! no anisotropy
  iflag_aniso = 0

  ! elastic material
  idomain_id = IDOMAIN_ELASTIC

  end subroutine model_external_values

