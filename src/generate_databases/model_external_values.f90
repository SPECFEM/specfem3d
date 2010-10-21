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

! generic model file
!
! note: the idea is to super-impose velocity model values on the GLL points,
!          additional to the ones assigned on the CUBIT mesh
!
! most of the routines here are place-holders, please add/implement your own routines
!

  module external_model

!---
!
! ADD YOUR MODEL HERE
!
!---
 
  ! only here to illustrate an example
  !  type model_external_variables
  !    sequence
  !    double precision dvs(0:dummy_size)
  !  end type model_external_variables
  !  type (model_external_variables) MEXT_V

  end module external_model

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_external_broadcast(myrank)

! standard routine to setup model 

  use external_model
  
  implicit none

  include "constants.h"
  ! standard include of the MPI library
  include 'mpif.h'

  integer :: myrank
  
  ! local parameters
  integer :: idummy

  ! dummy to ignore compiler warnings
  idummy = myrank

!---
!
! ADD YOUR MODEL HERE
!
!---

  ! the variables read are declared and stored in structure MEXT_V      
  !if(myrank == 0) call read_external_model()
      
  ! broadcast the information read on the master to the nodes
  !call MPI_BCAST(MEXT_V%dvs,size(MEXT_V%dvs),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

  end subroutine model_external_broadcast
  
!
!-------------------------------------------------------------------------------------------------
!

!
!  subroutine read_external_model()
!
!  use external_model
!  
!  implicit none
!
!  include "constants.h"
!---
!
! ADD YOUR MODEL HERE
!
!---
!
!  end subroutine read_external_model


!
!-------------------------------------------------------------------------------------------------
!


  subroutine model_external_values(i,j,k,ispec,idomain_id,imaterial_id,&
                            nspec,ibool, &
                            iflag_aniso,iflag_atten, &
                            rho,vp,vs, &
                            c11,c12,c13,c14,c15,c16, &
                            c22,c23,c24,c25,c26,c33, &
                            c34,c35,c36,c44,c45,c46, &
                            c55,c56,c66,ANISOTROPY)

! given a GLL point, returns super-imposed velocity model values

  use external_model
  use create_regions_mesh_ext_par
  
  implicit none

  ! GLL point indices
  integer :: i,j,k,ispec
  
  ! acoustic/elastic/.. domain flag ( 1 = acoustic / 2 = elastic / ... )
  integer :: idomain_id
  
  ! associated material flag (in cubit, this would be the volume id number)
  integer :: imaterial_id

  ! local-to-global index array
  integer :: nspec
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  
  ! anisotropy flag
  integer :: iflag_aniso
  
  ! attenuation flag
  integer :: iflag_atten
  
  ! density, Vp and Vs
  real(kind=CUSTOM_REAL) :: vp,vs,rho  
  
  ! all anisotropy coefficients
  real(kind=CUSTOM_REAL) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25, &
                        c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66
  logical :: ANISOTROPY

  ! local parameters
  real(kind=CUSTOM_REAL) :: x,y,z
  real(kind=CUSTOM_REAL) :: xmin,xmax,ymin,ymax,zmin,zmax
  real(kind=CUSTOM_REAL) :: depth
  integer :: iglob,idummy

!---
!
! ADD YOUR MODEL HERE
!
!---
  
  ! GLL point location
  iglob = ibool(i,j,k,ispec)
  x = xstore_dummy(iglob)
  y = ystore_dummy(iglob)
  z = zstore_dummy(iglob)

  ! model dimensions
  xmin = 0. ! minval(xstore_dummy)
  xmax = 134000. ! maxval(xstore_dummy)
  ymin = 0.  !minval(ystore_dummy)
  ymax = 134000. ! maxval(ystore_dummy)
  zmin = 0. ! minval(zstore_dummy)
  zmax = -60000. ! maxval(zstore_dummy)

  ! depth in Z-direction
  depth = zmax - z
  
  ! normalizes depth between 0 and 1
  if( abs( zmax - zmin ) > TINYVAL ) depth = depth / (zmax - zmin)


  ! super-imposes values
  !rho = 2.6910d0+0.6924d0*depth
  !vp = 4.1875d0+3.9382d0*depth
  !vs = 2.1519d0+2.3481d0*depth

  ! adds a velocity depth gradient 
  ! (e.g. from PREM mantle gradients: 
  !     vp : 3.9382*6371/5.5 
  !     vs : 2.3481*6371/5.5 
  !     rho : 0.6924*6371/5.5 )
  rho = rho + 802.d0 * depth
  vp = vp + 4562.d0 * depth
  vs = vs + 2720.d0 * depth  
  
  ! adds anisotropic velocity values
  if( ANISOTROPY ) &
    call model_aniso(iflag_aniso,rho,vp,vs,c11,c12,c13,c14,c15,c16, &
                     c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45, &
                     c46,c55,c56,c66) 

  ! to avoid compiler warnings
  idummy = imaterial_id
  idummy = idomain_id
  idummy = iflag_atten
      
  end subroutine model_external_values

  