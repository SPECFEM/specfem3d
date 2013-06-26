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
!
! United States and French Government Sponsorship Acknowledged.

subroutine pml_compute_accel_contribution_elastic(ispec,ispec_CPML,displ,veloc,rmemory_displ_elastic)

  ! calculates contribution from each C-PML element to update acceleration to the global mesh

  ! second-order accurate convolution term calculation from equation (21) of
  ! Shumin Wang, Robert Lee, and Fernando L. Teixeira,
  ! Anisotropic-Medium PML for Vector FETD With Modified Basis Functions,
  ! IEEE Transactions on Antennas and Propagation, vol. 54, no. 1, (2006)

  use specfem_par, only: NGLOB_AB,it,deltat,wgll_cube,jacobian,ibool,rhostore
  use pml_par, only: CPML_regions,d_store_x,d_store_y,d_store_z,K_store_x,K_store_y,K_store_z,&
                     alpha_store,NSPEC_CPML,accel_elastic_CPML                 
  use constants, only: CUSTOM_REAL,NDIM,NGLLX,NGLLY,NGLLZ,CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY, &
                       CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ

  implicit none

  integer, intent(in) :: ispec,ispec_CPML
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB), intent(in) :: displ,veloc
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ,NSPEC_CPML,3) :: rmemory_displ_elastic

  ! local parameters
  integer :: i,j,k,iglob
  real(kind=CUSTOM_REAL) :: wgllcube,rhol,jacobianl
  real(kind=CUSTOM_REAL) :: bb,coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2,coef0_3,coef1_3,coef2_3
  real(kind=CUSTOM_REAL) :: A0,A1,A2,A3,A4,A5,temp_A3! for convolution of acceleration

  do k=1,NGLLZ
     do j=1,NGLLY
        do i=1,NGLLX
           rhol = rhostore(i,j,k,ispec)
           jacobianl = jacobian(i,j,k,ispec)
           iglob = ibool(i,j,k,ispec)
           wgllcube = wgll_cube(i,j,k)

           if( CPML_regions(ispec_CPML) == CPML_X_ONLY ) then
              !------------------------------------------------------------------------------
              !---------------------------- X-surface C-PML ---------------------------------
              !------------------------------------------------------------------------------

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.0d0 - exp(-bb * deltat/2.0d0) ) / bb
                 coef2_1 = (1.0d0 - exp(-bb * deltat/2.0d0) ) * exp(-bb * deltat/2.0d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) = coef0_1 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) &
                   + (displ(1,iglob) + deltat * veloc(1,iglob)) * coef1_1 + displ(1,iglob) * coef2_1
              rmemory_displ_elastic(1,i,j,k,ispec_CPML,2) = 0.0
              rmemory_displ_elastic(1,i,j,k,ispec_CPML,3) = 0.0

              rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) = coef0_1 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) &
                   + (displ(2,iglob) + deltat * veloc(2,iglob)) * coef1_1 + displ(2,iglob) * coef2_1
              rmemory_displ_elastic(2,i,j,k,ispec_CPML,2) = 0.0
              rmemory_displ_elastic(2,i,j,k,ispec_CPML,3) = 0.0

              rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) = coef0_1 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) &
                   + (displ(3,iglob) + deltat * veloc(3,iglob)) * coef1_1 + displ(3,iglob) * coef2_1
              rmemory_displ_elastic(3,i,j,k,ispec_CPML,2) = 0.0
              rmemory_displ_elastic(3,i,j,k,ispec_CPML,3) = 0.0

              !---------------------- A0, A1, A2, A3, A4 and A5 --------------------------
              A0 = k_store_x(i,j,k,ispec_CPML)
              A1 = d_store_x(i,j,k,ispec_CPML)
              A2 = - alpha_store(i,j,k,ispec_CPML) * d_store_x(i,j,k,ispec_CPML)
              A3 = d_store_x(i,j,k,ispec_CPML) * alpha_store(i,j,k,ispec_CPML) ** 2
              A4 = 0.d0
              A5 = 0.d0

           else if( CPML_regions(ispec_CPML) == CPML_Y_ONLY ) then
              !------------------------------------------------------------------------------
              !---------------------------- Y-surface C-PML ---------------------------------
              !------------------------------------------------------------------------------

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.0d0 - exp(-bb * deltat/2.0d0) ) / bb
                 coef2_1 = (1.0d0 - exp(-bb * deltat/2.0d0) ) * exp(-bb * deltat/2.0d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) = coef0_1 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) &
                   + (displ(1,iglob) + deltat * veloc(1,iglob)) * coef1_1 + displ(1,iglob) * coef2_1
              rmemory_displ_elastic(1,i,j,k,ispec_CPML,2) = 0.0
              rmemory_displ_elastic(1,i,j,k,ispec_CPML,3) = 0.0

              rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) = coef0_1 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) &
                   + (displ(2,iglob) + deltat * veloc(2,iglob)) * coef1_1 + displ(2,iglob) * coef2_1
              rmemory_displ_elastic(2,i,j,k,ispec_CPML,2) = 0.0
              rmemory_displ_elastic(2,i,j,k,ispec_CPML,3) = 0.0

              rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) = coef0_1 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) &
                   + (displ(3,iglob) + deltat * veloc(3,iglob)) * coef1_1 + displ(3,iglob) * coef2_1
              rmemory_displ_elastic(3,i,j,k,ispec_CPML,2) = 0.0
              rmemory_displ_elastic(3,i,j,k,ispec_CPML,3) = 0.0

              !---------------------- A0, A1, A2, A3, A4 and A5 --------------------------
              A0 = k_store_y(i,j,k,ispec_CPML)
              A1 = d_store_y(i,j,k,ispec_CPML)
              A2 = - alpha_store(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML)
              A3 = d_store_y(i,j,k,ispec_CPML) * alpha_store(i,j,k,ispec_CPML) ** 2
              A4 = 0.d0
              A5 = 0.d0

           else if( CPML_regions(ispec_CPML) == CPML_Z_ONLY ) then
              !------------------------------------------------------------------------------
              !---------------------------- Z-surface C-PML ---------------------------------
              !------------------------------------------------------------------------------

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.0d0 - exp(-bb * deltat/2.0d0) ) / bb
                 coef2_1 = (1.0d0 - exp(-bb * deltat/2.0d0) ) * exp(-bb * deltat/2.0d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) = coef0_1 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) &
                   + (displ(1,iglob) + deltat * veloc(1,iglob)) * coef1_1 + displ(1,iglob) * coef2_1
              rmemory_displ_elastic(1,i,j,k,ispec_CPML,2) = 0.0
              rmemory_displ_elastic(1,i,j,k,ispec_CPML,3) = 0.0

              rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) = coef0_1 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) &
                   + (displ(2,iglob) + deltat * veloc(2,iglob)) * coef1_1 + displ(2,iglob) * coef2_1
              rmemory_displ_elastic(2,i,j,k,ispec_CPML,2) = 0.0
              rmemory_displ_elastic(2,i,j,k,ispec_CPML,3) = 0.0

              rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) = coef0_1 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) &
                   + (displ(3,iglob) + deltat * veloc(3,iglob)) * coef1_1 + displ(3,iglob) * coef2_1
              rmemory_displ_elastic(3,i,j,k,ispec_CPML,2) = 0.0
              rmemory_displ_elastic(3,i,j,k,ispec_CPML,3) = 0.0

              !---------------------- A0, A1, A2, A3, A4 and A5 --------------------------
              A0 = k_store_z(i,j,k,ispec_CPML)
              A1 = d_store_z(i,j,k,ispec_CPML)
              A2 = - alpha_store(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)
              A3 = d_store_z(i,j,k,ispec_CPML) * alpha_store(i,j,k,ispec_CPML) ** 2
              A4 = 0.d0
              A5 = 0.d0

           else if( CPML_regions(ispec_CPML) == CPML_XY_ONLY ) then
              !------------------------------------------------------------------------------
              !---------------------------- XY-edge C-PML -----------------------------------
              !------------------------------------------------------------------------------

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.0d0 - exp(-bb * deltat/2.0d0) ) / bb
                 coef2_1 = (1.0d0 - exp(-bb * deltat/2.0d0) ) * exp(-bb * deltat/2.0d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              coef0_2 = coef0_1
              coef1_2 = coef1_1
              coef2_2 = coef2_1

              rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) = coef0_1 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) &
                   + (displ(1,iglob) + deltat * veloc(1,iglob)) * coef1_1 + displ(1,iglob) * coef2_1
              rmemory_displ_elastic(1,i,j,k,ispec_CPML,2) = coef0_2 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,2) &
                   + (displ(1,iglob) + deltat * veloc(1,iglob)) * it*deltat * coef1_2 &
                   + displ(1,iglob) * it*deltat * coef2_2
              rmemory_displ_elastic(1,i,j,k,ispec_CPML,3) = 0.0

              rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) = coef0_1 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) &
                   + (displ(2,iglob) + deltat * veloc(2,iglob)) * coef1_1 + displ(2,iglob) * coef2_1
              rmemory_displ_elastic(2,i,j,k,ispec_CPML,2) = coef0_2 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,2) &
                   + (displ(2,iglob) + deltat * veloc(2,iglob)) * it*deltat * coef1_2 &
                   + displ(2,iglob) * it*deltat * coef2_2
              rmemory_displ_elastic(2,i,j,k,ispec_CPML,3) = 0.0

              rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) = coef0_1 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) &
                   + (displ(3,iglob) + deltat * veloc(3,iglob)) * coef1_1 + displ(3,iglob) * coef2_1
              rmemory_displ_elastic(3,i,j,k,ispec_CPML,2) = coef0_2 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,2) &
                   + (displ(3,iglob) + deltat * veloc(3,iglob)) * it*deltat * coef1_2 &
                   + displ(3,iglob) * it*deltat * coef2_2
              rmemory_displ_elastic(3,i,j,k,ispec_CPML,3) = 0.0

              !---------------------- A0, A1, A2, A3, A4 and A5 --------------------------
              A0 = k_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML)
              A1 = d_store_y(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) &
                   + d_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML)
              A2 = d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) &
                   - alpha_store(i,j,k,ispec_CPML) * d_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) &
                   - alpha_store(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML)
              A3 = -2.0 * alpha_store(i,j,k,ispec_CPML) * d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) &
                   + alpha_store(i,j,k,ispec_CPML)**2 * ( d_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) &
                   + d_store_y(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) ) + alpha_store(i,j,k,ispec_CPML)**2 &
                   * it*deltat * d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML)
              A4 = -alpha_store(i,j,k,ispec_CPML)**2 * d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML)
              A5 = 0.0

           else if( CPML_regions(ispec_CPML) == CPML_XZ_ONLY ) then
              !------------------------------------------------------------------------------
              !---------------------------- XZ-edge C-PML -----------------------------------
              !------------------------------------------------------------------------------

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.0d0 - exp(-bb * deltat/2.0d0) ) / bb
                 coef2_1 = (1.0d0 - exp(-bb * deltat/2.0d0) ) * exp(-bb * deltat/2.0d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              coef0_2 = coef0_1
              coef1_2 = coef1_1
              coef2_2 = coef2_1

              rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) = coef0_1 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) &
                   + (displ(1,iglob) + deltat * veloc(1,iglob)) * coef1_1 + displ(1,iglob) * coef2_1
              rmemory_displ_elastic(1,i,j,k,ispec_CPML,2) = coef0_2 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,2) &
                   + (displ(1,iglob) + deltat * veloc(1,iglob)) * it*deltat * coef1_2 &
                   + displ(1,iglob) * it*deltat * coef2_2
              rmemory_displ_elastic(1,i,j,k,ispec_CPML,3)= 0.0

              rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) = coef0_1 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) &
                   + (displ(2,iglob) + deltat * veloc(2,iglob)) * coef1_1 + displ(2,iglob) * coef2_1
              rmemory_displ_elastic(2,i,j,k,ispec_CPML,2) = coef0_2 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,2) &
                   + (displ(2,iglob) + deltat * veloc(2,iglob)) * it*deltat * coef1_2 &
                   + displ(2,iglob) * it*deltat * coef2_2
              rmemory_displ_elastic(2,i,j,k,ispec_CPML,3)= 0.0

              rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) = coef0_1 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) &
                   + (displ(3,iglob) + deltat * veloc(3,iglob)) * coef1_1 + displ(3,iglob) * coef2_1
              rmemory_displ_elastic(3,i,j,k,ispec_CPML,2) = coef0_2 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,2) &
                   + (displ(3,iglob) + deltat * veloc(3,iglob)) * it*deltat * coef1_2 &
                   + displ(3,iglob) * it*deltat * coef2_2
              rmemory_displ_elastic(3,i,j,k,ispec_CPML,3)= 0.0

              !---------------------- A0, A1, A2, A3, A4 and A5 --------------------------
              A0 = k_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML)
              A1 = d_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML)&
                   + d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML)
              A2 = d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) &
                   - alpha_store(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) &
                   - alpha_store(i,j,k,ispec_CPML) * d_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML)
              A3 = -2.0 * alpha_store(i,j,k,ispec_CPML) * d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) &
                   + alpha_store(i,j,k,ispec_CPML)**2 * ( d_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) &
                   + d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) ) + alpha_store(i,j,k,ispec_CPML)**2 &
                   * it*deltat * d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)
              A4 = -alpha_store(i,j,k,ispec_CPML)**2 * d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)
              A5 = 0.0

           else if( CPML_regions(ispec_CPML) == CPML_YZ_ONLY ) then
              !------------------------------------------------------------------------------
              !---------------------------- YZ-edge C-PML -----------------------------------
              !------------------------------------------------------------------------------

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.0d0 - exp(-bb * deltat/2.0d0) ) / bb
                 coef2_1 = (1.0d0 - exp(-bb * deltat/2.0d0) ) * exp(-bb * deltat/2.0d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              coef0_2 = coef0_1
              coef1_2 = coef1_1
              coef2_2 = coef2_1

              rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) = coef0_1 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) &
                   + (displ(1,iglob) + deltat * veloc(1,iglob)) * coef1_1 + displ(1,iglob) * coef2_1
              rmemory_displ_elastic(1,i,j,k,ispec_CPML,2) = coef0_2 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,2) &
                   + (displ(1,iglob) + deltat * veloc(1,iglob)) * it*deltat * coef1_2 &
                   + displ(1,iglob) * it*deltat * coef2_2
              rmemory_displ_elastic(1,i,j,k,ispec_CPML,3)=0.d0

              rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) = coef0_1 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) &
                   + (displ(2,iglob) + deltat * veloc(2,iglob)) * coef1_1 + displ(2,iglob) * coef2_1
              rmemory_displ_elastic(2,i,j,k,ispec_CPML,2) = coef0_2 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,2) &
                   + (displ(2,iglob) + deltat * veloc(2,iglob)) * it*deltat * coef1_2 &
                   + displ(2,iglob) * it*deltat * coef2_2
              rmemory_displ_elastic(2,i,j,k,ispec_CPML,3)=0.d0

              rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) = coef0_1 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) &
                   + (displ(3,iglob) + deltat * veloc(3,iglob)) * coef1_1 + displ(3,iglob) * coef2_1
              rmemory_displ_elastic(3,i,j,k,ispec_CPML,2) = coef0_2 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,2) &
                   + (displ(3,iglob) + deltat * veloc(3,iglob)) * it*deltat * coef1_2 &
                   + displ(3,iglob) * it*deltat * coef2_2
              rmemory_displ_elastic(3,i,j,k,ispec_CPML,3)=0.d0

              !---------------------- A0, A1, A2, A3, A4 and A5 --------------------------
              A0 = k_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML)
              A1 = d_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) &
                   + d_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML)
              A2 = d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) &
                   - alpha_store(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) &
                   - alpha_store(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML)
              A3 = -2.0 * alpha_store(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) &
                   + alpha_store(i,j,k,ispec_CPML)**2 * ( d_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) &
                   + d_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) ) + alpha_store(i,j,k,ispec_CPML)**2 &
                   * it*deltat * d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)
              A4 = -alpha_store(i,j,k,ispec_CPML)**2 * d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)
              A5 = 0.0

           else if( CPML_regions(ispec_CPML) == CPML_XYZ ) then
              !------------------------------------------------------------------------------
              !---------------------------- XYZ-corner C-PML --------------------------------
              !------------------------------------------------------------------------------

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.0d0 - exp(-bb * deltat/2.0d0) ) / bb
                 coef2_1 = (1.0d0 - exp(-bb * deltat/2.0d0) ) * exp(-bb * deltat/2.0d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              coef0_2 = coef0_1
              coef1_2 = coef1_1
              coef2_2 = coef2_1

              coef0_3 = coef0_1
              coef1_3 = coef1_1
              coef2_3 = coef2_1

              rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) = coef0_1 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) &
                   + (displ(1,iglob) + deltat * veloc(1,iglob)) * coef1_1 + displ(1,iglob) * coef2_1
              rmemory_displ_elastic(1,i,j,k,ispec_CPML,2) = coef0_2 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,2) &
                   + (displ(1,iglob) + deltat * veloc(1,iglob)) * it*deltat * coef1_2 &
                   + displ(1,iglob) * it*deltat * coef2_2
              rmemory_displ_elastic(1,i,j,k,ispec_CPML,3) = coef0_3 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,3) &
                   + (displ(1,iglob) + deltat * veloc(1,iglob)) * (it*deltat)**2 * coef1_3 &
                   + displ(1,iglob) * (it*deltat)**2 * coef2_3

              rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) = coef0_1 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) &
                   + (displ(2,iglob) + deltat * veloc(2,iglob)) * coef1_1 + displ(2,iglob) * coef2_1
              rmemory_displ_elastic(2,i,j,k,ispec_CPML,2) = coef0_2 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,2) &
                   + (displ(2,iglob) + deltat * veloc(2,iglob)) * it*deltat * coef1_2 &
                   + displ(2,iglob) * it*deltat * coef2_2
              rmemory_displ_elastic(2,i,j,k,ispec_CPML,3) = coef0_3 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,3) &
                   + (displ(2,iglob) + deltat * veloc(2,iglob)) * (it*deltat)**2 * coef1_3 &
                   + displ(2,iglob) * (it*deltat)**2 * coef2_3

              rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) = coef0_1 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) &
                   + (displ(3,iglob) + deltat * veloc(3,iglob)) * coef1_1 + displ(3,iglob) * coef2_1
              rmemory_displ_elastic(3,i,j,k,ispec_CPML,2) = coef0_2 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,2) &
                   + (displ(3,iglob) + deltat * veloc(3,iglob)) * it*deltat * coef1_2 &
                   + displ(3,iglob) * it*deltat * coef2_2
              rmemory_displ_elastic(3,i,j,k,ispec_CPML,3) = coef0_3 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,3) &
                   + (displ(3,iglob) + deltat * veloc(3,iglob)) * (it*deltat)**2 * coef1_3 &
                   + displ(3,iglob) * (it*deltat)**2 * coef2_3

              !---------------------- A0, A1, A2, A3, A4 and A5 --------------------------
              A0 = k_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML)
              A1 = d_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) + &
                   d_store_y(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) + &
                   d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML)
              A2 = k_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) + &
                   k_store_y(i,j,k,ispec_CPML) * d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) + &
                   k_store_z(i,j,k,ispec_CPML) * d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) - &
                   d_store_x(i,j,k,ispec_CPML) * alpha_store(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) * &
                   k_store_z(i,j,k,ispec_CPML) - d_store_y(i,j,k,ispec_CPML) * alpha_store(i,j,k,ispec_CPML) * &
                   k_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) - d_store_z(i,j,k,ispec_CPML) * &
                   alpha_store(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML)
              temp_A3 = d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) - &
                   2.0 * alpha_store(i,j,k,ispec_CPML) * ( &
                   d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) + &
                   d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) + &
                   d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) &
                   ) + alpha_store(i,j,k,ispec_CPML)**2 * ( &
                   d_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) + &
                   d_store_y(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) + &
                   d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) &
                   )
!             temp_A4 = -2.0 * alpha_store(i,j,k,ispec_CPML) * d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) * &
!                  d_store_z(i,j,k,ispec_CPML) + alpha_store(i,j,k,ispec_CPML)**2 * ( &
!                  d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) + &
!                  d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) + &
!                  d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) &
!                  )
!             temp_A5 = 0.5 * d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)
!             A3 = temp_A3 + (it+0.0)*deltat*temp_A4 + ((it+0.0)*deltat)**2*temp_A5
!             A4 = -temp_A4 -2.0*(it+0.0)*deltat*temp_A5
!             A5 = temp_A5
!!! the full experssion of A3,A4,A5 are given by above equation, here we use reduced 
!!! exprssion of A3,A4,A5 in order to stabilized the code. 

              A3 = temp_A3
              A4 = 0.0
              A5 = 0.0

           endif

           accel_elastic_CPML(1,i,j,k) =  wgllcube * rhol * jacobianl * &
                ( A1 * veloc(1,iglob) + A2 * displ(1,iglob) + &
                  A3 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) + &
                  A4 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,2) + &
                  A5 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,3)  &
                )

           accel_elastic_CPML(2,i,j,k) =  wgllcube * rhol * jacobianl * &
                ( A1 * veloc(2,iglob) + A2 * displ(2,iglob) + &
                  A3 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) + &
                  A4 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,2) + &
                  A5 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,3)  &
                )

           accel_elastic_CPML(3,i,j,k) =  wgllcube * rhol * jacobianl * &
                ( A1 * veloc(3,iglob) + A2 * displ(3,iglob) + &
                  A3 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) + &
                  A4 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,2) + &
                  A5 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,3)  &
                )
        enddo
     enddo
  enddo

end subroutine pml_compute_accel_contribution_elastic

!
!=====================================================================
!
! 

subroutine pml_compute_accel_contribution_acoustic(ispec,ispec_CPML,potential_acoustic,&
                                                   potential_dot_acoustic,rmemory_potential_acoustic)

  ! calculates contribution from each C-PML element to update acceleration to the global mesh

  ! second-order accurate convolution term calculation from equation (21) of
  ! Shumin Wang, Robert Lee, and Fernando L. Teixeira,
  ! Anisotropic-Medium PML for Vector FETD With Modified Basis Functions,
  ! IEEE Transactions on Antennas and Propagation, vol. 54, no. 1, (2006)

  use specfem_par, only: NGLOB_AB,it,deltat,wgll_cube,jacobian,ibool,kappastore
  use pml_par, only: CPML_regions,NSPEC_CPML,d_store_x,d_store_y,d_store_z,K_store_x,K_store_y,K_store_z,&
                     alpha_store, potential_dot_dot_acoustic_CPML
  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY, &
                       CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ

  implicit none

  integer, intent(in) :: ispec,ispec_CPML
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB), intent(in) :: potential_acoustic,potential_dot_acoustic 
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CPML,3) :: rmemory_potential_acoustic


  ! local parameters
  integer :: i,j,k,iglob
  real(kind=CUSTOM_REAL) :: wgllcube,kappal_inv,jacobianl
  real(kind=CUSTOM_REAL) :: bb,coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2,coef0_3,coef1_3,coef2_3
  real(kind=CUSTOM_REAL) :: A0,A1,A2,A3,A4,A5,temp_A3 ! for convolution of acceleration

  do k=1,NGLLZ
     do j=1,NGLLY
        do i=1,NGLLX
           kappal_inv = 1.d0 / kappastore(i,j,k,ispec)
           jacobianl = jacobian(i,j,k,ispec)
           iglob = ibool(i,j,k,ispec)
           wgllcube = wgll_cube(i,j,k)

           if( CPML_regions(ispec_CPML) == CPML_X_ONLY ) then
              !------------------------------------------------------------------------------
              !---------------------------- X-surface C-PML ---------------------------------
              !------------------------------------------------------------------------------

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.0d0 - exp(-bb * deltat/2.0d0) ) / bb
                 coef2_1 = (1.0d0 - exp(-bb * deltat/2.0d0) ) * exp(-bb * deltat/2.0d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              rmemory_potential_acoustic(i,j,k,ispec_CPML,1)=coef0_1 * rmemory_potential_acoustic(i,j,k,ispec_CPML,1) &
                   + (potential_acoustic(iglob) + deltat*potential_dot_acoustic(iglob)) * coef1_1 &
                   + potential_acoustic(iglob) * coef2_1
              rmemory_potential_acoustic(i,j,k,ispec_CPML,2)=0.0
              rmemory_potential_acoustic(i,j,k,ispec_CPML,3)=0.0

              !---------------------- A0, A1, A2, A3, A4 and A5 --------------------------
              A0 = k_store_x(i,j,k,ispec_CPML)
              A1 = d_store_x(i,j,k,ispec_CPML)
              A2 = - alpha_store(i,j,k,ispec_CPML) * d_store_x(i,j,k,ispec_CPML)
              A3 = d_store_x(i,j,k,ispec_CPML) * alpha_store(i,j,k,ispec_CPML) ** 2
              A4 = 0.d0
              A5 = 0.d0

           else if( CPML_regions(ispec_CPML) == CPML_Y_ONLY ) then
              !------------------------------------------------------------------------------
              !---------------------------- Y-surface C-PML ---------------------------------
              !------------------------------------------------------------------------------

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.0d0 - exp(-bb * deltat/2.0d0) ) / bb
                 coef2_1 = (1.0d0 - exp(-bb * deltat/2.0d0) ) * exp(-bb * deltat/2.0d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              rmemory_potential_acoustic(i,j,k,ispec_CPML,1)=coef0_1 * rmemory_potential_acoustic(i,j,k,ispec_CPML,1) &
                    + (potential_acoustic(iglob) + deltat*potential_dot_acoustic(iglob)) * coef1_1 &
                    + potential_acoustic(iglob) * coef2_1
              rmemory_potential_acoustic(i,j,k,ispec_CPML,2)=0.0
              rmemory_potential_acoustic(i,j,k,ispec_CPML,3)=0.0

              !---------------------- A0, A1, A2, A3, A4 and A5 --------------------------
              A0 = k_store_y(i,j,k,ispec_CPML)
              A1 = d_store_y(i,j,k,ispec_CPML)
              A2 = - alpha_store(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML)
              A3 = d_store_y(i,j,k,ispec_CPML) * alpha_store(i,j,k,ispec_CPML) ** 2
              A4 = 0.d0
              A5 = 0.d0

           else if( CPML_regions(ispec_CPML) == CPML_Z_ONLY ) then
              !------------------------------------------------------------------------------
              !---------------------------- Z-surface C-PML ---------------------------------
              !------------------------------------------------------------------------------

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.0d0 - exp(-bb * deltat/2.0d0) ) / bb
                 coef2_1 = (1.0d0 - exp(-bb * deltat/2.0d0) ) * exp(-bb * deltat/2.0d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              rmemory_potential_acoustic(i,j,k,ispec_CPML,1)=coef0_1 * rmemory_potential_acoustic(i,j,k,ispec_CPML,1) &
                    + (potential_acoustic(iglob) + deltat*potential_dot_acoustic(iglob)) * coef1_1 &
                    + potential_acoustic(iglob) * coef2_1
              rmemory_potential_acoustic(i,j,k,ispec_CPML,2)=0.0
              rmemory_potential_acoustic(i,j,k,ispec_CPML,3)=0.0

              !---------------------- A0, A1, A2, A3, A4 and A5 --------------------------
              A0 = k_store_z(i,j,k,ispec_CPML)
              A1 = d_store_z(i,j,k,ispec_CPML)
              A2 = - alpha_store(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)
              A3 = d_store_z(i,j,k,ispec_CPML) * alpha_store(i,j,k,ispec_CPML) ** 2
              A4 = 0.d0
              A5 = 0.d0

           else if( CPML_regions(ispec_CPML) == CPML_XY_ONLY ) then
              !------------------------------------------------------------------------------
              !---------------------------- XY-edge C-PML -----------------------------------
              !------------------------------------------------------------------------------

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.0d0 - exp(-bb * deltat/2.0d0) ) / bb
                 coef2_1 = (1.0d0 - exp(-bb * deltat/2.0d0) ) * exp(-bb * deltat/2.0d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              coef0_2 = coef0_1
              coef1_2 = coef1_1
              coef2_2 = coef2_1

              rmemory_potential_acoustic(i,j,k,ispec_CPML,1)=coef0_1 * rmemory_potential_acoustic(i,j,k,ispec_CPML,1) &
                    + (potential_acoustic(iglob) + deltat*potential_dot_acoustic(iglob)) * coef1_1 &
                    + potential_acoustic(iglob) * coef2_1
              rmemory_potential_acoustic(i,j,k,ispec_CPML,2)=coef0_2 * rmemory_potential_acoustic(i,j,k,ispec_CPML,2) &
                    + (potential_acoustic(iglob) + deltat*potential_dot_acoustic(iglob)) * it*deltat * coef1_2 &
                    + potential_acoustic(iglob) * it*deltat * coef2_2
              rmemory_potential_acoustic(i,j,k,ispec_CPML,3)=0.0

              !---------------------- A0, A1, A2, A3, A4 and A5 --------------------------
              A0 = k_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML)
              A1 = d_store_y(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) &
                   + d_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML)
              A2 = d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) &
                   - alpha_store(i,j,k,ispec_CPML) * d_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) &
                   - alpha_store(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML)
              A3 = -2.0 * alpha_store(i,j,k,ispec_CPML) * d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) + &
                   alpha_store(i,j,k,ispec_CPML)**2 * ( d_store_x(i,j,k,ispec_CPML)*k_store_y(i,j,k,ispec_CPML) &
                   + d_store_y(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) ) + alpha_store(i,j,k,ispec_CPML)**2 &
                   * it*deltat * d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML)
              A4 = -alpha_store(i,j,k,ispec_CPML)**2 * d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML)
              A5 = 0.0

           else if( CPML_regions(ispec_CPML) == CPML_XZ_ONLY ) then
              !------------------------------------------------------------------------------
              !---------------------------- XZ-edge C-PML -----------------------------------
              !------------------------------------------------------------------------------

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.0d0 - exp(-bb * deltat/2.0d0) ) / bb
                 coef2_1 = (1.0d0 - exp(-bb * deltat/2.0d0) ) * exp(-bb * deltat/2.0d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              coef0_2 = coef0_1
              coef1_2 = coef1_1
              coef2_2 = coef2_1

              rmemory_potential_acoustic(i,j,k,ispec_CPML,1)=coef0_1 * rmemory_potential_acoustic(i,j,k,ispec_CPML,1) &
                    + (potential_acoustic(iglob) + deltat*potential_dot_acoustic(iglob)) * coef1_1 &
                    + potential_acoustic(iglob) * coef2_1
              rmemory_potential_acoustic(i,j,k,ispec_CPML,2)=coef0_2 * rmemory_potential_acoustic(i,j,k,ispec_CPML,2) &
                    + (potential_acoustic(iglob) + deltat*potential_dot_acoustic(iglob)) * it*deltat * coef1_2 &
                    + potential_acoustic(iglob) * it*deltat * coef2_2
              rmemory_potential_acoustic(i,j,k,ispec_CPML,3)= 0.0

              !---------------------- A0, A1, A2, A3, A4 and A5 --------------------------
              A0 = k_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML)
              A1 = d_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML)&
                   + d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML)
              A2 = d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) &
                   - alpha_store(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) &
                   - alpha_store(i,j,k,ispec_CPML) * d_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML)
              A3 = -2.0 * alpha_store(i,j,k,ispec_CPML) * d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) &
                   + alpha_store(i,j,k,ispec_CPML)**2 * ( d_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) &
                   + d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) ) + alpha_store(i,j,k,ispec_CPML)**2 &
                   * it*deltat * d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)
              A4 = -alpha_store(i,j,k,ispec_CPML)**2 * d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)
              A5 = 0.0

           else if( CPML_regions(ispec_CPML) == CPML_YZ_ONLY ) then
              !------------------------------------------------------------------------------
              !---------------------------- YZ-edge C-PML -----------------------------------
              !------------------------------------------------------------------------------

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.0d0 - exp(-bb * deltat/2.0d0) ) / bb
                 coef2_1 = (1.0d0 - exp(-bb * deltat/2.0d0) ) * exp(-bb * deltat/2.0d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              coef0_2 = coef0_1
              coef1_2 = coef1_1
              coef2_2 = coef2_1

              rmemory_potential_acoustic(i,j,k,ispec_CPML,1)=coef0_1 * rmemory_potential_acoustic(i,j,k,ispec_CPML,1) &
                    + (potential_acoustic(iglob) + deltat*potential_dot_acoustic(iglob)) * coef1_1 &
                    + potential_acoustic(iglob) * coef2_1
              rmemory_potential_acoustic(i,j,k,ispec_CPML,2)=coef0_2 * rmemory_potential_acoustic(i,j,k,ispec_CPML,2) &
                    + (potential_acoustic(iglob) + deltat*potential_dot_acoustic(iglob)) * it*deltat * coef1_2 &
                    + potential_acoustic(iglob) * it*deltat * coef2_2
              rmemory_potential_acoustic(i,j,k,ispec_CPML,3)=0.d0

              !---------------------- A0, A1, A2, A3, A4 and A5 --------------------------
              A0 = k_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML)
              A1 = d_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) &
                   + d_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML)
              A2 = d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) &
                   - alpha_store(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) &
                   - alpha_store(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML)
              A3 = -2.0 * alpha_store(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) &
                   + alpha_store(i,j,k,ispec_CPML)**2 * ( d_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) &
                   + d_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) ) + alpha_store(i,j,k,ispec_CPML)**2 &
                   * it*deltat * d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)
              A4 = -alpha_store(i,j,k,ispec_CPML)**2 * d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)
              A5 = 0.0

           else if( CPML_regions(ispec_CPML) == CPML_XYZ ) then
              !------------------------------------------------------------------------------
              !---------------------------- XYZ-corner C-PML --------------------------------
              !------------------------------------------------------------------------------

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.0d0 - exp(-bb * deltat/2.0d0) ) / bb
                 coef2_1 = (1.0d0 - exp(-bb * deltat/2.0d0) ) * exp(-bb * deltat/2.0d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              coef0_2 = coef0_1
              coef1_2 = coef1_1
              coef2_2 = coef2_1

              coef0_3 = coef0_1
              coef1_3 = coef1_1
              coef2_3 = coef2_1

              rmemory_potential_acoustic(i,j,k,ispec_CPML,1)=coef0_1 * rmemory_potential_acoustic(i,j,k,ispec_CPML,1) &
                    + (potential_acoustic(iglob) + deltat*potential_dot_acoustic(iglob)) * coef1_1 &
                    + potential_acoustic(iglob) * coef2_1
              rmemory_potential_acoustic(i,j,k,ispec_CPML,2)=coef0_2 * rmemory_potential_acoustic(i,j,k,ispec_CPML,2) &
                    + (potential_acoustic(iglob) + deltat*potential_dot_acoustic(iglob)) * it*deltat * coef1_2 &
                    + potential_acoustic(iglob) * it*deltat * coef2_2
              rmemory_potential_acoustic(i,j,k,ispec_CPML,3)=coef0_3 * rmemory_potential_acoustic(i,j,k,ispec_CPML,3) &
                    + (potential_acoustic(iglob) + deltat*potential_dot_acoustic(iglob)) * (it*deltat)**2 * coef1_3 &
                    + potential_acoustic(iglob) * (it*deltat)**2 * coef2_3

              !---------------------- A0, A1, A2, A3, A4 and A5 --------------------------
              A0 = k_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML)
              A1 = d_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) + &
                   d_store_y(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) + &
                   d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML)
              A2 = k_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) + &
                   k_store_y(i,j,k,ispec_CPML) * d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) + &
                   k_store_z(i,j,k,ispec_CPML) * d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) - &
                   d_store_x(i,j,k,ispec_CPML) * alpha_store(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) * &
                   k_store_z(i,j,k,ispec_CPML) - d_store_y(i,j,k,ispec_CPML) * alpha_store(i,j,k,ispec_CPML) * &
                   k_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) - d_store_z(i,j,k,ispec_CPML) * &
                   alpha_store(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML)
              temp_A3 = d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) - &
                   2.0 * alpha_store(i,j,k,ispec_CPML) * ( &
                   d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) + &
                   d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) + &
                   d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) &
                   ) + alpha_store(i,j,k,ispec_CPML)**2 * ( &
                   d_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) + &
                   d_store_y(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) + &
                   d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) &
                   )
!             temp_A4 = -2.0 * alpha_store(i,j,k,ispec_CPML) * d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) * &
!                  d_store_z(i,j,k,ispec_CPML) + alpha_store(i,j,k,ispec_CPML)**2 * ( &
!                  d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) + &
!                  d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) + &
!                  d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) &
!                  )
!             temp_A5 = 0.5 * d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)
!             A3 = temp_A3 + it*deltat*temp_A4 + (it*deltat)**2*temp_A5
!             A4 = -temp_A4 -2.0*it*deltat*temp_A5
!             A5 = temp_A5

!!! the full experssion of A3,A4,A5 are given by above equation, here we use reduced 
!!! exprssion of A3,A4,A5 in order to stabilized the code. 

              A3 = temp_A3
              A4 = 0.0
              A5 = 0.0

           endif

           potential_dot_dot_acoustic_CPML(i,j,k) =  wgllcube * kappal_inv *jacobianl * &
                 ( A1 * potential_dot_acoustic(iglob) + A2 * potential_acoustic(iglob) + &
                   A3 * rmemory_potential_acoustic(i,j,k,ispec_CPML,1)+ &
                   A4 * rmemory_potential_acoustic(i,j,k,ispec_CPML,2)+ &
                   A5 * rmemory_potential_acoustic(i,j,k,ispec_CPML,3)  &
                 )
        enddo
     enddo
  enddo

end subroutine pml_compute_accel_contribution_acoustic
!
!=====================================================================
!
subroutine save_field_on_pml_interface(displ,veloc,accel,nglob_interface_PML_elastic,&
                                       b_PML_field,b_reclen_PML_field)

  use specfem_par, only: NGLOB_AB,it
  use constants, only: CUSTOM_REAL,NDIM
  implicit none

  integer, intent(in) :: nglob_interface_PML_elastic,b_reclen_PML_field
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB), intent(in) :: displ,veloc,accel
  real(kind=CUSTOM_REAL), dimension(9,nglob_interface_PML_elastic) :: b_PML_field

  integer :: iglob

  do iglob = 1, nglob_interface_PML_elastic
     b_PML_field(1,iglob) = displ(1,iglob) 
     b_PML_field(2,iglob) = displ(2,iglob) 
     b_PML_field(3,iglob) = displ(3,iglob)

     b_PML_field(4,iglob) = veloc(1,iglob) 
     b_PML_field(5,iglob) = veloc(2,iglob) 
     b_PML_field(6,iglob) = veloc(3,iglob) 

     b_PML_field(7,iglob) = accel(1,iglob)
     b_PML_field(8,iglob) = accel(2,iglob) 
     b_PML_field(9,iglob) = accel(3,iglob)   
  enddo

  call write_abs(0,b_PML_field,b_reclen_PML_field,it)

end subroutine save_field_on_pml_interface
!
!=====================================================================
!
subroutine read_field_on_pml_interface(b_accel,b_veloc,b_displ,nglob_interface_PML_elastic,&
                                       b_PML_field,b_reclen_PML_field)

  use specfem_par, only: NGLOB_AB,ibool,NSTEP,it
  use pml_par, only: NSPEC_CPML,CPML_to_spec
  use constants, only: CUSTOM_REAL,NDIM,NGLLX,NGLLY,NGLLZ
  implicit none

  integer, intent(in) :: nglob_interface_PML_elastic,b_reclen_PML_field
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: b_displ,b_veloc,b_accel
  real(kind=CUSTOM_REAL), dimension(9,nglob_interface_PML_elastic) :: b_PML_field

  integer :: iglob,ispec,ispec_pml,i,j,k

  do ispec_pml = 1, NSPEC_CPML
     ispec = CPML_to_spec(ispec_pml)
     do i = 1, NGLLX; do j = 1, NGLLY; do k = 1, NGLLZ
        iglob = ibool(i,j,k,ispec)
        b_displ(:,iglob) = 0._CUSTOM_REAL
        b_veloc(:,iglob) = 0._CUSTOM_REAL 
        b_accel(:,iglob) = 0._CUSTOM_REAL
     enddo; enddo; enddo
  enddo  

  call read_abs(0,b_PML_field,b_reclen_PML_field,NSTEP-it+1)

  do iglob = 1, nglob_interface_PML_elastic
     b_displ(1,iglob) = b_PML_field(1,iglob)  
     b_displ(2,iglob) = b_PML_field(2,iglob) 
     b_displ(3,iglob) = b_PML_field(3,iglob)

     b_veloc(1,iglob) = b_PML_field(4,iglob) 
     b_veloc(2,iglob) = b_PML_field(5,iglob) 
     b_veloc(3,iglob) = b_PML_field(6,iglob) 

     b_accel(1,iglob) = b_PML_field(7,iglob)
     b_accel(2,iglob) = b_PML_field(8,iglob) 
     b_accel(3,iglob) = b_PML_field(9,iglob)   
  enddo

end subroutine read_field_on_pml_interface
!
!=====================================================================
!
subroutine save_potential_on_pml_interface(potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic,&
                                           nglob_interface_PML_acoustic,b_PML_potential,b_reclen_PML_potential)

  use specfem_par, only: NGLOB_AB,it
  use constants, only: CUSTOM_REAL
  implicit none

  integer, intent(in) :: nglob_interface_PML_acoustic,b_reclen_PML_potential
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB), intent(in) :: potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic
  real(kind=CUSTOM_REAL), dimension(3,nglob_interface_PML_acoustic) :: b_PML_potential

  integer :: iglob
          
  do iglob = 1, nglob_interface_PML_acoustic
     b_PML_potential(1,iglob) = potential_acoustic(iglob)
     b_PML_potential(2,iglob) = potential_dot_acoustic(iglob) 
     b_PML_potential(3,iglob) = potential_dot_dot_acoustic(iglob)  
  enddo

  call write_abs(1,b_PML_potential,b_reclen_PML_potential,it)

end subroutine save_potential_on_pml_interface
!
!=====================================================================
!
subroutine read_potential_on_pml_interface(b_potential_dot_dot_acoustic,b_potential_dot_acoustic,b_potential_acoustic,&
                                           nglob_interface_PML_acoustic,b_PML_potential,b_reclen_PML_potential)

  use specfem_par, only: NGLOB_AB,ibool,NSTEP,it
  use pml_par, only: NSPEC_CPML,CPML_to_spec
  use constants, only: CUSTOM_REAL,NDIM,NGLLX,NGLLY,NGLLZ
  implicit none

  integer, intent(in) :: nglob_interface_PML_acoustic,b_reclen_PML_potential
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: b_potential_dot_dot_acoustic,b_potential_dot_acoustic,b_potential_acoustic
  real(kind=CUSTOM_REAL), dimension(3,nglob_interface_PML_acoustic) :: b_PML_potential

  integer :: iglob,ispec,ispec_pml,i,j,k

  do ispec_pml = 1, NSPEC_CPML
     ispec = CPML_to_spec(ispec_pml)
     do i = 1, NGLLX; do j = 1, NGLLY; do k = 1, NGLLZ
        iglob = ibool(i,j,k,ispec)
        b_potential_acoustic(iglob) = 0._CUSTOM_REAL
        b_potential_dot_acoustic(iglob) = 0._CUSTOM_REAL 
        b_potential_dot_dot_acoustic(iglob) = 0._CUSTOM_REAL      
     enddo; enddo; enddo
  enddo  

  call read_abs(1,b_PML_potential,b_reclen_PML_potential,NSTEP-it+1)

  do iglob = 1, nglob_interface_PML_acoustic
     b_potential_acoustic(iglob) = b_PML_potential(1,iglob)  
     b_potential_dot_acoustic(iglob) = b_PML_potential(2,iglob) 
     b_potential_dot_dot_acoustic(iglob) = b_PML_potential(3,iglob)   
  enddo

end subroutine read_potential_on_pml_interface
