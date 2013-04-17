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

subroutine pml_compute_accel_contribution(ispec,ispec_CPML,deltat,nspec_AB,jacobian,accel_elastic_CPML)

  ! calculates contribution from each C-PML element to update acceleration to the global mesh

  ! second-order accurate convolution term calculation from equation (21) of
  ! Shumin Wang, Robert Lee, and Fernando L. Teixeira,
  ! Anisotropic-Medium PML for Vector FETD With Modified Basis Functions,
  ! IEEE Transactions on Antennas and Propagation, vol. 54, no. 1, (2006)

  use specfem_par, only: ibool,wgllwgll_yz,wgllwgll_xz,wgllwgll_xy,it,kappastore,rhostore
  use specfem_par_elastic, only: displ,veloc,ispec_is_elastic
  use specfem_par_acoustic, only: potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic,ispec_is_acoustic
  use pml_par, only: NSPEC_CPML,rmemory_displ_elastic,rmemory_potential_acoustic,CPML_regions,spec_to_CPML,alpha_store, &
                     d_store_x,d_store_y,d_store_z,K_store_x,K_store_y,K_store_z,potential_dot_dot_acoustic_CPML
  use constants, only: CUSTOM_REAL,NDIM,NGLLX,NGLLY,NGLLZ, &
                       CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY,CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ

  implicit none

  integer, intent(in) :: ispec,ispec_CPML,nspec_AB

  real(kind=CUSTOM_REAL), intent(in) :: deltat

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ,NSPEC_CPML), intent(out) :: accel_elastic_CPML
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_AB), intent(in) :: jacobian

  ! local parameters
  integer :: i,j,k,iglob

  real(kind=CUSTOM_REAL) :: fac1,fac2,fac3,fac4,rhol,kappal,jacobianl
  real(kind=CUSTOM_REAL) :: bb,coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2,coef0_3,coef1_3,coef2_3
  real(kind=CUSTOM_REAL) :: A0,A1,A2,A3,A4,A5 ! for convolution of acceleration
  real(kind=CUSTOM_REAL) :: temp_A3,temp_A4,temp_A5

  do k=1,NGLLZ
     do j=1,NGLLY
        do i=1,NGLLX
           if( ispec_is_elastic(ispec) ) then
              rhol = rhostore(i,j,k,ispec)
              jacobianl = jacobian(i,j,k,ispec)
           else if( ispec_is_acoustic(ispec) ) then
              kappal = kappastore(i,j,k,ispec)
              jacobianl = jacobian(i,j,k,ispec)
           else
              stop 'CPML elements should be either acoustic or elastic; exiting...'
           endif

           iglob = ibool(i,j,k,ispec)

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

              if( ispec_is_elastic(ispec) ) then

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

              else if( ispec_is_acoustic(ispec) ) then

                 rmemory_potential_acoustic(i,j,k,ispec_CPML,1)=coef0_1 * rmemory_potential_acoustic(i,j,k,ispec_CPML,1) &
                      + (potential_acoustic(iglob) + deltat*potential_dot_acoustic(iglob)) * coef1_1 &
                      + potential_acoustic(iglob) * coef2_1
                 rmemory_potential_acoustic(i,j,k,ispec_CPML,2)=0.0
                 rmemory_potential_acoustic(i,j,k,ispec_CPML,3)=0.0
              endif

              !---------------------- A0, A1, A2, A3, A4 and A5 --------------------------
              A0 = k_store_x(i,j,k,ispec_CPML)
              A1 = d_store_x(i,j,k,ispec_CPML)
              A2 = - alpha_store(i,j,k,ispec_CPML) * d_store_x(i,j,k,ispec_CPML)
              A3 = d_store_x(i,j,k,ispec_CPML) * alpha_store(i,j,k,ispec_CPML) ** 2
              A4 = 0.d0
              A5 = 0.d0

              fac1 = wgllwgll_yz(j,k)
              fac2 = wgllwgll_xz(i,k)
              fac3 = wgllwgll_xy(i,j)
              fac4 = sqrt(fac1 * fac2 * fac3)

              if( ispec_is_elastic(ispec) ) then

                 accel_elastic_CPML(1,i,j,k,ispec_CPML) =  fac4 * rhol * jacobianl * &
                      ( A1 * veloc(1,iglob) + A2 * displ(1,iglob) + &
                      A3 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) + &
                      A4 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,2) + &
                      A5 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,3)  &
                      )

                 accel_elastic_CPML(2,i,j,k,ispec_CPML) =  fac4 * rhol * jacobianl * &
                      ( A1 * veloc(2,iglob) + A2 * displ(2,iglob) + &
                      A3 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) + &
                      A4 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,2) + &
                      A5 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,3)  &
                      )

                 accel_elastic_CPML(3,i,j,k,ispec_CPML) =  fac4 * rhol * jacobianl * &
                      ( A1 * veloc(3,iglob) + A2 * displ(3,iglob) + &
                      A3 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) + &
                      A4 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,2) + &
                      A5 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,3)  &
                      )

              else if( ispec_is_acoustic(ispec) ) then

                 potential_dot_dot_acoustic_CPML(i,j,k,ispec_CPML) =  fac4 * 1.0/kappal *jacobianl * &
                      ( A1 * potential_dot_acoustic(iglob) + A2 * potential_acoustic(iglob) + &
                      A3 * rmemory_potential_acoustic(i,j,k,ispec_CPML,1)+ &
                      A4 * rmemory_potential_acoustic(i,j,k,ispec_CPML,2)+ &
                      A5 * rmemory_potential_acoustic(i,j,k,ispec_CPML,3)  &
                      )
              endif

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

              if( ispec_is_elastic(ispec) ) then

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

              else if( ispec_is_acoustic(ispec) ) then

                 rmemory_potential_acoustic(i,j,k,ispec_CPML,1)=coef0_1 * rmemory_potential_acoustic(i,j,k,ispec_CPML,1) &
                      + (potential_acoustic(iglob) + deltat*potential_dot_acoustic(iglob)) * coef1_1 &
                      + potential_acoustic(iglob) * coef2_1
                 rmemory_potential_acoustic(i,j,k,ispec_CPML,2)=0.0
                 rmemory_potential_acoustic(i,j,k,ispec_CPML,3)=0.0
              endif

              !---------------------- A0, A1, A2, A3, A4 and A5 --------------------------
              A0 = k_store_y(i,j,k,ispec_CPML)
              A1 = d_store_y(i,j,k,ispec_CPML)
              A2 = - alpha_store(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML)
              A3 = d_store_y(i,j,k,ispec_CPML) * alpha_store(i,j,k,ispec_CPML) ** 2
              A4 = 0.d0
              A5 = 0.d0

              fac1 = wgllwgll_yz(j,k)
              fac2 = wgllwgll_xz(i,k)
              fac3 = wgllwgll_xy(i,j)
              fac4 = sqrt(fac1 * fac2 * fac3)

              if( ispec_is_elastic(ispec) ) then

                 accel_elastic_CPML(1,i,j,k,ispec_CPML) =  fac4 * rhol * jacobianl * &
                      ( A1 * veloc(1,iglob) + A2 * displ(1,iglob) + &
                      A3 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) + &
                      A4 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,2) + &
                      A5 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,3)  &
                      )

                 accel_elastic_CPML(2,i,j,k,ispec_CPML) =  fac4 * rhol * jacobianl * &
                      ( A1 * veloc(2,iglob) + A2 * displ(2,iglob) + &
                      A3 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) + &
                      A4 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,2) + &
                      A5 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,3)  &
                      )

                 accel_elastic_CPML(3,i,j,k,ispec_CPML) =  fac4 * rhol * jacobianl * &
                      ( A1 * veloc(3,iglob) + A2 * displ(3,iglob) + &
                      A3 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) + &
                      A4 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,2) + &
                      A5 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,3)  &
                      )

              else if( ispec_is_acoustic(ispec) ) then

                 potential_dot_dot_acoustic_CPML(i,j,k,ispec_CPML) =  fac4 * 1.0/kappal *jacobianl * &
                      ( A1 * potential_dot_acoustic(iglob) + A2 * potential_acoustic(iglob) + &
                      A3 * rmemory_potential_acoustic(i,j,k,ispec_CPML,1)+ &
                      A4 * rmemory_potential_acoustic(i,j,k,ispec_CPML,2)+ &
                      A5 * rmemory_potential_acoustic(i,j,k,ispec_CPML,3)  &
                      )
              endif

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

              if( ispec_is_elastic(ispec) ) then

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

              else if( ispec_is_acoustic(ispec) ) then

                 rmemory_potential_acoustic(i,j,k,ispec_CPML,1)=coef0_1 * rmemory_potential_acoustic(i,j,k,ispec_CPML,1) &
                      + (potential_acoustic(iglob) + deltat*potential_dot_acoustic(iglob)) * coef1_1 &
                      + potential_acoustic(iglob) * coef2_1
                 rmemory_potential_acoustic(i,j,k,ispec_CPML,2)=0.0
                 rmemory_potential_acoustic(i,j,k,ispec_CPML,3)=0.0
              endif

              !---------------------- A0, A1, A2, A3, A4 and A5 --------------------------
              A0 = k_store_z(i,j,k,ispec_CPML)
              A1 = d_store_z(i,j,k,ispec_CPML)
              A2 = - alpha_store(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)
              A3 = d_store_z(i,j,k,ispec_CPML) * alpha_store(i,j,k,ispec_CPML) ** 2
              A4 = 0.d0
              A5 = 0.d0

              fac1 = wgllwgll_yz(j,k)
              fac2 = wgllwgll_xz(i,k)
              fac3 = wgllwgll_xy(i,j)
              fac4 = sqrt(fac1 * fac2 * fac3)

              if( ispec_is_elastic(ispec) ) then

                 accel_elastic_CPML(1,i,j,k,ispec_CPML) =  fac4 * rhol * jacobianl * &
                      ( A1 * veloc(1,iglob) + A2 * displ(1,iglob) + &
                      A3 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) + &
                      A4 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,2) + &
                      A5 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,3)  &
                      )

                 accel_elastic_CPML(2,i,j,k,ispec_CPML) =  fac4 * rhol * jacobianl * &
                      ( A1 * veloc(2,iglob) + A2 * displ(2,iglob) + &
                      A3 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) + &
                      A4 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,2) + &
                      A5 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,3)  &
                      )

                 accel_elastic_CPML(3,i,j,k,ispec_CPML) =  fac4 * rhol * jacobianl * &
                      ( A1 * veloc(3,iglob) + A2 * displ(3,iglob) + &
                      A3 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) + &
                      A4 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,2) + &
                      A5 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,3)  &
                      )

              else if( ispec_is_acoustic(ispec) ) then

                 potential_dot_dot_acoustic_CPML(i,j,k,ispec_CPML) =  fac4 * 1.0/kappal *jacobianl * &
                      ( A1 * potential_dot_acoustic(iglob) + A2 * potential_acoustic(iglob) + &
                      A3 * rmemory_potential_acoustic(i,j,k,ispec_CPML,1)+ &
                      A4 * rmemory_potential_acoustic(i,j,k,ispec_CPML,2)+ &
                      A5 * rmemory_potential_acoustic(i,j,k,ispec_CPML,3)  &
                      )
              endif

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

              if( ispec_is_elastic(ispec) ) then

                 rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) = coef0_1 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) &
                      + (displ(1,iglob) + deltat * veloc(1,iglob)) * coef1_1 + displ(1,iglob) * coef2_1
                 rmemory_displ_elastic(1,i,j,k,ispec_CPML,2) = coef0_2 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,2) &
                      + (displ(1,iglob) + deltat * veloc(1,iglob)) * (it+0.0) * deltat * coef1_2 &
                      + displ(1,iglob) * (it-0.0) * deltat * coef2_2
                 rmemory_displ_elastic(1,i,j,k,ispec_CPML,3) = 0.0

                 rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) = coef0_1 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) &
                      + (displ(2,iglob) + deltat * veloc(2,iglob)) * coef1_1 + displ(2,iglob) * coef2_1
                 rmemory_displ_elastic(2,i,j,k,ispec_CPML,2) = coef0_2 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,2) &
                      + (displ(2,iglob) + deltat * veloc(2,iglob)) * (it+0.0) * deltat * coef1_2 &
                      + displ(2,iglob) * (it-0.0) * deltat * coef2_2
                 rmemory_displ_elastic(2,i,j,k,ispec_CPML,3) = 0.0

                 rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) = coef0_1 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) &
                      + (displ(3,iglob) + deltat * veloc(3,iglob)) * coef1_1 + displ(3,iglob) * coef2_1
                 rmemory_displ_elastic(3,i,j,k,ispec_CPML,2) = coef0_2 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,2) &
                      + (displ(3,iglob) + deltat * veloc(3,iglob)) * (it+0.0) * deltat * coef1_2 &
                      + displ(3,iglob) * (it-0.0) * deltat * coef2_2
                 rmemory_displ_elastic(3,i,j,k,ispec_CPML,3) = 0.0

              else if( ispec_is_acoustic(ispec) ) then

                 rmemory_potential_acoustic(i,j,k,ispec_CPML,1)=coef0_1 * rmemory_potential_acoustic(i,j,k,ispec_CPML,1) &
                      + (potential_acoustic(iglob) + deltat*potential_dot_acoustic(iglob)) * coef1_1 &
                      + potential_acoustic(iglob) * coef2_1
                 rmemory_potential_acoustic(i,j,k,ispec_CPML,2)=coef0_2 * rmemory_potential_acoustic(i,j,k,ispec_CPML,2) &
                      + (potential_acoustic(iglob) + deltat*potential_dot_acoustic(iglob)) * (it+0.5)*deltat * coef1_2 &
                      + potential_acoustic(iglob) * (it-0.5)*deltat * coef2_2
                 rmemory_potential_acoustic(i,j,k,ispec_CPML,3)=0.0
              endif

              !---------------------- A0, A1, A2, A3, A4 and A5 --------------------------
              A0 = k_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML)
              A1 = d_store_y(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) &
                   + d_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML)
              A2 = d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) &
                   - alpha_store(i,j,k,ispec_CPML) * d_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) &
                   - alpha_store(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML)
              if( ispec_is_elastic(ispec) ) then
                 A3 = -2.0 * alpha_store(i,j,k,ispec_CPML) * d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) &
                      + alpha_store(i,j,k,ispec_CPML)**2 * ( d_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) &
                      + d_store_y(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) ) + alpha_store(i,j,k,ispec_CPML)**2 &
                      * (it+0.0) * deltat * d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML)
              else if( ispec_is_acoustic(ispec) ) then
                 A3 = -2.0 * alpha_store(i,j,k,ispec_CPML) * d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) + &
                      alpha_store(i,j,k,ispec_CPML)**2 * ( d_store_x(i,j,k,ispec_CPML)*k_store_y(i,j,k,ispec_CPML) &
                      + d_store_y(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) ) + alpha_store(i,j,k,ispec_CPML)**2 &
                      * (it+0.5) * deltat * d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML)
              endif
              A4 = -alpha_store(i,j,k,ispec_CPML)**2 * d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML)
              A5 = 0.0

              fac1 = wgllwgll_yz(j,k)
              fac2 = wgllwgll_xz(i,k)
              fac3 = wgllwgll_xy(i,j)
              fac4 = sqrt(fac1 * fac2 * fac3)

              if( ispec_is_elastic(ispec) ) then

                 accel_elastic_CPML(1,i,j,k,ispec_CPML) =  fac4 * rhol * jacobianl * &
                      ( A1 * veloc(1,iglob) + A2 * displ(1,iglob) + &
                      A3 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) + &
                      A4 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,2) + &
                      A5 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,3)  &
                      )

                 accel_elastic_CPML(2,i,j,k,ispec_CPML) =  fac4 * rhol * jacobianl * &
                      ( A1 * veloc(2,iglob) + A2 * displ(2,iglob) + &
                      A3 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) + &
                      A4 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,2) + &
                      A5 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,3)  &
                      )

                 accel_elastic_CPML(3,i,j,k,ispec_CPML) =  fac4 * rhol * jacobianl * &
                      ( A1 * veloc(3,iglob) + A2 * displ(3,iglob) + &
                      A3 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) + &
                      A4 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,2) + &
                      A5 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,3)  &
                      )

              else if( ispec_is_acoustic(ispec) ) then

                 potential_dot_dot_acoustic_CPML(i,j,k,ispec_CPML) =  fac4 * 1.0/kappal *jacobianl * &
                      ( A1 * potential_dot_acoustic(iglob) + A2 * potential_acoustic(iglob) + &
                      A3 * rmemory_potential_acoustic(i,j,k,ispec_CPML,1)+ &
                      A4 * rmemory_potential_acoustic(i,j,k,ispec_CPML,2)+ &
                      A5 * rmemory_potential_acoustic(i,j,k,ispec_CPML,3)  &
                      )
              endif

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

              if( ispec_is_elastic(ispec) ) then

                 rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) = coef0_1 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) &
                      + (displ(1,iglob) + deltat * veloc(1,iglob)) * coef1_1 + displ(1,iglob) * coef2_1
                 rmemory_displ_elastic(1,i,j,k,ispec_CPML,2) = coef0_2 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,2) &
                      + (displ(1,iglob) + deltat * veloc(1,iglob)) * (it+0.0) * deltat * coef1_2 &
                      + displ(1,iglob) * (it-0.0) * deltat * coef2_2
                 rmemory_displ_elastic(1,i,j,k,ispec_CPML,3)= 0.0

                 rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) = coef0_1 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) &
                      + (displ(2,iglob) + deltat * veloc(2,iglob)) * coef1_1 + displ(2,iglob) * coef2_1
                 rmemory_displ_elastic(2,i,j,k,ispec_CPML,2) = coef0_2 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,2) &
                      + (displ(2,iglob) + deltat * veloc(2,iglob)) * (it+0.0) * deltat * coef1_2 &
                      + displ(2,iglob) * (it-0.0) * deltat * coef2_2
                 rmemory_displ_elastic(2,i,j,k,ispec_CPML,3)= 0.0

                 rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) = coef0_1 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) &
                      + (displ(3,iglob) + deltat * veloc(3,iglob)) * coef1_1 + displ(3,iglob) * coef2_1
                 rmemory_displ_elastic(3,i,j,k,ispec_CPML,2) = coef0_2 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,2) &
                      + (displ(3,iglob) + deltat * veloc(3,iglob)) * (it+0.0) * deltat * coef1_2 &
                      + displ(3,iglob) * (it-0.0) * deltat * coef2_2
                 rmemory_displ_elastic(3,i,j,k,ispec_CPML,3)= 0.0

              else if( ispec_is_acoustic(ispec) ) then

                 rmemory_potential_acoustic(i,j,k,ispec_CPML,1)=coef0_1 * rmemory_potential_acoustic(i,j,k,ispec_CPML,1) &
                      + (potential_acoustic(iglob) + deltat*potential_dot_acoustic(iglob)) * coef1_1 &
                      + potential_acoustic(iglob) * coef2_1
                 rmemory_potential_acoustic(i,j,k,ispec_CPML,2)=coef0_2 * rmemory_potential_acoustic(i,j,k,ispec_CPML,2) &
                      + (potential_acoustic(iglob) + deltat*potential_dot_acoustic(iglob)) * (it+0.5)*deltat * coef1_2 &
                      + potential_acoustic(iglob) * (it-0.5)*deltat * coef2_2
                 rmemory_potential_acoustic(i,j,k,ispec_CPML,3)= 0.0
              endif

              !---------------------- A0, A1, A2, A3, A4 and A5 --------------------------
              A0 = k_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML)
              A1 = d_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML)&
                   + d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML)
              A2 = d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) &
                   - alpha_store(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) &
                   - alpha_store(i,j,k,ispec_CPML) * d_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML)
              if( ispec_is_elastic(ispec) ) then
                 A3 = -2.0 * alpha_store(i,j,k,ispec_CPML) * d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) &
                      + alpha_store(i,j,k,ispec_CPML)**2 * ( d_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) &
                      + d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) ) + alpha_store(i,j,k,ispec_CPML)**2 &
                      * (it+0.0) * deltat * d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)
              else if( ispec_is_acoustic(ispec) ) then
                 A3 = -2.0 * alpha_store(i,j,k,ispec_CPML) * d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) &
                      + alpha_store(i,j,k,ispec_CPML)**2 * ( d_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) &
                      + d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) ) + alpha_store(i,j,k,ispec_CPML)**2 &
                      * (it+0.5) * deltat * d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)
              endif
              A4 = -alpha_store(i,j,k,ispec_CPML)**2 * d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)
              A5 = 0.0

              fac1 = wgllwgll_yz(j,k)
              fac2 = wgllwgll_xz(i,k)
              fac3 = wgllwgll_xy(i,j)
              fac4 = sqrt(fac1 * fac2 * fac3)

              if( ispec_is_elastic(ispec) ) then

                 accel_elastic_CPML(1,i,j,k,ispec_CPML) =  fac4 * rhol * jacobianl * &
                      ( A1 * veloc(1,iglob) + A2 * displ(1,iglob) + &
                      A3 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) + &
                      A4 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,2) + &
                      A5 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,3)  &
                      )

                 accel_elastic_CPML(2,i,j,k,ispec_CPML) =  fac4 * rhol * jacobianl * &
                      ( A1 * veloc(2,iglob) + A2 * displ(2,iglob) + &
                      A3 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) + &
                      A4 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,2) + &
                      A5 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,3)  &
                      )

                 accel_elastic_CPML(3,i,j,k,ispec_CPML) =  fac4 * rhol * jacobianl * &
                      ( A1 * veloc(3,iglob) + A2 * displ(3,iglob) + &
                      A3 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) + &
                      A4 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,2) + &
                      A5 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,3)  &
                      )

              else if( ispec_is_acoustic(ispec) ) then

                 potential_dot_dot_acoustic_CPML(i,j,k,ispec_CPML) =  fac4 * 1.0/kappal *jacobianl * &
                      ( A1 * potential_dot_acoustic(iglob) + A2 * potential_acoustic(iglob) + &
                      A3 * rmemory_potential_acoustic(i,j,k,ispec_CPML,1)+ &
                      A4 * rmemory_potential_acoustic(i,j,k,ispec_CPML,2)+ &
                      A5 * rmemory_potential_acoustic(i,j,k,ispec_CPML,3)  &
                      )
              endif

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

              if( ispec_is_elastic(ispec) ) then

                 rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) = coef0_1 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) &
                      + (displ(1,iglob) + deltat * veloc(1,iglob)) * coef1_1 + displ(1,iglob) * coef2_1
                 rmemory_displ_elastic(1,i,j,k,ispec_CPML,2) = coef0_2 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,2) &
                      + (displ(1,iglob) + deltat * veloc(1,iglob)) * (it+0.0) * deltat * coef1_2 &
                      + displ(1,iglob) * (it-0.0) * deltat * coef2_2
                 rmemory_displ_elastic(1,i,j,k,ispec_CPML,3)=0.d0

                 rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) = coef0_1 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) &
                      + (displ(2,iglob) + deltat * veloc(2,iglob)) * coef1_1 + displ(2,iglob) * coef2_1
                 rmemory_displ_elastic(2,i,j,k,ispec_CPML,2) = coef0_2 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,2) &
                      + (displ(2,iglob) + deltat * veloc(2,iglob)) * (it+0.0) * deltat * coef1_2 &
                      + displ(2,iglob) * (it-0.0) * deltat * coef2_2
                 rmemory_displ_elastic(2,i,j,k,ispec_CPML,3)=0.d0

                 rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) = coef0_1 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) &
                      + (displ(3,iglob) + deltat * veloc(3,iglob)) * coef1_1 + displ(3,iglob) * coef2_1
                 rmemory_displ_elastic(3,i,j,k,ispec_CPML,2) = coef0_2 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,2) &
                      + (displ(3,iglob) + deltat * veloc(3,iglob)) * (it+0.0) * deltat * coef1_2 &
                      + displ(3,iglob) * (it-0.0) * deltat * coef2_2
                 rmemory_displ_elastic(3,i,j,k,ispec_CPML,3)=0.d0

              else if( ispec_is_acoustic(ispec) ) then

                 rmemory_potential_acoustic(i,j,k,ispec_CPML,1)=coef0_1 * rmemory_potential_acoustic(i,j,k,ispec_CPML,1) &
                      + (potential_acoustic(iglob) + deltat*potential_dot_acoustic(iglob)) * coef1_1 &
                      + potential_acoustic(iglob) * coef2_1
                 rmemory_potential_acoustic(i,j,k,ispec_CPML,2)=coef0_2 * rmemory_potential_acoustic(i,j,k,ispec_CPML,2) &
                      + (potential_acoustic(iglob) + deltat*potential_dot_acoustic(iglob)) * (it+0.5)*deltat * coef1_2 &
                      + potential_acoustic(iglob) * (it-0.5)*deltat * coef2_2
                 rmemory_potential_acoustic(i,j,k,ispec_CPML,3)=0.d0
              endif

              !---------------------- A0, A1, A2, A3, A4 and A5 --------------------------
              A0 = k_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML)
              A1 = d_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) &
                   + d_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML)
              A2 = d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) &
                   - alpha_store(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) &
                   - alpha_store(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML)
              if( ispec_is_elastic(ispec) ) then
                 A3 = -2.0 * alpha_store(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) &
                      + alpha_store(i,j,k,ispec_CPML)**2 * ( d_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) &
                      + d_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) ) + alpha_store(i,j,k,ispec_CPML)**2 &
                      * (it+0.0) * deltat * d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)
              else if( ispec_is_acoustic(ispec) ) then
                 A3 = -2.0 * alpha_store(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) &
                      + alpha_store(i,j,k,ispec_CPML)**2 * ( d_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) &
                      + d_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) ) + alpha_store(i,j,k,ispec_CPML)**2 &
                      * (it+0.5) * deltat * d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)
              endif
              A4 = -alpha_store(i,j,k,ispec_CPML)**2 * d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)
              A5 = 0.0

              fac1 = wgllwgll_yz(j,k)
              fac2 = wgllwgll_xz(i,k)
              fac3 = wgllwgll_xy(i,j)
              fac4 = sqrt(fac1 * fac2 * fac3)

              if( ispec_is_elastic(ispec) ) then

                 accel_elastic_CPML(1,i,j,k,ispec_CPML) =  fac4 * rhol * jacobianl * &
                      ( A1 * veloc(1,iglob) + A2 * displ(1,iglob) + &
                      A3 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) + &
                      A4 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,2) + &
                      A5 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,3)  &
                      )

                 accel_elastic_CPML(2,i,j,k,ispec_CPML) =  fac4 * rhol * jacobianl * &
                      ( A1 * veloc(2,iglob) + A2 * displ(2,iglob) + &
                      A3 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) + &
                      A4 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,2) + &
                      A5 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,3)  &
                      )

                 accel_elastic_CPML(3,i,j,k,ispec_CPML) =  fac4 * rhol * jacobianl * &
                      ( A1 * veloc(3,iglob) + A2 * displ(3,iglob) + &
                      A3 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) + &
                      A4 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,2) + &
                      A5 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,3)  &
                      )

              else if( ispec_is_acoustic(ispec) ) then

                 potential_dot_dot_acoustic_CPML(i,j,k,ispec_CPML) =  fac4 * 1.0/kappal *jacobianl * &
                      ( A1 * potential_dot_acoustic(iglob) + A2 * potential_acoustic(iglob) + &
                      A3 * rmemory_potential_acoustic(i,j,k,ispec_CPML,1)+ &
                      A4 * rmemory_potential_acoustic(i,j,k,ispec_CPML,2)+ &
                      A5 * rmemory_potential_acoustic(i,j,k,ispec_CPML,3)  &
                      )
              endif

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

              if( ispec_is_elastic(ispec) ) then

                 rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) = coef0_1 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) &
                      + (displ(1,iglob) + deltat * veloc(1,iglob)) * coef1_1 + displ(1,iglob) * coef2_1
                 rmemory_displ_elastic(1,i,j,k,ispec_CPML,2) = coef0_2 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,2) &
                      + (displ(1,iglob) + deltat * veloc(1,iglob)) * (it+0.0) * deltat * coef1_2 &
                      + displ(1,iglob) * (it-0.0) * deltat * coef2_2
                 rmemory_displ_elastic(1,i,j,k,ispec_CPML,3) = coef0_3 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,3) &
                      + (displ(1,iglob) + deltat * veloc(1,iglob)) * ((it+0.0) * deltat)**2 * coef1_3 &
                      + displ(1,iglob) * ((it-0.0) * deltat)**2 * coef2_3

                 rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) = coef0_1 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) &
                      + (displ(2,iglob) + deltat * veloc(2,iglob)) * coef1_1 + displ(2,iglob) * coef2_1
                 rmemory_displ_elastic(2,i,j,k,ispec_CPML,2) = coef0_2 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,2) &
                      + (displ(2,iglob) + deltat * veloc(2,iglob)) * (it+0.0) * deltat * coef1_2 &
                      + displ(2,iglob) * (it-0.0) * deltat * coef2_2
                 rmemory_displ_elastic(2,i,j,k,ispec_CPML,3) = coef0_3 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,3) &
                      + (displ(2,iglob) + deltat * veloc(2,iglob)) * ((it+0.0) * deltat)**2 * coef1_3 &
                      + displ(2,iglob) * ((it-0.0) * deltat)**2 * coef2_3

                 rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) = coef0_1 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) &
                      + (displ(3,iglob) + deltat * veloc(3,iglob)) * coef1_1 + displ(3,iglob) * coef2_1
                 rmemory_displ_elastic(3,i,j,k,ispec_CPML,2) = coef0_2 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,2) &
                      + (displ(3,iglob) + deltat * veloc(3,iglob)) * (it+0.0) * deltat * coef1_2 &
                      + displ(3,iglob) * (it-0.0) * deltat * coef2_2
                 rmemory_displ_elastic(3,i,j,k,ispec_CPML,3) = coef0_3 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,3) &
                      + (displ(3,iglob) + deltat * veloc(3,iglob)) * ((it+0.0) * deltat)**2 * coef1_3 &
                      + displ(3,iglob) * ((it-0.0) * deltat)**2 * coef2_3

              else if( ispec_is_acoustic(ispec) ) then

                 rmemory_potential_acoustic(i,j,k,ispec_CPML,1)=coef0_1 * rmemory_potential_acoustic(i,j,k,ispec_CPML,1) &
                      + (potential_acoustic(iglob) + deltat*potential_dot_acoustic(iglob)) * coef1_1 &
                      + potential_acoustic(iglob) * coef2_1
                 rmemory_potential_acoustic(i,j,k,ispec_CPML,2)=coef0_2 * rmemory_potential_acoustic(i,j,k,ispec_CPML,2) &
                      + (potential_acoustic(iglob) + deltat*potential_dot_acoustic(iglob)) * (it+0.5)*deltat * coef1_2 &
                      + potential_acoustic(iglob) * (it-0.5)*deltat * coef2_2
                 rmemory_potential_acoustic(i,j,k,ispec_CPML,3)=coef0_3 * rmemory_potential_acoustic(i,j,k,ispec_CPML,3) &
                      + (potential_acoustic(iglob) + deltat*potential_dot_acoustic(iglob)) * ((it+0.5)*deltat)**2 * coef1_3 &
                      + potential_acoustic(iglob) * ((it-0.5)*deltat)**2 * coef2_3
              endif

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
              temp_A4 = -2.0 * alpha_store(i,j,k,ispec_CPML) * d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) * &
                   d_store_z(i,j,k,ispec_CPML) + alpha_store(i,j,k,ispec_CPML)**2 * ( &
                   d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) + &
                   d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) + &
                   d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) &
                   )
              temp_A5 = 0.5 * d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)

              !                  A3 = temp_A3 + (it+0.0) * deltat*temp_A4 + ((it+0.0) * deltat)**2*temp_A5
              !                  A4 = -temp_A4-2.0*(it+0.0) * deltat*temp_A5
              !                  A5 = temp_A5

              if( ispec_is_elastic(ispec) ) then
                 A3 = temp_A3 !+ (it+0.0) * deltat*temp_A4 !+ ((it+0.0) * deltat)**2*temp_A5
                 A4 = 0.0 !-temp_A4 ! -2.0*(it+0.0) * deltat*temp_A5
              else if( ispec_is_acoustic(ispec)) then
                 A3 = temp_A3 !+ (it+0.5)*deltat*temp_A4 !+ ((it+0.5)*deltat)**2*temp_A5
                 A4 = 0.0 !-temp_A4 !-2.0*(it+0.5)*deltat*temp_A5
              endif
              A5 = 0.0 ! temp_A5

              fac1 = wgllwgll_yz(j,k)
              fac2 = wgllwgll_xz(i,k)
              fac3 = wgllwgll_xy(i,j)
              fac4 = sqrt(fac1 * fac2 * fac3)

              if( ispec_is_elastic(ispec) ) then

                 accel_elastic_CPML(1,i,j,k,ispec_CPML) =  fac4 * rhol * jacobianl * &
                      ( A1 * veloc(1,iglob) + A2 * displ(1,iglob) + &
                      A3 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) + &
                      A4 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,2) + &
                      A5 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,3)  &
                      )

                 accel_elastic_CPML(2,i,j,k,ispec_CPML) =  fac4 * rhol * jacobianl * &
                      ( A1 * veloc(2,iglob) + A2 * displ(2,iglob) + &
                      A3 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) + &
                      A4 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,2) + &
                      A5 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,3)  &
                      )

                 accel_elastic_CPML(3,i,j,k,ispec_CPML) =  fac4 * rhol * jacobianl * &
                      ( A1 * veloc(3,iglob) + A2 * displ(3,iglob) + &
                      A3 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) + &
                      A4 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,2) + &
                      A5 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,3)  &
                      )

              else if( ispec_is_acoustic(ispec) ) then

                 potential_dot_dot_acoustic_CPML(i,j,k,ispec_CPML) =  fac4 * 1.0/kappal *jacobianl * &
                      ( A1 * potential_dot_acoustic(iglob) + A2 * potential_acoustic(iglob) + &
                      A3 * rmemory_potential_acoustic(i,j,k,ispec_CPML,1)+ &
                      A4 * rmemory_potential_acoustic(i,j,k,ispec_CPML,2)+ &
                      A5 * rmemory_potential_acoustic(i,j,k,ispec_CPML,3)  &
                      )
              endif
           endif
        enddo
     enddo
  enddo

end subroutine pml_compute_accel_contribution
