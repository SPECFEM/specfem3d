!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
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
  ! Anisotropic-medium PML for vector FETD with modified basis functions,
  ! IEEE Transactions on Antennas and Propagation, vol. 54, no. 1, (2006)

  use specfem_par, only: NGLOB_AB,deltat,wgll_cube,jacobian,ibool,rhostore
  use pml_par, only: CPML_regions,d_store_x,d_store_y,d_store_z,K_store_x,K_store_y,K_store_z, &
                     alpha_store_x, alpha_store_y, alpha_store_z, &
                     NSPEC_CPML,accel_elastic_CPML,PML_displ_old,PML_displ_new
  use constants, only: CUSTOM_REAL,NDIM,NGLLX,NGLLY,NGLLZ,CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY, &
                       CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ

  implicit none

  integer, intent(in) :: ispec,ispec_CPML
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB), intent(in) :: displ,veloc
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ,NSPEC_CPML,3) :: rmemory_displ_elastic

  ! local parameters
  integer :: i,j,k,iglob,CPML_region_local

  real(kind=CUSTOM_REAL) :: wgllcube,rhol,jacobianl
  real(kind=CUSTOM_REAL) :: alpha_x,alpha_y,alpha_z,d_x,d_y,d_z,kappa_x,kappa_y,kappa_z
  real(kind=CUSTOM_REAL) :: coef0_x,coef1_x,coef2_x,coef0_y,coef1_y,coef2_y,coef0_z,coef1_z,coef2_z
  real(kind=CUSTOM_REAL) :: A_0,A_1,A_2,A_3,A_4,A_5

  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        iglob = ibool(i,j,k,ispec)

        rhol = rhostore(i,j,k,ispec)
        jacobianl = jacobian(i,j,k,ispec)
        wgllcube = wgll_cube(i,j,k)

        ! PML coefficient values
        CPML_region_local = CPML_regions(ispec_CPML)

        kappa_x = k_store_x(i,j,k,ispec_CPML)
        kappa_y = k_store_y(i,j,k,ispec_CPML)
        kappa_z = k_store_z(i,j,k,ispec_CPML)

        d_x = d_store_x(i,j,k,ispec_CPML)
        d_y = d_store_y(i,j,k,ispec_CPML)
        d_z = d_store_z(i,j,k,ispec_CPML)

        alpha_x = alpha_store_x(i,j,k,ispec_CPML)
        alpha_y = alpha_store_y(i,j,k,ispec_CPML)
        alpha_z = alpha_store_z(i,j,k,ispec_CPML)

        call l_parameter_computation(deltat, &
               kappa_x, d_x, alpha_x, &
               kappa_y, d_y, alpha_y, &
               kappa_z, d_z, alpha_z, &
               CPML_region_local, &
               A_0, A_1, A_2, A_3, A_4, A_5, &
               coef0_x, coef1_x, coef2_x, &
               coef0_y, coef1_y, coef2_y, &
               coef0_z, coef1_z, coef2_z)

        ! updates memory variables
        rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) = coef0_x * rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) &
                + PML_displ_new(1,i,j,k,ispec_CPML) * coef1_x + PML_displ_old(1,i,j,k,ispec_CPML) * coef2_x
        rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) = coef0_x * rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) &
                + PML_displ_new(2,i,j,k,ispec_CPML) * coef1_x + PML_displ_old(2,i,j,k,ispec_CPML) * coef2_x
        rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) = coef0_x * rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) &
                + PML_displ_new(3,i,j,k,ispec_CPML) * coef1_x + PML_displ_old(3,i,j,k,ispec_CPML) * coef2_x

        rmemory_displ_elastic(1,i,j,k,ispec_CPML,2) = coef0_y * rmemory_displ_elastic(1,i,j,k,ispec_CPML,2) &
                + PML_displ_new(1,i,j,k,ispec_CPML) * coef1_y + PML_displ_old(1,i,j,k,ispec_CPML) * coef2_y
        rmemory_displ_elastic(2,i,j,k,ispec_CPML,2) = coef0_y * rmemory_displ_elastic(2,i,j,k,ispec_CPML,2) &
                + PML_displ_new(2,i,j,k,ispec_CPML) * coef1_y + PML_displ_old(2,i,j,k,ispec_CPML) * coef2_y
        rmemory_displ_elastic(3,i,j,k,ispec_CPML,2) = coef0_y * rmemory_displ_elastic(3,i,j,k,ispec_CPML,2) &
                + PML_displ_new(3,i,j,k,ispec_CPML) * coef1_y + PML_displ_old(3,i,j,k,ispec_CPML) * coef2_y

        rmemory_displ_elastic(1,i,j,k,ispec_CPML,3) = coef0_z * rmemory_displ_elastic(1,i,j,k,ispec_CPML,3) &
                + PML_displ_new(1,i,j,k,ispec_CPML) * coef1_z + PML_displ_old(1,i,j,k,ispec_CPML) * coef2_z
        rmemory_displ_elastic(2,i,j,k,ispec_CPML,3) = coef0_z * rmemory_displ_elastic(2,i,j,k,ispec_CPML,3) &
                + PML_displ_new(2,i,j,k,ispec_CPML) * coef1_z + PML_displ_old(2,i,j,k,ispec_CPML) * coef2_z
        rmemory_displ_elastic(3,i,j,k,ispec_CPML,3) = coef0_z * rmemory_displ_elastic(3,i,j,k,ispec_CPML,3) &
                + PML_displ_new(3,i,j,k,ispec_CPML) * coef1_z + PML_displ_old(3,i,j,k,ispec_CPML) * coef2_z

        ! updates PML acceleration
        accel_elastic_CPML(1,i,j,k) =  wgllcube * rhol * jacobianl * &
             ( A_1 * veloc(1,iglob) + A_2 * displ(1,iglob) + &
               A_3 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) + &
               A_4 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,2) + &
               A_5 * rmemory_displ_elastic(1,i,j,k,ispec_CPML,3)  &
             )

        accel_elastic_CPML(2,i,j,k) =  wgllcube * rhol * jacobianl * &
             ( A_1 * veloc(2,iglob) + A_2 * displ(2,iglob) + &
               A_3 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) + &
               A_4 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,2) + &
               A_5 * rmemory_displ_elastic(2,i,j,k,ispec_CPML,3)  &
             )

        accel_elastic_CPML(3,i,j,k) =  wgllcube * rhol * jacobianl * &
             ( A_1 * veloc(3,iglob) + A_2 * displ(3,iglob) + &
               A_3 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) + &
               A_4 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,2) + &
               A_5 * rmemory_displ_elastic(3,i,j,k,ispec_CPML,3)  &
             )
      enddo
    enddo
  enddo

end subroutine pml_compute_accel_contribution_elastic

!
!=====================================================================
!
subroutine pml_compute_accel_contribution_acoustic(ispec,ispec_CPML,potential_acoustic, &
                                                   potential_dot_acoustic,rmemory_potential_acoustic)

  ! calculates contribution from each C-PML element to update acceleration to the global mesh

  ! second-order accurate convolution term calculation from equation (21) of
  ! Shumin Wang, Robert Lee, and Fernando L. Teixeira,
  ! Anisotropic-medium PML for vector FETD with modified basis functions,
  ! IEEE Transactions on Antennas and Propagation, vol. 54, no. 1, (2006)

  use specfem_par, only: NGLOB_AB,deltat,wgll_cube,jacobian,ibool,kappastore
  use pml_par, only: CPML_regions,NSPEC_CPML,d_store_x,d_store_y,d_store_z,K_store_x,K_store_y,K_store_z, &
                     alpha_store_x, alpha_store_y, alpha_store_z, &
                     NSPEC_CPML,potential_dot_dot_acoustic_CPML, &
                     PML_potential_acoustic_old,PML_potential_acoustic_new
  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY, &
                       CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ

  implicit none

  integer, intent(in) :: ispec,ispec_CPML
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB), intent(in) :: potential_acoustic,potential_dot_acoustic
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CPML,3) :: rmemory_potential_acoustic


  ! local parameters
  integer :: i,j,k,iglob,CPML_region_local

  real(kind=CUSTOM_REAL) :: wgllcube,kappal_inv,jacobianl
  real(kind=CUSTOM_REAL) :: alpha_x,alpha_y,alpha_z,d_x,d_y,d_z,kappa_x,kappa_y,kappa_z
  real(kind=CUSTOM_REAL) :: coef0_x,coef1_x,coef2_x,coef0_y,coef1_y,coef2_y,coef0_z,coef1_z,coef2_z
  real(kind=CUSTOM_REAL) :: A_0,A_1,A_2,A_3,A_4,A_5

  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        iglob = ibool(i,j,k,ispec)

        wgllcube = wgll_cube(i,j,k)
        jacobianl = jacobian(i,j,k,ispec)
        kappal_inv = 1._CUSTOM_REAL / kappastore(i,j,k,ispec)

        ! PML coefficient values
        CPML_region_local = CPML_regions(ispec_CPML)

        kappa_x = k_store_x(i,j,k,ispec_CPML)
        kappa_y = k_store_y(i,j,k,ispec_CPML)
        kappa_z = k_store_z(i,j,k,ispec_CPML)

        d_x = d_store_x(i,j,k,ispec_CPML)
        d_y = d_store_y(i,j,k,ispec_CPML)
        d_z = d_store_z(i,j,k,ispec_CPML)

        alpha_x = alpha_store_x(i,j,k,ispec_CPML)
        alpha_y = alpha_store_y(i,j,k,ispec_CPML)
        alpha_z = alpha_store_z(i,j,k,ispec_CPML)

        call l_parameter_computation(deltat, &
               kappa_x, d_x, alpha_x, &
               kappa_y, d_y, alpha_y, &
               kappa_z, d_z, alpha_z, &
               CPML_region_local, &
               A_0, A_1, A_2, A_3, A_4, A_5, &
               coef0_x, coef1_x, coef2_x, &
               coef0_y, coef1_y, coef2_y, &
               coef0_z, coef1_z, coef2_z)

        ! updates memory variables
        rmemory_potential_acoustic(i,j,k,ispec_CPML,1) = coef0_x * rmemory_potential_acoustic(i,j,k,ispec_CPML,1) &
                + coef1_x * PML_potential_acoustic_new(i,j,k,ispec_CPML) &
                + coef2_x * PML_potential_acoustic_old(i,j,k,ispec_CPML)

        rmemory_potential_acoustic(i,j,k,ispec_CPML,2) = coef0_y * rmemory_potential_acoustic(i,j,k,ispec_CPML,2) &
                + coef1_y * PML_potential_acoustic_new(i,j,k,ispec_CPML) &
                + coef2_y * PML_potential_acoustic_old(i,j,k,ispec_CPML)

        rmemory_potential_acoustic(i,j,k,ispec_CPML,3) = coef0_z * rmemory_potential_acoustic(i,j,k,ispec_CPML,3) &
                + coef1_z * PML_potential_acoustic_new(i,j,k,ispec_CPML) &
                + coef2_z * PML_potential_acoustic_old(i,j,k,ispec_CPML)

        ! updates PML potential
        potential_dot_dot_acoustic_CPML(i,j,k) =  wgllcube * kappal_inv * jacobianl * &
                  ( A_1 * potential_dot_acoustic(iglob) + A_2 * potential_acoustic(iglob) &
                  + A_3 * rmemory_potential_acoustic(i,j,k,ispec_CPML,1) &
                  + A_4 * rmemory_potential_acoustic(i,j,k,ispec_CPML,2) &
                  + A_5 * rmemory_potential_acoustic(i,j,k,ispec_CPML,3) &
                  )
      enddo
    enddo
  enddo

end subroutine pml_compute_accel_contribution_acoustic
!
!=====================================================================
!
subroutine save_field_on_pml_interface(displ,veloc,accel,nglob_interface_PML_elastic, &
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
subroutine read_field_on_pml_interface(b_accel,b_veloc,b_displ,nglob_interface_PML_elastic, &
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
    do i = 1, NGLLX
      do j = 1, NGLLY
        do k = 1, NGLLZ
          iglob = ibool(i,j,k,ispec)
          b_displ(:,iglob) = 0._CUSTOM_REAL
          b_veloc(:,iglob) = 0._CUSTOM_REAL
          b_accel(:,iglob) = 0._CUSTOM_REAL
        enddo
      enddo
    enddo
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
subroutine save_potential_on_pml_interface(potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
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
subroutine read_potential_on_pml_interface(b_potential_dot_dot_acoustic,b_potential_dot_acoustic,b_potential_acoustic, &
                                           nglob_interface_PML_acoustic,b_PML_potential,b_reclen_PML_potential)

  use specfem_par, only: NGLOB_AB,ibool,NSTEP,it
  use pml_par, only: NSPEC_CPML,CPML_to_spec
  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ
  implicit none

  integer, intent(in) :: nglob_interface_PML_acoustic,b_reclen_PML_potential
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: b_potential_dot_dot_acoustic,b_potential_dot_acoustic,b_potential_acoustic
  real(kind=CUSTOM_REAL), dimension(3,nglob_interface_PML_acoustic) :: b_PML_potential

  integer :: iglob,ispec,ispec_pml,i,j,k

  do ispec_pml = 1, NSPEC_CPML
    ispec = CPML_to_spec(ispec_pml)
    do i = 1, NGLLX
      do j = 1, NGLLY
        do k = 1, NGLLZ
          iglob = ibool(i,j,k,ispec)
          b_potential_acoustic(iglob) = 0._CUSTOM_REAL
          b_potential_dot_acoustic(iglob) = 0._CUSTOM_REAL
          b_potential_dot_dot_acoustic(iglob) = 0._CUSTOM_REAL
        enddo
      enddo
    enddo
  enddo

  call read_abs(1,b_PML_potential,b_reclen_PML_potential,NSTEP-it+1)

  do iglob = 1, nglob_interface_PML_acoustic
    b_potential_acoustic(iglob) = b_PML_potential(1,iglob)
    b_potential_dot_acoustic(iglob) = b_PML_potential(2,iglob)
    b_potential_dot_dot_acoustic(iglob) = b_PML_potential(3,iglob)
  enddo

end subroutine read_potential_on_pml_interface
!
!=====================================================================
!
subroutine l_parameter_computation(deltat, &
               kappa_x, d_x, alpha_x, &
               kappa_y, d_y, alpha_y, &
               kappa_z, d_z, alpha_z, &
               CPML_region_local, &
               A_0, A_1, A_2, A_3, A_4, A_5, &
               coef0_x, coef1_x, coef2_x, &
               coef0_y, coef1_y, coef2_y, &
               coef0_z, coef1_z, coef2_z)

  use constants, only: CUSTOM_REAL, CPML_XYZ, &
                       CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY, &
                       CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY
  use pml_par, only: min_distance_between_CPML_parameter

  implicit none

  real(kind=CUSTOM_REAL), intent(in) :: deltat
  real(kind=CUSTOM_REAL) :: kappa_x,d_x,alpha_x,beta_x, &
                            kappa_y,d_y,alpha_y,beta_y, &
                            kappa_z,d_z,alpha_z,beta_z
  integer, intent(in) :: CPML_region_local

  real(kind=CUSTOM_REAL), intent(out) :: A_0, A_1, A_2, A_3, A_4, A_5
  real(kind=CUSTOM_REAL), intent(out) :: coef0_x, coef1_x, coef2_x, &
                                         coef0_y, coef1_y, coef2_y, &
                                         coef0_z, coef1_z, coef2_z

  !local variable
  real(kind=CUSTOM_REAL) :: bar_A_0, bar_A_1, bar_A_2, bar_A_3, bar_A_4, bar_A_5
  real(kind=CUSTOM_REAL) :: beta_xyz_1, beta_xyz_2, beta_xyz_3

  coef0_x = 0._CUSTOM_REAL
  coef1_x = 0._CUSTOM_REAL
  coef2_x = 0._CUSTOM_REAL
  coef0_y = 0._CUSTOM_REAL
  coef1_y = 0._CUSTOM_REAL
  coef2_y = 0._CUSTOM_REAL
  coef0_z = 0._CUSTOM_REAL
  coef1_z = 0._CUSTOM_REAL
  coef2_z = 0._CUSTOM_REAL

  if (CPML_region_local == CPML_XYZ) then

     if ( abs( alpha_x - alpha_y ) >= min_distance_between_CPML_parameter .and. &
          abs( alpha_x - alpha_z ) >= min_distance_between_CPML_parameter .and. &
          abs( alpha_y - alpha_z ) >= min_distance_between_CPML_parameter) then

       beta_x = alpha_x + d_x / kappa_x
       beta_y = alpha_y + d_y / kappa_y
       beta_z = alpha_z + d_z / kappa_z

       beta_xyz_1 = beta_x + beta_y + beta_z
       beta_xyz_2 = beta_x * beta_y + beta_x * beta_z + beta_y * beta_z
       beta_xyz_3 = beta_x * beta_y * beta_z

       bar_A_0 = kappa_x * kappa_y * kappa_z
       bar_A_1 = bar_A_0 * (beta_x + beta_y + beta_z - alpha_x - alpha_y - alpha_z)
       bar_A_2 = bar_A_0 * (beta_x - alpha_x) * (beta_y - alpha_y - alpha_x) &
               + bar_A_0 * (beta_y - alpha_y) * (beta_z - alpha_z - alpha_y) &
               + bar_A_0 * (beta_z - alpha_z) * (beta_x - alpha_x - alpha_z)

       bar_A_3 = bar_A_0 * alpha_x**2 &
               * (beta_x - alpha_x) * (beta_y - alpha_x) * (beta_z - alpha_x) &
               / ((alpha_y - alpha_x) * (alpha_z - alpha_x))
       bar_A_4 = bar_A_0 * alpha_y**2 &
               * (beta_x - alpha_y) * (beta_y - alpha_y) * (beta_z - alpha_y) &
               / ((alpha_x - alpha_y) * (alpha_z - alpha_y))
       bar_A_5 = bar_A_0 * alpha_z**2 &
               * (beta_x - alpha_z) * (beta_y - alpha_z) * (beta_z - alpha_z) &
               / ((alpha_y - alpha_z) * (alpha_x - alpha_z))

       A_0 = bar_A_0
       A_1 = bar_A_1
       A_2 = bar_A_2
       A_3 = bar_A_3
       A_4 = bar_A_4
       A_5 = bar_A_5

       call compute_convolution_coef(alpha_x, deltat, coef0_x, coef1_x, coef2_x)
       call compute_convolution_coef(alpha_y, deltat, coef0_y, coef1_y, coef2_y)
       call compute_convolution_coef(alpha_z, deltat, coef0_z, coef1_z, coef2_z)

     else
       stop 'error occured in l_parameter_computation in CPML_XYZ region'
     endif

  else if (CPML_region_local == CPML_XY_ONLY) then

    if (abs( alpha_x - alpha_y ) >= min_distance_between_CPML_parameter) then

      beta_x = alpha_x + d_x / kappa_x
      beta_y = alpha_y + d_y / kappa_y
      beta_z = alpha_z + d_z / kappa_z

      beta_xyz_1 = beta_x + beta_y
      beta_xyz_2 = beta_x * beta_y

      bar_A_0 = kappa_x * kappa_y
      bar_A_1 = bar_A_0 * (beta_x + beta_y - alpha_x - alpha_y)
      bar_A_2 = bar_A_0 * (beta_x - alpha_x) * (beta_y - alpha_y - alpha_x) &
              - bar_A_0 * (beta_y - alpha_y) * alpha_y
      bar_A_3 = bar_A_0 * alpha_x**2 &
              * (beta_x - alpha_x) * (beta_y - alpha_x) &
              / (alpha_y - alpha_x)
      bar_A_4 = bar_A_0 * alpha_y**2 &
              * (beta_x - alpha_y) * (beta_y - alpha_y)  &
              / (alpha_x - alpha_y)
      bar_A_5 = 0._CUSTOM_REAL

      A_0 = bar_A_0
      A_1 = bar_A_1
      A_2 = bar_A_2
      A_3 = bar_A_3
      A_4 = bar_A_4
      A_5 = bar_A_5

      call compute_convolution_coef(alpha_x, deltat, coef0_x, coef1_x, coef2_x)
      call compute_convolution_coef(alpha_y, deltat, coef0_y, coef1_y, coef2_y)
      call compute_convolution_coef(alpha_z, deltat, coef0_z, coef1_z, coef2_z)
    else
      stop 'error occured in l_parameter_computation in CPML_XY_ONLY region'
    endif

  else if (CPML_region_local == CPML_XZ_ONLY) then

    if (abs( alpha_x - alpha_z ) >= min_distance_between_CPML_parameter) then

      beta_x = alpha_x + d_x / kappa_x
      beta_y = alpha_y + d_y / kappa_y
      beta_z = alpha_z + d_z / kappa_z

      beta_xyz_1 = beta_x + beta_z
      beta_xyz_2 = beta_x * beta_z

      bar_A_0 = kappa_x * kappa_z
      bar_A_1 = bar_A_0 * (beta_x + beta_z - alpha_x - alpha_z)
      bar_A_2 = bar_A_0 * (beta_x - alpha_x) * (- alpha_x) &
              + bar_A_0 * (beta_z - alpha_z) * (beta_x - alpha_x - alpha_z)

      bar_A_3 = bar_A_0 * alpha_x**2 &
              * (beta_x - alpha_x) * (beta_z - alpha_x) &
              / (alpha_z - alpha_x)
      bar_A_4 = 0._CUSTOM_REAL
      bar_A_5 = bar_A_0 * alpha_z**2 &
              * (beta_x - alpha_z) * (beta_z - alpha_z)  &
              / (alpha_x - alpha_z)

      A_0 = bar_A_0
      A_1 = bar_A_1
      A_2 = bar_A_2
      A_3 = bar_A_3
      A_4 = bar_A_4
      A_5 = bar_A_5

      call compute_convolution_coef(alpha_x, deltat, coef0_x, coef1_x, coef2_x)
      call compute_convolution_coef(alpha_y, deltat, coef0_y, coef1_y, coef2_y)
      call compute_convolution_coef(alpha_z, deltat, coef0_z, coef1_z, coef2_z)

    else
      stop 'error occured in l_parameter_computation in CPML_XZ_ONLY region'
    endif

  else if (CPML_region_local == CPML_YZ_ONLY) then

    if (abs( alpha_y - alpha_z ) >= min_distance_between_CPML_parameter) then

      beta_x = alpha_x + d_x / kappa_x
      beta_y = alpha_y + d_y / kappa_y
      beta_z = alpha_z + d_z / kappa_z

      beta_xyz_1 = beta_y + beta_z
      beta_xyz_2 = beta_y * beta_z

      bar_A_0 = kappa_y * kappa_z
      bar_A_1 = bar_A_0 * (beta_y + beta_z - alpha_y - alpha_z)
      bar_A_2 = bar_A_0 * (beta_y - alpha_y) * (beta_z - alpha_z - alpha_y) &
              - bar_A_0 * (beta_z - alpha_z) * alpha_z
      bar_A_3 = 0._CUSTOM_REAL
      bar_A_4 = bar_A_0 * alpha_y**2 &
              * (beta_y - alpha_y) * (beta_z - alpha_y) &
              / (alpha_z - alpha_y)
      bar_A_5 = bar_A_0 * alpha_z**2 &
              * (beta_y - alpha_z) * (beta_z - alpha_z)  &
              / (alpha_y - alpha_z)

      A_0 = bar_A_0
      A_1 = bar_A_1
      A_2 = bar_A_2
      A_3 = bar_A_3
      A_4 = bar_A_4
      A_5 = bar_A_5

      call compute_convolution_coef(alpha_x, deltat, coef0_x, coef1_x, coef2_x)
      call compute_convolution_coef(alpha_y, deltat, coef0_y, coef1_y, coef2_y)
      call compute_convolution_coef(alpha_z, deltat, coef0_z, coef1_z, coef2_z)

    else
      stop 'error occured in l_parameter_computation in CPML_YZ_ONLY region'
    endif

  else if (CPML_region_local == CPML_X_ONLY) then

    beta_x = alpha_x + d_x / kappa_x
    beta_y = alpha_y + d_y / kappa_y
    beta_z = alpha_z + d_z / kappa_z

    bar_A_0 = kappa_x
    bar_A_1 = bar_A_0 * (beta_x - alpha_x)
    bar_A_2 = - bar_A_0 * alpha_x * (beta_x - alpha_x)
    bar_A_3 = bar_A_0 * alpha_x**2 * (beta_x - alpha_x)
    bar_A_4 = 0._CUSTOM_REAL
    bar_A_5 = 0._CUSTOM_REAL

    A_0 = bar_A_0
    A_1 = bar_A_1
    A_2 = bar_A_2
    A_3 = bar_A_3
    A_4 = bar_A_4
    A_5 = bar_A_5

    call compute_convolution_coef(alpha_x, deltat, coef0_x, coef1_x, coef2_x)
    call compute_convolution_coef(alpha_y, deltat, coef0_y, coef1_y, coef2_y)
    call compute_convolution_coef(alpha_z, deltat, coef0_z, coef1_z, coef2_z)

  else if (CPML_region_local == CPML_Y_ONLY) then

    beta_x = alpha_x + d_x / kappa_x
    beta_y = alpha_y + d_y / kappa_y
    beta_z = alpha_z + d_z / kappa_z

    bar_A_0 = kappa_y
    bar_A_1 = bar_A_0 * (beta_y - alpha_y)
    bar_A_2 = - bar_A_0 * alpha_y * (beta_y - alpha_y)
    bar_A_3 = 0._CUSTOM_REAL
    bar_A_4 = bar_A_0 * alpha_y**2 * (beta_y - alpha_y)
    bar_A_5 = 0._CUSTOM_REAL

    A_0 = bar_A_0
    A_1 = bar_A_1
    A_2 = bar_A_2
    A_3 = bar_A_3
    A_4 = bar_A_4
    A_5 = bar_A_5

    call compute_convolution_coef(alpha_x, deltat, coef0_x, coef1_x, coef2_x)
    call compute_convolution_coef(alpha_y, deltat, coef0_y, coef1_y, coef2_y)
    call compute_convolution_coef(alpha_z, deltat, coef0_z, coef1_z, coef2_z)

  else if (CPML_region_local == CPML_Z_ONLY) then

    beta_x = alpha_x + d_x / kappa_x
    beta_y = alpha_y + d_y / kappa_y
    beta_z = alpha_z + d_z / kappa_z

    bar_A_0 = kappa_z
    bar_A_1 = bar_A_0 * (beta_z - alpha_z)
    bar_A_2 = - bar_A_0 * alpha_z * (beta_z - alpha_z)
    bar_A_3 = 0._CUSTOM_REAL
    bar_A_4 = 0._CUSTOM_REAL
    bar_A_5 = bar_A_0 * alpha_z**2 * (beta_z - alpha_z)

    A_0 = bar_A_0
    A_1 = bar_A_1
    A_2 = bar_A_2
    A_3 = bar_A_3
    A_4 = bar_A_4
    A_5 = bar_A_5

    call compute_convolution_coef(alpha_x, deltat, coef0_x, coef1_x, coef2_x)
    call compute_convolution_coef(alpha_y, deltat, coef0_y, coef1_y, coef2_y)
    call compute_convolution_coef(alpha_z, deltat, coef0_z, coef1_z, coef2_z)

  endif

end subroutine l_parameter_computation

!
!=====================================================================
!
subroutine compute_convolution_coef(bb,deltat,coef0,coef1,coef2)

  use constants, only: CUSTOM_REAL
  use pml_par, only: min_distance_between_CPML_parameter

  implicit none

  real(kind=CUSTOM_REAL),intent(in) :: bb, deltat
  real(kind=CUSTOM_REAL),intent(out) :: coef0, coef1, coef2

  ! local parameters
  logical,parameter :: FIRST_ORDER_CONVOLUTION = .false.
  real(kind=CUSTOM_REAL) :: bbpow2,bbpow3
  real(kind=CUSTOM_REAL) :: deltatpow2,deltatpow3,deltatpow4,deltatpow5,deltatpow6,deltat_half
  real(kind=CUSTOM_REAL) :: prod1,prod1_half

  ! permanent factors (avoids divisions which are computationally expensive)
  real(kind=CUSTOM_REAL),parameter :: ONE_OVER_8 = 0.125_CUSTOM_REAL
  real(kind=CUSTOM_REAL),parameter :: ONE_OVER_12 = 1._CUSTOM_REAL / 12._CUSTOM_REAL
  real(kind=CUSTOM_REAL),parameter :: ONE_OVER_24 = 1._CUSTOM_REAL / 24._CUSTOM_REAL
  real(kind=CUSTOM_REAL),parameter :: ONE_OVER_48 = 1._CUSTOM_REAL / 48._CUSTOM_REAL
  real(kind=CUSTOM_REAL),parameter :: ONE_OVER_128 = 0.0078125_CUSTOM_REAL
  real(kind=CUSTOM_REAL),parameter :: ONE_OVER_384 = 1._CUSTOM_REAL / 384._CUSTOM_REAL
  real(kind=CUSTOM_REAL),parameter :: ONE_OVER_960 = 1._CUSTOM_REAL / 960._CUSTOM_REAL
  real(kind=CUSTOM_REAL),parameter :: ONE_OVER_1920 = 1._CUSTOM_REAL / 1920._CUSTOM_REAL
  real(kind=CUSTOM_REAL),parameter :: SEVEN_OVER_3840 = 7._CUSTOM_REAL / 3840._CUSTOM_REAL
  real(kind=CUSTOM_REAL),parameter :: FIVE_OVER_11520 = 5._CUSTOM_REAL/11520._CUSTOM_REAL

  ! helper variables
  bbpow2 = bb**2
  bbpow3 = bb**3

  deltatpow2 = deltat**2
  deltatpow3 = deltat**3
  deltatpow4 = deltat**4
  deltatpow5 = deltat**5
  deltatpow6 = deltat**6
  deltat_half = deltat * 0.5_CUSTOM_REAL

  prod1 = bb * deltat
  prod1_half = prod1 * 0.5_CUSTOM_REAL

  ! calculates coefficients
  coef0 = exp(-prod1)

  if (abs(bb) >= min_distance_between_CPML_parameter) then
    if (FIRST_ORDER_CONVOLUTION) then
      coef1 = (1._CUSTOM_REAL - exp(-prod1) ) / bb
      coef2 = 0._CUSTOM_REAL
    else
      coef1 = (1._CUSTOM_REAL - exp(-prod1_half) ) / bb
      coef2 = (1._CUSTOM_REAL - exp(-prod1_half) ) * exp(-prod1_half) / bb
    endif
  else
    if (FIRST_ORDER_CONVOLUTION) then
      coef1 = deltat
      coef2 = 0._CUSTOM_REAL
    else
      coef1 = deltat_half + &
              (- deltatpow2*bb*ONE_OVER_8 + &
               (deltatpow3*bbpow2*ONE_OVER_48 - &
                deltatpow4*bbpow3*ONE_OVER_384))
      coef2 = deltat_half + &
              (- 3._CUSTOM_REAL*deltatpow2*bb*ONE_OVER_8 + &
               (7._CUSTOM_REAL*deltatpow3*bbpow2*ONE_OVER_48 - &
                5._CUSTOM_REAL*deltatpow4*bbpow3*ONE_OVER_128))
    endif
  endif

end subroutine compute_convolution_coef
