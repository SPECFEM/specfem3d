!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
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
  ! Anisotropic-Medium PML for Vector FETD With Modified Basis Functions,
  ! IEEE Transactions on Antennas and Propagation, vol. 54, no. 1, (2006)

  use specfem_par, only: NGLOB_AB,it,deltat,wgll_cube,jacobian,ibool,rhostore
  use pml_par, only: CPML_regions,d_store_x,d_store_y,d_store_z,K_store_x,K_store_y,K_store_z,&
                     alpha_store_x, alpha_store_y, alpha_store_z, &
                     NSPEC_CPML,accel_elastic_CPML,displ_old,displ_new
  use constants, only: CUSTOM_REAL,NDIM,NGLLX,NGLLY,NGLLZ,CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY, &
                       CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ

  implicit none

  integer, intent(in) :: ispec,ispec_CPML
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB), intent(in) :: displ,veloc
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ,NSPEC_CPML,3) :: rmemory_displ_elastic

  ! local parameters
  integer :: i,j,k,iglob,CPML_region_local
  integer :: singularity_type_4, singularity_type_5
  real(kind=CUSTOM_REAL) :: wgllcube,rhol,jacobianl
  real(kind=CUSTOM_REAL) :: alpha_x,alpha_y,alpha_z,d_x,d_y,d_z,kappa_x,kappa_y,kappa_z
  real(kind=CUSTOM_REAL) :: coef0_x,coef1_x,coef2_x,coef0_y,coef1_y,coef2_y,coef0_z,coef1_z,coef2_z
  real(kind=CUSTOM_REAL) :: A_0,A_1,A_2,A_3,A_4,A_5
  real(kind=CUSTOM_REAL) :: time_nplus1

  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        rhol = rhostore(i,j,k,ispec)
        jacobianl = jacobian(i,j,k,ispec)
        iglob = ibool(i,j,k,ispec)
        wgllcube = wgll_cube(i,j,k)

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

        time_nplus1 = (it - 1._CUSTOM_REAL) * deltat

        call l_parameter_computation( &
               time_nplus1, deltat, &
               kappa_x, d_x, alpha_x, &
               kappa_y, d_y, alpha_y, &
               kappa_z, d_z, alpha_z, &
               CPML_region_local,  &
               A_0, A_1, A_2, A_3, A_4, A_5, &
               coef0_x, coef1_x, coef2_x, &
               coef0_y, coef1_y, coef2_y, &
               coef0_z, coef1_z, coef2_z, &
               singularity_type_4, singularity_type_5)

        rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) = coef0_x * rmemory_displ_elastic(1,i,j,k,ispec_CPML,1) &
                + displ_new(1,iglob) * coef1_x + displ_old(1,iglob) * coef2_x
        rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) = coef0_x * rmemory_displ_elastic(2,i,j,k,ispec_CPML,1) &
                + displ_new(2,iglob) * coef1_x + displ_old(2,iglob) * coef2_x
        rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) = coef0_x * rmemory_displ_elastic(3,i,j,k,ispec_CPML,1) &
                + displ_new(3,iglob) * coef1_x + displ_old(3,iglob) * coef2_x

        rmemory_displ_elastic(1,i,j,k,ispec_CPML,2) = coef0_y * rmemory_displ_elastic(1,i,j,k,ispec_CPML,2) &
                + displ_new(1,iglob) * coef1_y + displ_old(1,iglob) * coef2_y
        rmemory_displ_elastic(2,i,j,k,ispec_CPML,2) = coef0_y * rmemory_displ_elastic(2,i,j,k,ispec_CPML,2) &
                + displ_new(2,iglob) * coef1_y + displ_old(2,iglob) * coef2_y
        rmemory_displ_elastic(3,i,j,k,ispec_CPML,2) = coef0_y * rmemory_displ_elastic(3,i,j,k,ispec_CPML,2) &
                + displ_new(3,iglob) * coef1_y + displ_old(3,iglob) * coef2_y

        rmemory_displ_elastic(1,i,j,k,ispec_CPML,3) = coef0_z * rmemory_displ_elastic(1,i,j,k,ispec_CPML,3) &
                + displ_new(1,iglob) * coef1_z + displ_old(1,iglob) * coef2_z
        rmemory_displ_elastic(2,i,j,k,ispec_CPML,3) = coef0_z * rmemory_displ_elastic(2,i,j,k,ispec_CPML,3) &
                + displ_new(2,iglob) * coef1_z + displ_old(2,iglob) * coef2_z
        rmemory_displ_elastic(3,i,j,k,ispec_CPML,3) = coef0_z * rmemory_displ_elastic(3,i,j,k,ispec_CPML,3) &
                + displ_new(3,iglob) * coef1_z + displ_old(3,iglob) * coef2_z

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
subroutine pml_compute_accel_contribution_acoustic(ispec,ispec_CPML,potential_acoustic,&
                                                   potential_dot_acoustic,rmemory_potential_acoustic)

  ! calculates contribution from each C-PML element to update acceleration to the global mesh

  ! second-order accurate convolution term calculation from equation (21) of
  ! Shumin Wang, Robert Lee, and Fernando L. Teixeira,
  ! Anisotropic-Medium PML for Vector FETD With Modified Basis Functions,
  ! IEEE Transactions on Antennas and Propagation, vol. 54, no. 1, (2006)

  use specfem_par, only: NGLOB_AB,it,deltat,wgll_cube,jacobian,ibool,kappastore
  use pml_par, only: CPML_regions,NSPEC_CPML,d_store_x,d_store_y,d_store_z,K_store_x,K_store_y,K_store_z,&
                     alpha_store_x, alpha_store_y, alpha_store_z, &
                     NSPEC_CPML,potential_dot_dot_acoustic_CPML,potential_acoustic_old,potential_acoustic_new
  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY, &
                       CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ

  implicit none

  integer, intent(in) :: ispec,ispec_CPML
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB), intent(in) :: potential_acoustic,potential_dot_acoustic
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CPML,3) :: rmemory_potential_acoustic


  ! local parameters
  integer :: i,j,k,iglob,CPML_region_local
  integer :: singularity_type_4, singularity_type_5
  real(kind=CUSTOM_REAL) :: wgllcube,kappal_inv,jacobianl
  real(kind=CUSTOM_REAL) :: alpha_x,alpha_y,alpha_z,d_x,d_y,d_z,kappa_x,kappa_y,kappa_z
  real(kind=CUSTOM_REAL) :: coef0_x,coef1_x,coef2_x,coef0_y,coef1_y,coef2_y,coef0_z,coef1_z,coef2_z
  real(kind=CUSTOM_REAL) :: A_0,A_1,A_2,A_3,A_4,A_5
  real(kind=CUSTOM_REAL) :: time_nplus1

  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        kappal_inv = 1._CUSTOM_REAL / kappastore(i,j,k,ispec)
        jacobianl = jacobian(i,j,k,ispec)
        iglob = ibool(i,j,k,ispec)
        wgllcube = wgll_cube(i,j,k)

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

        time_nplus1 = (it - 1._CUSTOM_REAL) * deltat

        call l_parameter_computation( &
               time_nplus1, deltat, &
               kappa_x, d_x, alpha_x, &
               kappa_y, d_y, alpha_y, &
               kappa_z, d_z, alpha_z, &
               CPML_region_local,  &
               A_0, A_1, A_2, A_3, A_4, A_5, &
               coef0_x, coef1_x, coef2_x, &
               coef0_y, coef1_y, coef2_y, &
               coef0_z, coef1_z, coef2_z, &
               singularity_type_4, singularity_type_5)

        rmemory_potential_acoustic(i,j,k,ispec_CPML,1) = coef0_x * rmemory_potential_acoustic(i,j,k,ispec_CPML,1) &
                + coef1_x * potential_acoustic_new(iglob) + coef2_x * potential_acoustic_old(iglob)

        rmemory_potential_acoustic(i,j,k,ispec_CPML,2) = coef0_y * rmemory_potential_acoustic(i,j,k,ispec_CPML,2) &
                + coef1_y * potential_acoustic_new(iglob) + coef2_y * potential_acoustic_old(iglob)

        rmemory_potential_acoustic(i,j,k,ispec_CPML,3) = coef0_z * rmemory_potential_acoustic(i,j,k,ispec_CPML,3) &
                + coef1_z * potential_acoustic_new(iglob) + coef2_z * potential_acoustic_old(iglob)


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
subroutine l_parameter_computation( &
               timeval, deltat, &
               kappa_x, d_x, alpha_x, &
               kappa_y, d_y, alpha_y, &
               kappa_z, d_z, alpha_z, &
               CPML_region_local, &
               A_0, A_1, A_2, A_3, A_4, A_5, &
               coef0_x, coef1_x, coef2_x, &
               coef0_y, coef1_y, coef2_y, &
               coef0_z, coef1_z, coef2_z, &
               singularity_type_4, singularity_type_5)

  use constants, only: CUSTOM_REAL, CPML_XYZ, &
                       CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY, &
                       CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY

  implicit none

  real(kind=CUSTOM_REAL), intent(in) :: timeval,deltat
  real(kind=CUSTOM_REAL) :: kappa_x,d_x,alpha_x,beta_x, &
                            kappa_y,d_y,alpha_y,beta_y, &
                            kappa_z,d_z,alpha_z,beta_z
  integer, intent(in) :: CPML_region_local

  real(kind=CUSTOM_REAL), intent(out) :: A_0, A_1, A_2, A_3, A_4, A_5
  real(kind=CUSTOM_REAL), intent(out) :: coef0_x, coef1_x, coef2_x, &
                                         coef0_y, coef1_y, coef2_y, &
                                         coef0_z, coef1_z, coef2_z
  integer, intent(out) :: singularity_type_4, singularity_type_5

  !local variable
  real(kind=CUSTOM_REAL) :: time_nplus1, time_n
  real(kind=CUSTOM_REAL) :: bar_A_0, bar_A_1, bar_A_2, bar_A_3, bar_A_4, bar_A_5
  real(kind=CUSTOM_REAL) :: alpha_0, beta_xyz_1, beta_xyz_2, beta_xyz_3

  time_nplus1 = timeval
  time_n = timeval - deltat

  coef0_x = 0._CUSTOM_REAL
  coef1_x = 0._CUSTOM_REAL
  coef2_x = 0._CUSTOM_REAL
  coef0_y = 0._CUSTOM_REAL
  coef1_y = 0._CUSTOM_REAL
  coef2_y = 0._CUSTOM_REAL
  coef0_z = 0._CUSTOM_REAL
  coef1_z = 0._CUSTOM_REAL
  coef2_z = 0._CUSTOM_REAL

  if (CPML_region_local == CPML_XYZ ) then

     if ( &
          abs( alpha_x - alpha_y ) >= 1.e-5_CUSTOM_REAL .AND. &
          abs( alpha_x - alpha_z ) >= 1.e-5_CUSTOM_REAL .AND. &
          abs( alpha_y - alpha_z ) >= 1.e-5_CUSTOM_REAL &
        ) then

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

       singularity_type_4 = 0  ! 0 means no singularity, 1 means first order singularity, 2 means second order singularity
       singularity_type_5 = 0  ! 0 means no singularity, 1 means first order singularity, 2 means second order singularity

       call compute_convolution_coef(alpha_x, deltat, coef0_x, coef1_x, coef2_x, 0, time_nplus1, time_n)
       call compute_convolution_coef(alpha_y, deltat, coef0_y, coef1_y, coef2_y, singularity_type_4, time_nplus1, time_n)
       call compute_convolution_coef(alpha_z, deltat, coef0_z, coef1_z, coef2_z, singularity_type_5, time_nplus1, time_n)

     else if (&
             abs( alpha_x - alpha_y ) < 1.e-5_CUSTOM_REAL  .AND. &
             abs( alpha_x - alpha_z ) >= 1.e-5_CUSTOM_REAL .AND. &
             abs( alpha_y - alpha_z ) >= 1.e-5_CUSTOM_REAL &
            ) then

       alpha_0 = max(alpha_x,alpha_y)

       alpha_x = alpha_0
       alpha_y = alpha_0

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
       bar_A_3 = bar_A_0 * alpha_0 / (alpha_z - alpha_0)**2 * ( &
               - alpha_0**3 * (4._CUSTOM_REAL * alpha_0 - 5._CUSTOM_REAL * alpha_z) &
               + alpha_0**2 * (3._CUSTOM_REAL * alpha_0 - 4._CUSTOM_REAL * alpha_z) * beta_xyz_1 &
               - alpha_0 * (2._CUSTOM_REAL * alpha_0 - 3._CUSTOM_REAL * alpha_z) * beta_xyz_2 &
               + (alpha_0 - 2._CUSTOM_REAL * alpha_z) * beta_xyz_3 )
       bar_A_4 = bar_A_0 * alpha_0**2 &
               * (beta_x - alpha_0) * (beta_y - alpha_0) * (beta_z - alpha_0) &
               / (alpha_z - alpha_0)
       bar_A_5 = bar_A_0 * alpha_z**2 &
               * (beta_x - alpha_z) * (beta_y - alpha_z) * (beta_z - alpha_z) &
               / ((alpha_0 - alpha_z) * (alpha_0 - alpha_z))

       A_0 = bar_A_0
       A_1 = bar_A_1
       A_2 = bar_A_2
       A_3 = bar_A_3 + time_nplus1 * bar_A_4
       A_4 = - bar_A_4
       A_5 = bar_A_5

       singularity_type_4 = 1  ! 0 means no singularity, 1 means first order singularity, 2 means second order singularity
       singularity_type_5 = 0  ! 0 means no singularity, 1 means first order singularity, 2 means second order singularity

       call compute_convolution_coef(alpha_x, deltat, coef0_x, coef1_x, coef2_x, 0, time_nplus1, time_n)
       call compute_convolution_coef(alpha_y, deltat, coef0_y, coef1_y, coef2_y, singularity_type_4, time_nplus1, time_n)
       call compute_convolution_coef(alpha_z, deltat, coef0_z, coef1_z, coef2_z, singularity_type_5, time_nplus1, time_n)

     else if (&
             abs( alpha_x - alpha_y ) >= 1.e-5_CUSTOM_REAL .AND. &
             abs( alpha_x - alpha_z ) < 1.e-5_CUSTOM_REAL  .AND. &
             abs( alpha_y - alpha_z ) >= 1.e-5_CUSTOM_REAL &
            )then

       alpha_0 = max(alpha_x,alpha_z)

       alpha_x = alpha_0
       alpha_z = alpha_0

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

       bar_A_3 = bar_A_0 * alpha_0 / (alpha_y - alpha_0)**2 * ( &
               - alpha_0**3 * (4._CUSTOM_REAL * alpha_0 - 5._CUSTOM_REAL * alpha_y) &
               + alpha_0**2 * (3._CUSTOM_REAL * alpha_0 - 4._CUSTOM_REAL * alpha_y) * beta_xyz_1 &
               - alpha_0 * (2._CUSTOM_REAL * alpha_0 - 3._CUSTOM_REAL * alpha_y) * beta_xyz_2 &
               + (alpha_0 - 2._CUSTOM_REAL * alpha_y) * beta_xyz_3 )
       bar_A_4 = bar_A_0 * alpha_y**2 &
               * (beta_x - alpha_y) * (beta_y - alpha_y) * (beta_z - alpha_y) &
               / ((alpha_0 - alpha_y) * (alpha_0 - alpha_y))
       bar_A_5 = bar_A_0 * alpha_0**2 &
               * (beta_x - alpha_0) * (beta_y - alpha_0) * (beta_z - alpha_0) &
               / (alpha_y - alpha_0)

       A_0 = bar_A_0
       A_1 = bar_A_1
       A_2 = bar_A_2
       A_3 = bar_A_3 + time_nplus1 * bar_A_5
       A_4 = bar_A_4
       A_5 = - bar_A_5

       singularity_type_4 = 0
       singularity_type_5 = 1

       call compute_convolution_coef(alpha_x, deltat, coef0_x, coef1_x, coef2_x, 0, time_nplus1, time_n)
       call compute_convolution_coef(alpha_y, deltat, coef0_y, coef1_y, coef2_y, singularity_type_4, time_nplus1, time_n)
       call compute_convolution_coef(alpha_z, deltat, coef0_z, coef1_z, coef2_z, singularity_type_5, time_nplus1, time_n)

     else if (&
             abs( alpha_x - alpha_y ) >= 1.e-5_CUSTOM_REAL .AND. &
             abs( alpha_x - alpha_z ) >= 1.e-5_CUSTOM_REAL .AND. &
             abs( alpha_y - alpha_z ) < 1.e-5_CUSTOM_REAL  &
            )then

       alpha_0 = max(alpha_y,alpha_z)

       alpha_y = alpha_0
       alpha_z = alpha_0

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
               / ((alpha_0 - alpha_x) * (alpha_0 - alpha_x))
       bar_A_4 = bar_A_0 * alpha_0 / (alpha_x - alpha_0)**2 * ( &
               - alpha_0**3 * (4._CUSTOM_REAL * alpha_0 - 5._CUSTOM_REAL * alpha_x) &
               + alpha_0**2 * (3._CUSTOM_REAL * alpha_0 - 4._CUSTOM_REAL * alpha_x) * beta_xyz_1 &
               - alpha_0 * (2._CUSTOM_REAL * alpha_0 - 3._CUSTOM_REAL * alpha_x ) * beta_xyz_2 &
               + (alpha_0 - 2._CUSTOM_REAL * alpha_x) * beta_xyz_3 )
       bar_A_5 = bar_A_0 * alpha_0**2 &
               * (beta_x - alpha_0) * (beta_y - alpha_0) * (beta_z - alpha_0) &
               / (alpha_x - alpha_0)

       A_0 = bar_A_0
       A_1 = bar_A_1
       A_2 = bar_A_2
       A_3 = bar_A_3
       A_4 = bar_A_4 + time_nplus1 * bar_A_5
       A_5 = - bar_A_5

       singularity_type_4 = 0
       singularity_type_5 = 1

       call compute_convolution_coef(alpha_x, deltat, coef0_x, coef1_x, coef2_x, 0, time_nplus1, time_n)
       call compute_convolution_coef(alpha_y, deltat, coef0_y, coef1_y, coef2_y, singularity_type_4, time_nplus1, time_n)
       call compute_convolution_coef(alpha_z, deltat, coef0_z, coef1_z, coef2_z, singularity_type_5, time_nplus1, time_n)

     else if (&
             (abs( alpha_x - alpha_y ) < 1.e-5_CUSTOM_REAL .AND. &
              abs( alpha_x - alpha_z ) < 1.e-5_CUSTOM_REAL .AND. &
              abs( alpha_y - alpha_z ) < 1.e-5_CUSTOM_REAL).or.  &

             (abs( alpha_x - alpha_y ) >= 1.e-5_CUSTOM_REAL .AND. &
              abs( alpha_x - alpha_z ) < 1.e-5_CUSTOM_REAL  .AND. &
              abs( alpha_y - alpha_z ) < 1.e-5_CUSTOM_REAL) .or.  &

             (abs( alpha_x - alpha_y ) < 1.e-5_CUSTOM_REAL  .AND. &
              abs( alpha_x - alpha_z ) >= 1.e-5_CUSTOM_REAL .AND. &
              abs( alpha_y - alpha_z ) < 1.e-5_CUSTOM_REAL) .or.  &

             (abs( alpha_x - alpha_y ) < 1.e-5_CUSTOM_REAL  .AND. &
              abs( alpha_x - alpha_z ) < 1.e-5_CUSTOM_REAL  .AND. &
              abs( alpha_y - alpha_z ) >= 1.e-5_CUSTOM_REAL)) then

       alpha_0 = max(alpha_x,alpha_y,alpha_z)

       alpha_x = alpha_0
       alpha_y = alpha_0
       alpha_z = alpha_0

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
       bar_A_3 = bar_A_0 * ( &
               - 10._CUSTOM_REAL * alpha_0**3 &
               +  6._CUSTOM_REAL * alpha_0**2 * beta_xyz_1 &
               -  3._CUSTOM_REAL * alpha_0 * beta_xyz_2 &
               + beta_xyz_3 )
       bar_A_4 = bar_A_0 * alpha_0 * ( &
                 5._CUSTOM_REAL * alpha_0**3 &
               - 4._CUSTOM_REAL * alpha_0**2 * beta_xyz_1 &
               + 3._CUSTOM_REAL * alpha_0 * beta_xyz_2 &
               - 2._CUSTOM_REAL * beta_xyz_3 )
       bar_A_5 = bar_A_0 * alpha_0**2 / 2._CUSTOM_REAL &
               * (beta_x - alpha_0) * (beta_y - alpha_0) * (beta_z - alpha_0)

       A_0 = bar_A_0
       A_1 = bar_A_1
       A_2 = bar_A_2
       A_3 = bar_A_3 + time_nplus1 * bar_A_4 + time_nplus1**2 * bar_A_5
       A_4 = - bar_A_4 - 2._CUSTOM_REAL * time_nplus1 * bar_A_5
       A_5 = bar_A_5

       singularity_type_4 = 1
       singularity_type_5 = 2

       call compute_convolution_coef(alpha_x, deltat, coef0_x, coef1_x, coef2_x, 0, time_nplus1, time_n)
       call compute_convolution_coef(alpha_y, deltat, coef0_y, coef1_y, coef2_y, singularity_type_4, time_nplus1, time_n)
       call compute_convolution_coef(alpha_z, deltat, coef0_z, coef1_z, coef2_z, singularity_type_5, time_nplus1, time_n)
     else
       stop 'error occured in l_parameter_computation in CPML_XYZ region'
     endif

  else if ( CPML_region_local == CPML_XY_ONLY ) then

    if ( abs( alpha_x - alpha_y ) >= 1.e-5_CUSTOM_REAL ) then

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

      singularity_type_4 = 0
      singularity_type_5 = 0

      call compute_convolution_coef(alpha_x, deltat, coef0_x, coef1_x, coef2_x, 0, time_nplus1, time_n)
      call compute_convolution_coef(alpha_y, deltat, coef0_y, coef1_y, coef2_y, singularity_type_4, time_nplus1, time_n)
      call compute_convolution_coef(alpha_z, deltat, coef0_z, coef1_z, coef2_z, singularity_type_5, time_nplus1, time_n)

    else if ( abs( alpha_x - alpha_y ) < 1.e-5_CUSTOM_REAL ) then

      alpha_0 = max(alpha_x,alpha_y)

      alpha_x = alpha_0
      alpha_y = alpha_0

      beta_x = alpha_x + d_x / kappa_x
      beta_y = alpha_y + d_y / kappa_y
      beta_z = alpha_z + d_z / kappa_z

      beta_xyz_1 = beta_x + beta_y
      beta_xyz_2 = beta_x * beta_y

      bar_A_0 = kappa_x * kappa_y
      bar_A_1 = bar_A_0 * (beta_x + beta_y - alpha_x - alpha_y)
      bar_A_2 = bar_A_0 * (beta_x - alpha_x) * (beta_y - alpha_y - alpha_x) &
              - bar_A_0 * (beta_y - alpha_y) * alpha_y
      bar_A_3 = bar_A_0 * ( &
              - 4._CUSTOM_REAL * alpha_0**3  &
              + 3._CUSTOM_REAL * alpha_0**2 * beta_xyz_1 &
              - 2._CUSTOM_REAL * alpha_0 * beta_xyz_2 )
      bar_A_4 = bar_A_0 * alpha_0**2 * (beta_x - alpha_0) * (beta_y - alpha_0)
      bar_A_5 = 0._CUSTOM_REAL

      A_0 = bar_A_0
      A_1 = bar_A_1
      A_2 = bar_A_2
      A_3 = bar_A_3 + time_nplus1 * bar_A_4
      A_4 = -bar_A_4
      A_5 = bar_A_5

      singularity_type_4 = 1
      singularity_type_5 = 0

      call compute_convolution_coef(alpha_x, deltat, coef0_x, coef1_x, coef2_x, 0, time_nplus1, time_n)
      call compute_convolution_coef(alpha_y, deltat, coef0_y, coef1_y, coef2_y, singularity_type_4, time_nplus1, time_n)
      call compute_convolution_coef(alpha_z, deltat, coef0_z, coef1_z, coef2_z, singularity_type_5, time_nplus1, time_n)

    else
      stop 'error occured in l_parameter_computation in CPML_XY_ONLY region'
    endif

  else if ( CPML_region_local == CPML_XZ_ONLY ) then

    if ( abs( alpha_x - alpha_z ) >= 1.e-5_CUSTOM_REAL ) then

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

      singularity_type_4 = 0
      singularity_type_5 = 0

      call compute_convolution_coef(alpha_x, deltat, coef0_x, coef1_x, coef2_x, 0, time_nplus1, time_n)
      call compute_convolution_coef(alpha_y, deltat, coef0_y, coef1_y, coef2_y, singularity_type_4, time_nplus1, time_n)
      call compute_convolution_coef(alpha_z, deltat, coef0_z, coef1_z, coef2_z, singularity_type_5, time_nplus1, time_n)

    else if ( abs( alpha_x - alpha_z ) < 1.e-5_CUSTOM_REAL ) then

      alpha_0 = max(alpha_x,alpha_z)

      alpha_x = alpha_0
      alpha_z = alpha_0

      beta_x = alpha_x + d_x / kappa_x
      beta_y = alpha_y + d_y / kappa_y
      beta_z = alpha_z + d_z / kappa_z

      beta_xyz_1 = beta_x + beta_z
      beta_xyz_2 = beta_x * beta_z

      bar_A_0 = kappa_x * kappa_z
      bar_A_1 = bar_A_0 * (beta_x + beta_z - alpha_x - alpha_z)
      bar_A_2 = bar_A_0 * (beta_x - alpha_x) * (- alpha_x) &
              + bar_A_0 * (beta_z - alpha_z) * (beta_x - alpha_x - alpha_z)

      bar_A_3 = bar_A_0 * ( &
              - 4._CUSTOM_REAL * alpha_0**3  &
              + 3._CUSTOM_REAL * alpha_0**2 * beta_xyz_1 &
              - 2._CUSTOM_REAL * alpha_0 * beta_xyz_2 )
      bar_A_4 = 0._CUSTOM_REAL
      bar_A_5 = bar_A_0 * alpha_0**2 * (beta_x - alpha_0) * (beta_z - alpha_0)

      A_0 = bar_A_0
      A_1 = bar_A_1
      A_2 = bar_A_2
      A_3 = bar_A_3 + time_nplus1 * bar_A_5
      A_4 = bar_A_4
      A_5 = -bar_A_5

      singularity_type_4 = 0
      singularity_type_5 = 1

      call compute_convolution_coef(alpha_x, deltat, coef0_x, coef1_x, coef2_x, 0, time_nplus1, time_n)
      call compute_convolution_coef(alpha_y, deltat, coef0_y, coef1_y, coef2_y, singularity_type_4, time_nplus1, time_n)
      call compute_convolution_coef(alpha_z, deltat, coef0_z, coef1_z, coef2_z, singularity_type_5, time_nplus1, time_n)

    else
      stop 'error occured in l_parameter_computation in CPML_XZ_ONLY region'
    endif

  else if ( CPML_region_local == CPML_YZ_ONLY ) then

    if ( abs( alpha_y - alpha_z ) >= 1.e-5_CUSTOM_REAL ) then

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

      singularity_type_4 = 0
      singularity_type_5 = 0

      call compute_convolution_coef(alpha_x, deltat, coef0_x, coef1_x, coef2_x, 0, time_nplus1, time_n)
      call compute_convolution_coef(alpha_y, deltat, coef0_y, coef1_y, coef2_y, singularity_type_4, time_nplus1, time_n)
      call compute_convolution_coef(alpha_z, deltat, coef0_z, coef1_z, coef2_z, singularity_type_5, time_nplus1, time_n)

    else if ( abs( alpha_y - alpha_z ) < 1.e-5_CUSTOM_REAL ) then

      alpha_0 = max(alpha_y,alpha_z)

      alpha_y = alpha_0
      alpha_z = alpha_0

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
      bar_A_4 = bar_A_0 * ( &
              - 4._CUSTOM_REAL * alpha_0**3  &
              + 3._CUSTOM_REAL * alpha_0**2 * beta_xyz_1 &
              - 2._CUSTOM_REAL * alpha_0 * beta_xyz_2 )
      bar_A_5 = bar_A_0 * alpha_0**2 * (beta_y - alpha_0) * (beta_z - alpha_0)

      A_0 = bar_A_0
      A_1 = bar_A_1
      A_2 = bar_A_2
      A_3 = bar_A_3
      A_4 = bar_A_4 + time_nplus1 * bar_A_5
      A_5 = -bar_A_5

      singularity_type_4 = 0
      singularity_type_5 = 1

      call compute_convolution_coef(alpha_x, deltat, coef0_x, coef1_x, coef2_x, 0, time_nplus1, time_n)
      call compute_convolution_coef(alpha_y, deltat, coef0_y, coef1_y, coef2_y, singularity_type_4, time_nplus1, time_n)
      call compute_convolution_coef(alpha_z, deltat, coef0_z, coef1_z, coef2_z, singularity_type_5, time_nplus1, time_n)

    else
      stop 'error occured in l_parameter_computation in CPML_YZ_ONLY region'
    endif

  else if ( CPML_region_local == CPML_X_ONLY ) then

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

    singularity_type_4 = 0
    singularity_type_5 = 0

    call compute_convolution_coef(alpha_x, deltat, coef0_x, coef1_x, coef2_x, 0, time_nplus1, time_n)
    call compute_convolution_coef(alpha_y, deltat, coef0_y, coef1_y, coef2_y, singularity_type_4, time_nplus1, time_n)
    call compute_convolution_coef(alpha_z, deltat, coef0_z, coef1_z, coef2_z, singularity_type_5, time_nplus1, time_n)

  else if ( CPML_region_local == CPML_Y_ONLY ) then

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

    singularity_type_4 = 0
    singularity_type_5 = 0

    call compute_convolution_coef(alpha_x, deltat, coef0_x, coef1_x, coef2_x, 0, time_nplus1, time_n)
    call compute_convolution_coef(alpha_y, deltat, coef0_y, coef1_y, coef2_y, singularity_type_4, time_nplus1, time_n)
    call compute_convolution_coef(alpha_z, deltat, coef0_z, coef1_z, coef2_z, singularity_type_5, time_nplus1, time_n)

  else if ( CPML_region_local == CPML_Z_ONLY ) then

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

    singularity_type_4 = 0
    singularity_type_5 = 0

    call compute_convolution_coef(alpha_x, deltat, coef0_x, coef1_x, coef2_x, 0, time_nplus1, time_n)
    call compute_convolution_coef(alpha_y, deltat, coef0_y, coef1_y, coef2_y, singularity_type_4, time_nplus1, time_n)
    call compute_convolution_coef(alpha_z, deltat, coef0_z, coef1_z, coef2_z, singularity_type_5, time_nplus1, time_n)

  endif

end subroutine l_parameter_computation

!
!=====================================================================
!
Subroutine compute_convolution_coef(bb, deltat, coef0, coef1, coef2, singularity_type, time_nplus1, time_n)

  use constants, only: CUSTOM_REAL

  logical,parameter :: FIRST_ORDER_CONVOLUTION = .false.
  real(kind=CUSTOM_REAL) :: bb, deltat, coef0, coef1, coef2, time_nplus1, time_n
  integer :: singularity_type

  coef0 = exp(-bb * deltat)

  if (singularity_type == 0)then
    if ( abs(bb) >= 1.e-5_CUSTOM_REAL ) then
      if ( FIRST_ORDER_CONVOLUTION ) then
         coef1 = (1._CUSTOM_REAL - exp(-bb * deltat) ) / bb
         coef2 = 0._CUSTOM_REAL
      else
         coef1 = (1._CUSTOM_REAL - exp(-bb * deltat/2._CUSTOM_REAL) ) / bb
         coef2 = (1._CUSTOM_REAL - exp(-bb * deltat/2._CUSTOM_REAL) ) * exp(-bb * deltat/2._CUSTOM_REAL) / bb
      endif
    else
      if ( FIRST_ORDER_CONVOLUTION ) then
        coef1 = deltat
        coef2 = 0._CUSTOM_REAL
      else
        coef1 = deltat/2._CUSTOM_REAL + &
                (- deltat**2*bb/8._CUSTOM_REAL + &
                 (deltat**3*bb**2/48._CUSTOM_REAL - &
                  deltat**4*bb**3/384._CUSTOM_REAL))
        coef2 = deltat/2._CUSTOM_REAL + &
                (- 3._CUSTOM_REAL*deltat**2*bb/8._CUSTOM_REAL + &
                 (7._CUSTOM_REAL*deltat**3*bb**2/48._CUSTOM_REAL - &
                  5._CUSTOM_REAL*deltat**4*bb**3/128._CUSTOM_REAL))
      endif
    endif
  else if (singularity_type == 1)then
    if ( abs(bb) >= 1.e-5_CUSTOM_REAL ) then
      coef1 = (1._CUSTOM_REAL - exp(-bb * deltat/2._CUSTOM_REAL) ) / bb
      coef1 = time_nplus1 * coef1 + (deltat/2._CUSTOM_REAL*exp(-bb * deltat/2._CUSTOM_REAL) - coef1) / bb

      coef2 = (1._CUSTOM_REAL - exp(-bb * deltat/2._CUSTOM_REAL) ) * exp(-bb * deltat/2._CUSTOM_REAL) / bb
      coef2 = time_n * coef2 + ( deltat/2._CUSTOM_REAL*exp(-bb * deltat/2._CUSTOM_REAL) - coef2) / bb
    else
      coef1 = -deltat**2/8._CUSTOM_REAL + (deltat**3*bb/24._CUSTOM_REAL + &
              (- deltat**4*bb**2/128._CUSTOM_REAL + deltat**5*bb**3/960._CUSTOM_REAL)) + &
               time_nplus1* (deltat/2._CUSTOM_REAL + (- deltat**2*bb/8._CUSTOM_REAL + &
                     (deltat**3*bb**2/48._CUSTOM_REAL - deltat**4*bb**3/384._CUSTOM_REAL)))
      coef2 = deltat**2/8._CUSTOM_REAL + (-deltat**3*bb/12._CUSTOM_REAL + &
              (11._CUSTOM_REAL*deltat**4*bb**2/384._CUSTOM_REAL - 13._CUSTOM_REAL*deltat**5*bb**3/1920._CUSTOM_REAL)) + &
               time_n * (deltat/2._CUSTOM_REAL + (- 3._CUSTOM_REAL*deltat**2*bb/8._CUSTOM_REAL + &
                                    (7._CUSTOM_REAL*deltat**3*bb**2/48._CUSTOM_REAL - &
                                     5._CUSTOM_REAL*deltat**4*bb**3/128._CUSTOM_REAL)))
     endif
  else if (singularity_type == 2)then
    if ( abs(bb) >= 1.e-5_CUSTOM_REAL ) then
      coef1 = (1._CUSTOM_REAL - exp(-bb * deltat/2._CUSTOM_REAL) ) / bb
      coef1 = time_nplus1**2 * coef1 + (time_nplus1*(-2._CUSTOM_REAL/bb*coef1 + deltat/bb*exp(-bb * deltat/2._CUSTOM_REAL)) + &
              ((2._CUSTOM_REAL/bb**2* coef1 - exp(-bb * deltat/2._CUSTOM_REAL)*deltat/bb**2) - &
               (deltat/2._CUSTOM_REAL)**2/bb*exp(-bb * deltat/2._CUSTOM_REAL)))

      coef2 = (1._CUSTOM_REAL - exp(-bb * deltat/2._CUSTOM_REAL) ) * exp(-bb * deltat/2._CUSTOM_REAL) / bb
      coef2 = time_n**2 * coef2 + &
              (time_n*(-2._CUSTOM_REAL/bb*coef2 + deltat/bb*exp(-bb * deltat/2._CUSTOM_REAL)) + &
              ((2._CUSTOM_REAL/bb**2* coef2 - exp(-bb * deltat/2._CUSTOM_REAL)*deltat/bb**2) - &
               (deltat/2._CUSTOM_REAL)**2/bb*exp(-bb * deltat/2._CUSTOM_REAL)))

    else
      coef1 = deltat**3/24._CUSTOM_REAL + (- deltat**4*bb*3._CUSTOM_REAL/192._CUSTOM_REAL &
              + (deltat**5*bb**2*3._CUSTOM_REAL/960._CUSTOM_REAL - deltat**6*bb**3*5._CUSTOM_REAL/11520._CUSTOM_REAL)) + &
               time_nplus1*2._CUSTOM_REAL*(-deltat**2/8._CUSTOM_REAL + (deltat**3*bb/24._CUSTOM_REAL + &
                (- deltat**4*bb**2/128._CUSTOM_REAL + deltat**5*bb**3/960._CUSTOM_REAL))) + &
               time_nplus1**2*(deltat/2._CUSTOM_REAL + (- deltat**2*bb/8._CUSTOM_REAL + &
                       (deltat**3*bb**2/48._CUSTOM_REAL - deltat**4*bb**3/384._CUSTOM_REAL)))
      coef2 = deltat**3/24._CUSTOM_REAL + (- deltat**4*bb*5._CUSTOM_REAL/192._CUSTOM_REAL &
              + (deltat**5*bb**2*8._CUSTOM_REAL/960._CUSTOM_REAL - deltat**6*bb**3*7._CUSTOM_REAL/3840._CUSTOM_REAL)) + &
               time_n*2._CUSTOM_REAL*(deltat**2/8._CUSTOM_REAL + (-deltat**3*bb/12._CUSTOM_REAL + &
                (11._CUSTOM_REAL*deltat**4*bb**2/384._CUSTOM_REAL - 13._CUSTOM_REAL*deltat**5*bb**3/1920._CUSTOM_REAL))) + &
                time_n**2 * (deltat/2._CUSTOM_REAL + (- 3._CUSTOM_REAL*deltat**2*bb/8._CUSTOM_REAL + &
                                       (7._CUSTOM_REAL*deltat**3*bb**2/48._CUSTOM_REAL - &
                                        5._CUSTOM_REAL*deltat**4*bb**3/128._CUSTOM_REAL)))
    endif
  else
     stop "error in singularity_type in compute_convolution_coefficient"
  endif

end subroutine compute_convolution_coef
