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


  subroutine pml_compute_accel_contribution_elastic(ispec,ispec_CPML,displ,veloc, &
                                                    accel_elastic_CPML,rmemory_displ_elastic)

  ! calculates contribution from each C-PML element to update acceleration to the global mesh

  ! second-order accurate convolution term calculation from equation (21) of
  ! Shumin Wang, Robert Lee, and Fernando L. Teixeira,
  ! Anisotropic-medium PML for vector FETD with modified basis functions,
  ! IEEE Transactions on Antennas and Propagation, vol. 54, no. 1, (2006)

  use constants, only: CUSTOM_REAL,NDIM,NGLLX,NGLLY,NGLLZ

  use specfem_par, only: NGLOB_AB,wgll_cube,jacobianstore,ibool,rhostore,irregular_element_number,jacobian_regular

  use pml_par, only: pml_convolution_coef_alpha, &
                     pml_convolution_coef_abar, &
                     NSPEC_CPML, &
                     PML_displ_old,PML_displ_new

  implicit none

  integer, intent(in) :: ispec,ispec_CPML
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB), intent(in) :: displ,veloc
  ! stores C-PML contribution to update acceleration to the global mesh
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ), intent(out) :: accel_elastic_CPML

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ,NSPEC_CPML,3),intent(inout) :: rmemory_displ_elastic

  ! local parameters
  integer :: i,j,k,iglob,ispec_irreg

  real(kind=CUSTOM_REAL) :: wgllcube,rhol,jacobianl
  real(kind=CUSTOM_REAL) :: coef0_x,coef1_x,coef2_x,coef0_y,coef1_y,coef2_y,coef0_z,coef1_z,coef2_z
  real(kind=CUSTOM_REAL) :: A_1,A_2,A_3,A_4,A_5

  ispec_irreg = irregular_element_number(ispec)
  if (ispec_irreg == 0) jacobianl = jacobian_regular

  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        ! coefficients
        ! alpha_x
        coef0_x = pml_convolution_coef_alpha(1,i,j,k,ispec_CPML)
        coef1_x = pml_convolution_coef_alpha(2,i,j,k,ispec_CPML)
        coef2_x = pml_convolution_coef_alpha(3,i,j,k,ispec_CPML)

        ! alpha_y
        coef0_y = pml_convolution_coef_alpha(4,i,j,k,ispec_CPML)
        coef1_y = pml_convolution_coef_alpha(5,i,j,k,ispec_CPML)
        coef2_y = pml_convolution_coef_alpha(6,i,j,k,ispec_CPML)

        ! alpha_z
        coef0_z = pml_convolution_coef_alpha(7,i,j,k,ispec_CPML)
        coef1_z = pml_convolution_coef_alpha(8,i,j,k,ispec_CPML)
        coef2_z = pml_convolution_coef_alpha(9,i,j,k,ispec_CPML)

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

        iglob = ibool(i,j,k,ispec)

        rhol = rhostore(i,j,k,ispec)
        if (ispec_irreg /= 0) jacobianl = jacobianstore(i,j,k,ispec_irreg)
        wgllcube = wgll_cube(i,j,k)

        ! PML coefficient values
        A_1 = pml_convolution_coef_abar(1,i,j,k,ispec_CPML)
        A_2 = pml_convolution_coef_abar(2,i,j,k,ispec_CPML)
        A_3 = pml_convolution_coef_abar(3,i,j,k,ispec_CPML)
        A_4 = pml_convolution_coef_abar(4,i,j,k,ispec_CPML)
        A_5 = pml_convolution_coef_abar(5,i,j,k,ispec_CPML)

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
  subroutine pml_compute_accel_contribution_acoustic(ispec,ispec_CPML, &
                                                     potential_acoustic,potential_dot_acoustic, &
                                                     rmemory_potential_acoustic, &
                                                     potential_dot_dot_acoustic_CPML)

  ! calculates contribution from each C-PML element to update acceleration to the global mesh

  ! second-order accurate convolution term calculation from equation (21) of
  ! Shumin Wang, Robert Lee, and Fernando L. Teixeira,
  ! Anisotropic-medium PML for vector FETD with modified basis functions,
  ! IEEE Transactions on Antennas and Propagation, vol. 54, no. 1, (2006)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

  use specfem_par, only: NGLOB_AB,wgll_cube,jacobianstore,ibool,kappastore,irregular_element_number,jacobian_regular

  use pml_par, only: pml_convolution_coef_alpha, &
                     pml_convolution_coef_abar, &
                     NSPEC_CPML, &
                     PML_potential_acoustic_old,PML_potential_acoustic_new

  implicit none

  integer, intent(in) :: ispec,ispec_CPML
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB), intent(in) :: potential_acoustic,potential_dot_acoustic

  real(kind=CUSTOM_REAL), dimension(3,NGLLX,NGLLY,NGLLZ,NSPEC_CPML),intent(inout) :: rmemory_potential_acoustic
  ! stores C-PML contribution to update the second derivative of the potential to the global mesh
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ), intent(out) :: potential_dot_dot_acoustic_CPML

  ! local parameters
  integer :: i,j,k,iglob,ispec_irreg

  real(kind=CUSTOM_REAL) :: wgllcube,kappal_inv,jacobianl
  real(kind=CUSTOM_REAL) :: coef0_x,coef1_x,coef2_x,coef0_y,coef1_y,coef2_y,coef0_z,coef1_z,coef2_z
  real(kind=CUSTOM_REAL) :: A_1,A_2,A_3,A_4,A_5

  ! irregular element index
  ispec_irreg = irregular_element_number(ispec)
  if (ispec_irreg == 0) jacobianl = jacobian_regular

  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        ! PML coefficient values
        A_1 = pml_convolution_coef_abar(1,i,j,k,ispec_CPML)
        A_2 = pml_convolution_coef_abar(2,i,j,k,ispec_CPML)
        A_3 = pml_convolution_coef_abar(3,i,j,k,ispec_CPML)
        A_4 = pml_convolution_coef_abar(4,i,j,k,ispec_CPML)
        A_5 = pml_convolution_coef_abar(5,i,j,k,ispec_CPML)

        ! coefficients
        ! alpha_x
        coef0_x = pml_convolution_coef_alpha(1,i,j,k,ispec_CPML)
        coef1_x = pml_convolution_coef_alpha(2,i,j,k,ispec_CPML)
        coef2_x = pml_convolution_coef_alpha(3,i,j,k,ispec_CPML)

        ! alpha_y
        coef0_y = pml_convolution_coef_alpha(4,i,j,k,ispec_CPML)
        coef1_y = pml_convolution_coef_alpha(5,i,j,k,ispec_CPML)
        coef2_y = pml_convolution_coef_alpha(6,i,j,k,ispec_CPML)

        ! alpha_z
        coef0_z = pml_convolution_coef_alpha(7,i,j,k,ispec_CPML)
        coef1_z = pml_convolution_coef_alpha(8,i,j,k,ispec_CPML)
        coef2_z = pml_convolution_coef_alpha(9,i,j,k,ispec_CPML)

        ! updates memory variables
        rmemory_potential_acoustic(1,i,j,k,ispec_CPML) = coef0_x * rmemory_potential_acoustic(1,i,j,k,ispec_CPML) &
                + coef1_x * PML_potential_acoustic_new(i,j,k,ispec_CPML) &
                + coef2_x * PML_potential_acoustic_old(i,j,k,ispec_CPML)

        rmemory_potential_acoustic(2,i,j,k,ispec_CPML) = coef0_y * rmemory_potential_acoustic(2,i,j,k,ispec_CPML) &
                + coef1_y * PML_potential_acoustic_new(i,j,k,ispec_CPML) &
                + coef2_y * PML_potential_acoustic_old(i,j,k,ispec_CPML)

        rmemory_potential_acoustic(3,i,j,k,ispec_CPML) = coef0_z * rmemory_potential_acoustic(3,i,j,k,ispec_CPML) &
                + coef1_z * PML_potential_acoustic_new(i,j,k,ispec_CPML) &
                + coef2_z * PML_potential_acoustic_old(i,j,k,ispec_CPML)

        iglob = ibool(i,j,k,ispec)
        wgllcube = wgll_cube(i,j,k)
        if (ispec_irreg /= 0) jacobianl = jacobianstore(i,j,k,ispec_irreg)
        kappal_inv = 1._CUSTOM_REAL / kappastore(i,j,k,ispec)

        ! updates PML potential
        potential_dot_dot_acoustic_CPML(i,j,k) =  wgllcube * kappal_inv * jacobianl * &
                  ( A_1 * potential_dot_acoustic(iglob) + A_2 * potential_acoustic(iglob) &
                  + A_3 * rmemory_potential_acoustic(1,i,j,k,ispec_CPML) &
                  + A_4 * rmemory_potential_acoustic(2,i,j,k,ispec_CPML) &
                  + A_5 * rmemory_potential_acoustic(3,i,j,k,ispec_CPML) )
      enddo
    enddo
  enddo

  end subroutine pml_compute_accel_contribution_acoustic
!
!=====================================================================
!
  subroutine save_field_on_pml_interface(nglob_interface_PML_elastic,b_PML_field,b_reclen_PML_field)

  use constants, only: CUSTOM_REAL,NDIM
  use specfem_par, only: NGLOB_AB,it
  use specfem_par, only: GPU_MODE,Mesh_pointer
  use specfem_par_elastic, only: displ,veloc,accel

  use pml_par, only: points_interface_PML_elastic

  implicit none

  integer, intent(in) :: nglob_interface_PML_elastic,b_reclen_PML_field
  real(kind=CUSTOM_REAL), dimension(9,nglob_interface_PML_elastic) :: b_PML_field

  integer :: iglob_pml,iglob

  if (GPU_MODE) then
    ! needs to transfer forward displ,veloc,accel fields onto CPU for storing
    call transfer_fields_el_from_device(NDIM*NGLOB_AB,displ,veloc,accel,Mesh_pointer)
  endif

  do iglob_pml = 1, nglob_interface_PML_elastic
    iglob = points_interface_PML_elastic(iglob_pml)
    b_PML_field(1,iglob_pml) = displ(1,iglob)
    b_PML_field(2,iglob_pml) = displ(2,iglob)
    b_PML_field(3,iglob_pml) = displ(3,iglob)

    b_PML_field(4,iglob_pml) = veloc(1,iglob)
    b_PML_field(5,iglob_pml) = veloc(2,iglob)
    b_PML_field(6,iglob_pml) = veloc(3,iglob)

    b_PML_field(7,iglob_pml) = accel(1,iglob)
    b_PML_field(8,iglob_pml) = accel(2,iglob)
    b_PML_field(9,iglob_pml) = accel(3,iglob)
  enddo

  call write_abs(0,b_PML_field,b_reclen_PML_field,it)

  end subroutine save_field_on_pml_interface
!
!=====================================================================
!
  subroutine read_field_on_pml_interface(b_accel,b_veloc,b_displ,nglob_interface_PML_elastic, &
                                         b_PML_field,b_reclen_PML_field)

  use constants, only: CUSTOM_REAL,NDIM
  use specfem_par, only: NGLOB_AB,NSTEP,it,UNDO_ATTENUATION_AND_OR_PML ! ibool
  use specfem_par, only: GPU_MODE,Mesh_pointer
  use pml_par, only: points_interface_PML_elastic ! NSPEC_CPML,CPML_to_spec
  implicit none

  integer, intent(in) :: nglob_interface_PML_elastic,b_reclen_PML_field
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: b_displ,b_veloc,b_accel
  real(kind=CUSTOM_REAL), dimension(9,nglob_interface_PML_elastic) :: b_PML_field

  integer :: iglob,iglob_pml,it_temp ! i,j,k,ispec,ispec_pml

  !do ispec_pml = 1, NSPEC_CPML
  !  ispec = CPML_to_spec(ispec_pml)
  !  do i = 1, NGLLX
  !    do j = 1, NGLLY
  !      do k = 1, NGLLZ
  !        iglob = ibool(i,j,k,ispec)
  !        b_displ(:,iglob) = 0._CUSTOM_REAL
  !        b_veloc(:,iglob) = 0._CUSTOM_REAL
  !        b_accel(:,iglob) = 0._CUSTOM_REAL
  !      enddo
  !    enddo
  !  enddo
  !enddo

  ! time step index
  if (UNDO_ATTENUATION_AND_OR_PML) then
    ! forward in time
    it_temp = it
  else
    ! backward in time
    it_temp = NSTEP - it + 1
  endif

  call read_abs(0,b_PML_field,b_reclen_PML_field,it_temp)

  ! GPU routine
  ! not fully implemented yet, will need to transfer wavefields onto CPU and back again
  if (GPU_MODE) then
    ! transfers backward fields onto CPU for updating
    call transfer_b_fields_from_device(NDIM*NGLOB_AB,b_displ,b_veloc,b_accel,Mesh_pointer)
  endif

  ! adds PML contribution to backward/reconstructed fields
  do iglob_pml = 1, nglob_interface_PML_elastic
    iglob = points_interface_PML_elastic(iglob_pml)
    b_displ(1,iglob) = b_PML_field(1,iglob_pml)
    b_displ(2,iglob) = b_PML_field(2,iglob_pml)
    b_displ(3,iglob) = b_PML_field(3,iglob_pml)

    b_veloc(1,iglob) = b_PML_field(4,iglob_pml)
    b_veloc(2,iglob) = b_PML_field(5,iglob_pml)
    b_veloc(3,iglob) = b_PML_field(6,iglob_pml)

    b_accel(1,iglob) = b_PML_field(7,iglob_pml)
    b_accel(2,iglob) = b_PML_field(8,iglob_pml)
    b_accel(3,iglob) = b_PML_field(9,iglob_pml)
  enddo

  if (GPU_MODE) then
    ! transfers backward fields to GPU
    call transfer_b_fields_to_device(NDIM*NGLOB_AB,b_displ,b_veloc,b_accel,Mesh_pointer)
  endif

  end subroutine read_field_on_pml_interface
!
!=====================================================================
!
  subroutine save_potential_on_pml_interface(potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                                             nglob_interface_PML_acoustic,b_PML_potential,b_reclen_PML_potential)

  use specfem_par, only: NGLOB_AB,it
  use constants, only: CUSTOM_REAL
  use pml_par, only: points_interface_PML_acoustic
  implicit none

  integer, intent(in) :: nglob_interface_PML_acoustic,b_reclen_PML_potential
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB), intent(in) :: potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic
  real(kind=CUSTOM_REAL), dimension(3,nglob_interface_PML_acoustic) :: b_PML_potential

  integer :: iglob,iglob_pml

  do iglob_pml = 1, nglob_interface_PML_acoustic
    iglob = points_interface_PML_acoustic(iglob_pml)
    b_PML_potential(1,iglob_pml) = potential_acoustic(iglob)
    b_PML_potential(2,iglob_pml) = potential_dot_acoustic(iglob)
    b_PML_potential(3,iglob_pml) = potential_dot_dot_acoustic(iglob)
  enddo

  call write_abs(1,b_PML_potential,b_reclen_PML_potential,it)

  end subroutine save_potential_on_pml_interface
!
!=====================================================================
!
  subroutine read_potential_on_pml_interface(b_potential_dot_dot_acoustic,b_potential_dot_acoustic,b_potential_acoustic, &
                                             nglob_interface_PML_acoustic,b_PML_potential,b_reclen_PML_potential)

  use specfem_par, only: NGLOB_AB,NSTEP,it,UNDO_ATTENUATION_AND_OR_PML ! ibool
  use pml_par, only: points_interface_PML_acoustic ! NSPEC_CPML,CPML_to_spec
  use constants, only: CUSTOM_REAL

  implicit none

  integer, intent(in) :: nglob_interface_PML_acoustic,b_reclen_PML_potential
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: b_potential_dot_dot_acoustic,b_potential_dot_acoustic,b_potential_acoustic
  real(kind=CUSTOM_REAL), dimension(3,nglob_interface_PML_acoustic) :: b_PML_potential

  ! local parameter
  integer :: iglob,iglob_pml,it_temp ! i,j,k,ispec,ispec_pml

  !do ispec_pml = 1, NSPEC_CPML
  !  ispec = CPML_to_spec(ispec_pml)
  !  do i = 1, NGLLX
  !    do j = 1, NGLLY
  !      do k = 1, NGLLZ
  !        iglob = ibool(i,j,k,ispec)
  !        b_potential_acoustic(iglob) = 0._CUSTOM_REAL
  !        b_potential_dot_acoustic(iglob) = 0._CUSTOM_REAL
  !        b_potential_dot_dot_acoustic(iglob) = 0._CUSTOM_REAL
  !      enddo
  !    enddo
  !  enddo
  !enddo

  ! time step index
  if (UNDO_ATTENUATION_AND_OR_PML) then
    ! forward in time
    it_temp = it
  else
    ! backward in time
    it_temp = NSTEP - it + 1
  endif

  call read_abs(1,b_PML_potential,b_reclen_PML_potential,it_temp)

  do iglob_pml = 1, nglob_interface_PML_acoustic
    iglob = points_interface_PML_acoustic(iglob_pml)
    b_potential_acoustic(iglob) = b_PML_potential(1,iglob_pml)
    b_potential_dot_acoustic(iglob) = b_PML_potential(2,iglob_pml)
    b_potential_dot_dot_acoustic(iglob) = b_PML_potential(3,iglob_pml)
  enddo

  end subroutine read_potential_on_pml_interface

!
!=====================================================================
!

  subroutine l_parameter_computation(kappa_x, d_x, alpha_x, &
                                     kappa_y, d_y, alpha_y, &
                                     kappa_z, d_z, alpha_z, &
                                     CPML_region_local, &
                                     A_0, A_1, A_2, A_3, A_4, A_5)

! computes convolution coefficients
! see: Xie et al. 2014, appendix A1

  use constants, only: CUSTOM_REAL, CPML_XYZ, &
                       CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY, &
                       CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY
  use pml_par, only: min_distance_between_CPML_parameter

  implicit none

  real(kind=CUSTOM_REAL), intent(in) :: kappa_x,d_x,alpha_x, &
                                        kappa_y,d_y,alpha_y, &
                                        kappa_z,d_z,alpha_z
  integer, intent(in) :: CPML_region_local

  real(kind=CUSTOM_REAL), intent(out) :: A_0, A_1, A_2, A_3, A_4, A_5

  !local variable
  real(kind=CUSTOM_REAL) :: bar_A_0, bar_A_1, bar_A_2, bar_A_3, bar_A_4, bar_A_5
  real(kind=CUSTOM_REAL) :: beta_x,beta_y,beta_z

  beta_x = alpha_x + d_x / kappa_x
  beta_y = alpha_y + d_y / kappa_y
  beta_z = alpha_z + d_z / kappa_z

  !unused, now uses pre-computed arrays
  !call compute_convolution_coef(alpha_x, coef0_x, coef1_x, coef2_x)
  !call compute_convolution_coef(alpha_y, coef0_y, coef1_y, coef2_y)
  !call compute_convolution_coef(alpha_z, coef0_z, coef1_z, coef2_z)

  if (CPML_region_local == CPML_XYZ) then

     if ( abs( alpha_x - alpha_y ) >= min_distance_between_CPML_parameter .and. &
          abs( alpha_x - alpha_z ) >= min_distance_between_CPML_parameter .and. &
          abs( alpha_y - alpha_z ) >= min_distance_between_CPML_parameter) then

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
     else
       stop 'Error occured in l_parameter_computation in CPML_XYZ region'
     endif

  else if (CPML_region_local == CPML_XY_ONLY) then

    if (abs( alpha_x - alpha_y ) >= min_distance_between_CPML_parameter) then
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
    else
      stop 'Error occured in l_parameter_computation in CPML_XY_ONLY region'
    endif

  else if (CPML_region_local == CPML_XZ_ONLY) then

    if (abs( alpha_x - alpha_z ) >= min_distance_between_CPML_parameter) then
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
    else
      stop 'Error occured in l_parameter_computation in CPML_XZ_ONLY region'
    endif

  else if (CPML_region_local == CPML_YZ_ONLY) then

    if (abs( alpha_y - alpha_z ) >= min_distance_between_CPML_parameter) then
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
    else
      stop 'Error occured in l_parameter_computation in CPML_YZ_ONLY region'
    endif

  else if (CPML_region_local == CPML_X_ONLY) then

    bar_A_0 = kappa_x
    bar_A_1 = bar_A_0 * (beta_x - alpha_x)
    bar_A_2 = - bar_A_0 * alpha_x * (beta_x - alpha_x)
    bar_A_3 = bar_A_0 * alpha_x**2 * (beta_x - alpha_x)
    bar_A_4 = 0._CUSTOM_REAL
    bar_A_5 = 0._CUSTOM_REAL

  else if (CPML_region_local == CPML_Y_ONLY) then

    bar_A_0 = kappa_y
    bar_A_1 = bar_A_0 * (beta_y - alpha_y)
    bar_A_2 = - bar_A_0 * alpha_y * (beta_y - alpha_y)
    bar_A_3 = 0._CUSTOM_REAL
    bar_A_4 = bar_A_0 * alpha_y**2 * (beta_y - alpha_y)
    bar_A_5 = 0._CUSTOM_REAL

  else if (CPML_region_local == CPML_Z_ONLY) then

    bar_A_0 = kappa_z
    bar_A_1 = bar_A_0 * (beta_z - alpha_z)
    bar_A_2 = - bar_A_0 * alpha_z * (beta_z - alpha_z)
    bar_A_3 = 0._CUSTOM_REAL
    bar_A_4 = 0._CUSTOM_REAL
    bar_A_5 = bar_A_0 * alpha_z**2 * (beta_z - alpha_z)

  endif

  A_0 = bar_A_0
  A_1 = bar_A_1
  A_2 = bar_A_2
  A_3 = bar_A_3
  A_4 = bar_A_4
  A_5 = bar_A_5

  end subroutine l_parameter_computation

!
!=====================================================================
!

  subroutine pml_impose_boundary_condition_acoustic()

! impose Dirichlet conditions for the potential (i.e. Neumann for displacement) on the outer edges of the C-PML layers

  use specfem_par
  use specfem_par_acoustic
  use pml_par

  implicit none

  ! local parameters
  integer :: iface,ispec_CPML,ispec,igll,i,j,k,iglob

  if (.not. GPU_MODE) then
    ! on CPU
    do iface = 1,num_abs_boundary_faces
      ispec = abs_boundary_ispec(iface)
    !!! It is better to move this into do iphase=1,2 loop
        if (ispec_is_acoustic(ispec) .and. is_CPML(ispec)) then
          ! reference GLL points on boundary face
          ispec_CPML = spec_to_CPML(ispec)
          do igll = 1,NGLLSQUARE
            ! gets local indices for GLL point
            i = abs_boundary_ijk(1,igll,iface)
            j = abs_boundary_ijk(2,igll,iface)
            k = abs_boundary_ijk(3,igll,iface)

            iglob=ibool(i,j,k,ispec)

            potential_dot_dot_acoustic(iglob) = 0._CUSTOM_REAL
            potential_dot_acoustic(iglob) = 0._CUSTOM_REAL
            potential_acoustic(iglob) = 0._CUSTOM_REAL
            if (ELASTIC_SIMULATION) then
              PML_potential_acoustic_old(i,j,k,ispec_CPML) = 0._CUSTOM_REAL
              PML_potential_acoustic_new(i,j,k,ispec_CPML) = 0._CUSTOM_REAL
            endif
          enddo
        endif ! ispec_is_acoustic
    !!!   endif
    enddo
  else
    ! on GPU
    stop 'GPU_MODE for routine pml_impose_boundary_condition_acoustic() not implemented yet'
  endif

  end subroutine pml_impose_boundary_condition_acoustic

!
!=====================================================================
!

  subroutine pml_impose_boundary_condition_elastic()

! impose Dirichlet conditions on the outer edges of the C-PML layers

  use specfem_par
  use specfem_par_elastic
  use pml_par

  implicit none

  ! local parameters
  integer :: iface,ispec_CPML,ispec,igll,i,j,k,iglob

  if (.not. GPU_MODE) then
    ! on CPU
    do iface = 1,num_abs_boundary_faces
      ispec = abs_boundary_ispec(iface)
      !!! It is better to move this into do iphase=1,2 loop
      if (ispec_is_elastic(ispec) .and. is_CPML(ispec)) then
        ! reference GLL points on boundary face
        ispec_CPML = spec_to_CPML(ispec)
        do igll = 1,NGLLSQUARE
          ! gets local indices for GLL point
          i = abs_boundary_ijk(1,igll,iface)
          j = abs_boundary_ijk(2,igll,iface)
          k = abs_boundary_ijk(3,igll,iface)

          iglob = ibool(i,j,k,ispec)

          accel(:,iglob) = 0._CUSTOM_REAL
          veloc(:,iglob) = 0._CUSTOM_REAL
          displ(:,iglob) = 0._CUSTOM_REAL

          PML_displ_old(:,i,j,k,ispec_CPML) = 0._CUSTOM_REAL
          PML_displ_new(:,i,j,k,ispec_CPML) = 0._CUSTOM_REAL
        enddo
      endif ! ispec_is_elastic
      !!!        endif
    enddo
  else
    ! on GPU
    call pml_impose_boundary_condition_elastic_cuda(Mesh_pointer)
  endif

  end subroutine pml_impose_boundary_condition_elastic
