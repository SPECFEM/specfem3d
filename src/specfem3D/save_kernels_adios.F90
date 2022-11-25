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


!==============================================================================
!> \file save_kernels_adios.f90
!! \brief Save kernels arrays to file with the help of the ADIOS library.
!! \author MPBL
!==============================================================================
#include "config.fh"


!==============================================================================
!> Define all the kernels that will be written to the ADIOS file.
!! \note Everything is define in this single function, even the group size.

  subroutine define_kernel_adios_variables()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! Variables
  character(len=MAX_STRING_LEN) :: output_name, group_name
  integer(kind=8) :: group_size_inc
  integer(kind=8) :: local_dim

  !--- Variables to allreduce - wmax stands for world_max
  integer :: nspec_wmax

  ! Type inference for define_adios_global_array1D. Avoid additional args. (requires actual size of kernel arrays)
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: dummy_kernel

  ! initializes
  output_name = get_adios_filename(trim(LOCAL_PATH) // "/kernels")

  group_size_inc = 0

  ! initializes i/o group
  group_name = "SPECFEM3D_KERNELS"
  call init_adios_group(myadios_group,group_name)

  ! determines maximum values for nspec over all partition slices
  call max_allreduce_singlei(NSPEC_AB,nspec_wmax)

  ! defines variable entries
  call define_adios_scalar(myadios_group, group_size_inc, '', "nspec", NSPEC_AB)

  local_dim = NGLLX * NGLLY * NGLLZ * nspec_wmax

  if (SAVE_WEIGHTS) then
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "weights_kernel", dummy_kernel)
  endif

  if (ACOUSTIC_SIMULATION) then
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "rho_ac_kl", dummy_kernel)
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "kappa_ac_kl", dummy_kernel)
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "rhop_ac_kl", dummy_kernel)
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "alpha_ac_kl", dummy_kernel)
  endif

  if (ELASTIC_SIMULATION) then
    if (ANISOTROPIC_KL) then
      if (SAVE_TRANSVERSE_KL) then
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "alphav_kl", dummy_kernel)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "alphah_kl", dummy_kernel)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "betav_kl", dummy_kernel)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "betah_kl", dummy_kernel)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "eta_kl", dummy_kernel)
      else
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "rho_kl", dummy_kernel)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "c11_kl", dummy_kernel)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "c12_kl", dummy_kernel)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "c13_kl", dummy_kernel)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "c14_kl", dummy_kernel)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "c15_kl", dummy_kernel)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "c16_kl", dummy_kernel)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "c22_kl", dummy_kernel)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "c23_kl", dummy_kernel)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "c24_kl", dummy_kernel)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "c25_kl", dummy_kernel)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "c26_kl", dummy_kernel)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "c33_kl", dummy_kernel)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "c34_kl", dummy_kernel)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "c35_kl", dummy_kernel)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "c36_kl", dummy_kernel)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "c44_kl", dummy_kernel)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "c45_kl", dummy_kernel)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "c46_kl", dummy_kernel)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "c55_kl", dummy_kernel)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "c56_kl", dummy_kernel)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "c66_kl", dummy_kernel)
      endif
    else
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "rho_kl", dummy_kernel)
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "mu_kl", dummy_kernel)
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "kappa_kl", dummy_kernel)
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "rhop_kl", dummy_kernel)
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "beta_kl", dummy_kernel)
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "alpha_kl", dummy_kernel)
    endif
  endif

  if (POROELASTIC_SIMULATION) then
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "rhot_kl", dummy_kernel)
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "rhof_kl", dummy_kernel)
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "sm_kl", dummy_kernel)
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "eta_kl", dummy_kernel)
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "mufr_kl", dummy_kernel)
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "B_kl", dummy_kernel)
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "C_kl", dummy_kernel)
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "M_kl", dummy_kernel)
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "rhofb_kl", dummy_kernel)
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "phi_kl", dummy_kernel)
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "mufrb_kl", dummy_kernel)
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "Bb_kl", dummy_kernel)
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "Cb_kl", dummy_kernel)
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "Mb_kl", dummy_kernel)
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "rhofbb_kl", dummy_kernel)
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "phib_kl", dummy_kernel)
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "cs_kl", dummy_kernel)
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "cpI_kl", dummy_kernel)
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "cpII_kl", dummy_kernel)
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "ratio_kl", dummy_kernel)
  endif

  if (SAVE_MOHO_MESH) then
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "moho_kl", dummy_kernel)
  endif

  if (APPROXIMATE_HESS_KL) then
    if (ACOUSTIC_SIMULATION) then
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "hess_ac_kl", dummy_kernel)
    endif
    if (ELASTIC_SIMULATION) then
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', "hess_kl", dummy_kernel)
    endif
  endif

  !------------------------------------------------------------.
  ! Open the handle to file containing all the ADIOS variables |
  ! previously defined                                         |
  !------------------------------------------------------------'
  ! opens file for writing
  call open_file_adios_write(myadios_file,myadios_group,output_name,group_name)

  ! sets group size
  call set_adios_group_size(myadios_file,group_size_inc)

  ! writes out nspec
  call write_adios_scalar(myadios_file, myadios_group, "nspec",NSPEC_AB)

  end subroutine define_kernel_adios_variables

!==============================================================================
!> Perform the actual write of all the kernels variables to file.
!!
!! \note Obviously this is a general routine that should be extracted and used
!!       everywhere as the 'adios_handle' argument can be used for any kind of
!!       ADIOS file.
!!       The only reason such a routine is defined is to avoid using
!!       ADIOS modules in non ADIOS file, in case the ADIOS library is not
!!       available on the system.
  subroutine perform_write_adios_kernels()

  use manager_adios

  implicit none

  ! closes file
  call close_file_adios(myadios_file)

  end subroutine perform_write_adios_kernels


!==============================================================================
!> Save weights for volume integration,
!! in order to benchmark the kernels with analytical expressions.
!subroutine save_weights_kernel_adios(weights_kernel)
  !use specfem_par
  !use specfem_par_acoustic
  !use specfem_par_elastic
  !use specfem_par_poroelastic

  !implicit none
  !! local parameters
  !real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: weights_kernel

  !!allocate(weights_kernel(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)

  !return
!end subroutine save_weights_kernel_adios

!==============================================================================
!> Save acoustic related kernels
  subroutine save_kernels_acoustic_adios()

  use specfem_par
  use specfem_par_acoustic

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! Parameters
  ! local parameters
  integer(kind=8) :: local_dim
  !--- Variables to allreduce - wmax stands for world_max
  integer :: nspec_wmax

  ! determines maximum values for nspec over all partition slices
  call max_allreduce_singlei(NSPEC_AB,nspec_wmax)

  local_dim = NGLLX * NGLLY * NGLLZ * nspec_wmax

  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(rho_ac_kl))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(kappa_ac_kl))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(rhop_ac_kl))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(alpha_ac_kl))

  end subroutine save_kernels_acoustic_adios

!==============================================================================
!> Save elastic related kernels
  subroutine save_kernels_elastic_iso_adios(rhop_kl, alpha_kl, beta_kl)

  use specfem_par
  use specfem_par_elastic

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! Parameters
  ! isotropic kernels
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), intent(in) :: rhop_kl, alpha_kl, beta_kl

  ! local parameters
  integer(kind=8) :: local_dim
  !--- Variables to allreduce - wmax stands for world_max
  integer :: nspec_wmax

  ! determines maximum values for nspec over all partition slices
  call max_allreduce_singlei(NSPEC_AB,nspec_wmax)

  local_dim = NGLLX * NGLLY * NGLLZ * nspec_wmax

  ! save kernels to binary files
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(rho_kl))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(mu_kl))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(kappa_kl))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(rhop_kl))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(beta_kl))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(alpha_kl))

  end subroutine save_kernels_elastic_iso_adios

!==============================================================================
!> Save elastic related kernels
  subroutine save_kernels_elastic_aniso_adios(alphav_kl, alphah_kl, betav_kl, betah_kl, eta_kl, &
                                              c11_kl,c12_kl,c13_kl,c14_kl,c15_kl,c16_kl, &
                                              c22_kl,c23_kl,c24_kl,c25_kl,c26_kl, &
                                              c33_kl,c34_kl,c35_kl,c36_kl, &
                                              c44_kl,c45_kl,c46_kl, &
                                              c55_kl,c56_kl, &
                                              c66_kl)

  use specfem_par
  use specfem_par_elastic

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! Parameters
  ! Transverse isotropic kernels
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), intent(in) :: &
    alphav_kl,alphah_kl,betav_kl,betah_kl, eta_kl

  ! full Cijkl kernels
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), intent(in) :: &
    c11_kl,c12_kl,c13_kl,c14_kl,c15_kl,c16_kl, &
    c22_kl,c23_kl,c24_kl,c25_kl,c26_kl, &
    c33_kl,c34_kl,c35_kl,c36_kl, &
    c44_kl,c45_kl,c46_kl, &
    c55_kl,c56_kl, &
    c66_kl

  ! local parameters
  integer(kind=8) :: local_dim
  !--- Variables to allreduce - wmax stands for world_max
  integer :: nspec_wmax

  ! determines maximum values for nspec over all partition slices
  call max_allreduce_singlei(NSPEC_AB,nspec_wmax)

  local_dim = NGLLX * NGLLY * NGLLZ * nspec_wmax

  if (SAVE_TRANSVERSE_KL) then
    ! transverse isotropic kernels
    ! (alpha_v, alpha_h, beta_v, beta_h, eta, rho ) parameterization
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(alphav_kl))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(alphah_kl))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(betav_kl))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(betah_kl))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(eta_kl))
  else
    ! fully anisotropic kernels
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, "rho_kl", -rho_kl)
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(c11_kl))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(c12_kl))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(c13_kl))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(c14_kl))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(c15_kl))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(c16_kl))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(c22_kl))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(c23_kl))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(c24_kl))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(c25_kl))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(c26_kl))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(c33_kl))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(c34_kl))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(c35_kl))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(c36_kl))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(c44_kl))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(c45_kl))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(c46_kl))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(c55_kl))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(c56_kl))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(c66_kl))
  endif

  end subroutine save_kernels_elastic_aniso_adios

!==============================================================================
!> Save poroelastic related kernels
  subroutine save_kernels_poroelastic_adios()

  use specfem_par
  use specfem_par_poroelastic

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! Parameters
  ! local parameters
  integer(kind=8) :: local_dim
  !--- Variables to allreduce - wmax stands for world_max
  integer :: nspec_wmax

  ! determines maximum values for nspec over all partition slices
  call max_allreduce_singlei(NSPEC_AB,nspec_wmax)

  local_dim = NGLLX * NGLLY * NGLLZ * nspec_wmax

  ! primary kernels
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(rhot_kl))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(rhof_kl))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(sm_kl))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(eta_kl))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(mufr_kl))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(B_kl))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(C_kl))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(M_kl))

  ! density kernels
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(rhob_kl))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(rhofb_kl))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(phi_kl))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(mufrb_kl))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(Bb_kl))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(Cb_kl))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(Mb_kl))

  ! wavespeed kernels
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(rhobb_kl))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(rhofbb_kl))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(phib_kl))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(cs_kl))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(cpI_kl))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(cpII_kl))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(ratio_kl))

  end subroutine save_kernels_poroelastic_adios

!==============================================================================
!> Save Moho boundary kernels
  subroutine save_kernels_moho_adios()

  use specfem_par
  use specfem_par_elastic

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! Parameters
  ! local parameters
  integer(kind=8) :: local_dim
  !--- Variables to allreduce - wmax stands for world_max
  integer :: nspec_wmax

  ! safety check
  if (.not. SAVE_MOHO_MESH) return

  ! determines maximum values for nspec over all partition slices
  call max_allreduce_singlei(NSPEC_AB,nspec_wmax)

  ! save moho kernels to binary files
  local_dim = NGLLX * NGLLY * NGLLZ * nspec_wmax
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(moho_kl))

  end subroutine save_kernels_moho_adios


!==============================================================================
!> Save Hessians
  subroutine save_kernels_Hessian_adios()

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! Parameters
  integer(kind=8) :: local_dim
  !--- Variables to allreduce - wmax stands for world_max
  integer :: nspec_wmax

  ! determines maximum values for nspec over all partition slices
  call max_allreduce_singlei(NSPEC_AB,nspec_wmax)

  local_dim = NGLLX * NGLLY * NGLLZ * nspec_wmax

  ! acoustic domains
  if (ACOUSTIC_SIMULATION) then
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(hess_ac_kl))
  endif
  ! elastic domains
  if (ELASTIC_SIMULATION) then
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, STRINGIFY_VAR(hess_kl))
  endif

  end subroutine save_kernels_Hessian_adios

