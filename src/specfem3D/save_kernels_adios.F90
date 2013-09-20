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


!==============================================================================
!> \file save_kernels_adios.f90
!! \brief Save kernels arrays to file with the help of the ADIOS library.
!! \author MPBL
!==============================================================================

!> \def STRINGIFY_VAR(a)
!! Macro taking a variable and returning the stringified variable and
!! the variable itself.
!! STRINGIFY_VAR(x) expand as:
!!   "x", x
!! x being the variable name inside the code.
#ifdef __INTEL_COMPILER
#define STRINGIFY_VAR(a) #a, a
#else
#define STRINGIFY_VAR(a) "a", a
#endif

!==============================================================================
!> Define all the kernels that will be written to the ADIOS file.
!! \param[INOUT] adios_handle The handle pointing on the open ADIOS file
!!                            intended to store kernels data.
!! \note Everything is define in this single function, even the group size.
!!       It is the reason why this function require only an handle on an ADIOS
!!       file as an argument.
subroutine define_kernel_adios_variables(handle, SAVE_WEIGHTS)

  use mpi
  use adios_write_mod

  use adios_helpers_mod

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  implicit none

  ! Parameters
  integer(kind=8), intent(INOUT) :: handle
  logical, intent(IN) :: SAVE_WEIGHTS
  ! Variables
  character(len=256) :: output_name, group_name
  integer(kind=8) :: group, groupsize, adios_totalsize
  integer :: local_dim, comm, adios_err, ierr
  !--- Variables to allreduce - wmax stands for world_max
  integer :: nspec_wmax, ier
  integer, parameter :: num_vars = 1
  integer, dimension(num_vars) :: max_global_values
  ! Type inference for define_adios_global_array1D. Avoid additional args.
  real(kind=CUSTOM_REAL), dimension(1,1,1,1) :: dummy_kernel

  output_name = LOCAL_PATH(1:len_trim(LOCAL_PATH))// "/kernels.bp"
  group_name = "SPECFEM3D_KERNELS"
  call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)

  groupsize = 0
  call adios_declare_group(group, group_name, "", 0, adios_err)
  call adios_select_method(group, "MPI", "", "", adios_err)

  max_global_values(1) = NSPEC_AB

  call MPI_Allreduce(MPI_IN_PLACE, max_global_values, num_vars, &
                     MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ier)
  if( ier /= 0 ) call exit_MPI(myrank,'Allreduce to get max values failed.')

  nspec_wmax = max_global_values(1)

  call define_adios_scalar(group, groupsize, "", "nspec", NSPEC_AB)

  local_dim = NGLLX * NGLLY * NGLLZ * nspec_wmax

  if( SAVE_WEIGHTS ) then
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", "weights_kernel", dummy_kernel)
  endif

  if( ACOUSTIC_SIMULATION ) then
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", "rho_ac_kl", dummy_kernel)
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", "kappa_ac_kl", dummy_kernel)
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", "rhop_ac_kl", dummy_kernel)
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", "alpha_ac_kl", dummy_kernel)
  endif

  if( ELASTIC_SIMULATION ) then
    if (ANISOTROPIC_KL) then
      if (SAVE_TRANSVERSE_KL) then
        call define_adios_global_array1D(group, groupsize, local_dim, &
                                         "", "alphav_kl", dummy_kernel)
        call define_adios_global_array1D(group, groupsize, local_dim, &
                                         "", "alphah_kl", dummy_kernel)
        call define_adios_global_array1D(group, groupsize, local_dim, &
                                         "", "betav_kl", dummy_kernel)
        call define_adios_global_array1D(group, groupsize, local_dim, &
                                         "", "betah_kl", dummy_kernel)
        call define_adios_global_array1D(group, groupsize, local_dim, &
                                         "", "eta_kl", dummy_kernel)
        call define_adios_global_array1D(group, groupsize, local_dim, &
                                         "", "alpha_kl", dummy_kernel)
        call define_adios_global_array1D(group, groupsize, local_dim, &
                                         "", "beta_kl", dummy_kernel)
      else
        call define_adios_global_array1D(group, groupsize, local_dim, &
                                         "", "rho_kl", dummy_kernel)
        call define_adios_global_array1D(group, groupsize, local_dim, &
                                         "", "cijkl_kl", dummy_kernel)
      endif
    else
      call define_adios_global_array1D(group, groupsize, local_dim, &
                                       "", "rho_kl", dummy_kernel)
      call define_adios_global_array1D(group, groupsize, local_dim, &
                                       "", "mu_kl", dummy_kernel)
      call define_adios_global_array1D(group, groupsize, local_dim, &
                                       "", "kappa_kl", dummy_kernel)
      call define_adios_global_array1D(group, groupsize, local_dim, &
                                       "", "rhop_kl", dummy_kernel)
      call define_adios_global_array1D(group, groupsize, local_dim, &
                                       "", "beta_kl", dummy_kernel)
      call define_adios_global_array1D(group, groupsize, local_dim, &
                                       "", "alpha_kl", dummy_kernel)
    endif
    if (SAVE_MOHO_MESH) then
      call define_adios_global_array1D(group, groupsize, local_dim, &
                                       "", "moho_kl", dummy_kernel)
    endif
  endif

  if( POROELASTIC_SIMULATION ) then
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", "rhot_kl", dummy_kernel)
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", "rhof_kl", dummy_kernel)
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", "sm_kl", dummy_kernel)
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", "eta_kl", dummy_kernel)
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", "mufr_kl", dummy_kernel)
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", "B_kl", dummy_kernel)
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", "C_kl", dummy_kernel)
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", "M_kl", dummy_kernel)
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", "rhofb_kl", dummy_kernel)
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", "phi_kl", dummy_kernel)
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", "mufrb_kl", dummy_kernel)
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", "Bb_kl", dummy_kernel)
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", "Cb_kl", dummy_kernel)
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", "Mb_kl", dummy_kernel)
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", "rhofbb_kl", dummy_kernel)
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", "phib_kl", dummy_kernel)
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", "cs_kl", dummy_kernel)
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", "cpI_kl", dummy_kernel)
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", "cpII_kl", dummy_kernel)
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", "ratio_kl", dummy_kernel)
  endif

  if ( APPROXIMATE_HESS_KL ) then
    if( ACOUSTIC_SIMULATION ) then
      call define_adios_global_array1D(group, groupsize, local_dim, &
                                       "", "hess_ac_kl", dummy_kernel)
    endif
    if( ELASTIC_SIMULATION ) then
      call define_adios_global_array1D(group, groupsize, local_dim, &
                                       "", "hess_kl", dummy_kernel)
    endif
  endif

  !------------------------------------------------------------.
  ! Open the handle to file containing all the ADIOS variables |
  ! previously defined                                         |
  !------------------------------------------------------------'
  call adios_open (handle, group_name, output_name, "w", comm, adios_err)
  call adios_group_size (handle, groupsize, adios_totalsize, adios_err)

  call adios_write(handle, "nspec", NSPEC_AB, ier)
end subroutine define_kernel_adios_variables

!==============================================================================
!> Perform the actual write of all the kernels variables to file.
!! \param[IN] adios_handle The handle pointing on the open ADIOS file intended
!!                         to store kernels.
!!
!! \note Obviously this is a general routine that should be extracted and used
!!       everywhere as the 'adios_handle' argument can be used for any kind of
!!       ADIOS file.
!!       The only reason such a routine is defined is to avoid using
!!       ADIOS modules in non ADIOS file, in case the ADIOS library is not
!!       available on the system.
subroutine perform_write_adios_kernels(handle)

  use adios_write_mod

  implicit none

  ! Parameters
  integer(kind=8), intent(in) :: handle
  ! Variables
  integer :: adios_err

  call adios_close(handle, adios_err)
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
subroutine save_kernels_acoustic_adios(handle)

  use specfem_par
  use specfem_par_acoustic
  use adios_helpers_mod

  implicit none

  ! Parameters
  integer(kind=8), intent(in) :: handle
  ! local parameters
  integer:: local_dim

  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_AB

  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(rho_ac_kl))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(kappa_ac_kl))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(rhop_ac_kl))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(alpha_ac_kl))
end subroutine save_kernels_acoustic_adios

!==============================================================================
!> Save elastic related kernels
subroutine save_kernels_elastic_adios(handle, alphav_kl, alphah_kl, &
                                      betav_kl, betah_kl, eta_kl,   &
                                      rhop_kl, alpha_kl, beta_kl)

  use specfem_par
  use specfem_par_elastic
  use adios_helpers_mod

  implicit none

  ! Parameters
  integer(kind=8), intent(in) :: handle
  ! local parameters
  integer:: local_dim

  ! Transverse isotropic paramters
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    alphav_kl,alphah_kl,betav_kl,betah_kl, &
    eta_kl, rhop_kl, alpha_kl, beta_kl

  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_AB

  if (ANISOTROPIC_KL) then
    ! outputs transverse isotropic kernels only
    if (SAVE_TRANSVERSE_KL) then
      ! transverse isotropic kernels
      ! (alpha_v, alpha_h, beta_v, beta_h, eta, rho ) parameterization
      call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                       local_dim, STRINGIFY_VAR(alphav_kl))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                       local_dim, STRINGIFY_VAR(alphah_kl))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                       local_dim, STRINGIFY_VAR(betav_kl))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                       local_dim, STRINGIFY_VAR(betah_kl))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                       local_dim, STRINGIFY_VAR(eta_kl))

      ! transverse isotropic test kernels
      call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                       local_dim, STRINGIFY_VAR(alpha_kl))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                       local_dim, STRINGIFY_VAR(beta_kl))
    else
      ! fully anisotropic kernels
      call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                       local_dim, "rho_kl", -rho_kl)
      call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                       local_dim, "cijkl_kl", -cijkl_kl)
    endif
  else
    ! save kernels to binary files
    call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                     local_dim, STRINGIFY_VAR(rho_kl))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                     local_dim, STRINGIFY_VAR(mu_kl))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                     local_dim, STRINGIFY_VAR(kappa_kl))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                     local_dim, STRINGIFY_VAR(rhop_kl))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                     local_dim, STRINGIFY_VAR(beta_kl))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                     local_dim, STRINGIFY_VAR(alpha_kl))
  endif

  if (SAVE_MOHO_MESH) then
    call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                     local_dim, STRINGIFY_VAR(moho_kl))
  endif

end subroutine save_kernels_elastic_adios

!==============================================================================
!> Save poroelastic related kernels
subroutine save_kernels_poroelastic_adios(handle)

  use specfem_par
  use specfem_par_poroelastic
  use adios_helpers_mod

  implicit none

  ! Parameters
  integer(kind=8), intent(in) :: handle
  ! local parameters
  integer :: local_dim

  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_AB

  ! primary kernels
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(rhot_kl))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(rhof_kl))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(sm_kl))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(eta_kl))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(mufr_kl))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(B_kl))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(C_kl))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(M_kl))

  ! density kernels
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(rhob_kl))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(rhofb_kl))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(phi_kl))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(mufrb_kl))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(Bb_kl))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(Cb_kl))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(Mb_kl))

  ! wavespeed kernels
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(rhobb_kl))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(rhofbb_kl))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(phib_kl))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(cs_kl))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(cpI_kl))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(cpII_kl))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(ratio_kl))

end subroutine save_kernels_poroelastic_adios

!==============================================================================
!> Save hessians
subroutine save_kernels_hessian_adios(handle)

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use adios_helpers_mod

  implicit none

  ! Parameters
  integer(kind=8), intent(in) :: handle
  integer :: local_dim

  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_AB

  ! acoustic domains
  if( ACOUSTIC_SIMULATION ) then
    call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                     local_dim, STRINGIFY_VAR(hess_ac_kl))
  endif
  ! elastic domains
  if( ELASTIC_SIMULATION ) then
    call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                     local_dim, STRINGIFY_VAR(hess_kl))
  endif

end subroutine save_kernels_hessian_adios

