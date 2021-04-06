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

!==============================================================================
!> \file save_forward_arrays_adios.F90
!!
!! \author MPBL
!==============================================================================

#include "config.fh"


  subroutine save_forward_arrays_adios()

  use adios_helpers_mod
  use adios_manager_mod, only: comm_adios

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic
  use pml_par

  implicit none

  !--- Local parameters for ADIOS ---
  character(len=MAX_STRING_LEN) :: output_name
  character(len=*), parameter :: group_name = "SPECFEM3D_DATABASES"
  integer(kind=8) :: group, handle
  integer(kind=8) :: groupsize, totalsize
  integer(kind=8) :: local_dim

  !--- Variables to allreduce - wmax stands for world_max
  integer :: nglob_wmax, NSPEC_ATTENUATION_wmax, &
             NSPEC_STRAIN_wmax, N_SLS_wmax
  integer, parameter :: num_vars = 4
  integer, dimension(num_vars) :: max_global_values

  integer :: ier
  integer :: comm

  !-----------------------------------------------------------------.
  ! Get maximum value for each variable used to define a local_dim. |
  ! ADIOS write equally sized chunks for each processor.            |
  !-----------------------------------------------------------------'
  ! Filling a temporary array to avoid doing allreduces for each var.
  max_global_values(1) = NGLOB_AB
  max_global_values(2) =  NSPEC_ATTENUATION_AB
  max_global_values(3) =  NSPEC_STRAIN_ONLY
  max_global_values(4) =  N_SLS

  call max_allreduce_i(max_global_values,num_vars)

  nglob_wmax                   = max_global_values(1)
  NSPEC_ATTENUATION_wmax       = max_global_values(2)
  NSPEC_STRAIN_wmax            = max_global_values(3)
  N_SLS_wmax                   = max_global_values(4)

  !-----------------------------------.
  ! Setup ADIOS for the current group |
  !-----------------------------------'
  groupsize = 0
  output_name = LOCAL_PATH(1:len_trim(LOCAL_PATH))// "/forward_arrays.bp"

  call adios_declare_group(group, group_name, '', 0, ier)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  !call check_adios_err(ier,"Error declare group")

  call adios_select_method(group, ADIOS_TRANSPORT_METHOD, '', '', ier)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  !call check_adios_err(ier,"Error select method")

  !------------------------.
  ! Define ADIOS Variables |
  !------------------------'
  call define_adios_scalar(group, groupsize, '', STRINGIFY_VAR(ngllx))
  call define_adios_scalar(group, groupsize, '', STRINGIFY_VAR(nglly))
  call define_adios_scalar(group, groupsize, '', STRINGIFY_VAR(ngllz))

  call define_adios_scalar(group, groupsize, '', "nglob", NGLOB_AB)
  call define_adios_scalar(group, groupsize, '', STRINGIFY_VAR(NSPEC_ATTENUATION_AB))
  call define_adios_scalar(group, groupsize, '', STRINGIFY_VAR(NSPEC_STRAIN_ONLY))
  call define_adios_scalar(group, groupsize, "", STRINGIFY_VAR(N_SLS))

  if (ACOUSTIC_SIMULATION) then
    local_dim = nglob_wmax
    call define_adios_global_array1D(group, groupsize, local_dim, '', &
                                     STRINGIFY_VAR(potential_acoustic))
    call define_adios_global_array1D(group, groupsize, local_dim, '', &
                                     STRINGIFY_VAR(potential_dot_acoustic))
    call define_adios_global_array1D(group, groupsize, local_dim, '', &
                                     STRINGIFY_VAR(potential_dot_dot_acoustic))
  endif
  if (ELASTIC_SIMULATION) then
    local_dim = NDIM * nglob_wmax
    call define_adios_global_array1D(group, groupsize, local_dim, '', &
                                     STRINGIFY_VAR(displ))
    call define_adios_global_array1D(group, groupsize, local_dim, '', &
                                     STRINGIFY_VAR(veloc))
    call define_adios_global_array1D(group, groupsize, local_dim, '', &
                                     STRINGIFY_VAR(accel))
    if (ATTENUATION) then
      local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_ATTENUATION_wmax * N_SLS_wmax
      call define_adios_global_array1D(group, groupsize, local_dim, '', &
                                       STRINGIFY_VAR(R_xx))
      call define_adios_global_array1D(group, groupsize, local_dim, '', &
                                       STRINGIFY_VAR(R_yy))
      call define_adios_global_array1D(group, groupsize, local_dim, '', &
                                       STRINGIFY_VAR(R_xy))
      call define_adios_global_array1D(group, groupsize, local_dim, '', &
                                       STRINGIFY_VAR(R_xz))
      call define_adios_global_array1D(group, groupsize, local_dim, '', &
                                       STRINGIFY_VAR(R_yz))
      call define_adios_global_array1D(group, groupsize, local_dim, "", &
                                       STRINGIFY_VAR(R_trace))

      local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_STRAIN_wmax
      call define_adios_global_array1D(group, groupsize, local_dim, '', &
                                       STRINGIFY_VAR(epsilondev_xx))
      call define_adios_global_array1D(group, groupsize, local_dim, '', &
                                       STRINGIFY_VAR(epsilondev_yy))
      call define_adios_global_array1D(group, groupsize, local_dim, '', &
                                       STRINGIFY_VAR(epsilondev_xy))
      call define_adios_global_array1D(group, groupsize, local_dim, '', &
                                       STRINGIFY_VAR(epsilondev_xz))
      call define_adios_global_array1D(group, groupsize, local_dim, '', &
                                       STRINGIFY_VAR(epsilondev_yz))
      call define_adios_global_array1D(group, groupsize, local_dim, "", &
                                       STRINGIFY_VAR(epsilondev_trace))
    endif
  endif
  if (POROELASTIC_SIMULATION) then
    local_dim = NDIM * nglob_wmax
    call define_adios_global_array1D(group, groupsize, local_dim, '', &
                                     STRINGIFY_VAR(displs_poroelastic))
    call define_adios_global_array1D(group, groupsize, local_dim, '', &
                                     STRINGIFY_VAR(velocs_poroelastic))
    call define_adios_global_array1D(group, groupsize, local_dim, '', &
                                     STRINGIFY_VAR(accels_poroelastic))
    call define_adios_global_array1D(group, groupsize, local_dim, '', &
                                     STRINGIFY_VAR(displw_poroelastic))
    call define_adios_global_array1D(group, groupsize, local_dim, '', &
                                     STRINGIFY_VAR(velocw_poroelastic))
    call define_adios_global_array1D(group, groupsize, local_dim, '', &
                                     STRINGIFY_VAR(accelw_poroelastic))
  endif

  !------------------------------------------------------------.
  ! Open an handler to the ADIOS file and setup the group size |
  !------------------------------------------------------------'
  ! get MPI communicator
  comm = comm_adios

  call adios_open(handle, group_name, output_name, "w", comm, ier)
  call adios_group_size (handle, groupsize, totalsize, ier)

  !------------------------------------------.
  ! Write previously defined ADIOS variables |
  !------------------------------------------'
  call adios_write(handle, STRINGIFY_VAR(ngllx), ier)
  call adios_write(handle, STRINGIFY_VAR(nglly), ier)
  call adios_write(handle, STRINGIFY_VAR(ngllz), ier)

  call adios_write(handle, "nglob", NGLOB_AB, ier)
  call adios_write(handle, STRINGIFY_VAR(NSPEC_ATTENUATION_AB), ier)
  call adios_write(handle, STRINGIFY_VAR(NSPEC_STRAIN_ONLY), ier)
  call adios_write(handle, STRINGIFY_VAR(N_SLS), ier)

  if (ACOUSTIC_SIMULATION) then
    local_dim = nglob_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(potential_acoustic))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(potential_dot_acoustic))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(potential_dot_dot_acoustic))
  endif
  if (ELASTIC_SIMULATION) then
    local_dim = NDIM * nglob_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(displ))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(veloc))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(accel))
    if (ATTENUATION) then
      local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_ATTENUATION_wmax * N_SLS_wmax
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(R_xx))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(R_yy))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(R_xy))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(R_xz))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(R_yz))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(R_trace))

      local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_STRAIN_wmax
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(epsilondev_xx))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(epsilondev_yy))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(epsilondev_xy))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(epsilondev_xz))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(epsilondev_yz))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(epsilondev_trace))
    endif
  endif
  if (POROELASTIC_SIMULATION) then
    local_dim = NDIM * nglob_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(displs_poroelastic))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(velocs_poroelastic))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(accels_poroelastic))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(displw_poroelastic))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(velocw_poroelastic))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(accelw_poroelastic))
  endif

  !----------------------------------.
  ! Perform the actual write to disk |
  !----------------------------------'
  call adios_set_path(handle, '', ier)
  call adios_close(handle, ier)

  end subroutine save_forward_arrays_adios
