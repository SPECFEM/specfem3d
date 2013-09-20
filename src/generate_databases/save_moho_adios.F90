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

!==============================================================================
!> \file save_arrays_solver_adios.F90
!!
!! \author MPBL
!==============================================================================

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
!> Save Moho informtaion using ADIOS
subroutine crm_save_moho_adios()

  use mpi
  use adios_helpers_mod
  use generate_databases_par, only : myrank, sizeprocs, LOCAL_PATH, &
                                     NSPEC_AB
  use create_regions_mesh_ext_par

  implicit none

  ! local parameters
  integer :: ier

  !--- Local parameters for ADIOS ---
  character(len=256) :: output_name
  character(len=64), parameter :: group_name  = "SPECFEM3D_MOHO"
  integer(kind=8) :: group, handle
  integer(kind=8) :: groupsize, totalsize
  integer :: local_dim

  !--- Variables to allreduce - wmax stands for world_max
  integer :: nspec_wmax, nspec2d_moho_wmax

  integer, parameter :: num_vars = 2
  integer, dimension(num_vars) :: max_global_values

  !-----------------------------------------------------------------.
  ! Get maximum value for each variable used to define a local_dim. |
  ! ADIOS write equally sized chunks for each processor.            |
  !-----------------------------------------------------------------'
  ! Filling a temporary array to avoid doing allreduces for each var.
  max_global_values(1) = nspec_ab
  max_global_values(2) = nspec2d_moho

  call MPI_Allreduce(MPI_IN_PLACE, max_global_values, num_vars, &
                     MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ier)
  if( ier /= 0 ) call exit_MPI(myrank,'Allreduce to get max values failed.')

  nspec_wmax        = max_global_values(1)
  nspec2d_moho_wmax = max_global_values(2)

  !-----------------------------------.
  ! Setup ADIOS for the current group |
  !-----------------------------------'
  groupsize = 0
  output_name = LOCAL_PATH(1:len_trim(LOCAL_PATH)) // "/moho.bp"
  call adios_declare_group(group, group_name, "", 1, ier)
  call adios_select_method(group, ADIOS_TRANSPORT_METHOD, "", "", ier)

  !------------------------.
  ! Define ADIOS Variables |
  !------------------------'
  call define_adios_scalar(group, groupsize, "", "nspec", NSPEC_AB)
  call define_adios_scalar(group, groupsize, "", STRINGIFY_VAR(nspec2d_moho))

  local_dim = nspec2d_moho_wmax
  call define_adios_global_array1D(group, groupsize, &
                                   local_dim, "",    &
                                   STRINGIFY_VAR(ibelm_moho_top))
  call define_adios_global_array1D(group, groupsize, &
                                   local_dim, "",    &
                                   STRINGIFY_VAR(ibelm_moho_bot))

  local_dim = 3 * NGLLSQUARE * nspec2d_moho_wmax
  call define_adios_global_array1D(group, groupsize, &
                                   local_dim, "",    &
                                   STRINGIFY_VAR(ijk_moho_top))
  call define_adios_global_array1D(group, groupsize, &
                                   local_dim, "",    &
                                   STRINGIFY_VAR(ijk_moho_bot))

  local_dim = NDIM * NGLLSQUARE * nspec2d_moho_wmax
  call define_adios_global_array1D(group, groupsize, &
                                   local_dim, "",    &
                                   STRINGIFY_VAR(normal_moho_top))
  call define_adios_global_array1D(group, groupsize, &
                                   local_dim, "",    &
                                   STRINGIFY_VAR(normal_moho_bot))

  local_dim = nspec_wmax
  call define_adios_global_array1D(group, groupsize, &
                                   local_dim, "",    &
                                   STRINGIFY_VAR(is_moho_top))
  call define_adios_global_array1D(group, groupsize, &
                                   local_dim, "",    &
                                   STRINGIFY_VAR(is_moho_bot))

  !------------------------------------------------------------.
  ! Open an handler to the ADIOS file and setup the group size |
  !------------------------------------------------------------'
  call adios_open(handle, group_name, output_name, "w", &
                  MPI_COMM_WORLD, ier);
  call adios_group_size (handle, groupsize, totalsize, ier)

  !------------------------------------------.
  ! Write previously defined ADIOS variables |
  !------------------------------------------'
  call adios_write(handle, "nspec", NSPEC_AB, ier)
  call adios_write(handle, STRINGIFY_VAR(nspec2d_moho), ier)

  local_dim = nspec2d_moho_wmax
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(ibelm_moho_top))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(ibelm_moho_bot))

  local_dim = 3 * NGLLSQUARE * nspec2d_moho_wmax
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(ijk_moho_top))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(ijk_moho_bot))

  local_dim = NDIM * NGLLSQUARE * nspec2d_moho_wmax
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(normal_moho_top))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(normal_moho_bot))

  local_dim = nspec_wmax
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(is_moho_top))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(is_moho_bot))

  !----------------------------------.
  ! Perform the actual write to disk |
  !----------------------------------'
  call adios_set_path(handle, "", ier)
  call adios_close(handle, ier)
end subroutine crm_save_moho_adios

