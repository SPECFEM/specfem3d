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
!> \file save_arrays_solver_adios.F90
!!
!! \author MPBL
!==============================================================================
#include "config.fh"


!==============================================================================
!> Save Moho informtaion using ADIOS
  subroutine crm_save_moho_adios()

  use generate_databases_par, only: myrank,sizeprocs,LOCAL_PATH,NSPEC_AB,NDIM,NGLLSQUARE,IMAIN

  use create_regions_mesh_ext_par

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! local parameters
  !--- Local parameters for ADIOS ---
  character(len=MAX_STRING_LEN) :: output_name
  character(len=*), parameter :: group_name = "SPECFEM3D_MOHO"
  integer(kind=8) :: group_size_inc
  integer(kind=8) :: local_dim

  !--- Variables to allreduce - wmax stands for world_max
  integer :: nspec_wmax, nspec2d_moho_wmax

  !-----------------------------------------------------------------.
  ! Get maximum value for each variable used to define a local_dim. |
  ! ADIOS write equally sized chunks for each processor.            |
  !-----------------------------------------------------------------'
  ! determines maximum values for nspec over all partition slices
  call max_allreduce_singlei(NSPEC_AB,NSPEC_wmax)
  call max_allreduce_singlei(nspec2d_moho,nspec2d_moho_wmax)

  !-----------------------------------.
  ! Setup ADIOS for the current group |
  !-----------------------------------'
  group_size_inc = 0
  output_name = get_adios_filename(trim(LOCAL_PATH) // "/moho")

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  moho file: ',trim(output_name)
#if defined(USE_ADIOS)
    write(IMAIN,*) '  using ADIOS1 file format'
#elif defined(USE_ADIOS2)
    write(IMAIN,*) '  using ADIOS2 file format'
#endif
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! initializes i/o group
  call init_adios_group(myadios_group,group_name)

  !------------------------.
  ! Define ADIOS Variables |
  !------------------------'
  call define_adios_scalar(myadios_group, group_size_inc, '', "nspec", NSPEC_AB)
  call define_adios_scalar(myadios_group, group_size_inc, '', STRINGIFY_VAR(nspec2d_moho))

  local_dim = nspec2d_moho_wmax
  call define_adios_global_array1D(myadios_group, group_size_inc, &
                                   local_dim, '', &
                                   STRINGIFY_VAR(ibelm_moho_top))
  call define_adios_global_array1D(myadios_group, group_size_inc, &
                                   local_dim, '', &
                                   STRINGIFY_VAR(ibelm_moho_bot))

  local_dim = 3 * NGLLSQUARE * nspec2d_moho_wmax
  call define_adios_global_array1D(myadios_group, group_size_inc, &
                                   local_dim, '', &
                                   STRINGIFY_VAR(ijk_moho_top))
  call define_adios_global_array1D(myadios_group, group_size_inc, &
                                   local_dim, '', &
                                   STRINGIFY_VAR(ijk_moho_bot))

  local_dim = NDIM * NGLLSQUARE * nspec2d_moho_wmax
  call define_adios_global_array1D(myadios_group, group_size_inc, &
                                   local_dim, '', &
                                   STRINGIFY_VAR(normal_moho_top))
  call define_adios_global_array1D(myadios_group, group_size_inc, &
                                   local_dim, '', &
                                   STRINGIFY_VAR(normal_moho_bot))

  local_dim = NSPEC_wmax
  call define_adios_global_array1D(myadios_group, group_size_inc, &
                                   local_dim, '', &
                                   STRINGIFY_VAR(is_moho_top))
  call define_adios_global_array1D(myadios_group, group_size_inc, &
                                   local_dim, '', &
                                   STRINGIFY_VAR(is_moho_bot))

  !------------------------------------------------------------.
  ! Open an handler to the ADIOS file and setup the group size |
  !------------------------------------------------------------'
  ! opens file for writing
  call open_file_adios_write(myadios_file,myadios_group,output_name,group_name)

  ! sets group size
  call set_adios_group_size(myadios_file,group_size_inc)

  !------------------------------------------.
  ! Write previously defined ADIOS variables |
  !------------------------------------------'
  call write_adios_scalar(myadios_file, myadios_group, "nspec",NSPEC_AB)
  call write_adios_scalar(myadios_file, myadios_group, STRINGIFY_VAR(nspec2d_moho))

  local_dim = nspec2d_moho_wmax
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(ibelm_moho_top))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(ibelm_moho_bot))

  local_dim = 3 * NGLLSQUARE * nspec2d_moho_wmax
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(ijk_moho_top))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(ijk_moho_bot))

  local_dim = NDIM * NGLLSQUARE * nspec2d_moho_wmax
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(normal_moho_top))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(normal_moho_bot))

  local_dim = NSPEC_wmax
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(is_moho_top))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(is_moho_bot))

  !----------------------------------.
  ! Perform the actual write to disk |
  !----------------------------------'
  call write_adios_perform(myadios_file)

  ! closes file
  call close_file_adios(myadios_file)

  end subroutine crm_save_moho_adios

