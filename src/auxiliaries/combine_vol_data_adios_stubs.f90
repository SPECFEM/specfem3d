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

module combine_vol_data_adios_mod

  use adios_manager_mod

  implicit none

contains

!=============================================================================

  subroutine print_usage_adios()

  call no_adios_err()

  end subroutine print_usage_adios

!=============================================================================

  subroutine read_args_adios(arg, MAX_NUM_NODES, node_list, num_node, &
                             var_name, value_file_name, mesh_file_name, &
                             outdir, ires)
  ! Arguments
  character(len=*) :: arg(:)
  integer :: MAX_NUM_NODES
  integer :: node_list(:)
  integer :: num_node, ires
  character(len=*) :: var_name, value_file_name, mesh_file_name, outdir

  integer(kind=4) :: unused_i4

  unused_i4 = len_trim(arg(1))
  unused_i4 = MAX_NUM_NODES
  unused_i4 = node_list(1)
  unused_i4 = num_node
  unused_i4 = ires
  unused_i4 = len_trim(var_name)
  unused_i4 = len_trim(value_file_name)
  unused_i4 = len_trim(mesh_file_name)
  unused_i4 = len_trim(outdir)

  call no_adios_err()

  end subroutine read_args_adios

!=============================================================================

  subroutine init_adios(value_file_name, mesh_file_name, value_handle, mesh_handle)

  ! Parameters
  character(len=*) :: value_file_name, mesh_file_name
  integer(kind=8) :: value_handle, mesh_handle

  integer(kind=4) :: unused_i4
  integer(kind=8) :: unused_i8

  unused_i4 = len_trim(value_file_name)
  unused_i4 = len_trim(mesh_file_name)
  unused_i8 = value_handle
  unused_i8 = mesh_handle

  call no_adios_err()

  end subroutine init_adios

!=============================================================================

  subroutine clean_adios(value_handle, mesh_handle)

  ! Parameters
  integer(kind=8) :: value_handle, mesh_handle

  integer(kind=8) :: unused_i8

  unused_i8 = value_handle
  unused_i8 = mesh_handle

  call no_adios_err()

  end subroutine clean_adios

!=============================================================================

  subroutine read_scalars_adios_mesh(mesh_handle, iproc, NGLOB_AB, NSPEC_AB, ibool_offset, x_global_offset)

  ! Parameters
  integer(kind=8) :: mesh_handle
  integer :: iproc
  integer :: NGLOB_AB, NSPEC_AB
  integer(kind=8) :: ibool_offset, x_global_offset

  integer(kind=4) :: unused_i4
  integer(kind=8) :: unused_i8

  unused_i8 = mesh_handle
  unused_i4 = iproc
  unused_i4 = NGLOB_AB
  unused_i4 = NSPEC_AB
  unused_i8 = ibool_offset
  unused_i8 = x_global_offset

  call no_adios_err()

  end subroutine read_scalars_adios_mesh

!=============================================================================

  subroutine read_ibool_adios_mesh(mesh_handle, ibool_offset, NGLLX, NGLLY, NGLLZ, NSPEC_AB, ibool)

  ! Parameters
  integer(kind=8) :: mesh_handle,ibool_offset
  integer :: NGLLX, NGLLY, NGLLZ, NSPEC_AB
  integer, dimension(:,:,:,:) :: ibool

  integer(kind=4) :: unused_i4
  integer(kind=8) :: unused_i8

  unused_i8 = mesh_handle
  unused_i8 = ibool_offset
  unused_i4 = NGLLX
  unused_i4 = NGLLY
  unused_i4 = NGLLZ
  unused_i4 = NSPEC_AB
  unused_i4 = ibool(1,1,1,1)

  call no_adios_err()

  end subroutine read_ibool_adios_mesh

!=============================================================================

  subroutine read_coordinates_adios_mesh(mesh_handle, x_global_offset, NGLOB_AB, xstore, ystore, zstore)

  use constants

  ! Parameters
  integer(kind=8) :: mesh_handle,x_global_offset
  integer :: NGLOB_AB
  real(kind=CUSTOM_REAL),dimension(:) :: xstore, ystore, zstore

  integer(kind=4) :: unused_i4
  integer(kind=8) :: unused_i8
  real(kind=CUSTOM_REAL) :: unused_cr

  unused_i8 = mesh_handle
  unused_i8 = x_global_offset
  unused_i4 = NGLOB_AB
  unused_cr = xstore(1)
  unused_cr = ystore(1)
  unused_cr = zstore(1)

  call no_adios_err()

  end subroutine read_coordinates_adios_mesh

!=============================================================================

  subroutine read_double_values_adios(value_handle, var_name, ibool_offset, NSPEC_AB, dat)

  ! Parameters
  integer(kind=8) :: value_handle,ibool_offset
  character(len=*) :: var_name
  integer :: NSPEC_AB
  double precision, dimension(:,:,:,:) :: dat

  integer(kind=4) :: unused_i4
  integer(kind=8) :: unused_i8
  double precision :: unused_dp

  unused_i8 = value_handle
  unused_i4 = len_trim(var_name)
  unused_i8 = ibool_offset
  unused_i4 = NSPEC_AB
  unused_dp = dat(1,1,1,1)

  call no_adios_err()

  end subroutine read_double_values_adios

!=============================================================================

  subroutine read_float_values_adios(value_handle, var_name, ibool_offset, NSPEC_AB, dat)

  ! Parameters
  integer(kind=8) :: value_handle,ibool_offset
  character(len=*) :: var_name
  integer :: NSPEC_AB
  real, dimension(:,:,:,:) :: dat

  integer(kind=4) :: unused_i4
  integer(kind=8) :: unused_i8
  double precision :: unused_dp

  unused_i8 = value_handle
  unused_i4 = len_trim(var_name)
  unused_i8 = ibool_offset
  unused_i4 = NSPEC_AB
  unused_dp = dat(1,1,1,1)

  call no_adios_err()

  end subroutine read_float_values_adios

end module combine_vol_data_adios_mod
