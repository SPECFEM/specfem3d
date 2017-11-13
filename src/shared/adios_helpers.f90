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


!===============================================================================
!> Helpers to set up adios features.
!! * Error checking
!! * Scalar definition
!! * Global arrays definition
!!
!! \author MPBL
!-------------------------------------------------------------------------------
module adios_helpers_mod

  use adios_helpers_definitions_mod
  use adios_helpers_writers_mod

  implicit none

  private

  ! from adios_helpers_definitions_mod
  public :: define_adios_scalar
  public :: define_adios_global_real_1d_array
  public :: define_adios_global_double_1d_array
  public :: define_adios_global_integer_1d_array
  public :: define_adios_global_long_1d_array
  public :: define_adios_global_logical_1d_array
  public :: define_adios_global_string_1d_array
  public :: define_adios_local_string_1d_array
  public :: define_adios_global_array1D

  ! from adios_helpers_writers_mod
  public :: write_adios_global_real_1d_array
  public :: write_adios_global_double_1d_array
  public :: write_adios_global_integer_1d_array
  public :: write_adios_global_long_1d_array
  public :: write_adios_global_logical_1d_array
  public :: write_adios_global_string_1d_array
  public :: write_adios_global_1d_array

  public :: write_adios_global_1d_array_offset
  public :: write_adios_global_real_1d_array_offset
  public :: write_adios_global_integer_1d_array_offset

  ! new from adios_helpers_writers_mod
  public :: check_adios_err

end module adios_helpers_mod
