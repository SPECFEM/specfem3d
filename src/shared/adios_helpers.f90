!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            April 2011
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

  ! from this module. No 'imports'
  public :: check_adios_err

  ! from adios_helpers_definitions_mod
  public :: define_adios_scalar
  public :: define_adios_global_real_1d_array
  public :: define_adios_global_double_1d_array
  public :: define_adios_global_integer_1d_array
  public :: define_adios_global_long_1d_array
  public :: define_adios_global_logical_1d_array
  public :: define_adios_global_array1D

  ! from adios_helpers_writers_mod
  public :: write_adios_global_real_1d_array
  public :: write_adios_global_double_1d_array
  public :: write_adios_global_integer_1d_array
  public :: write_adios_global_long_1d_array
  public :: write_adios_global_logical_1d_array
  public :: write_adios_global_1d_array
contains

!===============================================================================
!> Get the ADIOS error message from an adios error number if there is an error.
!! \param adios_err The error code considered.
subroutine check_adios_err(myrank, adios_err)
  use adios_read_mod
  implicit none
  integer, intent(in) :: myrank, adios_err
  character(len=1024) :: msg

  if (adios_err /= 0) then
    call adios_errmsg(msg)
    print *, "process: ", myrank, ", error: ", msg
    stop
  endif
end subroutine check_adios_err

end module adios_helpers_mod
