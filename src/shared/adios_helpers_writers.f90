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
!! * Scalar definition
!! * Global arrays definition
!!
!! \note We do not define function to write scalars variables into adios
!!       since it is already a single function call.
!!
!! \author MPBL
!-------------------------------------------------------------------------------
module adios_helpers_writers_mod

  implicit none

  private

  public :: write_adios_global_real_1d_array
  public :: write_adios_global_double_1d_array
  public :: write_adios_global_integer_1d_array
  public :: write_adios_global_long_1d_array
  public :: write_adios_global_logical_1d_array
  public :: write_adios_global_string_1d_array
  public :: write_adios_global_1d_array

  interface write_adios_global_real_1d_array
    module procedure write_adios_global_1d_real_1d
    module procedure write_adios_global_1d_real_2d
    module procedure write_adios_global_1d_real_3d
    module procedure write_adios_global_1d_real_4d
    module procedure write_adios_global_1d_real_5d
  end interface write_adios_global_real_1d_array

  interface write_adios_global_double_1d_array
    module procedure write_adios_global_1d_double_1d
    module procedure write_adios_global_1d_double_2d
    module procedure write_adios_global_1d_double_3d
    module procedure write_adios_global_1d_double_4d
    module procedure write_adios_global_1d_double_5d
  end interface write_adios_global_double_1d_array

  interface write_adios_global_integer_1d_array
    module procedure write_adios_global_1d_integer_1d
    module procedure write_adios_global_1d_integer_2d
    module procedure write_adios_global_1d_integer_3d
    module procedure write_adios_global_1d_integer_4d
    module procedure write_adios_global_1d_integer_5d
  end interface write_adios_global_integer_1d_array

  interface write_adios_global_long_1d_array
    module procedure write_adios_global_1d_long_1d
    module procedure write_adios_global_1d_long_2d
    module procedure write_adios_global_1d_long_3d
    module procedure write_adios_global_1d_long_4d
    module procedure write_adios_global_1d_long_5d
  end interface write_adios_global_long_1d_array

  interface write_adios_global_logical_1d_array
    module procedure write_adios_global_1d_logical_1d
    module procedure write_adios_global_1d_logical_2d
    module procedure write_adios_global_1d_logical_3d
    module procedure write_adios_global_1d_logical_4d
    module procedure write_adios_global_1d_logical_5d
  end interface write_adios_global_logical_1d_array

  interface write_adios_global_string_1d_array
    module procedure write_adios_global_1d_string_1d
  end interface write_adios_global_string_1d_array

  interface write_adios_global_1d_array
    module procedure write_adios_global_1d_integer_1d
    module procedure write_adios_global_1d_integer_2d
    module procedure write_adios_global_1d_integer_3d
    module procedure write_adios_global_1d_integer_4d
    module procedure write_adios_global_1d_integer_5d

    module procedure write_adios_global_1d_long_1d
    module procedure write_adios_global_1d_long_2d
    module procedure write_adios_global_1d_long_3d
    module procedure write_adios_global_1d_long_4d
    module procedure write_adios_global_1d_long_5d

    module procedure write_adios_global_1d_logical_1d
    module procedure write_adios_global_1d_logical_2d
    module procedure write_adios_global_1d_logical_3d
    module procedure write_adios_global_1d_logical_4d
    module procedure write_adios_global_1d_logical_5d

    module procedure write_adios_global_1d_real_1d
    module procedure write_adios_global_1d_real_2d
    module procedure write_adios_global_1d_real_3d
    module procedure write_adios_global_1d_real_4d
    module procedure write_adios_global_1d_real_5d

    module procedure write_adios_global_1d_double_1d
    module procedure write_adios_global_1d_double_2d
    module procedure write_adios_global_1d_double_3d
    module procedure write_adios_global_1d_double_4d
    module procedure write_adios_global_1d_double_5d
  end interface write_adios_global_1d_array


  !------------------------------------
  !
  ! with additional offset info
  !
  !------------------------------------
  public :: write_adios_global_real_1d_array_offset
  public :: write_adios_global_integer_1d_array_offset
  public :: write_adios_global_1d_array_offset

  interface write_adios_global_real_1d_array_offset
    module procedure write_adios_global_1d_real_1d_offset
  end interface write_adios_global_real_1d_array_offset

  interface write_adios_global_integer_1d_array_offset
    module procedure write_adios_global_1d_integer_1d_offset
  end interface write_adios_global_integer_1d_array_offset

  interface write_adios_global_1d_array_offset
    module procedure write_adios_global_1d_integer_1d_offset
    module procedure write_adios_global_1d_long_1d_offset
    module procedure write_adios_global_1d_logical_1d_offset
    module procedure write_adios_global_1d_real_1d_offset
    module procedure write_adios_global_1d_double_1d_offset
  end interface write_adios_global_1d_array_offset

  ! error checking
  public :: check_adios_err

contains


!===============================================================================

  subroutine write_1D_global_array_adios_dims(adios_handle, myrank, local_dim, sizeprocs, path)

  implicit none

  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: sizeprocs, myrank
  integer(kind=8), intent(in) :: local_dim
  character(len=*), intent(in) :: path

  integer :: adios_err
  integer(kind=8) :: offset
  integer(kind=8) :: global_dim ! should be long for large cases.

  ! checks if file handle valid
  if (adios_handle == 0) stop 'Invalid ADIOS file handle in write_1D_global_array_adios_dims()'

  ! global dimension
  global_dim = local_dim * int(sizeprocs,kind=8)

  ! process offset
  ! note: assumes that myrank starts is within [0,sizeprocs-1]
  offset = local_dim * int(myrank,kind=8)

  call adios_write(adios_handle, trim(path) // "/local_dim", local_dim, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_write(adios_handle, trim(path) // "/global_dim", global_dim, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_write(adios_handle, trim(path) // "/offset", offset, adios_err)
  call check_adios_err(myrank,adios_err)

  end subroutine write_1D_global_array_adios_dims


!===============================================================================
!> Schedule an ADIOS single precision global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
  subroutine write_adios_global_1d_real_1d(adios_handle, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  character(len=*) :: array_name
  real, dimension(:), intent(in) :: array
  ! Variables
  integer :: adios_err

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in write_adios_global_1d'

  ! adds local_dim/global_dim/offset infos
  call write_1D_global_array_adios_dims(adios_handle, myrank, local_dim, sizeprocs, array_name)

  call adios_write(adios_handle, trim(array_name) // "/array", array, adios_err)
  call check_adios_err(myrank,adios_err)

  end subroutine write_adios_global_1d_real_1d


!===============================================================================
!> Schedule an ADIOS single precision global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
  subroutine write_adios_global_1d_real_2d(adios_handle, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  character(len=*) :: array_name
  real, dimension(:,:), intent(in) :: array
  ! Variables
  integer :: adios_err

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in write_adios_global_1d'

  ! adds local_dim/global_dim/offset infos
  call write_1D_global_array_adios_dims(adios_handle, myrank, local_dim, sizeprocs, array_name)

  call adios_write(adios_handle, trim(array_name) // "/array", array, adios_err)
  call check_adios_err(myrank,adios_err)

  end subroutine write_adios_global_1d_real_2d


!===============================================================================
!> Schedule an ADIOS single precision global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
  subroutine write_adios_global_1d_real_3d(adios_handle, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  character(len=*) :: array_name
  real, dimension(:,:,:), intent(in) :: array
  ! Variables
  integer :: adios_err
  !character(len=1024) :: msg

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in write_adios_global_1d'

  ! adds local_dim/global_dim/offset infos
  call write_1D_global_array_adios_dims(adios_handle, myrank, local_dim, sizeprocs, array_name)

  call adios_write(adios_handle, trim(array_name) // "/array", array, adios_err)
  call check_adios_err(myrank,adios_err)

  end subroutine write_adios_global_1d_real_3d


!===============================================================================
!> Schedule an ADIOS single precision global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
  subroutine write_adios_global_1d_real_4d(adios_handle, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  character(len=*) :: array_name
  real, dimension(:,:,:,:), intent(in) :: array
  ! Variables
  integer :: adios_err

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in write_adios_global_1d'

  ! adds local_dim/global_dim/offset infos
  call write_1D_global_array_adios_dims(adios_handle, myrank, local_dim, sizeprocs, array_name)

  call adios_write(adios_handle, trim(array_name) // "/array", array, adios_err)
  call check_adios_err(myrank,adios_err)

  end subroutine write_adios_global_1d_real_4d


!===============================================================================
!> Schedule an ADIOS single precision global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
  subroutine write_adios_global_1d_real_5d(adios_handle, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  character(len=*) :: array_name
  real, dimension(:,:,:,:,:), intent(in) :: array
  ! Variables
  integer :: adios_err

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in write_adios_global_1d'

  ! adds local_dim/global_dim/offset infos
  call write_1D_global_array_adios_dims(adios_handle, myrank, local_dim, sizeprocs, array_name)

  call adios_write(adios_handle, trim(array_name) // "/array", array, adios_err)
  call check_adios_err(myrank,adios_err)

  end subroutine write_adios_global_1d_real_5d


!===============================================================================
!> Schedule an ADIOS double precision global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
  subroutine write_adios_global_1d_double_1d(adios_handle, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  character(len=*) :: array_name
  real(kind=8), dimension(:), intent(in) :: array
  ! Variables
  integer :: adios_err

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in write_adios_global_1d'

  ! adds local_dim/global_dim/offset infos
  call write_1D_global_array_adios_dims(adios_handle, myrank, local_dim, sizeprocs, array_name)

  call adios_write(adios_handle, trim(array_name) // "/array", array, adios_err)
  call check_adios_err(myrank,adios_err)

  end subroutine write_adios_global_1d_double_1d


!===============================================================================
!> Schedule an ADIOS double precision global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
  subroutine write_adios_global_1d_double_2d(adios_handle, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  character(len=*) :: array_name
  real(kind=8), dimension(:,:), intent(in) :: array
  ! Variables
  integer :: adios_err

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in write_adios_global_1d'

  ! adds local_dim/global_dim/offset infos
  call write_1D_global_array_adios_dims(adios_handle, myrank, local_dim, sizeprocs, array_name)

  call adios_write(adios_handle, trim(array_name) // "/array", array, adios_err)
  call check_adios_err(myrank,adios_err)

  end subroutine write_adios_global_1d_double_2d


!===============================================================================
!> Schedule an ADIOS double precision global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
  subroutine write_adios_global_1d_double_3d(adios_handle, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  character(len=*) :: array_name
  real(kind=8), dimension(:,:,:), intent(in) :: array
  ! Variables
  integer :: adios_err

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in write_adios_global_1d'

  ! adds local_dim/global_dim/offset infos
  call write_1D_global_array_adios_dims(adios_handle, myrank, local_dim, sizeprocs, array_name)

  call adios_write(adios_handle, trim(array_name) // "/array", array, adios_err)
  call check_adios_err(myrank,adios_err)

  end subroutine write_adios_global_1d_double_3d


!===============================================================================
!> Schedule an ADIOS double precision global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
  subroutine write_adios_global_1d_double_4d(adios_handle, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  character(len=*) :: array_name
  real(kind=8), dimension(:,:,:,:), intent(in) :: array
  ! Variables
  integer :: adios_err

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in write_adios_global_1d'

  ! adds local_dim/global_dim/offset infos
  call write_1D_global_array_adios_dims(adios_handle, myrank, local_dim, sizeprocs, array_name)

  call adios_write(adios_handle, trim(array_name) // "/array", array, adios_err)
  call check_adios_err(myrank,adios_err)

  end subroutine write_adios_global_1d_double_4d


!===============================================================================
!> Schedule an ADIOS double precision global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
  subroutine write_adios_global_1d_double_5d(adios_handle, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  character(len=*) :: array_name
  real(kind=8), dimension(:,:,:,:,:), intent(in) :: array
  ! Variables
  integer :: adios_err

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in write_adios_global_1d'

  ! adds local_dim/global_dim/offset infos
  call write_1D_global_array_adios_dims(adios_handle, myrank, local_dim, sizeprocs, array_name)

  call adios_write(adios_handle, trim(array_name) // "/array", array, adios_err)
  call check_adios_err(myrank,adios_err)

  end subroutine write_adios_global_1d_double_5d


!===============================================================================
!> Schedule an ADIOS integer global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
  subroutine write_adios_global_1d_integer_1d(adios_handle, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  character(len=*) :: array_name
  integer(kind=4), dimension(:), intent(in) :: array
  ! Variables
  integer :: adios_err

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in write_adios_global_1d'

  ! adds local_dim/global_dim/offset infos
  call write_1D_global_array_adios_dims(adios_handle, myrank, local_dim, sizeprocs, array_name)

  call adios_write(adios_handle, trim(array_name) // "/array", array, adios_err)
  call check_adios_err(myrank,adios_err)

  end subroutine write_adios_global_1d_integer_1d


!===============================================================================
!> Schedule an ADIOS integer global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
  subroutine write_adios_global_1d_integer_2d(adios_handle, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  character(len=*) :: array_name
  integer(kind=4), dimension(:,:), intent(in) :: array
  ! Variables
  integer :: adios_err

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in write_adios_global_1d'

  ! adds local_dim/global_dim/offset infos
  call write_1D_global_array_adios_dims(adios_handle, myrank, local_dim, sizeprocs, array_name)

  call adios_write(adios_handle, trim(array_name) // "/array", array, adios_err)
  call check_adios_err(myrank,adios_err)

  end subroutine write_adios_global_1d_integer_2d


!===============================================================================
!> Schedule an ADIOS integer global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
  subroutine write_adios_global_1d_integer_3d(adios_handle, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  character(len=*) :: array_name
  integer(kind=4), dimension(:,:,:), intent(in) :: array
  ! Variables
  integer :: adios_err

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in write_adios_global_1d'

  ! adds local_dim/global_dim/offset infos
  call write_1D_global_array_adios_dims(adios_handle, myrank, local_dim, sizeprocs, array_name)

  call adios_write(adios_handle, trim(array_name) // "/array", array, adios_err)
  call check_adios_err(myrank,adios_err)

  end subroutine write_adios_global_1d_integer_3d


!===============================================================================
!> Schedule an ADIOS integer global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
  subroutine write_adios_global_1d_integer_4d(adios_handle, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  character(len=*) :: array_name
  integer(kind=4), dimension(:,:,:,:), intent(in) :: array
  ! Variables
  integer :: adios_err

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in write_adios_global_1d'

  ! adds local_dim/global_dim/offset infos
  call write_1D_global_array_adios_dims(adios_handle, myrank, local_dim, sizeprocs, array_name)

  call adios_write(adios_handle, trim(array_name) // "/array", array, adios_err)
  call check_adios_err(myrank,adios_err)

  end subroutine write_adios_global_1d_integer_4d


!===============================================================================
!> Schedule an ADIOS integer global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
  subroutine write_adios_global_1d_integer_5d(adios_handle, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  character(len=*) :: array_name
  integer(kind=4), dimension(:,:,:,:,:), intent(in) :: array
  ! Variables
  integer :: adios_err

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in write_adios_global_1d'

  ! adds local_dim/global_dim/offset infos
  call write_1D_global_array_adios_dims(adios_handle, myrank, local_dim, sizeprocs, array_name)

  call adios_write(adios_handle, trim(array_name) // "/array", array, adios_err)
  call check_adios_err(myrank,adios_err)

  end subroutine write_adios_global_1d_integer_5d


!===============================================================================
!> Schedule an ADIOS long global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
  subroutine write_adios_global_1d_long_1d(adios_handle, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  character(len=*) :: array_name
  integer(kind=8), dimension(:), intent(in) :: array
  ! Variables
  integer :: adios_err

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in write_adios_global_1d'

  ! adds local_dim/global_dim/offset infos
  call write_1D_global_array_adios_dims(adios_handle, myrank, local_dim, sizeprocs, array_name)

  call adios_write(adios_handle, trim(array_name) // "/array", array, adios_err)
  call check_adios_err(myrank,adios_err)

  end subroutine write_adios_global_1d_long_1d


!===============================================================================
!> Schedule an ADIOS long global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
  subroutine write_adios_global_1d_long_2d(adios_handle, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  character(len=*) :: array_name
  integer(kind=8), dimension(:,:), intent(in) :: array
  ! Variables
  integer :: adios_err

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in write_adios_global_1d'

  ! adds local_dim/global_dim/offset infos
  call write_1D_global_array_adios_dims(adios_handle, myrank, local_dim, sizeprocs, array_name)

  call adios_write(adios_handle, trim(array_name) // "/array", array, adios_err)
  call check_adios_err(myrank,adios_err)

  end subroutine write_adios_global_1d_long_2d


!===============================================================================
!> Schedule an ADIOS long global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
  subroutine write_adios_global_1d_long_3d(adios_handle, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  character(len=*) :: array_name
  integer(kind=8), dimension(:,:,:), intent(in) :: array
  ! Variables
  integer :: adios_err

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in write_adios_global_1d'

  ! adds local_dim/global_dim/offset infos
  call write_1D_global_array_adios_dims(adios_handle, myrank, local_dim, sizeprocs, array_name)

  call adios_write(adios_handle, trim(array_name) // "/array", array, adios_err)
  call check_adios_err(myrank,adios_err)

  end subroutine write_adios_global_1d_long_3d


!===============================================================================
!> Schedule an ADIOS integer global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
  subroutine write_adios_global_1d_long_4d(adios_handle, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  character(len=*) :: array_name
  integer(kind=8), dimension(:,:,:,:), intent(in) :: array
  ! Variables
  integer :: adios_err

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in write_adios_global_1d'

  ! adds local_dim/global_dim/offset infos
  call write_1D_global_array_adios_dims(adios_handle, myrank, local_dim, sizeprocs, array_name)

  call adios_write(adios_handle, trim(array_name) // "/array", array, adios_err)
  call check_adios_err(myrank,adios_err)

  end subroutine write_adios_global_1d_long_4d


!===============================================================================
!> Schedule an ADIOS long global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
  subroutine write_adios_global_1d_long_5d(adios_handle, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  character(len=*) :: array_name
  integer(kind=8), dimension(:,:,:,:,:), intent(in) :: array
  ! Variables
  integer :: adios_err

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in write_adios_global_1d'

  ! adds local_dim/global_dim/offset infos
  call write_1D_global_array_adios_dims(adios_handle, myrank, local_dim, sizeprocs, array_name)

  call adios_write(adios_handle, trim(array_name) // "/array", array, adios_err)
  call check_adios_err(myrank,adios_err)

  end subroutine write_adios_global_1d_long_5d


!===============================================================================
!> Schedule an ADIOS logical global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
  subroutine write_adios_global_1d_logical_1d(adios_handle, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  character(len=*) :: array_name
  logical, dimension(:), intent(in) :: array
  ! Variables
  integer :: adios_err

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in write_adios_global_1d'

  ! adds local_dim/global_dim/offset infos
  call write_1D_global_array_adios_dims(adios_handle, myrank, local_dim, sizeprocs, array_name)

  call adios_write(adios_handle, trim(array_name) // "/array", array, adios_err)
  call check_adios_err(myrank,adios_err)

  end subroutine write_adios_global_1d_logical_1d


!===============================================================================
!> Schedule an ADIOS integer global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
  subroutine write_adios_global_1d_logical_2d(adios_handle, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  character(len=*) :: array_name
  logical, dimension(:,:), intent(in) :: array
  ! Variables
  integer :: adios_err

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in write_adios_global_1d'

  ! adds local_dim/global_dim/offset infos
  call write_1D_global_array_adios_dims(adios_handle, myrank, local_dim, sizeprocs, array_name)

  call adios_write(adios_handle, trim(array_name) // "/array", array, adios_err)
  call check_adios_err(myrank,adios_err)

  end subroutine write_adios_global_1d_logical_2d


!===============================================================================
!> Schedule an ADIOS integer global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
  subroutine write_adios_global_1d_logical_3d(adios_handle, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  character(len=*) :: array_name
  logical, dimension(:,:,:), intent(in) :: array
  ! Variables
  integer :: adios_err

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in write_adios_global_1d'

  ! adds local_dim/global_dim/offset infos
  call write_1D_global_array_adios_dims(adios_handle, myrank, local_dim, sizeprocs, array_name)

  call adios_write(adios_handle, trim(array_name) // "/array", array, adios_err)
  call check_adios_err(myrank,adios_err)

  end subroutine write_adios_global_1d_logical_3d


!===============================================================================
!> Schedule an ADIOS logical global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
  subroutine write_adios_global_1d_logical_4d(adios_handle, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  character(len=*) :: array_name
  logical, dimension(:,:,:,:), intent(in) :: array
  ! Variables
  integer :: adios_err

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in write_adios_global_1d'

  ! adds local_dim/global_dim/offset infos
  call write_1D_global_array_adios_dims(adios_handle, myrank, local_dim, sizeprocs, array_name)

  call adios_write(adios_handle, trim(array_name) // "/array", array, adios_err)
  call check_adios_err(myrank,adios_err)

  end subroutine write_adios_global_1d_logical_4d


!===============================================================================
!> Schedule an ADIOS logical global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
  subroutine write_adios_global_1d_logical_5d(adios_handle, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  character(len=*) :: array_name
  logical, dimension(:,:,:,:,:), intent(in) :: array
  ! Variables
  integer :: adios_err

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in write_adios_global_1d'

  ! adds local_dim/global_dim/offset infos
  call write_1D_global_array_adios_dims(adios_handle, myrank, local_dim, sizeprocs, array_name)

  call adios_write(adios_handle, trim(array_name) // "/array", array, adios_err)
  call check_adios_err(myrank,adios_err)

  end subroutine write_adios_global_1d_logical_5d

!===============================================================================

  subroutine write_1D_string_array_adios_dims(adios_handle, myrank, local_dim, global_dim, offset, sizeprocs, path)

  implicit none

  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: sizeprocs, myrank
  integer(kind=8), intent(in) :: local_dim
  integer(kind=8), intent(in) :: global_dim, offset
  character(len=*), intent(in) :: path

  integer :: adios_err
  integer :: idummy

  ! checks if file handle valid
  if (adios_handle == 0) stop 'Invalid ADIOS file handle in write_1D_string_array_adios_dims()'

  call adios_write(adios_handle, trim(path)// "/local_dim",local_dim, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_write(adios_handle, trim(path)// "/global_dim",global_dim, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_write(adios_handle, trim(path)// "/offset", offset, adios_err)
  call check_adios_err(myrank,adios_err)

  ! to avoid compiler warnings
  idummy = myrank
  idummy = sizeprocs

  end subroutine write_1D_string_array_adios_dims

!===============================================================================

!string subroutine added

  subroutine write_adios_global_1d_string_1d(adios_handle, myrank, sizeprocs, local_dim, global_dim, offset, &
                                             array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  integer(kind=8), intent(in) :: global_dim, offset
  character(len=*) :: array_name
  character(len=*), intent(in) :: array
  ! Variables
  integer :: adios_err

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in write_adios_global_1d_string_1d()'

  !debug
  !print *,"tag1: ",trim(array_name)," local_dim/global_dim/offset: ",local_dim,global_dim,offset

  ! adds local_dim/global_dim/offset infos
  call write_1D_string_array_adios_dims(adios_handle, myrank, local_dim, global_dim, offset, sizeprocs, array_name)

  !call write_1D_global_array_adios_dims(adios_handle, myrank, local_dim, sizeprocs, array_name)

  !debug
  !print *,"tag2: ",trim(array)

  call adios_write(adios_handle, trim(array_name)// "/array", array(1:local_dim), adios_err)
  call check_adios_err(myrank,adios_err)

  end subroutine write_adios_global_1d_string_1d

!===============================================================================
!
! with offset infos
!
!===============================================================================
  subroutine write_1D_global_array_adios_dims_offset(adios_handle, myrank, local_dim, global_dim, offset, sizeprocs, path)

  implicit none

  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: sizeprocs, myrank
  integer(kind=8), intent(in) :: local_dim
  integer(kind=8), intent(in) :: global_dim, offset
  character(len=*), intent(in) :: path

  ! local parameters
  integer :: adios_err
  integer :: idummy

  ! checks if file handle valid
  if (adios_handle == 0) stop 'Invalid ADIOS file handle in write_1D_string_array_adios_dims()'

  call adios_write(adios_handle, trim(path)// "/local_dim", local_dim, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_write(adios_handle, trim(path)// "/global_dim", global_dim, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_write(adios_handle, trim(path)// "/offset", offset, adios_err)
  call check_adios_err(myrank,adios_err)

  ! to avoid compiler warnings
  idummy = myrank
  idummy = sizeprocs

  end subroutine write_1D_global_array_adios_dims_offset

!===============================================================================
!> Schedule an ADIOS single precision global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
  subroutine write_adios_global_1d_real_1d_offset(adios_handle, myrank, sizeprocs, &
                                                  local_dim, global_dim, offset, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  integer(kind=8), intent(in) :: global_dim, offset
  character(len=*) :: array_name
  real, dimension(:), intent(in) :: array
  ! Variables
  integer :: adios_err

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in write_adios_global_1d_offset'

  ! adds local_dim/global_dim/offset infos
  call write_1D_global_array_adios_dims_offset(adios_handle, myrank,local_dim, global_dim, offset, sizeprocs, array_name)

  call adios_write(adios_handle, trim(array_name)// "/array", array, adios_err)
  call check_adios_err(myrank,adios_err)

  end subroutine write_adios_global_1d_real_1d_offset


!===============================================================================
!> Schedule an ADIOS double precision global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
  subroutine write_adios_global_1d_double_1d_offset(adios_handle, myrank, sizeprocs, &
                                                    local_dim, global_dim, offset, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  integer(kind=8), intent(in) :: global_dim, offset
  character(len=*) :: array_name
  real(kind=8), dimension(:), intent(in) :: array
  ! Variables
  integer :: adios_err

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in write_adios_global_1d_offset'

  ! adds local_dim/global_dim/offset infos
  call write_1D_global_array_adios_dims_offset(adios_handle, myrank, local_dim, global_dim, offset, sizeprocs, array_name)

  call adios_write(adios_handle, trim(array_name)// "/array", array, adios_err)
  call check_adios_err(myrank,adios_err)

  end subroutine write_adios_global_1d_double_1d_offset


!===============================================================================
!> Schedule an ADIOS integer global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
  subroutine write_adios_global_1d_integer_1d_offset(adios_handle, myrank, sizeprocs, &
                                                     local_dim, global_dim, offset, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  integer(kind=8), intent(in) :: global_dim, offset
  character(len=*) :: array_name
  integer(kind=4), dimension(:), intent(in) :: array
  ! Variables
  integer :: adios_err

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in write_adios_global_1d_offset'

  ! adds local_dim/global_dim/offset infos
  call write_1D_global_array_adios_dims_offset(adios_handle, myrank, local_dim, global_dim, offset, sizeprocs, array_name)

  call adios_write(adios_handle, trim(array_name)// "/array", array, adios_err)
  call check_adios_err(myrank,adios_err)

  end subroutine write_adios_global_1d_integer_1d_offset


!===============================================================================
!> Schedule an ADIOS long global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
  subroutine write_adios_global_1d_long_1d_offset(adios_handle, myrank, sizeprocs, &
                                                  local_dim, global_dim, offset, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  integer(kind=8), intent(in) :: global_dim, offset
  character(len=*) :: array_name
  integer(kind=8), dimension(:), intent(in) :: array
  ! Variables
  integer :: adios_err

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in write_adios_global_1d_offset'

  ! adds local_dim/global_dim/offset infos
  call write_1D_global_array_adios_dims_offset(adios_handle, myrank, local_dim, global_dim, offset, sizeprocs, array_name)

  call adios_write(adios_handle, trim(array_name)// "/array", array, adios_err)
  call check_adios_err(myrank,adios_err)

  end subroutine write_adios_global_1d_long_1d_offset


!===============================================================================
!> Schedule an ADIOS logical global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
  subroutine write_adios_global_1d_logical_1d_offset(adios_handle, myrank, sizeprocs, &
                                                     local_dim, global_dim, offset, array_name, array)
  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  integer(kind=8), intent(in) :: global_dim, offset
  character(len=*) :: array_name
  logical, dimension(:), intent(in) :: array
  ! Variables
  integer :: adios_err

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in write_adios_global_1d_offset'

  ! adds local_dim/global_dim/offset infos
  call write_1D_global_array_adios_dims_offset(adios_handle, myrank, local_dim, global_dim, offset, sizeprocs, array_name)

  call adios_write(adios_handle, trim(array_name)// "/array", array, adios_err)
  call check_adios_err(myrank,adios_err)

  end subroutine write_adios_global_1d_logical_1d_offset


!===============================================================================
!
! error checking
!
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
    print *, "process ", myrank, "has ADIOS error ",adios_err,' ',trim(msg)
    stop 'adios error'
  endif

  end subroutine check_adios_err

end module adios_helpers_writers_mod
