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

!-------------------------------------------------------------------------------
!> Helpers to set up adios features.
!  * writing out
!-------------------------------------------------------------------------------

#include "config.fh"


module adios_helpers_writers_mod

#if defined(USE_ADIOS2)
  use adios2
#endif

  use manager_adios, only: check_adios_err

  implicit none

  ! writing
  public :: write_adios_perform
  public :: write_adios_begin_step
  public :: write_adios_end_step

  public :: write_adios_scalar
  interface write_adios_scalar
    module procedure write_adios_scalar_int
    module procedure write_adios_scalar_real
    module procedure write_adios_scalar_double
  end interface write_adios_scalar

  public :: write_adios_array_gll

  public :: write_adios_global_real_1d_array
  public :: write_adios_global_double_1d_array
  public :: write_adios_global_integer_1d_array
  public :: write_adios_global_long_1d_array
  public :: write_adios_global_string_1d_array
  public :: write_adios_global_1d_array
  public :: write_adios_global_logical_1d_array

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
    module procedure write_adios_global_1d_logical_1d       ! only used for ispec_is_tiso
    !module procedure write_adios_global_1d_logical_2d
    !module procedure write_adios_global_1d_logical_3d
    !module procedure write_adios_global_1d_logical_4d
    !module procedure write_adios_global_1d_logical_5d
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

    module procedure write_adios_global_1d_logical_1d   ! only used for ispec_is_tiso
    !module procedure write_adios_global_1d_logical_2d
    !module procedure write_adios_global_1d_logical_3d
    !module procedure write_adios_global_1d_logical_4d
    !module procedure write_adios_global_1d_logical_5d

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
  ! not used yet...
  public :: write_adios_global_1d_array_offset
  interface write_adios_global_1d_array_offset
    module procedure write_adios_global_1d_integer_1d_offset
    module procedure write_adios_global_1d_long_1d_offset
    !module procedure write_adios_global_1d_logical_1d_offset
    module procedure write_adios_global_1d_real_1d_offset
    module procedure write_adios_global_1d_double_1d_offset
  end interface write_adios_global_1d_array_offset

  private

#if defined(USE_ADIOS2)
  ! default adios2 write mode
  ! sync
  !integer, parameter :: myadios_writer_mode = adios2_mode_sync
  ! deferred (default)
  integer, parameter :: myadios_writer_mode = adios2_mode_deferred
#endif

contains


!---------------------------------------------------------------------------------
!
! wrappers
!
!---------------------------------------------------------------------------------

  subroutine write_adios_perform(adios_handle)

! actual performs writing out to file

  implicit none

#if defined(USE_ADIOS)
  integer(kind=8), intent(in) :: adios_handle
#elif defined(USE_ADIOS2)
  type(adios2_engine), intent(in) :: adios_handle
#endif
  ! local parameters
  integer :: ier

  TRACE_ADIOS_L2('write_adios_perform')

! note: for adios 1, the adios_write() routine will copy the data to buffer
!       (if large enough, otherwise just does the i/o directly).
!       the memory pointed to by the variable can be re-used when the function returns.
!       if the buffer is large enough, the actual i/o is performed when closing the file.
!
!       for adios 2, the adios2_put() routine has two basic modes, either deferred or sync.
!       in deferred mode, data is collected until perform/close/end_step and memory pointers should not be
!       modified. in sync mode, data is reusable immediately after the call.
!
!       here, we enforce the i/o writing which in most cases should not be needed.

  ! Perform the actual write to disk
#if defined(USE_ADIOS)
  ! note: adios1 has no perform_puts() command, it will perform the writing only when closing a file or
  !       when we reset the path here.
  !
  !       thus, careful when performing with this command.
  !       this will close the path and subsequent writes would have
  !       lost that information which can cause problems as we sometimes use the path to designate a variable name.
  !       only call at the end of a group write, for example when re-setting the path before creating a new group.
  call adios_set_path(adios_handle, '', ier)
  call check_adios_err(ier,"Error adios could not perform write with set path to zero")

#elif defined(USE_ADIOS2)
  call adios2_perform_puts(adios_handle, ier)
  call check_adios_err(ier,"Error adios2 perform puts in write_adios_perform() routine")

#endif

  end subroutine write_adios_perform

!
!---------------------------------------------------------------------------------
!

  subroutine write_adios_begin_step(adios_handle)

! actual performs writing out to file

  implicit none

#if defined(USE_ADIOS)
  integer(kind=8), intent(in) :: adios_handle
#elif defined(USE_ADIOS2)
  type(adios2_engine), intent(in) :: adios_handle
#endif
  ! local parameters
  integer :: ier

  TRACE_ADIOS_L2('write_adios_begin_step')

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! version 1 has no begin_step routine, we close and append data to the file
  ! just check to avoid compiler warning
  if (adios_handle == 0) stop 'Error invalid adios handle in adios_begin_step'

  ! to avoid compiler warning
  ier = 0

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! moves to next step, starts at 0
  call adios2_begin_step(adios_handle, ier)
  call check_adios_err(ier,"Error adios2 in write_adios_begin_step() routine")

#endif

  end subroutine write_adios_begin_step

!
!---------------------------------------------------------------------------------
!

  subroutine write_adios_end_step(adios_handle)

! actual performs writing out to file

  implicit none

#if defined(USE_ADIOS)
  integer(kind=8), intent(in) :: adios_handle
#elif defined(USE_ADIOS2)
  type(adios2_engine), intent(in) :: adios_handle
#endif
  ! local parameters
  integer :: ier

  TRACE_ADIOS_L2('write_adios_end_step')

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! version 1 has no end_step routine, we close and append data to the file

  ! note: do not close path (as in write_adios_perform()) as it would corrupt the group structure when appending
  !       the next step
  !       this would lead to an error: call adios_set_path(adios_handle, '', ier)
  !       instead, just close the file. closing will perform the write.

  ! check
  if (adios_handle == 0) stop 'Invalid adios file handle in write_adios_end_step'

  ! to avoid compiler warning
  ier = 0

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! end step to indicate output is completed. ADIOS2 can do I/O
  call adios2_end_step(adios_handle, ier)
  call check_adios_err(ier,"Error adios2 in write_adios_end_step() routine")

#endif

  end subroutine write_adios_end_step


!---------------------------------------------------------------------------------
!
! scalars
!
!---------------------------------------------------------------------------------

  subroutine write_adios_scalar_int(adios_handle,adios_group,scalar_name,scalar)

! writes a single scalar

  implicit none

#if defined(USE_ADIOS)
  integer(kind=8), intent(in) :: adios_handle
  integer(kind=8), intent(in) :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_engine), intent(in) :: adios_handle
  type(adios2_io), intent(in) :: adios_group
  ! variable
  type(adios2_variable) :: v
#endif
  integer, intent(in) :: scalar
  character(len=*), intent(in) :: scalar_name
  ! local parameters
  integer :: ier

  TRACE_ADIOS_L2_ARG('write_adios_scalar_int: ',trim(scalar_name))

  ! checks name
  if (len_trim(scalar_name) == 0) stop 'Error: scalar_name has zero length in write_adios_scalar_int()'

  ! writes scalar (either buffered or directly to disk)
#if defined(USE_ADIOS)
  ! ADIOS 1
  ! checks if file handle valid
  if (adios_handle == 0) stop 'Invalid ADIOS file handle in write_adios_scalar_int()'

  call adios_write(adios_handle, trim(scalar_name), scalar, ier)
  call check_adios_err(ier,"Error adios could not write parameter: "//trim(scalar_name))

  ! to avoid compiler warning
  ier = adios_group

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! checks if file handle valid
  if (.not. adios_handle%valid) stop 'Invalid ADIOS2 file handle in write_adios_scalar_int()'

  ! gets associated variable for array
  call adios2_inquire_variable(v, adios_group, trim(scalar_name),ier)
  call check_adios_err(ier,"Error adios2 write_adios_scalar_int(): inquire variable '"//trim(scalar_name)//"' failed")
  if (.not. v%valid) stop 'Error adios2 variable invalid'

  call adios2_put(adios_handle, v, scalar, adios2_mode_sync, ier)
  call check_adios_err(ier,"Error adios could not write parameter: "//trim(scalar_name))

#endif

  end subroutine write_adios_scalar_int

!===============================================================================

  subroutine write_adios_scalar_real(adios_handle,adios_group,scalar_name,scalar)

! writes a real-value scalar

  implicit none

#if defined(USE_ADIOS)
  integer(kind=8), intent(in) :: adios_handle
  integer(kind=8), intent(in) :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_engine), intent(in) :: adios_handle
  type(adios2_io), intent(in) :: adios_group
  ! variable
  type(adios2_variable) :: v
#endif
  real, intent(in) :: scalar
  character(len=*), intent(in) :: scalar_name
  ! local parameters
  integer :: ier

  TRACE_ADIOS_L2_ARG('write_adios_scalar_real: ',trim(scalar_name))

  ! checks name
  if (len_trim(scalar_name) == 0) stop 'Error: scalar_name has zero length in write_adios_scalar_real()'

  ! writes scalar (either buffered or directly to disk)
#if defined(USE_ADIOS)
  ! ADIOS 1
  ! checks if file handle valid
  if (adios_handle == 0) stop 'Invalid ADIOS file handle in write_adios_scalar_real()'

  call adios_write(adios_handle, trim(scalar_name), scalar, ier)
  call check_adios_err(ier,"Error adios could not write parameter: "//trim(scalar_name))

  ! to avoid compiler warning
  ier = adios_group

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! checks if file handle valid
  if (.not. adios_handle%valid) stop 'Invalid ADIOS2 file handle in write_adios_scalar_real()'

  ! gets associated variable for array
  call adios2_inquire_variable(v, adios_group, trim(scalar_name),ier)
  call check_adios_err(ier,"Error adios2 write_adios_scalar_real(): inquire variable '"//trim(scalar_name)//"' failed")
  if (.not. v%valid) stop 'Error adios2 variable invalid'

  call adios2_put(adios_handle, v, scalar, adios2_mode_sync, ier)
  call check_adios_err(ier,"Error adios could not write parameter: "//trim(scalar_name))

#endif

  end subroutine write_adios_scalar_real

!===============================================================================

  subroutine write_adios_scalar_double(adios_handle,adios_group,scalar_name,scalar)

! writes a single scalar

  implicit none

#if defined(USE_ADIOS)
  integer(kind=8), intent(in) :: adios_handle
  integer(kind=8), intent(in) :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_engine), intent(in) :: adios_handle
  type(adios2_io), intent(in) :: adios_group
  ! variable
  type(adios2_variable) :: v
#endif
  double precision, intent(in) :: scalar
  character(len=*), intent(in) :: scalar_name

  ! local parameters
  integer :: ier

  TRACE_ADIOS_L2_ARG('write_adios_scalar_double: ',trim(scalar_name))

  ! checks name
  if (len_trim(scalar_name) == 0) stop 'Error: scalar_name has zero length in write_adios_scalar_double()'

  ! writes scalar (either buffered or directly to disk)
#if defined(USE_ADIOS)
  ! ADIOS 1
  ! checks if file handle valid
  if (adios_handle == 0) stop 'Invalid ADIOS file handle in write_adios_scalar_double()'

  call adios_write(adios_handle, trim(scalar_name), scalar, ier)
  call check_adios_err(ier,"Error adios could not write parameter: "//trim(scalar_name))

  ! to avoid compiler warning
  ier = adios_group

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! checks if file handle valid
  if (.not. adios_handle%valid) stop 'Invalid ADIOS2 file handle in write_adios_scalar_double()'

  ! gets associated variable for array
  call adios2_inquire_variable(v, adios_group, trim(scalar_name), ier)
  call check_adios_err(ier,"Error adios2 write_adios_scalar_double(): inquire variable "//trim(scalar_name)//" failed")
  if (.not. v%valid) stop 'Error adios2 variable invalid'

  call adios2_put(adios_handle, v, scalar, adios2_mode_sync, ier)
  call check_adios_err(ier,"Error adios could not write parameter: "//trim(scalar_name))

#endif

  end subroutine write_adios_scalar_double


!===============================================================================
!
! GLL array writing
!
!===============================================================================

  subroutine write_adios_array_gll(adios_handle,adios_group,myrank,sizeprocs,nspec,array_name,array_gll)

! writes a GLL array

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

  implicit none

#if defined(USE_ADIOS)
  integer(kind=8), intent(in) :: adios_handle
  integer(kind=8), intent(in) :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_engine), intent(in) :: adios_handle
  type(adios2_io), intent(in) :: adios_group
#endif
  integer, intent(in) :: myrank,sizeprocs
  integer, intent(in) :: nspec

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec), intent(in) :: array_gll
  character(len=*), intent(in) :: array_name

  ! local parameters
  integer(kind=8) :: local_dim

  TRACE_ADIOS_L2_ARG('write_adios_array_gll: ',trim(array_name))

  ! checks if file handle valid
#if defined(USE_ADIOS)
  if (adios_handle == 0) stop 'Invalid ADIOS file handle in write_adios_array_gll()'
#elif defined(USE_ADIOS2)
  if (.not. adios_handle%valid) stop 'Invalid ADIOS2 file handle in write_adios_array_gll()'
#endif

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in write_adios_array_gll()'

  ! array size
  local_dim = NGLLX * NGLLY * NGLLZ * nspec

  ! writes array
  call write_adios_global_1d_array(adios_handle,adios_group,myrank,sizeprocs,local_dim,trim(array_name),array_gll)

  end subroutine write_adios_array_gll


!===============================================================================
!
! generic arrays
!
!===============================================================================

  subroutine write_1D_global_array_adios_dims(adios_handle, adios_group, myrank, local_dim, sizeprocs, path, array_size)

  implicit none

#if defined(USE_ADIOS)
  integer(kind=8), intent(in) :: adios_handle
  integer(kind=8), intent(in) :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_engine), intent(in) :: adios_handle
  type(adios2_io), intent(in) :: adios_group
  type(adios2_variable) :: v
#endif
  integer, intent(in) :: sizeprocs, myrank
  integer(kind=8), intent(in) :: local_dim
  character(len=*), intent(in) :: path
  integer(kind=8), intent(in) :: array_size

  integer :: adios_err
  integer(kind=8) :: offset
  integer(kind=8) :: global_dim ! should be long for large cases.

  TRACE_ADIOS_L2_ARG('write_1D_global_array_adios_dims: ',trim(path))

  ! global dimension
  global_dim = local_dim * int(sizeprocs,kind=8)

  ! process offset
  ! note: assumes that myrank starts is within [0,sizeprocs-1]
  offset = local_dim * int(myrank,kind=8)

#if defined(USE_ADIOS)
  ! ADIOS 1
  call adios_write(adios_handle, trim(path) // "/local_dim", local_dim, adios_err)
  call check_adios_err(adios_err,"Error writing "//trim(path) // "/local_dim")

  call adios_write(adios_handle, trim(path) // "/global_dim", global_dim, adios_err)
  call check_adios_err(adios_err,"Error writing "//trim(path) // "/global_dim")

  call adios_write(adios_handle, trim(path) // "/offset", offset, adios_err)
  call check_adios_err(adios_err,"Error writing "//trim(path) // "/offset")

  call adios_write(adios_handle, trim(path) // "/size", array_size, adios_err)
  call check_adios_err(adios_err,"Error writing "//trim(path) // "/size")

  ! to avoid compiler warning
  offset = adios_group

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! note: we add these additional values for storing an array. we will use use offset and global_dim info
  !       to specify array selections and/or estimate sizes.

  ! synchronizations:
  ! sync mode
  ! > call adios2_put(adios_handle, v, local_dim, adios2_mode_sync, adios_err)
  ! deferred mode
  ! > call adios2_put(adios_handle, v, local_dim, adios_err)
  !
  ! we will try to use deferred mode to group together different calls.
  ! however, this requires the data and data pointer to be still valid and unmodified
  ! until a perform_puts/close/end_step call is done.

  ! local_dim
  ! gets associated variable for array
  call adios2_inquire_variable(v, adios_group, trim(path) // "/local_dim", adios_err)
  call check_adios_err(adios_err,"Error adios2 write global dims: inquire variable "//trim(path) // "/local_dim"//" failed")
  if (.not. v%valid) stop 'Error adios2 variable invalid'

  call adios2_put(adios_handle, v, local_dim, adios_err)
  call check_adios_err(adios_err,"Error adios could not write parameter: "//trim(path) // "/local_dim")

  ! global_dim
  call adios2_inquire_variable(v, adios_group, trim(path) // "/global_dim", adios_err)
  call check_adios_err(adios_err,"Error adios2 write global dims: inquire variable "//trim(path) // "/global_dim"//" failed")
  if (.not. v%valid) stop 'Error adios2 variable invalid'

  call adios2_put(adios_handle, v, global_dim, adios_err)
  call check_adios_err(adios_err,"Error adios could not write parameter: "//trim(path) // "/global_dim")

  ! offset
  call adios2_inquire_variable(v, adios_group, trim(path) // "/offset", adios_err)
  call check_adios_err(adios_err,"Error adios2 write global dims: inquire variable "//trim(path) // "/offset"//" failed")
  if (.not. v%valid) stop 'Error adios2 variable invalid'

  call adios2_put(adios_handle, v, offset, adios_err)
  call check_adios_err(adios_err,"Error adios could not write parameter: "//trim(path) // "/offset")

  ! size
  call adios2_inquire_variable(v, adios_group, trim(path) // "/size", adios_err)
  call check_adios_err(adios_err,"Error adios2 write global dims: inquire variable "//trim(path) // "/size"//" failed")
  if (.not. v%valid) stop 'Error adios2 variable invalid'

  call adios2_put(adios_handle, v, array_size, adios_err)
  call check_adios_err(adios_err,"Error adios could not write parameter: "//trim(path) // "/size")

  ! sync write
  call adios2_perform_puts(adios_handle, adios_err)
  call check_adios_err(adios_err,"Error adios2 perform puts in write_adios_perform() routine")

#endif

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
  subroutine write_adios_global_1d_real_1d(adios_handle, adios_group, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  real, dimension(:), intent(in) :: array

#include "adios_helpers_writers_1d_generic.inc"

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
  subroutine write_adios_global_1d_real_2d(adios_handle, adios_group, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  real, dimension(:,:), intent(in) :: array

#include "adios_helpers_writers_1d_generic.inc"

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
  subroutine write_adios_global_1d_real_3d(adios_handle, adios_group, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  real, dimension(:,:,:), intent(in) :: array

#include "adios_helpers_writers_1d_generic.inc"

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
  subroutine write_adios_global_1d_real_4d(adios_handle, adios_group, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  real, dimension(:,:,:,:), intent(in) :: array

#include "adios_helpers_writers_1d_generic.inc"

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
  subroutine write_adios_global_1d_real_5d(adios_handle, adios_group, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  real, dimension(:,:,:,:,:), intent(in) :: array

#include "adios_helpers_writers_1d_generic.inc"

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
  subroutine write_adios_global_1d_double_1d(adios_handle, adios_group, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  real(kind=8), dimension(:), intent(in) :: array

#include "adios_helpers_writers_1d_generic.inc"

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
  subroutine write_adios_global_1d_double_2d(adios_handle, adios_group, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  real(kind=8), dimension(:,:), intent(in) :: array

#include "adios_helpers_writers_1d_generic.inc"

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
  subroutine write_adios_global_1d_double_3d(adios_handle, adios_group, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  real(kind=8), dimension(:,:,:), intent(in) :: array

#include "adios_helpers_writers_1d_generic.inc"

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
  subroutine write_adios_global_1d_double_4d(adios_handle, adios_group, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  real(kind=8), dimension(:,:,:,:), intent(in) :: array

#include "adios_helpers_writers_1d_generic.inc"

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
  subroutine write_adios_global_1d_double_5d(adios_handle, adios_group, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  real(kind=8), dimension(:,:,:,:,:), intent(in) :: array

#include "adios_helpers_writers_1d_generic.inc"

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
  subroutine write_adios_global_1d_integer_1d(adios_handle, adios_group, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=4), dimension(:), intent(in) :: array

#include "adios_helpers_writers_1d_generic.inc"

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
  subroutine write_adios_global_1d_integer_2d(adios_handle, adios_group, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=4), dimension(:,:), intent(in) :: array

#include "adios_helpers_writers_1d_generic.inc"

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
  subroutine write_adios_global_1d_integer_3d(adios_handle, adios_group, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=4), dimension(:,:,:), intent(in) :: array

#include "adios_helpers_writers_1d_generic.inc"

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
  subroutine write_adios_global_1d_integer_4d(adios_handle, adios_group, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=4), dimension(:,:,:,:), intent(in) :: array

#include "adios_helpers_writers_1d_generic.inc"

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
  subroutine write_adios_global_1d_integer_5d(adios_handle, adios_group, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=4), dimension(:,:,:,:,:), intent(in) :: array

#include "adios_helpers_writers_1d_generic.inc"

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
  subroutine write_adios_global_1d_long_1d(adios_handle, adios_group, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), dimension(:), intent(in) :: array

#include "adios_helpers_writers_1d_generic.inc"

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
  subroutine write_adios_global_1d_long_2d(adios_handle, adios_group, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), dimension(:,:), intent(in) :: array

#include "adios_helpers_writers_1d_generic.inc"

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
  subroutine write_adios_global_1d_long_3d(adios_handle, adios_group, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), dimension(:,:,:), intent(in) :: array

#include "adios_helpers_writers_1d_generic.inc"

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
  subroutine write_adios_global_1d_long_4d(adios_handle, adios_group, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), dimension(:,:,:,:), intent(in) :: array

#include "adios_helpers_writers_1d_generic.inc"

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
  subroutine write_adios_global_1d_long_5d(adios_handle, adios_group, myrank, sizeprocs, local_dim, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), dimension(:,:,:,:,:), intent(in) :: array

#include "adios_helpers_writers_1d_generic.inc"

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
!
! note: adios1 seems to handle the adios_write routine a bit different than adios2 which checks the array type.
!       both don't support logical, but store it as integer*4 size.
!
! adios1 would do this implicitly and we could call:
!  subroutine write_adios_global_1d_logical_1d(adios_handle, adios_group, myrank, sizeprocs, local_dim, array_name, array)
!
!  implicit none
!  ! Parameters
!  logical, dimension(:), intent(in) :: array
!
!#include "adios_helpers_writers_1d_generic.inc"
!
!  end subroutine write_adios_global_1d_logical_1d
!
! adios2 has no adios_put() routine for logicals.
! we use explicit integer array instead for storing/reading:

  subroutine write_adios_global_1d_logical_1d(adios_handle, adios_group, myrank, sizeprocs, local_dim, array_name, array_l)

  implicit none
  ! Parameters
  logical, dimension(:), intent(in) :: array_l
  integer, dimension(size(array_l)) :: array

!modified from:
!#include "adios_helpers_writers_1d_generic.inc"

  ! common parameters
#if defined(USE_ADIOS)
  integer(kind=8), intent(in) :: adios_handle
  integer(kind=8), intent(in) :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_engine), intent(in) :: adios_handle
  type(adios2_io), intent(in) :: adios_group
  type(adios2_variable) :: v
#endif
  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  character(len=*) :: array_name
  ! Variables
  integer :: adios_err
  integer(kind=8) :: array_size

  TRACE_ADIOS_L2_ARG('write_adios_global_1d_logical_1d: ',trim(array_name))

  ! explicitly converts logical to integer array
  array(:) = 0
  where(array_l(:) .eqv. .true.) array(:) = 1

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in adios_helpers_writers_1d_generic'

  ! gets array size
  array_size = size(array,kind=8)

  ! adds local_dim/global_dim/offset infos
  call write_1D_global_array_adios_dims(adios_handle, adios_group, myrank, local_dim, sizeprocs, array_name, array_size)

#if defined(USE_ADIOS)
  ! ADIOS 1
  call adios_write(adios_handle, trim(array_name) // "/array", array, adios_err)
  call check_adios_err(adios_err,"Error writing adios array "//trim(array_name) // "/array")

! note: performing the write will be done when closing file. do not call set_path to close path as it will corrupt
!       the group structure in subsequent write calls.

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! gets associated variable for array from current group
  call adios2_inquire_variable(v, adios_group, trim(array_name)// "/array", adios_err)
  call check_adios_err(adios_err,"Error adios2 write logical_1d: inquire variable "//trim(array_name)// "/array" //" failed")

  if (.not. v%valid) stop 'Error adios2 variable invalid'

  ! uses synchronized put
  ! since we might loose the temporary scope of the array before we do a perform/close/end_step call
  call adios2_put(adios_handle, v, array, adios2_mode_sync, adios_err)
  call check_adios_err(adios_err,"Error writing adios2 array "//trim(array_name) // "/array")

#endif

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
! unused so far...
!  subroutine write_adios_global_1d_logical_2d(adios_handle, adios_group, myrank, sizeprocs, local_dim, array_name, array)
!
!  implicit none
!  ! Parameters
!  logical, dimension(:,:), intent(in) :: array
!
!#include "adios_helpers_writers_1d_generic.inc"
!
!  end subroutine write_adios_global_1d_logical_2d
!
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
! unused so far...
!  subroutine write_adios_global_1d_logical_3d(adios_handle, adios_group, myrank, sizeprocs, local_dim, array_name, array)
!
!  implicit none
!  ! Parameters
!  logical, dimension(:,:,:), intent(in) :: array
!
!#include "adios_helpers_writers_1d_generic.inc"
!
!  end subroutine write_adios_global_1d_logical_3d

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
! unused so far...
!  subroutine write_adios_global_1d_logical_4d(adios_handle, adios_group, myrank, sizeprocs, local_dim, array_name, array)
!
!  implicit none
!  ! Parameters
!  logical, dimension(:,:,:,:), intent(in) :: array
!
!#include "adios_helpers_writers_1d_generic.inc"
!
!  end subroutine write_adios_global_1d_logical_4d

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
! unused so far...
!  subroutine write_adios_global_1d_logical_5d(adios_handle, adios_group, myrank, sizeprocs, local_dim, array_name, array)
!
!  implicit none
!  ! Parameters
!  logical, dimension(:,:,:,:,:), intent(in) :: array
!
!#include "adios_helpers_writers_1d_generic.inc"
!
!  end subroutine write_adios_global_1d_logical_5d

!===============================================================================

  subroutine write_1D_string_array_adios_dims(adios_handle, adios_group, myrank, local_dim, global_dim, offset, sizeprocs, &
                                              path, array_size)

  implicit none

#if defined(USE_ADIOS)
  integer(kind=8), intent(in) :: adios_handle
  integer(kind=8), intent(in) :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_engine), intent(in) :: adios_handle
  type(adios2_io), intent(in) :: adios_group
  type(adios2_variable) :: v
#endif

  integer, intent(in) :: sizeprocs, myrank
  integer(kind=8), intent(in) :: local_dim
  integer(kind=8), intent(in) :: global_dim, offset
  character(len=*), intent(in) :: path
  integer(kind=8), intent(in) :: array_size

  integer :: adios_err
  integer :: idummy

  TRACE_ADIOS_L2_ARG('write_1D_string_array_adios_dims: ',trim(path))

#if defined(USE_ADIOS)
  ! ADIOS 1
  call adios_write(adios_handle, trim(path)// "/local_dim",local_dim, adios_err)
  call check_adios_err(adios_err,"Error writing "//trim(path) // "/local_dim")

  call adios_write(adios_handle, trim(path)// "/global_dim",global_dim, adios_err)
  call check_adios_err(adios_err,"Error writing "//trim(path) // "/global_dim")

  call adios_write(adios_handle, trim(path)// "/offset", offset, adios_err)
  call check_adios_err(adios_err,"Error writing "//trim(path) // "/offset")

  call adios_write(adios_handle, trim(path)// "/size", array_size, adios_err)
  call check_adios_err(adios_err,"Error writing "//trim(path) // "/size")

  ! to avoid compiler warning
  idummy = adios_group

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! note: not sure if we really need these additional values for storing an array

  ! local_dim
  ! gets associated variable for array
  call adios2_inquire_variable(v, adios_group, trim(path) // "/local_dim", adios_err)
  call check_adios_err(adios_err,"Error adios2 write string dims: inquire variable "//trim(path) // "/local_dim"//" failed")
  if (.not. v%valid) stop 'Error adios2 variable invalid'

  call adios2_put(adios_handle, v, local_dim, adios_err)
  call check_adios_err(adios_err,"Error adios could not write parameter: "//trim(path) // "/local_dim")

  ! global_dim
  call adios2_inquire_variable(v, adios_group, trim(path) // "/global_dim", adios_err)
  call check_adios_err(adios_err,"Error adios2 write string dims: inquire variable "//trim(path) // "/global_dim"//" failed")
  if (.not. v%valid) stop 'Error adios2 variable invalid'

  call adios2_put(adios_handle, v, global_dim, adios_err)
  call check_adios_err(adios_err,"Error adios could not write parameter: "//trim(path) // "/global_dim")

  ! offset
  call adios2_inquire_variable(v, adios_group, trim(path) // "/offset", adios_err)
  call check_adios_err(adios_err,"Error adios2 write string dims: inquire variable "//trim(path) // "/offset"//" failed")
  if (.not. v%valid) stop 'Error adios2 variable invalid'

  call adios2_put(adios_handle, v, offset, adios_err)
  call check_adios_err(adios_err,"Error adios could not write parameter: "//trim(path) // "/offset")

  ! size
  call adios2_inquire_variable(v, adios_group, trim(path) // "/size", adios_err)
  call check_adios_err(adios_err,"Error adios2 write string dims: inquire variable "//trim(path) // "/size"//" failed")
  if (.not. v%valid) stop 'Error adios2 variable invalid'

  call adios2_put(adios_handle, v, array_size, adios_err)
  call check_adios_err(adios_err,"Error adios could not write parameter: "//trim(path) // "/size")

  ! sync write
  call adios2_perform_puts(adios_handle, adios_err)
  call check_adios_err(adios_err,"Error adios2 perform puts in write_adios_perform() routine")

#endif

  ! to avoid compiler warnings
  idummy = myrank
  idummy = sizeprocs

  end subroutine write_1D_string_array_adios_dims

!===============================================================================

!string subroutine added

  subroutine write_adios_global_1d_string_1d(adios_handle, adios_group, myrank, sizeprocs, &
                                             local_dim, global_dim, offset, array_name, array)

  implicit none
  ! Parameters
#if defined(USE_ADIOS)
  integer(kind=8), intent(in) :: adios_handle
  integer(kind=8), intent(in) :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_engine), intent(in) :: adios_handle
  type(adios2_io), intent(in) :: adios_group
  type(adios2_variable) :: v
#endif

  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  integer(kind=8), intent(in) :: global_dim, offset
  character(len=*) :: array_name
  character(len=*), intent(in) :: array
  ! Variables
  integer :: adios_err
  integer(kind=8) :: array_size

  TRACE_ADIOS_L2_ARG('write_adios_global_1d_string_1d: ',trim(array_name))

  ! debug
  !print *,"tag2:",trim(array_name)
  !print *,"tag2:",trim(array)

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in write_adios_global_1d_string_1d()'

  ! gets size
  array_size = len(array)

  ! adds local_dim/global_dim/offset infos
  call write_1D_string_array_adios_dims(adios_handle, adios_group, myrank, local_dim, global_dim, offset, sizeprocs, &
                                        array_name, array_size)

#if defined(USE_ADIOS)
  ! ADIOS 1
  call adios_write(adios_handle, trim(array_name)// "/array", array(1:local_dim), adios_err)
  call check_adios_err(adios_err,"Error writing "//trim(array_name) // "/array")

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! note: not sure if we really need these additional values for storing an array

  ! local_dim
  ! gets associated variable for array
  call adios2_inquire_variable(v, adios_group, trim(array_name)// "/array", adios_err)
  call check_adios_err(adios_err,"Error adios2 write 1d string: inquire variable "//trim(array_name)// "/array" //" failed")

  if (.not. v%valid) stop 'Error adios2 variable invalid'

  call adios2_put(adios_handle, v, array(1:local_dim), adios_err)
  call check_adios_err(adios_err,"Error adios could not write parameter: "//trim(array_name)// "/array")

#endif

  end subroutine write_adios_global_1d_string_1d

!===============================================================================
!
! with offset infos
!
!===============================================================================
  subroutine write_1D_global_array_adios_dims_offset(adios_handle, adios_group, myrank, &
                                                     local_dim, global_dim, offset, sizeprocs, path, array_size)

  implicit none

#if defined(USE_ADIOS)
  integer(kind=8), intent(in) :: adios_handle
  integer(kind=8), intent(in) :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_engine), intent(in) :: adios_handle
  type(adios2_io), intent(in) :: adios_group
  type(adios2_variable) :: v
#endif
  integer, intent(in) :: sizeprocs, myrank
  integer(kind=8), intent(in) :: local_dim
  integer(kind=8), intent(in) :: global_dim, offset
  character(len=*), intent(in) :: path
  integer(kind=8), intent(in) :: array_size

  ! local parameters
  integer :: adios_err
  integer :: idummy

  TRACE_ADIOS_L2_ARG('write_1D_global_array_adios_dims_offset: ',trim(path))

#if defined(USE_ADIOS)
  ! ADIOS 1
  call adios_write(adios_handle, trim(path)// "/local_dim", local_dim, adios_err)
  call check_adios_err(adios_err,"Error writing "//trim(path) // "/local_dim")

  call adios_write(adios_handle, trim(path)// "/global_dim", global_dim, adios_err)
  call check_adios_err(adios_err,"Error writing "//trim(path) // "/global_dim")

  call adios_write(adios_handle, trim(path)// "/offset", offset, adios_err)
  call check_adios_err(adios_err,"Error writing "//trim(path) // "/offset")

  call adios_write(adios_handle, trim(path)// "/size", array_size, adios_err)
  call check_adios_err(adios_err,"Error writing "//trim(path) // "/size")

  ! to avoid compiler warning
  idummy = adios_group

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! note: not sure if we really need these additional values for storing an array

  ! local_dim
  ! gets associated variable for array
  call adios2_inquire_variable(v, adios_group, trim(path) // "/local_dim", adios_err)
  call check_adios_err(adios_err,"Error adios2 write dims offset: inquire variable "//trim(path) // "/local_dim"//" failed")
  if (.not. v%valid) stop 'Error adios2 variable invalid'

  call adios2_put(adios_handle, v, local_dim, adios_err)
  call check_adios_err(adios_err,"Error adios could not write parameter: "//trim(path) // "/local_dim")

  ! global_dim
  call adios2_inquire_variable(v, adios_group, trim(path) // "/global_dim", adios_err)
  call check_adios_err(adios_err,"Error adios2 write dims offset: inquire variable "//trim(path) // "/global_dim"//" failed")
  if (.not. v%valid) stop 'Error adios2 variable invalid'

  call adios2_put(adios_handle, v, global_dim, adios_err)
  call check_adios_err(adios_err,"Error adios could not write parameter: "//trim(path) // "/global_dim")

  ! offset
  call adios2_inquire_variable(v, adios_group, trim(path) // "/offset", adios_err)
  call check_adios_err(adios_err,"Error adios2 write dims offset: inquire variable "//trim(path) // "/offset"//" failed")
  if (.not. v%valid) stop 'Error adios2 variable invalid'

  call adios2_put(adios_handle, v, offset, adios_err)
  call check_adios_err(adios_err,"Error adios could not write parameter: "//trim(path) // "/offset")

  ! size
  call adios2_inquire_variable(v, adios_group, trim(path) // "/size", adios_err)
  call check_adios_err(adios_err,"Error adios2 write dims offset: inquire variable "//trim(path) // "/size"//" failed")
  if (.not. v%valid) stop 'Error adios2 variable invalid'

  call adios2_put(adios_handle, v, array_size, adios_err)
  call check_adios_err(adios_err,"Error adios could not write parameter: "//trim(path) // "/size")

  ! sync write
  call adios2_perform_puts(adios_handle, adios_err)
  call check_adios_err(adios_err,"Error adios2 perform puts in write_adios_perform() routine")

#endif

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
  subroutine write_adios_global_1d_real_1d_offset(adios_handle, adios_group, myrank, sizeprocs, &
                                                  local_dim, global_dim, offset, array_name, array)

  implicit none
  ! Parameters
  real, dimension(:), intent(in) :: array

#include "adios_helpers_writers_1d_generic_offset.inc"

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
  subroutine write_adios_global_1d_double_1d_offset(adios_handle, adios_group, myrank, sizeprocs, &
                                                    local_dim, global_dim, offset, array_name, array)

  implicit none
  ! Parameters
  real(kind=8), dimension(:), intent(in) :: array

#include "adios_helpers_writers_1d_generic_offset.inc"

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
  subroutine write_adios_global_1d_integer_1d_offset(adios_handle, adios_group, myrank, sizeprocs, &
                                                     local_dim, global_dim, offset, array_name, array)

  implicit none
  ! Parameters
  integer(kind=4), dimension(:), intent(in) :: array

#include "adios_helpers_writers_1d_generic_offset.inc"

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
  subroutine write_adios_global_1d_long_1d_offset(adios_handle, adios_group, myrank, sizeprocs, &
                                                  local_dim, global_dim, offset, array_name, array)

  implicit none
  ! Parameters
  integer(kind=8), dimension(:), intent(in) :: array

#include "adios_helpers_writers_1d_generic_offset.inc"

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
! unused so far...
!  subroutine write_adios_global_1d_logical_1d_offset(adios_handle, adios_group, myrank, sizeprocs, &
!                                                     local_dim, global_dim, offset, array_name, array)
!
!  implicit none
!  ! Parameters
!  logical, dimension(:), intent(in) :: array
!
!#include "adios_helpers_writers_1d_generic_offset.inc"
!
!  end subroutine write_adios_global_1d_logical_1d_offset
!
end module adios_helpers_writers_mod
