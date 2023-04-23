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
!  * reading in
!-------------------------------------------------------------------------------

#include "config.fh"


module adios_helpers_readers_mod

#if defined(USE_ADIOS2)
  use adios2
#endif

  use manager_adios, only: check_adios_err !,myrank_adios

  implicit none

  ! reading
  public :: read_adios_perform
  public :: read_adios_begin_step
  public :: read_adios_end_step

  public :: read_adios_array
  interface read_adios_array
    module procedure read_adios_array_gll
    module procedure read_adios_array_gll_int
    module procedure read_adios_array_1d
    module procedure read_adios_array_1d_int
  end interface read_adios_array

  public :: read_adios_array_gll_check

  public :: read_adios_scalar
  interface read_adios_scalar
    module procedure read_adios_scalar_int
    module procedure read_adios_scalar_long
    module procedure read_adios_scalar_real
    module procedure read_adios_scalar_double
  end interface read_adios_scalar

  public :: read_adios_scalar_local_dim

  public :: read_adios_schedule_array
  interface read_adios_schedule_array
    module procedure read_adios_schedule_array_global_1d_integer_1d
    module procedure read_adios_schedule_array_global_1d_integer_2d
    module procedure read_adios_schedule_array_global_1d_integer_3d
    module procedure read_adios_schedule_array_global_1d_integer_4d
    module procedure read_adios_schedule_array_global_1d_real_1d
    module procedure read_adios_schedule_array_global_1d_customreal_2d
    module procedure read_adios_schedule_array_global_1d_customreal_3d
    module procedure read_adios_schedule_array_global_1d_customreal_4d
    module procedure read_adios_schedule_array_global_1d_customreal_5d
    module procedure read_adios_schedule_array_global_1d_double_1d
    module procedure read_adios_schedule_array_global_1d_double_2d
    module procedure read_adios_schedule_array_global_1d_logical_1d
    module procedure read_adios_schedule_array_global_1d_string_1d
  end interface read_adios_schedule_array

  private

#if defined(USE_ADIOS2)
  ! default adios2 read mode:
  ! sync
  !integer, parameter :: myadios_reader_mode = adios2_mode_sync
  ! deferred (default)
  integer, parameter :: myadios_reader_mode = adios2_mode_deferred
#endif

contains

!---------------------------------------------------------------------------------
!
! additionals
!
!---------------------------------------------------------------------------------


  subroutine read_adios_perform(adios_handle)

  implicit none

#if defined(USE_ADIOS)
  ! adios
  integer(kind=8), intent(in) :: adios_handle
#elif defined(USE_ADIOS2)
  ! adios2
  type(adios2_engine), intent(in) :: adios_handle
#endif

  ! local parameters
  integer :: ier

  TRACE_ADIOS_L2('read_adios_perform')

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! performs read
  call adios_perform_reads(adios_handle, ier)
  if (ier /= 0 ) then
    print *,'Error adios: performing read failed'
    stop 'Error adios helper perform read failed'
  endif

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! performs read
  call adios2_perform_gets(adios_handle,ier)
  call check_adios_err(ier,"Error adios2 perform get failed")

#endif

  end subroutine read_adios_perform


!
!---------------------------------------------------------------------------------
!

  subroutine read_adios_begin_step(adios_handle)

  implicit none

#if defined(USE_ADIOS)
  integer(kind=8), intent(in) :: adios_handle
#elif defined(USE_ADIOS2)
  type(adios2_engine), intent(in) :: adios_handle
#endif
  ! local parameters
  integer :: ier

  TRACE_ADIOS_L2('read_adios_begin_step')

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! version 1 has no begin_step routine

  ! just check to avoid compiler warning
  if (adios_handle == 0) stop 'Error invalid adios handle in read_adios_begin_step'

  ! to avoid compiler warning
  ier = 0

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! moves to next step, starts at 0
  call adios2_begin_step(adios_handle, adios2_step_mode_read, ier)
  call check_adios_err(ier,"Error adios2 in read_adios_begin_step() routine")

#endif

  end subroutine read_adios_begin_step

!
!---------------------------------------------------------------------------------
!

  subroutine read_adios_end_step(adios_handle)

  implicit none

#if defined(USE_ADIOS)
  integer(kind=8), intent(in) :: adios_handle
#elif defined(USE_ADIOS2)
  type(adios2_engine), intent(in) :: adios_handle
#endif
  ! local parameters
  integer :: ier

  TRACE_ADIOS_L2('read_adios_end_step')

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! version 1 has no end_step routine

  ! just check
  if (adios_handle == 0) stop 'Invalid adios file handle in read_adios_end_step'

  ! to avoid compiler warning
  ier = 0

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! end step to indicate output is completed. ADIOS2 can do I/O
  call adios2_end_step(adios_handle, ier)
  call check_adios_err(ier,"Error adios2 in read_adios_end_step() routine")

#endif

  end subroutine read_adios_end_step

!---------------------------------------------------------------------------------
!
! synchronized scalar retrievals
!
!---------------------------------------------------------------------------------

  subroutine read_adios_scalar_int(adios_handle,adios_group,rank,scalar_name,scalar,step)

! reads in a single integer value

  implicit none

  integer, intent(out) :: scalar

#if defined(USE_ADIOS)
  integer(kind=8), intent(in) :: adios_handle
  integer(kind=8), intent(in) :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_engine), intent(in) :: adios_handle
  type(adios2_io), intent(in) :: adios_group
#endif
  integer, intent(in) :: rank
  character(len=*), intent(in) :: scalar_name
  integer(kind=8), intent(in), optional :: step

  ! local parameters
  integer :: ier
#if defined(USE_ADIOS)
  integer(kind=8) :: sel
  integer :: istep
#elif defined(USE_ADIOS2)
  type(adios2_variable) :: v
  integer(kind=8) :: step_start !,nsteps
  integer(kind=8) :: start(1),count(1)
  !debug
  logical, parameter :: DEBUG = .false.
  integer(kind=8) :: nsteps
#endif

  TRACE_ADIOS_L2_ARG('read_adios_scalar_int: ',trim(scalar_name))

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! selects data block
  call adios_selection_writeblock(sel, rank)

  ! reads array
  if (present(step)) then
    istep = int(step)
    call adios_schedule_read(adios_handle, sel, trim(scalar_name), istep, 1, scalar, ier)
  else
    call adios_schedule_read(adios_handle, sel, trim(scalar_name), 0, 1, scalar, ier)
  endif
  call check_adios_err(ier,"Error could not read parameter: "//trim(scalar_name))

  call adios_perform_reads(adios_handle, ier)
  call check_adios_err(ier,"Error could not perform read parameter: "//trim(scalar_name))

  ! frees selection
  call adios_selection_delete(sel)

  ! to avoid compiler warning
  ier = adios_group

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! gets associated variable for array
  call adios2_inquire_variable(v, adios_group, trim(scalar_name), ier)
  call check_adios_err(ier,"Error adios2 read_adios_scalar_int(): inquire variable "//trim(scalar_name)//" failed")

  ! checks
  if (.not. v%valid) stop 'Error adios2 variable invalid'

  if (v%type /= adios2_type_integer4) then
    print *,'Error: adios2 variable type mismatch: ',v%type,' instead of ',adios2_type_integer4
    call check_adios_err(1,"Error adios2 variable type mismatch for "//trim(scalar_name))
  endif

  !debug steps
  if (DEBUG) then
    ! for engine
    call adios2_steps(nsteps,adios_handle,ier)
    call check_adios_err(ier, "Error adios2 get engine steps for "//trim(scalar_name)//" failed")

    call adios2_current_step(step_start,adios_handle,ier)
    call check_adios_err(ier, "Error adios2 get current engine steps for "//trim(scalar_name)//" failed")

    print *,'debug adios: ',rank,trim(scalar_name),' engine nsteps = ',nsteps,' current = ',step_start; flush(6)

    ! for variable
    call adios2_variable_steps(nsteps, v, ier)
    call check_adios_err(ier, "Error adios2 get steps for "//trim(scalar_name)//" failed")

    print *,'debug adios: ',rank,trim(scalar_name),' variable nsteps = ',nsteps; flush(6)

    !only available in C:
    !> call adios2_variable_steps_start(step_start, v, ier)
    !> call check_adios_err(ier, "Error adios2 get start step for "//trim(scalar_name)//" failed")
    !print *,'debug adios: ',trim(scalar_name),' nsteps = ',nsteps,' ndims = ',v%ndims
    !only available in C:
    !> call adios2_selection_size(nsteps, v, ier)
    !> call check_adios_err(ier, "Error adios2 get start step for "//trim(scalar_name)//" failed")
    !print *,'debug adios: ',trim(scalar_name),' selection size = ',nsteps

    ! selection for local scalar variables
    ! this will fail for variables appended to the file (e.g., reg2/nspec and reg3/nspec variables in solver_data.bp etc.)
    ! maybe in future, adios2 will work with this.
    !
    !call adios2_set_block_selection(v, int(rank,kind=8), ier)
    !call check_adios_err(ier, "Error adios2 set block selection for "//trim(scalar_name)//" failed")

    call synchronize_all()
  endif

  ! selection for scalar as 1-D array single entry
  start(1) = 1 * int(rank,kind=8)
  count(1) = 1
  call adios2_set_selection(v, 1, start, count, ier)
  call check_adios_err(ier,"Error adios2 set selection for "//trim(scalar_name)//" failed")

  ! selection
  if (present(step)) then
    step_start = step
  else
    step_start = 0
  endif
  ! relative to first step available of variable
  call adios2_set_step_selection(v, step_start, int(1,kind=8), ier)
  call check_adios_err(ier, "Error adios2 set step variable for "//trim(scalar_name)//" failed")

  ! reads array data
  call adios2_get(adios_handle, v, scalar, adios2_mode_sync, ier)
  call check_adios_err(ier,"Error adios2 get for array "//trim(scalar_name)//" failed")

  ! debug
  if (DEBUG) then
    print *,'debug adios: ',rank,trim(scalar_name),' scalar value = ',scalar; flush(6)
    call synchronize_all()
  endif

#endif

  end subroutine read_adios_scalar_int

!
!---------------------------------------------------------------------------------
!

  subroutine read_adios_scalar_long(adios_handle,adios_group,rank,scalar_name,scalar,step)

! reads in a single integer value

  implicit none

  integer(kind=8), intent(out) :: scalar

#if defined(USE_ADIOS)
  integer(kind=8), intent(in) :: adios_handle
  integer(kind=8), intent(in) :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_engine), intent(in) :: adios_handle
  type(adios2_io), intent(in) :: adios_group
#endif
  integer, intent(in) :: rank
  character(len=*), intent(in) :: scalar_name
  integer(kind=8), intent(in), optional :: step

  ! local parameters
  integer :: ier
#if defined(USE_ADIOS)
  integer(kind=8) :: sel
  integer :: istep
#elif defined(USE_ADIOS2)
  type(adios2_variable) :: v
  integer(kind=8) :: step_start
  integer(kind=8) :: start(1),count(1)
#endif

  TRACE_ADIOS_L2_ARG('read_adios_scalar_long: ',trim(scalar_name))

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! selects data block
  call adios_selection_writeblock(sel, rank)

  ! reads array
  if (present(step)) then
    istep = int(step)
    call adios_schedule_read(adios_handle, sel, trim(scalar_name), istep, 1, scalar, ier)
  else
    call adios_schedule_read(adios_handle, sel, trim(scalar_name), 0, 1, scalar, ier)
  endif
  call check_adios_err(ier,"Error could not read parameter: "//trim(scalar_name))

  call adios_perform_reads(adios_handle, ier)
  call check_adios_err(ier,"Error could not perform read parameter: "//trim(scalar_name))

  ! frees selection
  call adios_selection_delete(sel)

  ! to avoid compiler warning
  ier = adios_group

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! gets associated variable for array
  call adios2_inquire_variable(v, adios_group, trim(scalar_name), ier)
  call check_adios_err(ier,"Error adios2 read_adios_scalar_long(): inquire variable "//trim(scalar_name)//" failed")

  ! checks
  if (.not. v%valid) stop 'Error adios2 variable invalid'

  if (v%type /= adios2_type_integer8) then
    print *,'Error: adios2 variable type mismatch: ',v%type,' instead of ',adios2_type_integer8
    call check_adios_err(1,"Error adios2 variable type mismatch for "//trim(scalar_name))
  endif

  ! selection for local scalar variables
  ! this will fail for variables appended to the file (e.g., reg2/nspec and reg3/nspec variables in solver_data.bp etc.)
  ! maybe in future, adios2 will work with this.
  !
  !call adios2_set_block_selection(v, int(rank,kind=8), ier)
  !call check_adios_err(ier, "Error adios2 set block selection for "//trim(scalar_name)//" failed")
  !
  ! selection for scalar as 1-D array single entry
  start(1) = 1 * int(rank,kind=8)
  count(1) = 1
  call adios2_set_selection(v, 1, start, count, ier)
  call check_adios_err(ier,"Error adios2 set selection for "//trim(scalar_name)//" failed")

  ! selection
  if (present(step)) then
    step_start = step
  else
    step_start = 0
  endif
  call adios2_set_step_selection(v, step_start, int(1,kind=8), ier)
  call check_adios_err(ier, "Error adios2 set step variable for "//trim(scalar_name)//" failed")

  ! reads array data
  call adios2_get(adios_handle, v, scalar, adios2_mode_sync, ier)
  call check_adios_err(ier,"Error adios2 get for array "//trim(scalar_name)//" failed")

#endif

  end subroutine read_adios_scalar_long

!
!---------------------------------------------------------------------------------
!

  subroutine read_adios_scalar_real(adios_handle,adios_group,rank,scalar_name,scalar,step)

! reads in a single precision real value

  implicit none

  real, intent(out) :: scalar

#if defined(USE_ADIOS)
  integer(kind=8), intent(in) :: adios_handle
  integer(kind=8), intent(in) :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_engine), intent(in) :: adios_handle
  type(adios2_io), intent(in) :: adios_group
#endif
  integer, intent(in) :: rank
  character(len=*), intent(in) :: scalar_name
  integer(kind=8), intent(in), optional :: step

  ! local parameters
  integer :: ier
#if defined(USE_ADIOS)
  integer(kind=8) :: sel
  integer :: istep
#elif defined(USE_ADIOS2)
  type(adios2_variable) :: v
  integer(kind=8) :: step_start
  integer(kind=8) :: start(1),count(1)
#endif

  TRACE_ADIOS_L2_ARG('read_adios_scalar_real: ',trim(scalar_name))

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! selects data block
  call adios_selection_writeblock(sel, rank)

  ! reads array
  if (present(step)) then
    istep = int(step)
    call adios_schedule_read(adios_handle, sel, trim(scalar_name), istep, 1, scalar, ier)
  else
    call adios_schedule_read(adios_handle, sel, trim(scalar_name), 0, 1, scalar, ier)
  endif
  call check_adios_err(ier,"Error could not read parameter: "//trim(scalar_name))

  call adios_perform_reads(adios_handle, ier)
  call check_adios_err(ier,"Error could not perform read parameter: "//trim(scalar_name))

  ! frees selection
  call adios_selection_delete(sel)

  ! to avoid compiler warning
  ier = adios_group

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! gets associated variable for array
  call adios2_inquire_variable(v, adios_group, trim(scalar_name), ier)
  call check_adios_err(ier,"Error adios2 read_adios_scalar_real(): inquire variable "//trim(scalar_name)//" failed")

  if (.not. v%valid) stop 'Error adios2 variable invalid'

  if (v%type /= adios2_type_real4) then
    print *,'Error: adios2 variable type mismatch: ',v%type,' instead of ',adios2_type_real4
    call check_adios_err(1,"Error adios2 variable type mismatch for "//trim(scalar_name))
  endif

  ! selection for local scalar variables
  ! this will fail for variables appended to the file (e.g., reg2/nspec and reg3/nspec variables in solver_data.bp etc.)
  ! maybe in future, adios2 will work with this.
  !
  !call adios2_set_block_selection(v, int(rank,kind=8), ier)
  !call check_adios_err(ier, "Error adios2 set block selection for "//trim(scalar_name)//" failed")
  !
  ! selection for scalar as 1-D array single entry
  start(1) = 1 * int(rank,kind=8)
  count(1) = 1
  call adios2_set_selection(v, 1, start, count, ier)
  call check_adios_err(ier,"Error adios2 set selection for "//trim(scalar_name)//" failed")

  ! selection
  if (present(step)) then
    step_start = step
  else
    step_start = 0
  endif
  call adios2_set_step_selection(v, step_start, int(1,kind=8), ier)
  call check_adios_err(ier, "Error adios2 set step variable for "//trim(scalar_name)//" failed")

  ! reads array data
  call adios2_get(adios_handle, v, scalar, adios2_mode_sync, ier)
  call check_adios_err(ier,"Error adios2 get for array "//trim(scalar_name)//" failed")

#endif

  end subroutine read_adios_scalar_real

!
!---------------------------------------------------------------------------------
!

  subroutine read_adios_scalar_double(adios_handle,adios_group,rank,scalar_name,scalar,step)

! reads in a single integer value

  implicit none

  double precision, intent(out) :: scalar

#if defined(USE_ADIOS)
  integer(kind=8), intent(in) :: adios_handle
  integer(kind=8), intent(in) :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_engine), intent(in) :: adios_handle
  type(adios2_io), intent(in) :: adios_group
#endif
  integer, intent(in) :: rank
  character(len=*), intent(in) :: scalar_name
  integer(kind=8), intent(in), optional :: step

  ! local parameters
  integer :: ier
#if defined(USE_ADIOS)
  integer(kind=8) :: sel
  integer :: istep
#elif defined(USE_ADIOS2)
  type(adios2_variable) :: v
  integer(kind=8) :: step_start
  integer(kind=8) :: start(1),count(1)
#endif

  TRACE_ADIOS_L2_ARG('read_adios_scalar_double: ',trim(scalar_name))

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! selects data block
  call adios_selection_writeblock(sel, rank)

  ! reads array
  if (present(step)) then
    istep = int(step)
    call adios_schedule_read(adios_handle, sel, trim(scalar_name), istep, 1, scalar, ier)
  else
    call adios_schedule_read(adios_handle, sel, trim(scalar_name), 0, 1, scalar, ier)
  endif
  call check_adios_err(ier,"Error could not read parameter: "//trim(scalar_name))

  call adios_perform_reads(adios_handle, ier)
  call check_adios_err(ier,"Error could not perform read parameter: "//trim(scalar_name))

  ! frees selection
  call adios_selection_delete(sel)

  ! to avoid compiler warning
  ier = adios_group

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! gets associated variable for array
  call adios2_inquire_variable(v, adios_group, trim(scalar_name), ier)
  call check_adios_err(ier,"Error adios2 read_adios_scalar_double(): inquire variable "//trim(scalar_name)//" failed")

  if (.not. v%valid) stop 'Error adios2 variable invalid'

  if (v%type /= adios2_type_real8) then
    print *,'Error: adios2 variable type mismatch: ',v%type,' instead of ',adios2_type_real8
    call check_adios_err(1,"Error adios2 variable type mismatch for "//trim(scalar_name))
  endif

  ! selection for local scalar variables
  ! this will fail for variables appended to the file (e.g., reg2/nspec and reg3/nspec variables in solver_data.bp etc.)
  ! maybe in future, adios2 will work with this.
  !
  !call adios2_set_block_selection(v, int(rank,kind=8), ier)
  !call check_adios_err(ier, "Error adios2 set block selection for "//trim(scalar_name)//" failed")
  !
  ! selection for scalar as 1-D array single entry
  start(1) = 1 * int(rank,kind=8)
  count(1) = 1
  call adios2_set_selection(v, 1, start, count, ier)
  call check_adios_err(ier,"Error adios2 set selection for "//trim(scalar_name)//" failed")

  ! selection
  if (present(step)) then
    step_start = step
  else
    step_start = 0
  endif
  call adios2_set_step_selection(v, step_start, int(1,kind=8), ier)
  call check_adios_err(ier, "Error adios2 set step variable for "//trim(scalar_name)//" failed")

  ! reads array data
  call adios2_get(adios_handle, v, scalar, adios2_mode_sync, ier)
  call check_adios_err(ier,"Error adios2 get for array "//trim(scalar_name)//" failed")

#endif

  end subroutine read_adios_scalar_double

!
!---------------------------------------------------------------------------------
!

  subroutine read_adios_scalar_local_dim(adios_handle, adios_group, rank, array_name, local_dim)

  implicit none

#if defined(USE_ADIOS)
  integer(kind=8), intent(in) :: adios_handle
  integer(kind=8), intent(in) :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_engine), intent(in) :: adios_handle
  type(adios2_io), intent(in) :: adios_group
#endif
  integer, intent(in) :: rank
  character(len=*), intent(in) :: array_name
  integer(kind=8), intent(out) :: local_dim

  ! local parameters
  integer :: ier
#if defined(USE_ADIOS2)
  integer(kind=8) :: start(1),count(1)
  type(adios2_variable) :: v
#endif
  character(len=256) :: full_name

  TRACE_ADIOS_L2_ARG('read_adios_scalar_local_dim: ',trim(array_name))

  ! checks
  if (len_trim(array_name) == 0) stop 'Error adios invalid array name in read_adios_scalar_local_dim()'

  ! initializes
  local_dim = 0

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! note: adios_get_scalar here retrieves the same local_dim for everyone (from writer rank 0)
  full_name = trim(array_name) // "/local_dim"
  call adios_get_scalar(adios_handle, trim(full_name),local_dim, ier)
  if (ier /= 0) then
    print *,'Error: reading adios local_dim for array: ',trim(full_name)
    print *,'Please check if ADIOS file contains this array: ',trim(full_name)
    call check_adios_err(ier,"Error adios get scalar "//trim(full_name)//" failed")
  endif

  ! to avoid compiler warning
  ier = adios_group
  ier = rank

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! gets dimension associated to array
  full_name = trim(array_name) // "/local_dim"

  ! note: might need to check in future if adding a .. // C_NULL_CHAR helps for avoiding problems with passing strings.
  !       since we call the adio2_** Fortran wrappers, this should deal with such issues.

  ! one could try to use the following:
  !
  ! reads in local dimension
  ! call read_adios_scalar(adios_handle,adios_group,rank,trim(full_name),local_dim)
  !
  ! or more explicitely, to get the local dimension:
  call adios2_inquire_variable(v, adios_group, trim(full_name), ier)
  call check_adios_err(ier,"Error adios2 read_adios_scalar_local_dim(): inquire variable "//trim(full_name)//" failed")

  if (.not. v%valid) stop 'Error adios2 variable invalid'

  ! selection for scalar as 1-D array single entry
  start(1) = 1 * int(rank,kind=8)
  count(1) = 1
  call adios2_set_selection(v, 1, start, count, ier)
  call check_adios_err(ier,"Error adios2 set selection for "//trim(full_name)//" failed")

  ! step selection
  call adios2_set_step_selection(v, int(0,kind=8), int(1,kind=8), ier)
  call check_adios_err(ier, "Error adios2 set step variable for "//trim(full_name)//" failed")

  ! note: adios2_get here retrieves the same local_dim for everyone (from writer rank 0)
  call adios2_get(adios_handle, v, local_dim, adios2_mode_sync, ier)
  call check_adios_err(ier,"Error adios2 get for array "//trim(full_name)//" failed")

#endif

  end subroutine read_adios_scalar_local_dim

!---------------------------------------------------------------------------------
!
! synchronized array retrieval functions
!
!---------------------------------------------------------------------------------

  subroutine read_adios_array_gll(adios_handle, adios_group, rank, nspec, array_name, array_gll)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

  implicit none

#if defined(USE_ADIOS)
  integer(kind=8), intent(in) :: adios_handle
  integer(kind=8), intent(in) :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_engine), intent(in) :: adios_handle
  type(adios2_io), intent(in) :: adios_group
#endif
  integer, intent(in) :: rank
  integer, intent(in) :: nspec

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec), intent(out) :: array_gll
  character(len=*), intent(in) :: array_name

  ! local parameters
  integer :: ier
  integer(kind=8) :: local_dim
  integer(kind=8) :: start(1),count(1)
#if defined(USE_ADIOS)
  integer(kind=8) :: sel
#elif defined(USE_ADIOS2)
  type(adios2_variable) :: v
#endif
  character(len=256) :: full_name

  TRACE_ADIOS_L2_ARG('read_adios_array_gll: ',trim(array_name))

  ! checks
  if (len_trim(array_name) == 0) stop 'Error adios invalid array name in read_adios_array_gll()'

  ! initializes
  array_gll(:,:,:,:) = 0.0_CUSTOM_REAL
  local_dim = 0

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! note: adios_get_scalar here retrieves the same local_dim for everyone (from writer rank 0)
  full_name = trim(array_name) // "/local_dim"
  call adios_get_scalar(adios_handle, trim(full_name),local_dim, ier)
  if (ier /= 0) then
    print *,'Error: reading adios GLL array: ',trim(full_name)
    print *,'Please check if ADIOS file contains this array: ',trim(full_name)
    call check_adios_err(ier,"Error adios get scalar "//trim(full_name)//" failed")
  endif

  ! allow any rank to read from another rank-segment
  !! checks rank
  !!if (myrank_adios /= rank) &
  !!  stop 'Error invalid rank for reading adios GLL array'

  start(1) = local_dim * int(rank,kind=8)
  count(1) = NGLLX * NGLLY * NGLLZ * nspec
  call adios_selection_boundingbox(sel, 1, start, count)

  ! reads selected array
  full_name = trim(array_name) // "/array"
  call adios_schedule_read(adios_handle, sel, trim(full_name), 0, 1, array_gll, ier)
  if (ier /= 0 ) then
    print *,'Error adios: scheduling read of array ',trim(full_name),' failed'
    stop 'Error adios helper schedule read array'
  endif

  call adios_perform_reads(adios_handle, ier)
  if (ier /= 0 ) then
    print *,'Error adios: performing read of array ',trim(full_name),' failed'
    stop 'Error adios helper perform read array'
  endif

  ! to avoid compiler warning
  ier = adios_group

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! gets dimension associated to array
  full_name = trim(array_name) // "/local_dim"

  ! note: might need to check in future if adding a .. // C_NULL_CHAR helps for avoiding problems with passing strings.
  !       since we call the adio2_** Fortran wrappers, this should deal with such issues.

  ! one could try to use the following:
  !
  ! reads in local dimension
  ! call read_adios_scalar(adios_handle,adios_group,rank,trim(full_name),local_dim)
  !
  ! or more explicitely, to get the local dimension:
  call adios2_inquire_variable(v, adios_group, trim(full_name), ier)
  call check_adios_err(ier,"Error adios2 read_adios_array_gll(): inquire variable "//trim(full_name)//" failed")

  if (.not. v%valid) stop 'Error adios2 variable invalid'

  ! selection for scalar as 1-D array single entry
  start(1) = 1 * int(rank,kind=8)
  count(1) = 1
  call adios2_set_selection(v, 1, start, count, ier)
  call check_adios_err(ier,"Error adios2 set selection for "//trim(full_name)//" failed")

  ! step selection
  call adios2_set_step_selection(v, int(0,kind=8), int(1,kind=8), ier)
  call check_adios_err(ier, "Error adios2 set step variable for "//trim(full_name)//" failed")

  ! note: adios2_get here retrieves the same local_dim for everyone (from writer rank 0)
  call adios2_get(adios_handle, v, local_dim, adios2_mode_sync, ier)
  call check_adios_err(ier,"Error adios2 get for array "//trim(full_name)//" failed")

  ! allow any rank to read from another rank-segment
  !! checks rank
  !!if (myrank_adios /= rank) &
  !!  stop 'Error invalid rank for reading adios GLL array'

  ! array data
  ! gets associated variable for array
  full_name = trim(array_name) // "/array"
  call adios2_inquire_variable(v, adios_group, trim(full_name), ier)
  call check_adios_err(ier,"Error adios2 read_adios_array_gll(): inquire variable "//trim(full_name)//" failed")

  if (.not. v%valid) stop 'Error adios2 variable invalid'

  ! selection
  start(1) = local_dim * int(rank,kind=8)
  count(1) = NGLLX * NGLLY * NGLLZ * nspec
  call adios2_set_selection(v, 1, start, count, ier)
  call check_adios_err(ier,"Error adios2 set selection for "//trim(full_name)//" failed")

  call adios2_set_step_selection(v, int(0,kind=8), int(1,kind=8), ier)
  call check_adios_err(ier, "Error adios2 set step variable for "//trim(full_name)//" failed")

  ! reads array data
  call adios2_get(adios_handle, v, array_gll, adios2_mode_sync, ier)
  call check_adios_err(ier,"Error adios2 get for array "//trim(full_name)//" failed")

#endif

  end subroutine read_adios_array_gll

!
!---------------------------------------------------------------------------------
!

  subroutine read_adios_array_gll_check(adios_handle,adios_group,rank,nspec,array_name,array_gll,iexist)

! checks if array_name is found in adios file

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

  implicit none

#if defined(USE_ADIOS)
  integer(kind=8), intent(in) :: adios_handle
  integer(kind=8), intent(in) :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_engine), intent(in) :: adios_handle
  type(adios2_io), intent(in) :: adios_group
#endif
  integer, intent(in) :: rank
  integer, intent(in) :: nspec

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec), intent(out) :: array_gll
  character(len=*), intent(in) :: array_name

  integer, intent(out) :: iexist

  ! local parameters
  integer(kind=8) :: local_dim
  integer(kind=8) :: start(1),count(1)
  integer :: ier
#if defined(USE_ADIOS)
  integer(kind=8) :: sel
  ! inquiry
  integer :: i,variable_count, attribute_count
  integer :: timestep_first, timestep_last
  character (len=128), dimension(:), allocatable :: fnamelist
  logical :: found_par
#elif defined(USE_ADIOS2)
  type(adios2_variable) :: v
#endif
  character(len=256) :: full_name

  TRACE_ADIOS_L2_ARG('read_adios_array_gll_check: ',trim(array_name))

  ! checks
  if (len_trim(array_name) == 0) stop 'Error adios invalid array name in read_adios_array_gll_check()'

  ! initializes
  array_gll(:,:,:,:) = 0.0_CUSTOM_REAL
  local_dim = 0
  iexist = 0  ! 0 == does not exist; 1 == success, exists in file

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! file inquiry
  call adios_inq_file(adios_handle,variable_count,attribute_count,timestep_first,timestep_last,ier)
  if (ier /= 0) stop 'Error inquiring adios file for reading'

  ! variable names
  found_par = .false.
  if (variable_count > 0) then
    allocate (fnamelist(variable_count),stat=ier)
    if (ier /= 0) stop 'Error allocating namelist array'

    ! gets variable names
    call adios_inq_varnames(adios_handle, fnamelist, ier)
    if (ier /= 0) stop 'Error inquiring variable names'

    ! user output
    !print *,'variables: ',variable_count

    ! checks if a variable name matches the array_name
    do i = 1,variable_count
      !debug
      !print *,'  ',trim(fnamelist(i)),' compare to ',trim(array_name)//"/local_dim"

      ! compares with input name
      ! (adds /local_dims to name to have full name comparison, e.g., reg1/rho/local_dim)
      if (trim(fnamelist(i)) == trim(array_name)//"/local_dim") then
        ! debug
        !if (rank == 0) print *,'  found ',trim(array_name)

        found_par = .true.
        ! exit do-loop
        exit
      endif
    enddo
    deallocate(fnamelist)
  else
    print *,'ADIOS file contains no variables'
    return
  endif

  ! check if parameter found
  if (.not. found_par) then
    iexist = 0 ! returns zero if not found
    return
  endif

  ! gets dimension
  ! note: adios_get_scalar here retrieves the same local_dim for everyone (from writer rank 0)
  full_name = trim(array_name) // "/local_dim"
  call adios_get_scalar(adios_handle, trim(full_name),local_dim, ier)
  if (ier == 0) then
    ! allow any rank to read from another rank-segment
    !! checks rank
    !!if (myrank_adios /= rank) &
    !!  stop 'Error invalid rank for reading adios GLL array'

    start(1) = local_dim * int(rank,kind=8)
    count(1) = NGLLX * NGLLY * NGLLZ * nspec
    call adios_selection_boundingbox(sel, 1, start, count)

    ! reads selected array
    full_name = trim(array_name) // "/array"
    call adios_schedule_read(adios_handle, sel, trim(full_name),0, 1, array_gll, ier)
    if (ier /= 0 ) then
      print *,'Error adios: scheduling read of array ',trim(full_name),' failed'
      stop 'Error adios helper schedule read array'
    endif

    call adios_perform_reads(adios_handle, ier)
    if (ier /= 0 ) then
      print *,'Error adios: performing read of array ',trim(full_name),' failed'
      stop 'Error adios helper perform read array'
    endif

    ! found and read
    iexist = 1
  endif

  ! to avoid compiler warning
  ier = adios_group

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! gets dimension associated to array
  full_name = trim(array_name) // "/local_dim"
  call adios2_inquire_variable(v, adios_group, trim(full_name), ier)

  !debug
  !print *,'debug: ADIOS2 check variable ',rank,trim(full_name),ier

  if (ier == 0) then
    ! checks variable flag
    if (.not. v%valid) stop 'Error adios2 variable invalid'

    ! selection for scalar as 1-D array single entry
    start(1) = 1 * int(rank,kind=8)
    count(1) = 1
    call adios2_set_selection(v, 1, start, count, ier)
    call check_adios_err(ier,"Error adios2 set selection for "//trim(full_name)//" failed")

    ! step selection
    call adios2_set_step_selection(v, int(0,kind=8), int(1,kind=8), ier)
    call check_adios_err(ier, "Error adios2 set step variable for "//trim(full_name)//" failed")

    ! note: adios2_get here retrieves the same local_dim for everyone (from writer rank 0)
    call adios2_get(adios_handle, v, local_dim, adios2_mode_sync, ier)
    call check_adios_err(ier,"Error adios2 get for array "//trim(full_name)//" failed")

    ! allow any rank to read from another rank-segment
    !! checks rank
    !!if (myrank_adios /= rank) &
    !!  stop 'Error invalid rank for reading adios GLL array'

    ! array data
    ! gets associated variable for array
    full_name = trim(array_name) // "/array"
    call adios2_inquire_variable(v, adios_group, trim(full_name), ier)
    call check_adios_err(ier,"Error adios2 read_adios_array_gll_check(): inquire variable "//trim(full_name)//" failed")

    ! checks variable flag
    if (.not. v%valid) stop 'Error adios2 variable invalid'

    ! selection
    start(1) = local_dim * int(rank,kind=8)
    count(1) = NGLLX * NGLLY * NGLLZ * nspec
    call adios2_set_selection(v, 1, start, count, ier)
    call check_adios_err(ier,"Error adios2 set selection for "//trim(full_name)//" failed")

    call adios2_set_step_selection(v, int(0,kind=8), int(1,kind=8), ier)
    call check_adios_err(ier, "Error adios2 set step variable for "//trim(full_name)//" failed")

    ! reads array data
    call adios2_get(adios_handle, v, array_gll, adios2_mode_sync, ier)
    call check_adios_err(ier,"Error adios2 get for array "//trim(full_name)//" failed")

    ! found and read
    iexist = 1
  endif

#endif

  end subroutine read_adios_array_gll_check

!
!---------------------------------------------------------------------------------
!

  subroutine read_adios_array_gll_int(adios_handle,adios_group,rank,nspec,array_name,array_gll)

  use constants, only: NGLLX,NGLLY,NGLLZ

  implicit none

#if defined(USE_ADIOS)
  integer(kind=8), intent(in) :: adios_handle
  integer(kind=8), intent(in) :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_engine), intent(in) :: adios_handle
  type(adios2_io), intent(in) :: adios_group
#endif
  integer, intent(in) :: rank
  integer, intent(in) :: nspec

  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec), intent(out) :: array_gll
  character(len=*), intent(in) :: array_name

  ! local parameters
  integer :: ier
  integer(kind=8) :: local_dim
  integer(kind=8) :: start(1),count(1)
#if defined(USE_ADIOS)
  integer(kind=8) :: sel
#elif defined(USE_ADIOS2)
  type(adios2_variable) :: v
#endif
  character(len=256) :: full_name

  TRACE_ADIOS_L2_ARG('read_adios_array_gll_int: ',trim(array_name))

  ! checks
  if (len_trim(array_name) == 0) stop 'Error adios invalid array name in read_adios_array_gll_int()'

  ! initializes
  array_gll(:,:,:,:) = 0
  local_dim = 0

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! note: adios_get_scalar here retrieves the same local_dim for everyone (from writer rank 0)
  full_name = trim(array_name) // "/local_dim"
  call adios_get_scalar(adios_handle, trim(full_name),local_dim, ier)
  if (ier /= 0 ) then
    print *,'Error adios: reading array ',trim(full_name),' failed'
    stop 'Error adios helper read array'
  endif

  ! allow any rank to read from another rank-segment
  !! checks rank
  !!if (myrank_adios /= rank) &
  !!  stop 'Error invalid rank for reading adios GLL array'

  start(1) = local_dim * int(rank,kind=8)
  count(1) = NGLLX * NGLLY * NGLLZ * nspec

  call adios_selection_boundingbox(sel, 1, start, count)

  ! reads selected array
  full_name = trim(array_name) // "/array"
  call adios_schedule_read(adios_handle, sel, trim(full_name), 0, 1, array_gll, ier)
  if (ier /= 0 ) then
    print *,'Error adios: scheduling read of array ',trim(full_name),' failed'
    stop 'Error adios helper schedule read array'
  endif

  call adios_perform_reads(adios_handle, ier)
  if (ier /= 0 ) then
    print *,'Error adios: performing read of array ',trim(full_name),' failed'
    stop 'Error adios helper perform read array'
  endif

  ! to avoid compiler warning
  ier = adios_group

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! gets dimension associated to array
  full_name = trim(array_name) // "/local_dim"
  call adios2_inquire_variable(v, adios_group, trim(full_name), ier)
  call check_adios_err(ier,"Error adios2 read_adios_array_gll_int(): inquire variable "//trim(full_name)//" failed")

  if (.not. v%valid) stop 'Error adios2 variable invalid'

  ! selection for scalar as 1-D array single entry
  start(1) = 1 * int(rank,kind=8)
  count(1) = 1
  call adios2_set_selection(v, 1, start, count, ier)
  call check_adios_err(ier,"Error adios2 set selection for "//trim(full_name)//" failed")

  ! step selection
  call adios2_set_step_selection(v, int(0,kind=8), int(1,kind=8), ier)
  call check_adios_err(ier, "Error adios2 set step variable for "//trim(full_name)//" failed")

  ! note: adios2_get here retrieves the same local_dim for everyone (from writer rank 0)
  call adios2_get(adios_handle, v, local_dim, adios2_mode_sync, ier)
  call check_adios_err(ier,"Error adios2 get for array "//trim(full_name)//" failed")

  ! allow any rank to read from another rank-segment
  !! checks rank
  !!if (myrank_adios /= rank) &
  !!  stop 'Error invalid rank for reading adios GLL array'

  ! array data
  ! gets associated variable for array
  full_name = trim(array_name) // "/array"
  call adios2_inquire_variable(v, adios_group, trim(full_name), ier)
  call check_adios_err(ier,"Error adios2 read_adios_array_gll_int(): inquire variable "//trim(full_name)//" failed")

  if (.not. v%valid) stop 'Error adios2 variable invalid'

  ! selection
  start(1) = local_dim * int(rank,kind=8)
  count(1) = NGLLX * NGLLY * NGLLZ * nspec
  call adios2_set_selection(v, 1, start, count, ier)
  call check_adios_err(ier,"Error adios2 set selection for "//trim(full_name)//" failed")

  call adios2_set_step_selection(v, int(0,kind=8), int(1,kind=8), ier)
  call check_adios_err(ier, "Error adios2 set step variable for "//trim(full_name)//" failed")

  ! reads array data
  call adios2_get(adios_handle, v, array_gll, adios2_mode_sync, ier)
  call check_adios_err(ier,"Error adios2 get for array "//trim(full_name)//" failed")

#endif

  end subroutine read_adios_array_gll_int

!
!---------------------------------------------------------------------------------
!

  subroutine read_adios_array_1d(adios_handle,adios_group,rank,nsize,array_name,array_1d)

  use constants, only: CUSTOM_REAL

  implicit none

#if defined(USE_ADIOS)
  integer(kind=8), intent(in) :: adios_handle
  integer(kind=8), intent(in) :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_engine), intent(in) :: adios_handle
  type(adios2_io), intent(in) :: adios_group
#endif
  integer, intent(in) :: rank
  integer, intent(in) :: nsize

  real(kind=CUSTOM_REAL),dimension(nsize), intent(out) :: array_1d
  character(len=*), intent(in) :: array_name

  ! local parameters
  integer :: ier
  integer(kind=8) :: local_dim
  integer(kind=8) :: start(1),count(1)
#if defined(USE_ADIOS)
  integer(kind=8) :: sel
  !integer(kind=8) :: offset
#elif defined(USE_ADIOS2)
  type(adios2_variable) :: v
#endif
  character(len=256) :: full_name

  TRACE_ADIOS_L2_ARG('read_adios_array_1d: ',trim(array_name))

  ! checks
  if (len_trim(array_name) == 0) stop 'Error adios invalid array name in read_adios_array_1d()'

  ! allow any rank to read from another rank-segment
  !! checks rank
  !!if (myrank_adios /= rank) &
  !!  stop 'Error invalid rank for reading adios GLL array'

  ! initializes
  array_1d(:) = 0._CUSTOM_REAL
  local_dim = 0

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! gets local_dim metadata
  full_name = trim(array_name) // "/local_dim"
  call adios_get_scalar(adios_handle, trim(full_name),local_dim, ier)
  if (ier /= 0 ) then
    print *,'Error adios: reading array ',trim(full_name),' failed'
    stop 'Error adios helper read array'
  endif

  start(1) = local_dim * int(rank,kind=8)
  count(1) = int(nsize,kind=8)

  ! gets offset (offset is the same as: local_dim * rank)
  !call adios_selection_writeblock(sel,rank)
  !call adios_schedule_read(adios_handle, sel, trim(array_name)//"/offset", 0, 1, offset, ier)
  !if (ier /= 0 ) stop 'Error adios: reading offset'
  !call adios_perform_reads(adios_handle, ier)
  !if (ier /= 0 ) stop 'Error adios: perform reading mesh file offsets failed'
  !start(1) = offset
  !count(1) = nsize

  call adios_selection_boundingbox(sel, 1, start, count)

  ! reads selected array
  full_name = trim(array_name) // "/array"
  call adios_schedule_read(adios_handle, sel, trim(full_name), 0, 1, array_1d, ier)
  if (ier /= 0 ) then
    print *,'Error adios: scheduling read of array ',trim(full_name),' failed'
    stop 'Error adios helper schedule read array'
  endif

  call adios_perform_reads(adios_handle, ier)
  if (ier /= 0 ) then
    print *,'Error adios: performing read of array ',trim(full_name),' failed'
    stop 'Error adios helper perform read array'
  endif

  ! to avoid compiler warning
  ier = adios_group

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! gets dimension associated to array
  full_name = trim(array_name) // "/local_dim"
  call adios2_inquire_variable(v, adios_group, trim(full_name), ier)
  call check_adios_err(ier,"Error adios2 read_adios_array_1d(): inquire variable "//trim(full_name)//" failed")

  if (.not. v%valid) stop 'Error adios2 variable invalid'

  ! selection for scalar as 1-D array single entry
  start(1) = 1 * int(rank,kind=8)
  count(1) = 1
  call adios2_set_selection(v, 1, start, count, ier)
  call check_adios_err(ier,"Error adios2 set selection for "//trim(full_name)//" failed")

  ! step selection
  call adios2_set_step_selection(v, int(0,kind=8), int(1,kind=8), ier)
  call check_adios_err(ier, "Error adios2 set step variable for "//trim(full_name)//" failed")

  ! note: adios2_get here retrieves the same local_dim for everyone (from writer rank 0)
  call adios2_get(adios_handle, v, local_dim, adios2_mode_sync, ier)
  call check_adios_err(ier,"Error adios2 get for array "//trim(full_name)//" failed")

  ! allow any rank to read from another rank-segment
  !! checks rank
  !!if (myrank_adios /= rank) &
  !!  stop 'Error invalid rank for reading adios GLL array'

  ! array data
  ! gets associated variable for array
  full_name = trim(array_name) // "/array"
  call adios2_inquire_variable(v, adios_group, trim(full_name), ier)
  call check_adios_err(ier,"Error adios2 read_adios_array_1d(): inquire variable "//trim(full_name)//" failed")

  if (.not. v%valid) stop 'Error adios2 variable invalid'

  ! selection
  start(1) = local_dim * int(rank,kind=8)
  count(1) = int(nsize,kind=8)
  call adios2_set_selection(v, 1, start, count, ier)
  call check_adios_err(ier,"Error adios2 set selection for "//trim(full_name)//" failed")

  call adios2_set_step_selection(v, int(0,kind=8), int(1,kind=8), ier)
  call check_adios_err(ier, "Error adios2 set step variable for "//trim(full_name)//" failed")

  ! reads array data
  call adios2_get(adios_handle, v, array_1d, adios2_mode_sync, ier)
  call check_adios_err(ier,"Error adios2 get for array "//trim(full_name)//" failed")

#endif

  end subroutine read_adios_array_1d

!
!---------------------------------------------------------------------------------
!

  subroutine read_adios_array_1d_int(adios_handle,adios_group,rank,nsize,array_name,array_1d)

  implicit none

#if defined(USE_ADIOS)
  integer(kind=8), intent(in) :: adios_handle
  integer(kind=8), intent(in) :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_engine), intent(in) :: adios_handle
  type(adios2_io), intent(in) :: adios_group
#endif
  integer, intent(in) :: rank
  integer, intent(in) :: nsize

  integer, dimension(nsize), intent(out) :: array_1d
  character(len=*), intent(in) :: array_name

  ! local parameters
  integer :: ier
  integer(kind=8) :: local_dim
  integer(kind=8) :: start(1),count(1)
#if defined(USE_ADIOS)
  integer(kind=8) :: sel
  !integer(kind=8) :: offset
#elif defined(USE_ADIOS2)
  type(adios2_variable) :: v
#endif
  character(len=256) :: full_name

  TRACE_ADIOS_L2_ARG('read_adios_array_1d_int: ',trim(array_name))

  ! checks
  if (len_trim(array_name) == 0) stop 'Error adios invalid array name in read_adios_array_1d_int()'

  ! initializes
  array_1d(:) = 0
  local_dim = 0

#if defined(USE_ADIOS)
  ! ADIOS 1
  full_name = trim(array_name) // "/local_dim"
  call adios_get_scalar(adios_handle, trim(full_name),local_dim, ier)
  if (ier /= 0 ) then
    print *,'Error adios: reading array ',trim(full_name),' failed'
    stop 'Error adios helper read array'
  endif

  start(1) = local_dim * int(rank,kind=8)
  count(1) = int(nsize,kind=8)

  call adios_selection_boundingbox(sel, 1, start, count)

  ! reads selected array
  full_name = trim(array_name) // "/array"
  call adios_schedule_read(adios_handle, sel, trim(full_name), 0, 1, array_1d, ier)
  if (ier /= 0 ) then
    print *,'Error adios: scheduling read of array ',trim(full_name),' failed'
    stop 'Error adios helper schedule read array'
  endif

  call adios_perform_reads(adios_handle, ier)
  if (ier /= 0 ) then
    print *,'Error adios: performing read of array ',trim(full_name),' failed'
    stop 'Error adios helper perform read array'
  endif

  ! to avoid compiler warning
  ier = adios_group

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! gets dimension associated to array
  full_name = trim(array_name) // "/local_dim"
  call adios2_inquire_variable(v, adios_group, trim(full_name), ier)
  call check_adios_err(ier,"Error adios2 read_adios_array_1d_int(): inquire variable "//trim(full_name)//" failed")

  if (.not. v%valid) stop 'Error adios2 variable invalid'

  ! selection for scalar as 1-D array single entry
  start(1) = 1 * int(rank,kind=8)
  count(1) = 1
  call adios2_set_selection(v, 1, start, count, ier)
  call check_adios_err(ier,"Error adios2 set selection for "//trim(full_name)//" failed")

  ! step selection
  call adios2_set_step_selection(v, int(0,kind=8), int(1,kind=8), ier)
  call check_adios_err(ier, "Error adios2 set step variable for "//trim(full_name)//" failed")

  ! note: adios2_get here retrieves the same local_dim for everyone (from writer rank 0)
  call adios2_get(adios_handle, v, local_dim, adios2_mode_sync, ier)
  call check_adios_err(ier,"Error adios2 get for array "//trim(full_name)//" failed")

  ! array data
  ! gets associated variable for array
  full_name = trim(array_name) // "/array"
  call adios2_inquire_variable(v, adios_group, trim(full_name), ier)
  call check_adios_err(ier,"Error adios2 read_adios_array_1d_int(): inquire variable "//trim(full_name)//" failed")

  if (.not. v%valid) stop 'Error adios2 variable invalid'

  ! selection
  start(1) = local_dim * int(rank,kind=8)
  count(1) = int(nsize,kind=8)
  call adios2_set_selection(v, 1, start, count, ier)
  call check_adios_err(ier,"Error adios2 set selection for "//trim(full_name)//" failed")

  call adios2_set_step_selection(v, int(0,kind=8), int(1,kind=8), ier)
  call check_adios_err(ier, "Error adios2 set step variable for "//trim(full_name)//" failed")

  ! reads array data
  call adios2_get(adios_handle, v, array_1d, adios2_mode_sync, ier)
  call check_adios_err(ier,"Error adios2 get for array "//trim(full_name)//" failed")

#endif

  end subroutine read_adios_array_1d_int


!---------------------------------------------------------------------------------
!
! wrappers for schedule reads
!
!---------------------------------------------------------------------------------

  subroutine read_adios_schedule_array_global_1d_integer_1d(adios_handle,adios_group,sel,start,count,array_name,array,step)

  implicit none

  integer, dimension(:), intent(out) :: array

#include "adios_helpers_readers_1d_generic.inc"

  end subroutine read_adios_schedule_array_global_1d_integer_1d

!
!---------------------------------------------------------------------------------
!

  subroutine read_adios_schedule_array_global_1d_integer_2d(adios_handle,adios_group,sel,start,count,array_name,array,step)

  implicit none

  integer, dimension(:,:), intent(out) :: array

#include "adios_helpers_readers_1d_generic.inc"

  end subroutine read_adios_schedule_array_global_1d_integer_2d

!
!---------------------------------------------------------------------------------
!

  subroutine read_adios_schedule_array_global_1d_integer_3d(adios_handle,adios_group,sel,start,count,array_name,array,step)

  implicit none

  integer, dimension(:,:,:), intent(out) :: array

#include "adios_helpers_readers_1d_generic.inc"

  end subroutine read_adios_schedule_array_global_1d_integer_3d

!
!---------------------------------------------------------------------------------
!

  subroutine read_adios_schedule_array_global_1d_integer_4d(adios_handle,adios_group,sel,start,count,array_name,array,step)

  implicit none

  integer, dimension(:,:,:,:), intent(out) :: array

#include "adios_helpers_readers_1d_generic.inc"

  end subroutine read_adios_schedule_array_global_1d_integer_4d


!
!---------------------------------------------------------------------------------
!

  subroutine read_adios_schedule_array_global_1d_real_1d(adios_handle,adios_group,sel,start,count,array_name,array,step)

  implicit none

  real, dimension(:), intent(out) :: array

#include "adios_helpers_readers_1d_generic.inc"

  end subroutine read_adios_schedule_array_global_1d_real_1d

!
!---------------------------------------------------------------------------------
!

  subroutine read_adios_schedule_array_global_1d_customreal_2d(adios_handle,adios_group,sel,start,count,array_name,array,step)

  use constants, only: CUSTOM_REAL

  implicit none

  real(kind=CUSTOM_REAL), dimension(:,:), intent(out) :: array

#include "adios_helpers_readers_1d_generic.inc"

  end subroutine read_adios_schedule_array_global_1d_customreal_2d

!
!---------------------------------------------------------------------------------
!

  subroutine read_adios_schedule_array_global_1d_customreal_3d(adios_handle,adios_group,sel,start,count,array_name,array,step)

  use constants, only: CUSTOM_REAL

  implicit none

  real(kind=CUSTOM_REAL), dimension(:,:,:), intent(out) :: array

#include "adios_helpers_readers_1d_generic.inc"

  end subroutine read_adios_schedule_array_global_1d_customreal_3d

!
!---------------------------------------------------------------------------------
!

  subroutine read_adios_schedule_array_global_1d_customreal_4d(adios_handle,adios_group,sel,start,count,array_name,array,step)

  use constants, only: CUSTOM_REAL

  implicit none

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), intent(out) :: array

#include "adios_helpers_readers_1d_generic.inc"

  end subroutine read_adios_schedule_array_global_1d_customreal_4d

!
!---------------------------------------------------------------------------------
!

  subroutine read_adios_schedule_array_global_1d_customreal_5d(adios_handle,adios_group,sel,start,count,array_name,array,step)

  use constants, only: CUSTOM_REAL

  implicit none

  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), intent(out) :: array

#include "adios_helpers_readers_1d_generic.inc"

  end subroutine read_adios_schedule_array_global_1d_customreal_5d

!
!---------------------------------------------------------------------------------
!

  subroutine read_adios_schedule_array_global_1d_double_1d(adios_handle,adios_group,sel,start,count,array_name,array,step)

  implicit none

  double precision, dimension(:), intent(out) :: array

#include "adios_helpers_readers_1d_generic.inc"

  end subroutine read_adios_schedule_array_global_1d_double_1d

!
!---------------------------------------------------------------------------------
!

  subroutine read_adios_schedule_array_global_1d_double_2d(adios_handle,adios_group,sel,start,count,array_name,array,step)

  implicit none

  double precision, dimension(:,:), intent(out) :: array

#include "adios_helpers_readers_1d_generic.inc"

  end subroutine read_adios_schedule_array_global_1d_double_2d

!
!---------------------------------------------------------------------------------
!

! note: adios1 seems to handle the schedule_read routine a bit different than adios2 which checks the array type.
!       both don't support logical, but store it as integer*4 size.
!
! adios1: would assign store the logical array implicitly as integer*4
! thus we could just call:
!
!  subroutine read_adios_schedule_array_global_1d_logical_1d(adios_handle,adios_group,sel,start,count,array_name,array,step)
!
!  implicit none
!
!  logical, dimension(:), intent(out) :: array
!
!#include "adios_helpers_readers_1d_generic.inc"
!
!  end subroutine read_adios_schedule_array_global_1d_logical_1d
!
!
! adios2 has no adios_get() routine for logicals.
! we use an explicit integer array instead for storing/reading:

  subroutine read_adios_schedule_array_global_1d_logical_1d(adios_handle,adios_group,sel,start,count,array_name,array_l,step)

  implicit none

  logical, dimension(:), intent(out) :: array_l
  integer, dimension(size(array_l)) :: array

#include "adios_helpers_readers_1d_generic.inc"

  ! performs actual reads before temporary array goes out of scope
  call read_adios_perform(adios_handle)

  ! explicitly converts from integer to logical array
  array_l(:) = .false.
  where(array(:) == 1) array_l(:) = .true.

  end subroutine read_adios_schedule_array_global_1d_logical_1d

!
!---------------------------------------------------------------------------------
!

  subroutine read_adios_schedule_array_global_1d_string_1d(adios_handle,adios_group,sel,start,count,array_name,array,step)

  implicit none

  character(len=*), intent(out) :: array

#include "adios_helpers_readers_1d_generic.inc"

  end subroutine read_adios_schedule_array_global_1d_string_1d

!===============================================================================

end module adios_helpers_readers_mod
