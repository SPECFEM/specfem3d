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
!> Tools for ADIOS/ADIOS2,
!> contains wrapper subroutines for common adios calls
!
! note: adios library calls use a format like "adios_do_something***()"
!       our own wrapper functions thus will rather use something like "do_something_adios***()"
!       to better distinguish between library functions and wrappers.
!-------------------------------------------------------------------------------

#include "config.fh"


module manager_adios

#if defined(USE_ADIOS)
  use adios_write_mod
  use adios_read_mod
#elif defined(USE_ADIOS2)
  use adios2
#endif

  implicit none

  private

  ! MPI copies of communicator and rank
  integer :: comm_adios
  integer, public :: myrank_adios
  integer, public :: sizeprocs_adios

  ! initialized flag
  logical :: is_adios_initialized
  logical, public :: is_adios_version1
  logical, public :: is_adios_version2

  ! compression flag
  logical, public :: use_adios_compression

#if defined(USE_ADIOS)
  ! adios
  character(len=*),parameter :: ADIOS_VERBOSITY = "verbose=1" ! lowest level: verbose=1, .., debug=4
  ! default file handle for read/write
  integer(kind=8), public :: myadios_file
  ! IO group
  integer(kind=8), public :: myadios_group

  ! additional file handle for read/write value file
  integer(kind=8), public :: myadios_val_file
  integer(kind=8), public :: myadios_val_group

  ! for undo att
  integer(kind=8), public :: myadios_fwd_group
  integer(kind=8), public :: myadios_fwd_file
  logical, public :: is_initialized_fwd_group

#elif defined(USE_ADIOS2)
  ! note: we're using save attribute to be able to compile with a flag like -std=f2003
  !       without it, a compilation error with gfortran (v7.5.0) would occur:
  !       ..
  !         Error: Fortran 2008: Implied SAVE for module variable 'myadios2_obj' at (1), needed due to the default initialization
  !       ..
  !       this seems to be needed only for type(..) variables.
  !
  ! adios2 main object
  type(adios2_adios), public, save :: myadios2_obj
  ! default file handle for read/write
  type(adios2_engine), public, save :: myadios_file
  ! IO group
  type(adios2_io), public, save :: myadios_group

  ! additional file handle for read/write value file
  type(adios2_engine), public, save :: myadios_val_file
  type(adios2_io), public, save :: myadios_val_group

  ! for undo_att
  type(adios2_io), public, save :: myadios_fwd_group
  type(adios2_engine), public, save :: myadios_fwd_file
  logical, public :: is_initialized_fwd_group

  ! debugging mode
  ! note: adios2 still in development stage, let's keep debug mode on
  logical,parameter :: USE_ADIOS2_DEBUG_MODE = .true.

#endif

  ! public accessibility
  public :: initialize_adios
  public :: finalize_adios

! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools
#if defined(USE_ADIOS) || defined(USE_ADIOS2)
  ! file opening
  public :: open_file_adios_read_and_init_method
  public :: open_file_adios_read_only_rank
  public :: open_file_adios_read
  public :: open_file_adios_write
  public :: open_file_adios_write_append

  ! file closing
  public :: close_file_adios
  public :: close_file_adios_read_and_finalize_method
  public :: close_file_adios_read_and_finalize_method_only_rank
  public :: close_file_adios_read

  ! groups
  public :: init_adios_group
  public :: init_adios_group_undo_att
  public :: set_adios_group_size
  public :: set_selection_boundingbox
  public :: delete_adios_selection
  public :: delete_adios_group
  public :: get_adios_group
  public :: flush_adios_group_all

  ! check
  public :: check_adios_err
  public :: show_adios_file_variables
  public :: get_adios_filename

#endif  /* USE_ADIOS or USE_ADIOS2 */

contains

!-------------------------------------------------------------------------------
!
! public ADIOS wrapper routines (also available without adios compilation support)
!
!-------------------------------------------------------------------------------

  subroutine initialize_adios()

!> Initialize ADIOS and setup the xml output file

#if defined(USE_ADIOS)
  use constants, only: ADIOS_BUFFER_SIZE_IN_MB
!#ifdef ADIOS_VERSION_OLD
  ! ADIOS versions <= 1.9
  ! adios_set_max_buffer_size not defined yet
!#else
  ! ADIOS versions >= 1.10
  !use adios_write_mod, only: adios_set_max_buffer_size
!#endif
#endif

  implicit none

  ! local parameters
#if defined(USE_ADIOS) || defined(USE_ADIOS2)
  integer :: ier
#endif

  TRACE_ADIOS('initialize_adios')

  ! initializes
  is_adios_initialized = .false.
  is_adios_version1 = .false.
  is_adios_version2 = .false.
  comm_adios = 0
  myrank_adios = -1
  sizeprocs_adios = 0

#if defined(USE_ADIOS) || defined(USE_ADIOS2)
  ! gets MPI communicator for adios calls
  call world_duplicate(comm_adios)

  ! gets rank from (duplicate) adios communicator
  call world_rank_comm(myrank_adios,comm_adios)

  ! number of MPI processes
  call world_size_comm(sizeprocs_adios,comm_adios)

  ! checks
  if (sizeprocs_adios == 0) &
    stop 'Error adios initialization got zero processes'

#if defined(USE_ADIOS)
  ! ADIOS 1
  is_adios_version1 = .true.

  call adios_init_noxml (comm_adios, ier)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  !       e.g., version 1.5.0 returns 1 here
  !print *,'adios init return: ',ier
  if (ier /= 0) stop 'Error setting up ADIOS: calling adios_init_noxml() routine failed. Please use an ADIOS version >= 1.6'

! ask/check at configuration step for adios version 1.10 or higher?
#ifdef ADIOS_VERSION_OLD
  ! ADIOS versions <= 1.9
  ! note: for newer versions ( >= 1.10), this will produce a warning might not be supported anymore
  call adios_allocate_buffer(ADIOS_BUFFER_SIZE_IN_MB, ier)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  !       e.g., version 1.5.0 returns 1 if called first time, 0 if already called
  !print *,'adios allocate buffer return: ',ier
  !call check_adios_err(ier,"Error allocate buffer")
#else
  ! ADIOS versions >= 1.10
  call adios_set_max_buffer_size(ADIOS_BUFFER_SIZE_IN_MB)
#endif

  ! initializes file handles
  myadios_file = 0
  myadios_group = 0
  myadios_fwd_file = 0
  myadios_fwd_group = 0
  myadios_val_file = 0
  myadios_val_group = 0

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  is_adios_version2 = .true.

  ! Create adios handler passing the communicator, debug mode and error flag
  ! adios2 duplicates the communicator for its internal use
  if (USE_ADIOS2_DEBUG_MODE) then
    call adios2_init(myadios2_obj, comm_adios, adios2_debug_mode_on, ier)
  else
    call adios2_init(myadios2_obj, comm_adios, adios2_debug_mode_off, ier)
  endif
  if (ier /= 0) stop 'Error setting up ADIOS2: calling adios2_init() routine failed'

#endif

  ! sets flag
  is_adios_initialized = .true.

  ! for undo att
  is_initialized_fwd_group = .false.

#else
  ! no adios compilation support
  ! gets rank
  call world_rank(myrank_adios)

  ! compilation without ADIOS support
  if (myrank_adios == 0) then
    print *, "Error: ADIOS enabled without ADIOS Support."
    print *, "To enable ADIOS support, reconfigure with --with-adios or --with-adios2 flag."
  endif
  ! safety stop
  call exit_MPI(myrank_adios,"Error ADIOS manager: intitialize called without compilation support")

#endif

  end subroutine initialize_adios

!
!-------------------------------------------------------------------------------
!

  subroutine finalize_adios()

!> Finalize ADIOS. Must be called once everything is written down.

  implicit none

  ! local parameters
#if defined(USE_ADIOS) || defined(USE_ADIOS2)
  integer :: ier
  logical, external :: is_valid_comm
#endif

  TRACE_ADIOS('finalize_adios')

  ! synchronizes all first
  call synchronize_all_comm(comm_adios)

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! checks to close file at the end of run
  if (is_initialized_fwd_group .and. myadios_fwd_file /= 0) then
    ! debug
    print *,'Warning: rank ',myrank_adios,' has myadios_fwd_file still open, please check reading adios forward snapshots.'
    ! this might crash with a seg-fault if still open for kernel runs...
    call adios_close(myadios_fwd_file, ier)
    call check_adios_err(ier,"Error closing adios fwd file")
  endif

  ! double-check file handles
  if (myadios_file /= 0) stop 'Error adios file myadios_file still open, please check.'
  if (myadios_val_file /= 0) stop 'Error adios file myadios_val_file still open, please check.'
  if (myadios_fwd_file /= 0) stop 'Error adios file myadios_fwd_file still open, please check.'

  ! wait until finalized, synchronizes all using adios communicator
  call synchronize_all_comm(comm_adios)

  ! finalize
  call adios_finalize(myrank_adios, ier)
  if (ier /= 0 ) stop 'Error cleaning up ADIOS: calling adios_finalize() routine failed'

  ! frees (duplicate) MPI communicator
  if (is_valid_comm(comm_adios)) call world_comm_free(comm_adios)

#elif defined(USE_ADIOS2)
  ! checks to close file at the end of run
  if (is_initialized_fwd_group .and. myadios_fwd_file%valid) then
    call adios2_close(myadios_fwd_file, ier)
    call check_adios_err(ier,"Error closing adios fwd file")
  endif

  ! double-check file handles
  if (myadios_file%valid) stop 'Error adios2 file myadios_file still open, please check.'
  if (myadios_val_file%valid) stop 'Error adios2 file myadios_val_file still open, please check.'
  if (myadios_fwd_file%valid) stop 'Error adios2 file myadios_fwd_file still open, please check.'

  ! wait until finalized, synchronizes all using adios communicator
  call synchronize_all_comm(comm_adios)

  ! finalize
  call adios2_finalize(myadios2_obj, ier)
  if (ier /= 0 ) stop 'Error cleaning up ADIOS2: calling adios2_finalize() routine failed'

  ! adios2 internally calls MPI_Comm_dup and frees (duplicate) MPI communicator when finalizing.
  ! since we called an explicit duplicator, we also free it.
  if (is_valid_comm(comm_adios)) call world_comm_free(comm_adios)

#else
  ! safety stop
  call exit_MPI(myrank_adios,"Error ADIOS manager: finalize called without compilation support")

#endif

  end subroutine finalize_adios


!-------------------------------------------------------------------------------
!
! ADIOS wrapper routines (only available with adios compilation support)
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!
! file opening
!
!-------------------------------------------------------------------------------

#if defined(USE_ADIOS) || defined(USE_ADIOS2)
! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools

  subroutine open_file_adios_read_and_init_method(adios_handle,adios_group,filename)

! useful to read in data from the same number of processors
! as the data was written from

#if defined(USE_ADIOS2)
  use constants, only: ADIOS2_ENGINE_DEFAULT,ADIOS2_ENGINE_PARAMS_DEFAULT
#endif

  implicit none

#if defined(USE_ADIOS)
  ! adios
  integer(kind=8), intent(out) :: adios_handle
  integer(kind=8), intent(in) :: adios_group
#elif defined(USE_ADIOS2)
  ! adios2
  type(adios2_engine), intent(out) :: adios_handle
  type(adios2_io), intent(inout) :: adios_group
  ! debug
  !integer(kind=8) :: steps
#endif
  character(len=*), intent(in) :: filename

  ! local parameters
  integer :: ier

  TRACE_ADIOS_ARG('open_file_adios_read_and_init_method: file '//trim(filename)//' - rank ',myrank_adios)

  ! check
  if (.not. is_adios_initialized) &
    stop 'Error intializing for adios read: please initialize adios first using intitialize_adios() routine'

  ! initializes read method
#if defined(USE_ADIOS)
  ! ADIOS 1
  call adios_read_init_method(ADIOS_READ_METHOD_BP, comm_adios, ADIOS_VERBOSITY, ier)
  call check_adios_err(ier,"Error initializing read adios for file: "//trim(filename))

  ! opens file
  call adios_read_open_file(adios_handle, trim(filename), 0, comm_adios, ier)
  call check_adios_err(ier,"Error opening adios file for reading: "//trim(filename))

  ! to avoid compiler warning
  ier = adios_group

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! in case no io_group given, will need one for inquiry
  if (.not. adios_group%valid) then
    call adios2_declare_io(adios_group, myadios2_obj, "Reader", ier)
    call check_adios_err(ier,"Error declaring an ADIOS2 IO group in open_file_adios_read_and_init_method()")

    ! Set engine and parameters
    call adios2_set_engine(adios_group, ADIOS2_ENGINE_DEFAULT, ier)
    call check_adios_err(ier,"Error setting engine for ADIOS2 IO group in open_file_adios_read_and_init_method()")

    ! Set default parameters
    call adios2_set_parameters(adios_group, ADIOS2_ENGINE_PARAMS_DEFAULT, ier)
    call check_adios_err(ier,"Error setting parameters for ADIOS2 IO group in open_file_adios_read_and_init_method()")
  else
    ! debug
    !print *,'debug adios: open_file_adios_read_and_init_method() has valid adios group'
  endif

  ! synchronizes all processes to make sure engine & parameters has been set for all procs
  call synchronize_all_comm(comm_adios)

  ! Open the handle to file containing all the ADIOS variables for the current io group
  call adios2_open(adios_handle, adios_group, trim(filename), adios2_mode_read, comm_adios, ier)
  call check_adios_err(ier,"Error opening adios2 file "//trim(filename))

  ! debug
  !call adios2_steps(steps,adios_handle,ier)
  !call check_adios_err(ier,"Error adios2 getting steps for file "//trim(filename))
  !print *,'debug adios: file steps ',steps

#endif

  ! synchronizes all processes
  call synchronize_all_comm(comm_adios)

  end subroutine open_file_adios_read_and_init_method

#endif
!
!---------------------------------------------------------------------------------
!
#if defined(USE_ADIOS) || defined(USE_ADIOS2)
! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools

  subroutine open_file_adios_read_only_rank(adios_handle,adios_group,rank,filename)

! only single process is reading, useful for file inquiry

  implicit none

#if defined(USE_ADIOS)
  ! adios
  integer(kind=8), intent(out) :: adios_handle
  integer(kind=8), intent(in) :: adios_group
#elif defined(USE_ADIOS2)
  ! adios2
  type(adios2_engine), intent(out) :: adios_handle
  type(adios2_io), intent(in) :: adios_group
#endif
  integer, intent(in) :: rank
  character(len=*), intent(in) :: filename

  ! local parameters
  integer :: ier
  integer :: comm_dummy
  character(len=128) :: name !uses a string copy, trying to prevent a memory corruption issue somewhere in adios...
  logical, parameter :: DEBUG = .false.

  TRACE_ADIOS_ARG('open_file_adios_read_only_rank: file '//trim(filename)//' - rank ',myrank_adios)

  ! only specified rank proceeds
  if (myrank_adios /= rank) return

  ! gets MPI communicator for only a single process
  call world_get_comm_self(comm_dummy)

  ! copies name
  if (len_trim(filename) > 128) then
    stop 'Error adios filename provided is too long'
  else
    name = trim(filename)
  endif

  ! initializes read method
#if defined(USE_ADIOS)
  call adios_read_init_method(ADIOS_READ_METHOD_BP, comm_dummy, ADIOS_VERBOSITY, ier)
  call check_adios_err(ier,"Error initializing read adios by main for file: "//trim(name))

  ! opens file
  call adios_read_open_file(adios_handle, trim(name), 0, comm_dummy, ier)
  call check_adios_err(ier,"Error opening adios file for reading: "//trim(name))

  ! to avoid compiler warning
  ier = adios_group

#elif defined(USE_ADIOS2)
  ! MPI version only to pass a communicator other than the one from adios_init
  call adios2_open(adios_handle, adios_group, trim(name), adios2_mode_read, comm_dummy, ier)
  call check_adios_err(ier,"Error opening adios2 file with MPI comm for reading: "//trim(name))

#endif

  ! shows file contents
  if (DEBUG) then
    call show_adios_file_variables(adios_handle,adios_group,name)
  endif

  ! do not synchronize across adios processes
  ! this is only called by a single process (using MPI_COMM_SELF for adios i/o).

  end subroutine open_file_adios_read_only_rank

#endif
!
!---------------------------------------------------------------------------------
!
#if defined(USE_ADIOS) || defined(USE_ADIOS2)
! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools

  subroutine open_file_adios_read(adios_handle,adios_group,filename)

! useful to read in data from the same number of processors as the data was written from;
! assumes read init method has been set already, will only open file

  implicit none

#if defined(USE_ADIOS)
  ! adios
  integer(kind=8), intent(out) :: adios_handle
  integer(kind=8), intent(in) :: adios_group
#elif defined(USE_ADIOS2)
  ! adios2
  type(adios2_engine), intent(out) :: adios_handle
  type(adios2_io), intent(in) :: adios_group
#endif
  character(len=*), intent(in) :: filename

  ! local parameters
  integer :: ier

  TRACE_ADIOS_ARG('open_file_adios_read: file '//trim(filename)//' - rank ',myrank_adios)

  ! check
  if (.not. is_adios_initialized) &
    stop 'Error intializing for adios read: please initialize adios first using intitialize_adios() routine'

  ! initializes read method
#if defined(USE_ADIOS)
  ! ADIOS 1

  ! opens file
  call adios_read_open_file(adios_handle, trim(filename), 0, comm_adios, ier)
  call check_adios_err(ier,"Error opening adios file for reading: "//trim(filename))

  ! to avoid compiler warning
  ier = adios_group

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! Open the handle to file containing all the ADIOS variables for the current io group
  call adios2_open(adios_handle, adios_group, trim(filename), adios2_mode_read, comm_adios, ier)
  call check_adios_err(ier,"Error opening adios2 file "//trim(filename))
#endif

  ! synchronizes all processes
  call synchronize_all_comm(comm_adios)

  end subroutine open_file_adios_read

#endif
!
!-------------------------------------------------------------------------------
!
#if defined(USE_ADIOS) || defined(USE_ADIOS2)
! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools

  subroutine open_file_adios_write(adios_handle,adios_group,filename,group_name)

! opens adios file for writing

  implicit none

#if defined(USE_ADIOS)
  ! adios
  integer(kind=8), intent(out) :: adios_handle
  integer(kind=8), intent(in) :: adios_group
#elif defined(USE_ADIOS2)
  ! adios2
  type(adios2_engine), intent(out) :: adios_handle
  type(adios2_io), intent(in) :: adios_group
#endif
  character(len=*), intent(in) :: filename
  character(len=*) :: group_name

  ! local parameters
  integer :: ier

  TRACE_ADIOS_ARG('open_file_adios_write: file '//trim(filename)//' (for writing) - rank ',myrank_adios)

  ! check
  if (.not. is_adios_initialized) &
    stop 'Error intializing for adios write: please initialize adios first using initialize_adios() routine'

  ! opens file
#if defined(USE_ADIOS)
  ! ADIOS 1
  ! checks group
  if (len_trim(group_name) == 0) stop 'Error: group name has zero length in open_file_adios_write() routine'

  call adios_open(adios_handle, group_name, trim(filename), "w", comm_adios, ier)
  call check_adios_err(ier,"Error opening adios file for writing: "//trim(filename)//" group: "//trim(group_name))

  ! to avoid compiler warning
  ier = adios_group

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! checks group
  if (.not. adios_group%valid) stop 'Invalid adios2 group in open_file_adios_write() routine'

  ! Open the handle to file containing all the ADIOS variables
  call adios2_open(adios_handle, adios_group, trim(filename), adios2_mode_write, comm_adios, ier)
  call check_adios_err(ier,"Error opening adios file for writing: "//trim(filename)//" group: "//trim(group_name))

#endif

  ! synchronizes all processes
  call synchronize_all_comm(comm_adios)

  end subroutine open_file_adios_write

#endif
!
!-------------------------------------------------------------------------------
!
#if defined(USE_ADIOS) || defined(USE_ADIOS2)
! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools

  subroutine open_file_adios_write_append(adios_handle,adios_group,filename,group_name)

! open adios file for appending data

  implicit none

#if defined(USE_ADIOS)
  ! adios
  integer(kind=8), intent(out) :: adios_handle
  integer(kind=8), intent(in) :: adios_group
#elif defined(USE_ADIOS2)
  ! adios2
  type(adios2_engine), intent(out) :: adios_handle
  type(adios2_io), intent(in) :: adios_group
  !integer :: step_status
#endif
  character(len=*), intent(in) :: filename
  character(len=*) :: group_name

  ! local parameters
  integer :: ier

  TRACE_ADIOS_ARG('open_file_adios_write_append: file '//trim(filename)//' (for appending) - rank ',myrank_adios)

  ! check
  if (.not. is_adios_initialized) &
    stop 'Error intializing for adios write: please initialize adios first using intitialize_adios() routine'

  ! opens file
#if defined(USE_ADIOS)
  ! ADIOS 1
  ! checks group
  if (len_trim(group_name) == 0) &
    stop 'Error: group name has zero length in open_file_adios_write_append() routine'

  ! opens file in append mode
  call adios_open(adios_handle, group_name, trim(filename), "a", comm_adios, ier)
  call check_adios_err(ier,"Error opening adios file for appending: "//trim(filename)//" group: "//trim(group_name))

  ! to avoid compiler warning
  ier = adios_group

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! checks group
  if (.not. adios_group%valid) &
    stop 'Invalid adios2 group in open_file_adios_write_append() routine'

  ! opens file in append mode
  call adios2_open(adios_handle, adios_group, trim(filename), adios2_mode_append, comm_adios, ier)
  call check_adios_err(ier,"Error opening adios2 file for appending: "//trim(filename)//" group: "//trim(group_name))

  ! note: opening the same file will increase the step count for all subsequent variable writes.
  !       this can cause problems when reading back in and selecting a corresponding variable selection.
  ! this has no effect though...
  !call adios2_begin_step(adios_handle, adios2_step_mode_append, -1.0, step_status, ier)
  !call check_adios_err(ier,"Error begin step for appending: "//trim(filename)//" group: "//trim(group_name))
  !if (step_status /= adios2_step_status_ok) stop 'Error adios2 begin step for appending to file'

#endif

  ! synchronizes all processes
  call synchronize_all_comm(comm_adios)

  end subroutine open_file_adios_write_append

#endif


!-------------------------------------------------------------------------------
!
! file closing
!
!---------------------------------------------------------------------------------

#if defined(USE_ADIOS) || defined(USE_ADIOS2)
! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools

  subroutine close_file_adios(adios_handle)

  implicit none

#if defined(USE_ADIOS)
  ! adios
  integer(kind=8), intent(inout) :: adios_handle
#elif defined(USE_ADIOS2)
  ! adios2
  type(adios2_engine), intent(inout) :: adios_handle
#endif

  ! local parameters
  integer :: ier

  TRACE_ADIOS('close_file_adios')

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! closes file
  call adios_close(adios_handle, ier)
  call check_adios_err(ier,"Error closing adios file")

  ! sets explicitly to zero
  adios_handle = 0

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! close file
  call adios2_close(adios_handle, ier)
  call check_adios_err(ier,"Error closing adios2 file")

#endif

  ! synchronizes all processes
  call synchronize_all_comm(comm_adios)

  end subroutine close_file_adios

#endif
!
!---------------------------------------------------------------------------------
!
#if defined(USE_ADIOS) || defined(USE_ADIOS2)
! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools

  subroutine close_file_adios_read_and_finalize_method(adios_handle)

  implicit none

#if defined(USE_ADIOS)
  ! adios
  integer(kind=8), intent(inout) :: adios_handle
#elif defined(USE_ADIOS2)
  ! adios2
  type(adios2_engine), intent(inout) :: adios_handle
#endif

  ! local parameters
  integer :: ier

  TRACE_ADIOS('close_file_adios_read_and_finalize_method')

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! closes file
  call adios_read_close(adios_handle,ier)
  if (ier /= 0) stop 'Error helper adios read close in close_file_adios_read_and_finalize() routine'

  ! sets explicitly to zero
  adios_handle = 0

  ! finalizes file opened with adios_read_init_method()
  call adios_read_finalize_method(ADIOS_READ_METHOD_BP, ier)
  if (ier /= 0) stop 'Error helper adios read finalize failed'

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! no special case, file just has been opened with adios2_mode_read flag
  call adios2_close(adios_handle, ier)
  call check_adios_err(ier,"Error closing adios file close_file_adios_read_and_finalize() routine")

  ! note: closing the file will not delete the adios io yet.
  !       thus, the adios group handlers could still be valid when opening a file with the init method next time.
  !       this leads to issues when used together with begin/end steps, producing adios2 errors like:
  !         ..
  !         ERROR: variable reg2/rhostore/offset exists in IO object Reader, in call to DefineVariable
  !         ..
  !
  !       we will use explicit calls to delete_adios_group() routine as a work-around.

#endif

  ! synchronizes all processes
  call synchronize_all_comm(comm_adios)

  end subroutine close_file_adios_read_and_finalize_method

#endif
!
!---------------------------------------------------------------------------------
!
#if defined(USE_ADIOS) || defined(USE_ADIOS2)
! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools

  subroutine close_file_adios_read_and_finalize_method_only_rank(adios_handle,rank)

! only single process closes

  implicit none

#if defined(USE_ADIOS)
  ! adios
  integer(kind=8), intent(inout) :: adios_handle
#elif defined(USE_ADIOS2)
  ! adios2
  type(adios2_engine), intent(inout) :: adios_handle
#endif
  integer, intent(in) :: rank
  ! local parameters
  integer :: ier

  TRACE_ADIOS('close_file_adios_read_and_finalize')

  ! only specified rank proceeds
  if (myrank_adios /= rank) return

#if defined(USE_ADIOS)
  ! ADIOS 1
  call adios_read_close(adios_handle,ier)
  if (ier /= 0 ) stop 'Error helper adios read close in close_file_adios_read_and_finalize() routine'

  ! sets explicitly to zero
  adios_handle = 0

  ! finalizes file opened with adios_read_init_method()
  call adios_read_finalize_method(ADIOS_READ_METHOD_BP, ier)
  if (ier /= 0 ) stop 'Error helper adios read finalize'

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! no special case, file just has been opened with adios2_mode_read flag
  call adios2_close(adios_handle, ier)
  call check_adios_err(ier,"Error closing adios file close_file_adios_read_and_finalize() routine")

#endif

  ! do not synchronize after closing, as it will be called also by a single process in interpolate_model.F90
  ! and thus would stall the program execution.
  !
  ! do not: synchronizes all processes
  ! >call synchronize_all_comm(comm_adios)

  end subroutine close_file_adios_read_and_finalize_method_only_rank

#endif

!
!---------------------------------------------------------------------------------
!
#if defined(USE_ADIOS) || defined(USE_ADIOS2)
! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools

  subroutine close_file_adios_read(adios_handle)

  implicit none

#if defined(USE_ADIOS)
  ! adios
  integer(kind=8), intent(inout) :: adios_handle
#elif defined(USE_ADIOS2)
  ! adios2
  type(adios2_engine), intent(inout) :: adios_handle
#endif

  ! local parameters
  integer :: ier

  TRACE_ADIOS('close_file_adios_read')

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! check
  if (adios_handle == 0) stop 'Error: adios handle invalid in close_file_adios_read()'

  ! closes
  call adios_read_close(adios_handle,ier)
  if (ier /= 0 ) stop 'Error helper adios read close with file handle'

  ! sets explicitly to zero
  adios_handle = 0

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  if (.not. adios_handle%valid) stop 'Error: adios2 handle invalid in close_file_adios_read()'

  ! no special case, file just has been opened with adios2_mode_read flag
  call adios2_close(adios_handle, ier)
  call check_adios_err(ier,"Error closing adios2 file with handle")

#endif

  ! synchronizes all processes
  call synchronize_all_comm(comm_adios)

  end subroutine close_file_adios_read

#endif


!-------------------------------------------------------------------------------
!
! I/O groups
!
!---------------------------------------------------------------------------------

#if defined(USE_ADIOS) || defined(USE_ADIOS2)
! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools

  subroutine init_adios_group(adios_group,group_name)

! useful to read in data from the same number of processors
! as the data was written from

#if defined(USE_ADIOS)
  use constants, only: ADIOS_TRANSPORT_METHOD,ADIOS_METHOD_PARAMS
#elif defined(USE_ADIOS2)
  use constants, only: ADIOS2_ENGINE_DEFAULT,ADIOS2_ENGINE_PARAMS_DEFAULT
#endif

  implicit none

#if defined(USE_ADIOS)
  integer(kind=8), intent(inout) :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_io), intent(inout) :: adios_group
#endif

  character(len=*), intent(in) :: group_name

  ! local parameters
  integer :: ier

  TRACE_ADIOS_ARG('init_adios_group: group '//trim(group_name)//' - rank ',myrank_adios)

  ! check
  if (.not. is_adios_initialized) &
    stop 'Error intializing for adios group init: please initialize adios first using initialize_adios() routine'

  ! initializes adios group
#if defined(USE_ADIOS)
  ! ADIOS 1
  call adios_declare_group(adios_group, group_name, '', 0, ier)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  call check_adios_err(ier,"Error declare group")

  ! We set the transport method to 'MPI'. This seems to be the correct choice
  ! for now. We might want to move this to the constant.h file later on.
  call adios_select_method(adios_group, ADIOS_TRANSPORT_METHOD, ADIOS_METHOD_PARAMS, '', ier)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  call check_adios_err(ier,"Error select method")

#elif defined(USE_ADIOS2)
  ! Create the ADIOS IO group which will contain all variables and attributes
  call adios2_declare_io(adios_group, myadios2_obj, group_name, ier)
  call check_adios_err(ier,"Error declaring an ADIOS2 IO group in init_adios_group()")

  ! Set engine and parameters
  call adios2_set_engine(adios_group, ADIOS2_ENGINE_DEFAULT, ier)
  call check_adios_err(ier,"Error setting engine for ADIOS2 IO group in init_adios_group()")

  ! Set default parameters
  call adios2_set_parameters(adios_group, ADIOS2_ENGINE_PARAMS_DEFAULT, ier)
  call check_adios_err(ier,"Error setting parameters for ADIOS2 IO group in init_adios_group()")

  ! sets current group in case we have no group yet set
  if (.not. myadios_group%valid) myadios_group = adios_group
#endif

  end subroutine init_adios_group

#endif
!
!---------------------------------------------------------------------------------
!
#if defined(USE_ADIOS) || defined(USE_ADIOS2)
! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools

  subroutine init_adios_group_undo_att(adios_group,group_name)

#if defined(USE_ADIOS)
  use constants, only: ADIOS_TRANSPORT_METHOD_UNDO_ATT,ADIOS_METHOD_PARAMS_UNDO_ATT
#elif defined(USE_ADIOS2)
  use constants, only: ADIOS2_ENGINE_UNDO_ATT,ADIOS2_ENGINE_PARAMS_UNDO_ATT
#endif

  implicit none

#if defined(USE_ADIOS)
  integer(kind=8), intent(inout) :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_io), intent(inout) :: adios_group
#endif
  character(len=*), intent(in) :: group_name

  ! local parameters
  integer :: ier

  TRACE_ADIOS_ARG('init_adios_group_undo_att: group '//trim(group_name)//' - rank ',myrank_adios)

  ! initializes adios group
#if defined(USE_ADIOS)
  ! ADIOS 1
  call adios_declare_group(adios_group, group_name, "iter", 0, ier)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  call check_adios_err(ier,"Error declare group")

  ! sets transport method
  call adios_select_method(adios_group, ADIOS_TRANSPORT_METHOD_UNDO_ATT, ADIOS_METHOD_PARAMS_UNDO_ATT, '', ier)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  call check_adios_err(ier,"Error select method")

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! note: no special case, same engine & parameters as for "normal" io groups.
  !       we could try out different transport methods & parameters for performance

  ! Create the ADIOS IO group which will contain all variables and attributes
  call adios2_declare_io(adios_group, myadios2_obj, group_name, ier)
  call check_adios_err(ier,"Error declaring an ADIOS2 IO group in init_adios_group_undo_att()")

  ! Set engine and parameters: ADIOS2_ENGINE_UNDO_ATT or ADIOS2_ENGINE_DEFAULT
  call adios2_set_engine(adios_group, ADIOS2_ENGINE_UNDO_ATT, ier)
  call check_adios_err(ier,"Error setting engine for ADIOS2 IO group in init_adios_group_undo_att()")

  ! Set parameters: ADIOS2_ENGINE_PARAMS_UNDO_ATT or ADIOS2_ENGINE_PARAMS_DEFAULT
  call adios2_set_parameters(adios_group, ADIOS2_ENGINE_PARAMS_UNDO_ATT, ier)
  call check_adios_err(ier,"Error setting parameters for ADIOS2 IO group in init_adios_group_undo_att()")

#endif

  end subroutine init_adios_group_undo_att

#endif
!
!---------------------------------------------------------------------------------
!
#if defined(USE_ADIOS) || defined(USE_ADIOS2)
! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools

  subroutine set_adios_group_size(adios_handle,group_size_inc)

  implicit none

#if defined(USE_ADIOS)
  ! adios
  integer(kind=8), intent(inout) :: adios_handle
#elif defined(USE_ADIOS2)
  ! adios2
  type(adios2_engine), intent(inout) :: adios_handle
#endif
  integer(kind=8), intent(in) :: group_size_inc

  ! local parameters
  integer :: ier
  integer(kind=8) :: totalsize

  TRACE_ADIOS_L2('set_adios_group_size')

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! checks if file handle valid
  if (adios_handle == 0) stop 'Invalid ADIOS file handle argument in set_adios_group_size()'

  call adios_group_size(adios_handle, group_size_inc, totalsize, ier)
  if (ier /= 0 ) stop 'Error calling adios_group_size() routine failed'

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! no need to explicitly specify io group size (adios_group_size has become optional since adios 1.10)
  ! check to avoid compiler warning
  if (.not. adios_handle%valid) stop 'Invalid ADIOS2 file handle argument in set_adios_group_size()'

  ! to avoid compiler warning
  totalsize = group_size_inc
  ier = 0
#endif

  end subroutine set_adios_group_size

#endif
!
!---------------------------------------------------------------------------------
!
#if defined(USE_ADIOS) || defined(USE_ADIOS2)
! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools

  subroutine get_adios_group(adios_group,group_name,adios_handle)

  implicit none

#if defined(USE_ADIOS)
  ! adios
  integer(kind=8), intent(inout) :: adios_group
  integer(kind=8), intent(in) :: adios_handle
  ! groups
  character (len=128), dimension(:), allocatable :: fnamelist
  integer :: group_count,i
#elif defined(USE_ADIOS2)
  ! adios2
  type(adios2_io), intent(inout) :: adios_group
  type(adios2_engine), intent(in) :: adios_handle
#endif
  character(len=*), intent(in) :: group_name

  ! local parameters
  integer :: ier

  TRACE_ADIOS('get_adios_group')

  ! check
  if (.not. is_adios_initialized) &
    stop 'Error: ADIOS is not intialized for get_adios_group() call'

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! check
  if (adios_handle == 0) stop 'Invalid ADIOS file handle in get_adios_group()'

  ! gets groups
  call adios_inq_ngroups(adios_handle, group_count, ier)
  if (ier /= 0) stop 'Error calling adios_inq_ngroups()'

  if (group_count > 0) then
    allocate (fnamelist(group_count),stat=ier)
    if (ier /= 0) stop 'Error allocating namelist array'

    ! gets group names
    call adios_inq_groupnames(adios_handle, fnamelist, ier)
    if (ier /= 0) stop 'Error calling adios_inq_groupnames()'

    ! selects group
    do i = 1,group_count
      if (trim(group_name) == trim(fnamelist(i))) then
        ! adios group ids start from 0 to N-1
        adios_group = i - 1
        exit
      endif
    enddo
    deallocate (fnamelist)
  else
    ! group not found
    adios_group = 0
    print *, 'Error: no group found in get_adios_group() by name :',trim(group_name)
    stop 'adios group not found in get_adios_group()'
  endif

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! check
  if (.not. myadios2_obj%valid) &
    stop 'Invalid ADIOS2 object component in get_adios_group()'
  if (.not. adios_handle%valid) &
    stop 'Invalid ADIOS2 file handle in get_adios_group()'

  ! gets io group by name
  call adios2_at_io(adios_group, myadios2_obj, group_name, ier)
  call check_adios_err(ier,"Error calling adios2_at_io() for group "//trim(group_name))

  if (.not. adios_group%valid) then
    print *,'Error: adios2 group not found for group name: ',trim(group_name)
    stop 'Invalid ADIOS2 group in get_adios_group()'
  endif

#endif

  end subroutine get_adios_group

#endif
!
!---------------------------------------------------------------------------------
!
#if defined(USE_ADIOS) || defined(USE_ADIOS2)
! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools

  subroutine flush_adios_group_all(adios_group)

  implicit none

#if defined(USE_ADIOS)
  integer(kind=8), intent(inout) :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_io), intent(inout) :: adios_group
#endif

  ! local parameters
  integer :: ier

  TRACE_ADIOS('flush_adios_group_all')

  ! initializes adios group
#if defined(USE_ADIOS)
  ! ADIOS 1
  ! no flush all
  ! to avoid compiler warning
  ier = adios_group

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! flush all engines
  if (adios_group%valid) then
    call adios2_flush_all_engines(adios_group,ier)
    if (ier /= 0 ) stop 'Error cleaning up ADIOS2: calling adios2_flush_all_engines() failed'
  endif

#endif

  end subroutine flush_adios_group_all

#endif

!-------------------------------------------------------------------------------
!
! data selections
!
!---------------------------------------------------------------------------------

#if defined(USE_ADIOS) || defined(USE_ADIOS2)
! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools

  subroutine set_selection_boundingbox(sel,start,count)

  implicit none

  integer(kind=8), intent(inout) :: sel
  integer(kind=8), dimension(1), intent(in) :: start, count

  TRACE_ADIOS_L2('set_selection_boundingbox')

#if defined(USE_ADIOS)
  ! ADIOS 1
  call adios_selection_boundingbox(sel , 1, start, count)

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! no need to explicitly specify, will work with start/count directly

  ! to avoid compiler warning
  sel = count(1)
  sel = start(1)
  sel = 0

#endif

  end subroutine set_selection_boundingbox

#endif
!
!---------------------------------------------------------------------------------
!
#if defined(USE_ADIOS) || defined(USE_ADIOS2)
! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools

  subroutine delete_adios_selection(sel)

  implicit none

  integer(kind=8), intent(inout) :: sel

  TRACE_ADIOS_L2('delete_adios_selection')

#if defined(USE_ADIOS)
  ! ADIOS 1
  call adios_selection_delete(sel)

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! unused, just set it to zero
  sel = 0
#endif

  end subroutine delete_adios_selection

#endif
!
!---------------------------------------------------------------------------------
!
#if defined(USE_ADIOS) || defined(USE_ADIOS2)
! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools

  subroutine delete_adios_group(adios_group,group_name)

! removes all variables from group and deletes it

  implicit none

#if defined(USE_ADIOS)
  integer(kind=8), intent(inout) :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_io), intent(inout) :: adios_group
#endif
  character(len=*), intent(in) :: group_name

  ! local parameters
  integer :: ier
  logical :: result

  TRACE_ADIOS_ARG('delete_adios_group: group '//trim(group_name)//' - rank ',myrank_adios)

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! deletes all variable definitions from a group
  ! we can define a new set of variables for the next output step
  call adios_delete_vardefs(adios_group, ier)
  if (ier /= 0) then
    print *,'Error: adios delete group with group name ',trim(group_name),' failed'
    stop 'Error helper adios delete group variables failed in delete_adios_group() routine'
  endif

  ! to avoid compiler warnings
  result = .true.
  ! nothing left to do, no explicit delete routine for the group in ADIOS1

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! note: closing a file will not delete the adios io yet.
  !       thus, the adios group handlers could still be valid when opening a file with the init method next time.
  !       this leads to issues when used together with begin/end steps, producing adios2 errors like:
  !         ..
  !         ERROR: variable reg2/rhostore/offset exists in IO object Reader, in call to DefineVariable
  !         ..
  !
  ! cleans out the group
  ! deletes all variable definitions from a group
  call adios2_remove_all_variables(adios_group, ier)
  call check_adios_err(ier,"Error removing group variables in delete_adios_group() routine")

  ! removes io
  ! (will recreate a new group in a next new open-and-init call)
  call adios2_remove_io(result, myadios2_obj, group_name, ier)
  ! checks result flag
  if (result .neqv. .true.) stop 'Error failed removing io in delete_adios_group() routine'
  call check_adios_err(ier,"Error removing io in delete_adios_group() routine")

  ! explicitly resets group
  if (adios_group%valid) adios_group%valid = .false.

  ! this would fail, as io object is no more defined...
  !call adios2_at_io(adios_group, myadios2_obj, group_name, ier)
  !call check_adios_err(ier,"Error getting io in delete_adios_group() routine")

  !! This below doesn't work as it will lead to issues with xsum_kernels_adios and possibly other tools...
  !!
  !! Avoid removing all ios: this will lead to problems when one io was opened for reading, while a second one
  !!                         should stay valid for writing out (see e.g. sum_kernels.F90)
  !!
  !! removes io handlers created with adios2_declare_io()
  !!if (myadios_group%valid .or. myadios_fwd_group%valid .or. myadios_val_group%valid) then
  !!  ! removes ios
  !!  call adios2_remove_all_ios(myadios2_obj, ier)
  !!  call check_adios_err(ier,"Error removing all ios in delete_adios_group() routine")
  !!  ! reset groups
  !!  myadios_group%valid = .false.
  !!  myadios_fwd_group%valid = .false.
  !!  myadios_val_group%valid = .false.
  !!endif

#endif

  end subroutine delete_adios_group

#endif


!-------------------------------------------------------------------------------
!
! inquiry / user outputs
!
!---------------------------------------------------------------------------------

#if defined(USE_ADIOS) || defined(USE_ADIOS2)
! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools

  subroutine show_adios_file_variables(adios_handle,adios_group,filename)

! file inquiry showing all adios file variables

  implicit none

#if defined(USE_ADIOS)
  ! adios
  integer(kind=8), intent(in) :: adios_handle
  integer(kind=8), intent(in) :: adios_group
#elif defined(USE_ADIOS2)
  ! adios2
  type(adios2_engine), intent(in) :: adios_handle
  type(adios2_io), intent(in) :: adios_group
#endif
  character(len=*), intent(in) :: filename

  ! local parameters
  integer :: ier
  integer :: variable_count, attribute_count

  ! file inquiry
#if defined(USE_ADIOS)
  ! ADIOS 1
  integer :: timestep_first, timestep_last
  character (len=128), dimension(:), allocatable :: fnamelist
  character (len=128) :: vname
  integer :: vtype, vnsteps, vndim
  integer(kind=8), dimension(3) :: dims
  integer :: group_count,i

  ! user output
  if (myrank_adios == 0) then
    print *
    print *,'show ADIOS file variables: ',trim(filename)
    print *
  endif

  ! file inquiry
  call adios_inq_file(adios_handle, variable_count, attribute_count, timestep_first, timestep_last, ier)
  call check_adios_err(ier,"Error inquiring adios file for reading: "//trim(filename))

  ! debug
  !print *,'debug: rank ',myrank_adios,' adios_inq_file: ',variable_count,attribute_count,timestep_first,timestep_last,ier

  ! user output
  print *,'file name: ',trim(filename)
  ! variables
  if (variable_count > 0) then
    ! user output
    print *,'variables: ',variable_count

    allocate (fnamelist(variable_count),stat=ier)
    if (ier /= 0) stop 'Error allocating namelist array'

    ! gets variable names
    call adios_inq_varnames(adios_handle, fnamelist, ier)
    if (ier /= 0) stop 'Error in adios_inq_varnames call'

    ! user output
    do i = 1,variable_count
      vname = trim(fnamelist(i)) ! variable name
      print *,'  id:',i-1,' ',trim(vname)
      dims(:) = 0
      ! vtype
      ! 0 byte
      ! 1 short
      ! 2 integer
      ! 4 long
      ! 5 real
      ! 6 double
      ! 7 long double
      ! 9 string
      ! 10 complex
      ! 11 double_complex
      ! 50 unsigned_byte
      ! 51 unsigned_short
      ! 52 unsigned_integer
      ! 54 unsigned_long
      call adios_inq_var (adios_handle, trim(vname), vtype, vnsteps, vndim, dims, ier)
      call check_adios_err(ier,"Error inquiring adios variable failed: "//trim(vname))
      if (vndim == 0) then
        ! scalar variable
        print *,'    scalar type/nsteps/ndim ',vtype,vnsteps,vndim
      else
        ! array variable
        print *,'    array  type/nsteps/ndim ',vtype,vnsteps,vndim,' dims ',dims(1:vndim)
      endif
    enddo
    deallocate(fnamelist)
  else
    print *,'  no variables'
  endif
  print *

  ! attributes
  if (attribute_count > 0) then
    ! user output
    print *,'attributes: ',attribute_count

    allocate (fnamelist(attribute_count),stat=ier)
    if (ier /= 0) stop 'Error allocating namelist array'

    ! gets attribute names
    call adios_inq_attrnames(adios_handle, fnamelist, ier)
    if (ier /= 0) stop 'Error in adios_inq_attrnames call'

    ! user output
    do i = 1,attribute_count
      print *,'  id:',i-1,' ',trim(fnamelist(i))
    enddo
    deallocate(fnamelist)
  else
    print *,'  no attributes'
  endif
  print *
  print *,'timesteps: first/last',timestep_first,'/',timestep_last

  ! groups
  call adios_inq_ngroups(adios_handle, group_count, ier)

  if (group_count > 0) then
    ! user output
    print *,'groups: ',group_count

    allocate (fnamelist(group_count),stat=ier)
    if (ier /= 0) stop 'Error allocating namelist array'

    ! gets group names
    call adios_inq_groupnames(adios_handle, fnamelist, ier)

    ! user output
    do i = 1,group_count
      print *,'  ',trim(fnamelist(i))
    enddo
    deallocate (fnamelist)
  else
    print *,'  no groups'
  endif
  print *

  ! to avoid compiler warning
  ier = adios_group

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! note: inquiry functionality in adios 2 quite different to version 1
  !       for Fortran bindings, unfortunately there is no adios2_inquire_all_*** functionality
  !       thus, we have no way to inquire and list all variables in a file using Fortran commands.
  !
  ! user output
  if (myrank_adios == 0) then
    print *
    print *,'show ADIOS2 file variables: ',trim(filename)
    print *
  endif

  ! only checks
  if (.not. myadios2_obj%valid) then
    call check_adios_err(ier,"Error inquiring adios2 file for reading: "//trim(filename))
    stop 'Error adios2 object not valid yet'
  endif
  if (.not. adios_handle%valid) stop 'Error adios2 invalid file handle in show_adios_file_variables()'
  if (.not. adios_group%valid) stop 'Error adios2 invalid group in show_adios_file_variables()'

  ! in case no io_group given, will need one for inquiry
  !if (.not. adios_group%valid) then
  !  call adios2_declare_io(adios_group, myadios2_obj, "Reader", ier)
  !  call check_adios_err(ier,"Error declaring an ADIOS2 IO group in show_adios_file_variables()")
  !endif

  ! user output
  if (myrank_adios == 0) then
    print *,'file name: ',trim(filename)

    ! inquire all variables
    ! Fortran binding not supported yet by adios2
    !call adios2_inquire_all_variables(vars,variable_count,adios_group)
    ! C-wrapper
    call get_adios2_all_variables_count(variable_count,adios_group)
    print *,'number of variables: ',variable_count
    ! variables
    if (variable_count > 0) then
      call show_adios2_all_variables(adios_group)
    else
      print *,'  no variables'
    endif
    print *

    ! inquire all attributes
    ! Fortran binding not supported yet by adios2
    !call adios2_inquire_all_attributes(atts,attribute_count,adios_group)
    ! C-wrapper
    call get_adios2_all_variables_count(attribute_count,adios_group)
    ! user output
    print *,'number of attributes: ',attribute_count
    ! attributes
    if (attribute_count > 0) then
      if (myrank_adios == 0) call show_adios2_all_variables(adios_group)
    else
      print *,'  no attributes'
    endif
    print *
  endif

#endif

  end subroutine show_adios_file_variables

#endif
!
!---------------------------------------------------------------------------------
!
#if defined(USE_ADIOS) || defined(USE_ADIOS2)
! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools

!> get the ADIOS filename
!! \param name file name (without extension)

  function get_adios_filename(name,ENGINE)

  use constants, only: MAX_STRING_LEN,ADIOS2_ENGINE_DEFAULT

  implicit none
  character(len=*), intent(in) :: name
  character(len=*), intent(in), optional :: ENGINE
  character(len=MAX_STRING_LEN) :: get_adios_filename
  ! local parameters
  character(len=MAX_STRING_LEN) :: tmp_engine
  character(len=4) :: ext
  integer :: irange,i

  ! selects engine type
  if (present(ENGINE)) then
    ! given by function call
    tmp_engine = trim(ENGINE)
  else
    ! default
    tmp_engine = trim(ADIOS2_ENGINE_DEFAULT)
  endif

  ! converts all string characters to lowercase (to make user input case-insensitive)
  irange = iachar('a') - iachar('A')
  do i = 1,len_trim(tmp_engine)
    if (lge(tmp_engine(i:i),'A') .and. lle(tmp_engine(i:i),'Z')) then
      tmp_engine(i:i) = achar(iachar(tmp_engine(i:i)) + irange)
    endif
  enddo

  ! debug
  !print *,'debug: get_adios_filename ',trim(name),' - engine: ',trim(tmp_engine)

  ! sets ending according to engine
  ! https://adios2.readthedocs.io/en/latest/engines/engines.html
  select case(trim(tmp_engine))
  case ("hdf5")
    ext = ".h5"
  case ("bp4")
    ext = ".bp"
  case default
    ext = ".bp"
  end select

  ! engine chooses the filename extension (e.g., **.bp, **.h5, ..)
  get_adios_filename = trim(name) // trim(ext)

  end function get_adios_filename

#endif


!-------------------------------------------------------------------------------
!
! error checking
!
!-------------------------------------------------------------------------------

#if defined(USE_ADIOS) || defined(USE_ADIOS2)
! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools

!> Get the ADIOS error message from an adios error number if there is an error.
!! \param ier The error code considered.
!! \param msg additional error message.

  subroutine check_adios_err(ier,msg)

  implicit none
  integer, intent(in) :: ier
  character(len=*), intent(in) :: msg

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! error message
  character(len=1024) :: err_message

  if (ier /= 0) then
    print *,'Error: process ',myrank_adios,' has ADIOS error ',ier
    print *,'Error message: ',trim(msg)
    call adios_errmsg(err_message)
    print *,'Error message ADIOS: ',trim(err_message)
    call exit_mpi(myrank_adios,'ADIOS error')
  endif

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  if (ier /= adios2_error_none) then
    print *,'Error: process ',myrank_adios,' has ADIOS2 error ',ier
    print *,'Error message: ',trim(msg)
    call exit_mpi(myrank_adios,'ADIOS2 error')
  endif
#endif

  end subroutine check_adios_err

#endif

end module manager_adios

