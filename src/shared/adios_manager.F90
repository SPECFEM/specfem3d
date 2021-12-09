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
!> Tools to setup and cleanup ADIOS
!------------------------------------------------------------------------------
module adios_manager_mod

  implicit none

  private

  ! MPI copies of communicator and rank
  integer,public :: comm_adios
  integer :: myrank_adios
  integer :: sizeprocs_adios

#if defined(USE_ADIOS)
  ! adios read
  character(len=*),parameter,public :: ADIOS_VERBOSITY = "verbose=1" ! lowest level: verbose=0,error_only=1, .., debug=4
#endif

  ! public accessibility
  public :: adios_setup
  public :: adios_cleanup
  public :: no_adios_err

contains

!==============================================================================
!> Initialize ADIOS and setup the xml output file
  subroutine adios_setup()

#if defined(USE_ADIOS)
  use constants, only: ADIOS_BUFFER_SIZE_IN_MB
#ifdef ADIOS_VERSION_OLD
! ADIOS versions <= 1.9
! adios_set_max_buffer_size not defined yet
#else
! ADIOS versions >= 1.10
  use adios_write_mod, only: adios_set_max_buffer_size
#endif
#endif

  implicit none

  ! local parameters
#if defined(USE_ADIOS)
  integer :: ier
#endif

  ! initializes
  comm_adios = 0
  myrank_adios = -1
  sizeprocs_adios = 0

#if defined(USE_ADIOS)
  ! ADIOS
  ! gets MPI communicator for adios calls
  call world_duplicate(comm_adios)

  ! gets rank from (duplicate) adios communicator
  call world_rank_comm(myrank_adios,comm_adios)

  ! number of MPI processes
  call world_size_comm(sizeprocs_adios,comm_adios)

  ! checks
  if (sizeprocs_adios == 0) &
    stop 'Error adios initialization got zero processes'

  call adios_init_noxml (comm_adios, ier)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  !       e.g., version 1.5.0 returns 1 here
  !print *,'adios init return: ',adios_err
  !if (adios_err /= 0) stop 'Error setting up ADIOS: calling adios_init_noxml() routine failed'

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

#elif defined(USE_ADIOS2)
  ! ADIOS2
  ! not supported yet
  ! gets rank
  call world_rank(myrank_adios)

  ! error info
  if (myrank_adios == 0) then
    print *, "Error: ADIOS2 is not fully supported yet. Please use with ADIOS version 1.x instead."
    print *, "To enable ADIOS support, reconfigure with --with-adios flag."
  endif
  ! safety stop
  call exit_MPI(myrank_adios,"Error ADIOS manager: ADIOS2 not supported yet")

#else
  ! no adios compilation support
  ! gets rank
  call world_rank(myrank_adios)

  ! compilation without ADIOS support
  if (myrank_adios == 0) then
    print *, "Error: ADIOS enabled without ADIOS Support."
    print *, "To enable ADIOS support, reconfigure with --with-adios flag."
  endif
  ! safety stop
  call exit_MPI(myrank_adios,"Error ADIOS manager: intitialize called without compilation support")

#endif

  end subroutine adios_setup

!==============================================================================
!> Finalize ADIOS. Must be called once everything is written down.
  subroutine adios_cleanup()

#if defined(USE_ADIOS)
  use adios_write_mod, only: adios_finalize
#endif

  implicit none

  ! local parameters
  integer :: myrank
#if defined(USE_ADIOS)
  integer :: ier
#endif

  ! synchronizes main local communicator
  call world_rank(myrank)
  call synchronize_all()

  ! synchronizes from comm_adios communicator
  call synchronize_all_comm(comm_adios)

#if defined(USE_ADIOS)
  ! finalize
  call adios_finalize (myrank_adios, ier)
  if (ier /= 0) stop 'Error cleaning up ADIOS: calling adios_finalize() routine failed'
#endif

  end subroutine adios_cleanup

!==============================================================================
!> safety routine to stop execution with a warning to configure w/ ADIOS.
  subroutine no_adios_err()

  implicit none

  integer :: myrank

  call world_rank(myrank)
  if (myrank == 0) then
    print *, "----------------------------------------------------"
    print *, "Not configured to compile or run with ADIOS."
    print *, "Check your Par_file and set ADIOS_ENABLED to .false."
    print *, "or reconfigure using --with-adios."
    print *, "----------------------------------------------------"
  endif
  call abort_mpi()

  end subroutine

end module adios_manager_mod
