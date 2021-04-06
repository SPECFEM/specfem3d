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

  ! MPI copies of communicator and rank
  integer :: comm_adios
  integer,public :: myrank_adios
  integer,public :: sizeprocs_adios

contains

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

!==============================================================================
!> Initialize ADIOS and setup the xml output file
  subroutine adios_setup()

  call no_adios_err()

  end subroutine adios_setup

!==============================================================================
!> Finalize ADIOS. Must be called once everything is written down.
  subroutine adios_cleanup()

  call no_adios_err()

  end subroutine adios_cleanup

end module adios_manager_mod
