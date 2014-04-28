!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
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

!==============================================================================
!> Tools to setup and cleanup ADIOS
!------------------------------------------------------------------------------
module adios_manager_mod

contains

!==============================================================================
!> Initialize ADIOS and setup the xml output file
subroutine adios_setup()

  use adios_write_mod, only: adios_init

  implicit none

  include 'constants.h'

  integer :: adios_err
  integer :: comm

  call world_get_comm(comm)

  call adios_init_noxml (comm, adios_err);
  call adios_allocate_buffer (ADIOS_BUFFER_SIZE_IN_MB, adios_err)

end subroutine adios_setup

!==============================================================================
!> Finalize ADIOS. Must be called once everything is written down.
subroutine adios_cleanup()

  use adios_write_mod, only: adios_finalize

  implicit none
  integer :: myrank
  integer :: adios_err

  call world_rank(myrank)
  call synchronize_all()

  call adios_finalize (myrank, adios_err)

end subroutine adios_cleanup

end module adios_manager_mod
