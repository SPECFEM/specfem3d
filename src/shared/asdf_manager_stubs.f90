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

!==============================================================================
!> Tools to setup and cleanup ASDF
!------------------------------------------------------------------------------
module asdf_manager_mod

contains

subroutine no_asdf_err()
  implicit none

  integer :: myrank

  call world_rank(myrank)
  if (myrank == 0) then
    print *, "----------------------------------------------------"
    print *, "Not configured to compile or run with ASDF."
    print *, "Check your Par_file and set ASDF_FORMAT to .false."
    print *, "or reconfigure using --with-asdf."
    print *, "----------------------------------------------------"
  endif
  call abort_mpi()
end subroutine

!==============================================================================
!> Initialize ASDF
subroutine asdf_setup()
  call no_asdf_err()
end subroutine asdf_setup

!==============================================================================
!> Finalize ASDF
subroutine asdf_cleanup()
  call no_asdf_err()
end subroutine asdf_cleanup

end module asdf_manager_mod

