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
!> Initialize ASDF for reading the adjoint sources
subroutine asdf_setup(file_id)

  use iso_c_binding, only: C_NULL_CHAR
  implicit none

  ! asdf file handle
  integer,intent(inout) :: file_id

  ! local parameters
  character(len=512) :: filename
  integer :: comm
  integer :: ier

  call world_duplicate(comm)
  call ASDF_initialize_hdf5_f(ier)

  filename = 'SEM/adjoint.h5' // C_NULL_CHAR

  call ASDF_open_read_only_f(filename, comm, file_id)

end subroutine asdf_setup

!==============================================================================
!> Finalize ASDF. Called once all adjoint sources have been read from the file
subroutine asdf_cleanup()

  implicit none

  ! local parameters
  integer :: myrank
  integer :: ier

  call world_rank(myrank)
  call synchronize_all()

  call ASDF_finalize_hdf5_f(ier);

  if (ier /= 0 ) stop 'Error cleaning up ASDF: calling ASDF_finalize_hdf5_f() routine failed'

end subroutine asdf_cleanup
