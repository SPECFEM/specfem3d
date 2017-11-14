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

  subroutine create_name_database(prname,iproc,LOCAL_PATH)

! create the name of the database for the mesher and the solver

  use constants, only: MAX_STRING_LEN

  implicit none

  integer iproc

! name of the database file
  character(len=MAX_STRING_LEN) :: prname,procname,LOCAL_PATH

! create the name for the database of the current slide and region
  write(procname,"('/proc',i6.6,'_')") iproc

! create full name with path
  prname = LOCAL_PATH(1:len_trim(LOCAL_PATH)) // procname

  end subroutine create_name_database

