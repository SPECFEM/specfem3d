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

! read values from parameter file, ignoring white lines and comments

  subroutine read_value_integer(value_to_read, name, ierr)

  use constants, only: MAX_STRING_LEN
  implicit none

  integer value_to_read
  character(len=*) :: name
  character(len=MAX_STRING_LEN) :: string_read
  integer ierr

  call param_read(string_read, len(string_read), name, len(name), ierr)
  if (ierr /= 0) return
  read(string_read,*,iostat=ierr) value_to_read

  end subroutine read_value_integer

!--------------------

  subroutine read_value_double_precision(value_to_read, name, ierr)

  use constants, only: MAX_STRING_LEN
  implicit none

  double precision value_to_read
  character(len=*) :: name
  character(len=MAX_STRING_LEN) :: string_read
  integer ierr

  call param_read(string_read, len(string_read), name, len(name), ierr)
  if (ierr /= 0) return
  read(string_read,*,iostat=ierr) value_to_read

  end subroutine read_value_double_precision

!--------------------

  subroutine read_value_logical(value_to_read, name, ierr)

  use constants, only: MAX_STRING_LEN
  implicit none

  logical value_to_read
  character(len=*) :: name
  character(len=MAX_STRING_LEN) :: string_read
  integer ierr

  call param_read(string_read, len(string_read), name, len(name), ierr)
  if (ierr /= 0) return
  read(string_read,*,iostat=ierr) value_to_read

  end subroutine read_value_logical

!--------------------

  subroutine read_value_string(value_to_read, name, ierr)

  use constants, only: MAX_STRING_LEN
  implicit none

  character(len=*) :: value_to_read
  character(len=*) :: name
  character(len=MAX_STRING_LEN) :: string_read
  integer ierr

  call param_read(string_read, len(string_read), name, len(name), ierr)
  if (ierr /= 0) return
  value_to_read = string_read

  end subroutine read_value_string

!--------------------

  subroutine open_parameter_file(ierr)

  use constants, only: MAX_STRING_LEN,IN_DATA_FILES_PATH

  implicit none

  integer ierr
  character(len=MAX_STRING_LEN) :: filename
  filename = IN_DATA_FILES_PATH(1:len_trim(IN_DATA_FILES_PATH))//'Par_file'

  call param_open(filename, len(filename), ierr)
  if (ierr /= 0) then
    print*
    print*,'opening file failed, please check your file path and run-directory.'
    stop 'error opening Par_file'
  endif

  end subroutine open_parameter_file

!--------------------

  subroutine close_parameter_file

  call param_close()

  end subroutine close_parameter_file
