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

! read values from parameter file, ignoring white lines and comments

  subroutine read_value_integer(value_to_read, name, ier)

  use constants, only: MAX_STRING_LEN
  implicit none

  integer :: value_to_read
  character(len=*) :: name
  character(len=MAX_STRING_LEN) :: string_read
  integer :: ier

  call param_read(string_read, len(string_read), name, len(name), ier)
  if (ier /= 0) return
  read(string_read,*,iostat=ier) value_to_read

  end subroutine read_value_integer

!--------------------

  subroutine read_value_double_precision(value_to_read, name, ier)

  use constants, only: MAX_STRING_LEN
  implicit none

  double precision :: value_to_read
  character(len=*) :: name
  character(len=MAX_STRING_LEN) :: string_read
  integer :: ier

  call param_read(string_read, len(string_read), name, len(name), ier)
  if (ier /= 0) return
  read(string_read,*,iostat=ier) value_to_read

  end subroutine read_value_double_precision

!--------------------

  subroutine read_value_logical(value_to_read, name, ier)

  use constants, only: MAX_STRING_LEN
  implicit none

  logical :: value_to_read
  character(len=*) :: name
  character(len=MAX_STRING_LEN) :: string_read
  integer :: ier

  call param_read(string_read, len(string_read), name, len(name), ier)
  if (ier /= 0) return
  read(string_read,*,iostat=ier) value_to_read

  end subroutine read_value_logical

!--------------------

  subroutine read_value_string(value_to_read, name, ier)

  use constants, only: MAX_STRING_LEN
  implicit none

  character(len=*) :: value_to_read
  character(len=*) :: name
  character(len=MAX_STRING_LEN) :: string_read
  integer :: ier

  call param_read(string_read, len(string_read), name, len(name), ier)
  if (ier /= 0) return
  value_to_read = string_read

  end subroutine read_value_string

!--------------------

  subroutine open_parameter_file_from_main_only(ier)

  use constants, only: MAX_STRING_LEN,IN_DATA_FILES

  implicit none

  integer :: ier
  character(len=MAX_STRING_LEN) :: filename_main,filename_run0001
  logical :: exists_main_Par_file,exists_run0001_Par_file

  filename_main = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'Par_file'

! also see if we are running several independent runs in parallel
! to do so, add the right directory for that run for the main process only here
  filename_run0001 = 'run0001/'//filename_main(1:len_trim(filename_main))

  call param_open(filename_main, len(filename_main), ier)
  if (ier == 0) then
    exists_main_Par_file = .true.
    call close_parameter_file()
  else
    exists_main_Par_file    = .false.
  endif

  call param_open(filename_run0001, len(filename_run0001), ier)
  if (ier == 0) then
    exists_run0001_Par_file = .true.
    call close_parameter_file()
  else
    exists_run0001_Par_file = .false.
  endif

  if (exists_main_Par_file .and. exists_run0001_Par_file) then
    print *
    print *,'Cannot have both DATA/Par_file and run0001/DATA/Par_file present, please remove one of them.'
    print *
    print *,'In case you want to run simultaneous runs with: NUMBER_OF_SIMULTANEOUS_RUNS > 1 '
    print *,'make sure to remove DATA/Par_file and only have run0***/DATA/Par_file available'
    print *
    stop 'Error: two different copies of the Par_file'
  endif

  call param_open(filename_main, len(filename_main), ier)
  if (ier /= 0) then
    call param_open(filename_run0001, len(filename_run0001), ier)
    if (ier /= 0) then
      print *
      print *,'opening file failed, please check your file path and run-directory.'
      print *,'checked: ',trim(filename_main)
      print *,'         ',trim(filename_run0001)
      stop 'error opening Par_file'
    endif
  endif

  end subroutine open_parameter_file_from_main_only

!--------------------

  subroutine open_parameter_file(ier)

  use constants, only: MAX_STRING_LEN,IN_DATA_FILES,mygroup
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS

  implicit none

  integer :: ier
  character(len=MAX_STRING_LEN) :: filename,path_to_add

  filename = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'Par_file'
! see if we are running several independent runs in parallel
! if so, add the right directory for that run
! (group numbers start at zero, but directory names start at run0001, thus we add one)
! a negative value for "mygroup" is a convention that indicates that groups (i.e. sub-communicators, one per run) are off
  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
    write(path_to_add,"('run',i4.4,'/')") mygroup + 1
    filename = path_to_add(1:len_trim(path_to_add))//filename(1:len_trim(filename))
  endif

  call param_open(filename, len(filename), ier)
  if (ier /= 0) then
    print *
    print *,'opening file failed: ',trim(filename)
    print *,'Please check your file path and run-directory.'
    stop 'Error opening Par_file'
  endif

  end subroutine open_parameter_file

!--------------------

  subroutine close_parameter_file

  call param_close()

  end subroutine close_parameter_file
