!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 4
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

! read values from parameter file, ignoring white lines and comments

  subroutine read_value_integer(value_to_read, name)

  implicit none

  integer value_to_read
  character(len=*) name
  character(len=100) string_read

  call unused_string(name)

  call read_next_line(string_read)
  read(string_read,*) value_to_read

  end subroutine read_value_integer

!--------------------

  subroutine read_value_double_precision(value_to_read, name)

  implicit none

  double precision value_to_read
  character(len=*) name
  character(len=100) string_read

  call unused_string(name)

  call read_next_line(string_read)
  read(string_read,*) value_to_read

  end subroutine read_value_double_precision

!--------------------

  subroutine read_value_logical(value_to_read, name)

  implicit none

  logical value_to_read
  character(len=*) name
  character(len=100) string_read

  call unused_string(name)

  call read_next_line(string_read)
  read(string_read,*) value_to_read

  end subroutine read_value_logical

!--------------------

  subroutine read_value_string(value_to_read, name)

  implicit none

  character(len=*) value_to_read
  character(len=*) name
  character(len=100) string_read

  call unused_string(name)

  call read_next_line(string_read)
  value_to_read = string_read

  end subroutine read_value_string

!--------------------

  subroutine read_next_line(string_read)

  implicit none

  include "constants.h"

  character(len=100) string_read

  integer index_equal_sign,ios

  do
    read(unit=IIN,fmt="(a100)",iostat=ios) string_read
    if(ios /= 0) stop 'error while reading parameter file'

! suppress leading white spaces, if any
    string_read = adjustl(string_read)

! suppress trailing carriage return (ASCII code 13) if any (e.g. if input text file coming from Windows/DOS)
    if(index(string_read,achar(13)) > 0) string_read = string_read(1:index(string_read,achar(13))-1)

! exit loop when we find the first line that is not a comment or a white line
    if(len_trim(string_read) == 0) cycle
    if(string_read(1:1) /= '#') exit

  enddo

! suppress trailing white spaces, if any
  string_read = string_read(1:len_trim(string_read))

! suppress trailing comments, if any
  if(index(string_read,'#') > 0) string_read = string_read(1:index(string_read,'#')-1)

! suppress leading junk (up to the first equal sign, included)
  index_equal_sign = index(string_read,'=')
  if(index_equal_sign <= 1 .or. index_equal_sign == len_trim(string_read)) stop 'incorrect syntax detected in DATA/Par_file'
  string_read = string_read(index_equal_sign + 1:len_trim(string_read))

! suppress leading and trailing white spaces again, if any, after having suppressed the leading junk
  string_read = adjustl(string_read)
  string_read = string_read(1:len_trim(string_read))

  end subroutine read_next_line

!--------------------

  subroutine open_parameter_file

  include "constants.h"

  open(unit=IIN,file='DATA/Par_file',status='old',action='read')

  end subroutine open_parameter_file

!--------------------

  subroutine close_parameter_file

  include "constants.h"

  close(IIN)

  end subroutine close_parameter_file

!--------------------

  integer function err_occurred()

  err_occurred = 0

  end function err_occurred

!--------------------

! dummy subroutine to avoid warnings about variable not used in other subroutines
  subroutine unused_string(s)

  character(len=*) s

  if (len(s) == 1) continue

  end subroutine unused_string

