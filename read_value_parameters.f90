!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 2
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology July 2004
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

! read values from parameter file, ignoring white lines and comments

  subroutine read_value_integer(value_to_read)

  implicit none

  integer value_to_read
  character(len=100) string_read

  call read_next_line(string_read)
  read(string_read,*) value_to_read

  end subroutine read_value_integer

!--------------------

  subroutine read_value_double_precision(value_to_read)

  implicit none

  double precision value_to_read
  character(len=100) string_read

  call read_next_line(string_read)
  read(string_read,*) value_to_read

  end subroutine read_value_double_precision

!--------------------

  subroutine read_value_logical(value_to_read)

  implicit none

  logical value_to_read
  character(len=100) string_read

  call read_next_line(string_read)
  read(string_read,*) value_to_read

  end subroutine read_value_logical

!--------------------

  subroutine read_value_string(value_to_read)

  implicit none

  character(len=*) value_to_read
  character(len=100) string_read

  call read_next_line(string_read)
  read(string_read,"(a)") value_to_read

  end subroutine read_value_string

!--------------------

  subroutine read_next_line(string_read)

  implicit none

  include "constants.h"

  character(len=100) string_read

  integer ios

  do
    read(unit=IIN,fmt=200,iostat=ios) string_read
    if(ios /= 0) stop 'error while reading input file'

! suppress leading white spaces, if any
    string_read = adjustl(string_read)

! exit loop when we find the first line that is not a comment or a white line
    if(len_trim(string_read) == 0) cycle
    if(string_read(1:1) /= '#') exit

  enddo

! suppress trailing white spaces, if any
  string_read = string_read(1:len_trim(string_read))

! suppress trailing comments, if any
  if(index(string_read,'#') > 0) string_read = string_read(1:index(string_read,'#')-1)

! suppress leading junk (first 34 characters)
  string_read = string_read(35:len_trim(string_read))

! format
 200 format(a100)

  end subroutine read_next_line

