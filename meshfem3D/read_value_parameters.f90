!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  1 . 4
!               ---------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
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

  subroutine read_value_integer(ignore_junk,value_to_read, name)

  implicit none

  logical ignore_junk
  integer value_to_read
  character(len=*) name
  character(len=100) string_read

  call unused_string(name)

  call read_next_line(ignore_junk,string_read)
  read(string_read,*) value_to_read

  end subroutine read_value_integer

!--------------------

  subroutine read_value_double_precision(ignore_junk,value_to_read, name)

  implicit none

  logical ignore_junk
  double precision value_to_read
  character(len=*) name
  character(len=100) string_read

  call unused_string(name)

  call read_next_line(ignore_junk,string_read)
  read(string_read,*) value_to_read

  end subroutine read_value_double_precision

!--------------------

  subroutine read_value_logical(ignore_junk,value_to_read, name)

  implicit none

  logical ignore_junk
  logical value_to_read
  character(len=*) name
  character(len=100) string_read

  call unused_string(name)

  call read_next_line(ignore_junk,string_read)
  read(string_read,*) value_to_read

  end subroutine read_value_logical

!--------------------

  subroutine read_value_string(ignore_junk,value_to_read, name)

  implicit none

  logical ignore_junk
  character(len=*) value_to_read
  character(len=*) name
  character(len=100) string_read

  call unused_string(name)

  call read_next_line(ignore_junk,string_read)
  value_to_read = string_read

  end subroutine read_value_string

!--------------------

  subroutine read_interface_parameters(SUPPRESS_UTM_PROJECTION,npx_interface,npy_interface,&
             orig_x_interface,orig_y_interface,spacing_x_interface,spacing_y_interface)

  implicit none

  include "constants.h"

  logical SUPPRESS_UTM_PROJECTION
  integer npx_interface,npy_interface
  double precision orig_x_interface,orig_y_interface
  double precision spacing_x_interface,spacing_y_interface
  character(len=100) string_read

  call read_next_line(DONT_IGNORE_JUNK,string_read)
  read(string_read,*) SUPPRESS_UTM_PROJECTION,npx_interface,npy_interface,&
             orig_x_interface,orig_y_interface,spacing_x_interface,spacing_y_interface 

  end subroutine read_interface_parameters

!--------------------

  subroutine read_material_parameters(i,rho,vp,vs,Q_flag,anisotropy_flag,domain_id)

  implicit none

  include "constants.h"

  integer i
  double precision rho,vp,vs,Q_flag,anisotropy_flag,domain_id
  character(len=100) string_read

  call read_next_line(DONT_IGNORE_JUNK,string_read)
  read(string_read,*)  i,rho,vp,vs,Q_flag,anisotropy_flag,domain_id

  end subroutine read_material_parameters

!--------------------

  subroutine read_region_parameters(ix_beg_region,ix_end_region,iy_beg_region,iy_end_region,&
          iz_beg_region,iz_end_region,imaterial_number)

  implicit none

  include "constants.h"

  integer ix_beg_region,ix_end_region,iy_beg_region,iy_end_region
  integer iz_beg_region,iz_end_region,imaterial_number
  character(len=100) string_read

  call read_next_line(DONT_IGNORE_JUNK,string_read)
  read(string_read,*) ix_beg_region,ix_end_region,iy_beg_region,iy_end_region,&
          iz_beg_region,iz_end_region,imaterial_number

  end subroutine read_region_parameters
  
!--------------------

  subroutine read_next_line(suppress_junk,string_read)

  implicit none

  include "constants.h"


  logical suppress_junk
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

  if(suppress_junk) then
! suppress leading junk (up to the first equal sign, included)
     index_equal_sign = index(string_read,'=')
     if(index_equal_sign <= 1 .or. index_equal_sign == len_trim(string_read)) stop 'incorrect syntax detected in DATA/Par_file'
     string_read = string_read(index_equal_sign + 1:len_trim(string_read))
  end if

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

