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

  subroutine read_value_integer_mesh(iunit,ignore_junk,value_to_read, name, ierr)

  implicit none

  logical ignore_junk
  integer iunit
  integer value_to_read
  integer ierr
  character(len=*) name
  character(len=512) string_read

  call unused_string(name)
  ierr = 0

  call read_next_line(iunit,ignore_junk,string_read,ierr)
  if (ierr /= 0) return
  read(string_read,*,iostat=ierr) value_to_read

  end subroutine read_value_integer_mesh

!--------------------

  subroutine read_value_dble_precision_mesh(iunit,ignore_junk,value_to_read, name, ierr)

  implicit none

  logical ignore_junk
  integer iunit
  double precision value_to_read
  integer ierr
  character(len=*) name
  character(len=512) string_read

  call unused_string(name)
  ierr = 0

  call read_next_line(iunit,ignore_junk,string_read,ierr)
  if (ierr /= 0) return
  read(string_read,*,iostat=ierr) value_to_read

  end subroutine read_value_dble_precision_mesh

!--------------------

  subroutine read_value_logical_mesh(iunit,ignore_junk,value_to_read, name, ierr)

  implicit none

  logical ignore_junk
  logical value_to_read
  integer iunit
  integer ierr
  character(len=*) name
  character(len=512) string_read

  call unused_string(name)
  ierr = 0

  call read_next_line(iunit,ignore_junk,string_read,ierr)
  if (ierr /= 0) return
  read(string_read,*,iostat=ierr) value_to_read

  end subroutine read_value_logical_mesh

!--------------------

  subroutine read_value_string_mesh(iunit,ignore_junk,value_to_read, name, ierr)

  implicit none

  logical ignore_junk
  integer iunit
  character(len=*) value_to_read
  integer ierr
  character(len=*) name
  character(len=512) string_read

  call unused_string(name)
  ierr = 0

  call read_next_line(iunit,ignore_junk,string_read,ierr)
  if (ierr /= 0) return
  value_to_read = string_read

  end subroutine read_value_string_mesh

!--------------------

  subroutine read_interface_parameters(iunit,SUPPRESS_UTM_PROJECTION,interface_top_file, &
       npx_interface,npy_interface,&
       orig_x_interface,orig_y_interface,spacing_x_interface,spacing_y_interface,ierr)

  implicit none

  include "constants.h"

  logical SUPPRESS_UTM_PROJECTION
  integer iunit
  integer ierr
  integer npx_interface,npy_interface
  double precision orig_x_interface,orig_y_interface
  double precision spacing_x_interface,spacing_y_interface
  character(len=50)  interface_top_file
  character(len=512) string_read

  ierr = 0
  call read_next_line(iunit,DONT_IGNORE_JUNK,string_read,ierr)
  if (ierr /= 0) return
  read(string_read,*,iostat=ierr) SUPPRESS_UTM_PROJECTION,npx_interface,npy_interface,&
             orig_x_interface,orig_y_interface,spacing_x_interface,spacing_y_interface
  call read_value_string_mesh(iunit,DONT_IGNORE_JUNK,interface_top_file,'INTERFACE_TOP',ierr)
  end subroutine read_interface_parameters

!--------------------

  subroutine read_material_parameters(iunit,i,rho,vp,vs,Q_flag,anisotropy_flag,domain_id, ierr)

  implicit none

  include "constants.h"

  integer iunit
  integer i
  integer ierr
  double precision rho,vp,vs,Q_flag,anisotropy_flag,domain_id
  character(len=512) string_read

  ierr = 0
  call read_next_line(iunit,DONT_IGNORE_JUNK,string_read,ierr)
  if (ierr /= 0) return
  read(string_read,*,iostat=ierr) i,rho,vp,vs,Q_flag,anisotropy_flag,domain_id

  end subroutine read_material_parameters

!--------------------

  subroutine read_region_parameters(iunit,ix_beg_region,ix_end_region,iy_beg_region,iy_end_region,&
          iz_beg_region,iz_end_region,imaterial_number, ierr)

  implicit none

  include "constants.h"

  integer iunit
  integer ierr
  integer ix_beg_region,ix_end_region,iy_beg_region,iy_end_region
  integer iz_beg_region,iz_end_region,imaterial_number
  character(len=512) string_read

  ierr = 0
  call read_next_line(iunit,DONT_IGNORE_JUNK,string_read,ierr)
  if (ierr /= 0) return
  read(string_read,*,iostat=ierr) ix_beg_region,ix_end_region,iy_beg_region,iy_end_region,&
          iz_beg_region,iz_end_region,imaterial_number

  end subroutine read_region_parameters

!--------------------

  subroutine read_next_line(iunit,suppress_junk,string_read, ierr)

  implicit none

  include "constants.h"

  integer ierr
  logical suppress_junk
  character(len=512) string_read
  integer index_equal_sign,ios,iunit

  ierr = 0
  do
    read(unit=iunit,fmt="(a512)",iostat=ios) string_read
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
     if(index_equal_sign <= 1 .or. index_equal_sign == len_trim(string_read)) stop 'incorrect syntax detected in Mesh_Par_file'
     string_read = string_read(index_equal_sign + 1:len_trim(string_read))
  endif

! suppress leading and trailing white spaces again, if any, after having suppressed the leading junk
  string_read = adjustl(string_read)
  string_read = string_read(1:len_trim(string_read))

  end subroutine read_next_line

!--------------------

  subroutine open_parameter_file_mesh

  implicit none
  include "constants.h"

  integer ierr

  open(unit=IIN,file=MF_IN_DATA_FILES_PATH(1:len_trim(MF_IN_DATA_FILES_PATH)) &
       //'Mesh_Par_file',status='old',action='read',iostat=ierr)
  if (ierr /= 0) then
    print*,'error opening file: ',MF_IN_DATA_FILES_PATH(1:len_trim(MF_IN_DATA_FILES_PATH))
    print*
    print*,'please check your file path and run-directory.'
    stop 'error opening Mesh_Par_file'
  endif

  end subroutine open_parameter_file_mesh

!--------------------

  subroutine close_parameter_file_mesh

  include "constants.h"

  close(IIN)

  end subroutine close_parameter_file_mesh

!--------------------

! dummy subroutine to avoid warnings about variable not used in other subroutines
  subroutine unused_string(s)

  character(len=*) s

  if (len(s) == 1) continue

  end subroutine unused_string

