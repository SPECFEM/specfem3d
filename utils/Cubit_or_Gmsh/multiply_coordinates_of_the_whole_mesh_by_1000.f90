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

! read an external mesh file and multiply its coordinates by 1000,
! for instance when it has been created in kilometers but the user wants it in meters

! Dimitri Komatitsch, CNRS, Marseille, France, June 2015.

  program multiply_coordinates_by_1000

  implicit none

!
! work in single or in double precision (4 or 8 bytes)
!
  integer, parameter :: CUSTOM_REAL = 4 ! 8

!------------------------------------------------------------------------------------------------

  integer, parameter :: NGNOD = 8                        ! number of control nodes for hexahedral elements (can only be 8 or 27)

  character(len=*), parameter :: nodes_coords_file     = 'MESH/nodes_coords_file'
  character(len=*), parameter :: nodes_coords_file_new = 'MESH/nodes_coords_file_new'

!------------------------------------------------------------------------------------------------

  integer :: NGLOB                    ! number of nodes

  integer :: i,iread,ier

  real(kind=CUSTOM_REAL) :: xtmp,ytmp,ztmp

  if (NGNOD /= 8) then
    print *,'error: multiply_coordinates_by_1000 only supports NGNOD == 8 for now'
    stop 'error in multiply_coordinates_by_1000'
  endif

! read the mesh
  print *
  print *,'start reading the existing node coordinate file: ',nodes_coords_file(1:len_trim(nodes_coords_file))
  print *,'and writing the new one multiplied by 1000: ',nodes_coords_file_new(1:len_trim(nodes_coords_file_new))

  open(unit=10,file=nodes_coords_file,status='old',action='read')
  open(unit=11,file=nodes_coords_file_new,status='unknown',action='write')

  read(10,*) NGLOB
  print *,'  number of points: ',NGLOB
  write(11,*) NGLOB

  do i = 1,NGLOB

    ! gets node ID and position
    read(10,*,iostat=ier) iread,xtmp,ytmp,ztmp

    ! check
    if (ier /= 0) then
      print *,'error point read:',i,iread,xtmp,ytmp,ztmp
      stop 'error while reading points'
    endif

    ! checks if out-of-range
    if (iread < 1 .or. iread > NGLOB) then
      print *,'error at i,iread = ',i,iread
      stop 'wrong ID input for a point'
    endif

    ! write the new values multiplied by 1000
    write(11,*) iread,xtmp*1000._CUSTOM_REAL,ytmp*1000._CUSTOM_REAL,ztmp*1000._CUSTOM_REAL

  enddo

  close(10)
  close(11)

  end program multiply_coordinates_by_1000

