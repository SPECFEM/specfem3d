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

  program convert_CUBIT_SPECFEM_to_DX

! Dimitri Komatitsch, CNRS Marseille, France, April 2013

! convert CUBIT files that are in SPECFEM3D_Cartesian format to OpenDX format for visualization

  implicit none

! this is for HEX8; the code below probably works for HEX27 as well, it will just display the first 8 points
! i.e. it will display HEX27 elements as if they were HEX8
  integer, parameter :: NGNOD = 8

  integer :: nspec,npoin
  integer :: ispec,ipoin
  integer :: ipoin_read,ispec_read,imat_read
  integer :: i1,i2,i3,i4,i5,i6,i7,i8

  real, dimension(:), allocatable :: x,y,z

  integer, dimension(:), allocatable :: imat

  integer, dimension(:,:), allocatable :: ibool

  real :: xread,yread,zread
  real :: val_color

  print *,'using material file to color the mesh elements'
  print *

! open SPECFEM3D_Cartesian mesh file to read the points
    open(unit=23,file='nodes_coords_file',status='old',action='read')
    read(23,*) npoin
    allocate(x(npoin))
    allocate(y(npoin))
    allocate(z(npoin))
    do ipoin = 1,npoin
      read(23,*) ipoin_read,xread,yread,zread
      x(ipoin_read) = xread
      y(ipoin_read) = yread
      z(ipoin_read) = zread
    enddo
    close(23)

! write DX header with element data
    open(unit=11,file='DX_fullmesh.dx',status='unknown',action='write')
    write(11,*) 'object 1 class array type float rank 1 shape 3 items ',npoin,' data follows'

! read local points in this slice and output global DX points
  do ipoin=1,npoin
    write(11,*) x(ipoin),y(ipoin),z(ipoin)
  enddo

! ************* generate elements ******************

! open SPECFEM3D_Cartesian topology file to read the mesh elements
  open(unit=23,file='mesh_file',status='old',action='read')
  read(23,*) nspec

  allocate(imat(nspec))
  allocate(ibool(NGNOD,nspec))

  write(11,*) 'object 2 class array type int rank 1 shape 8 items ',nspec,' data follows'

! read local elements in this slice and output global DX elements
  do ispec=1,nspec
    read(23,*) ispec_read,i1,i2,i3,i4,i5,i6,i7,i8
    ibool(1,ispec_read) = i1
    ibool(2,ispec_read) = i2
    ibool(3,ispec_read) = i3
    ibool(4,ispec_read) = i4
    ibool(5,ispec_read) = i5
    ibool(6,ispec_read) = i6
    ibool(7,ispec_read) = i7
    ibool(8,ispec_read) = i8
  enddo

! point order in OpenDX is 4,1,8,5,3,2,7,6, *not* 1,2,3,4,5,6,7,8 as in AVS or SPECFEM
! and point numbers start at 0 rather than 1
  do ispec=1,nspec
    write(11,"(i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9)") &
                 ibool(4,ispec)-1,ibool(1,ispec)-1,ibool(8,ispec)-1,ibool(5,ispec)-1, &
                 ibool(3,ispec)-1,ibool(2,ispec)-1,ibool(7,ispec)-1,ibool(6,ispec)-1
  enddo

  close(23)

! ************* generate element data values ******************

! output DX header for data
    write(11,*) 'attribute "element type" string "cubes"'
    write(11,*) 'attribute "ref" string "positions"'
    write(11,*) 'object 3 class array type float rank 0 items ',nspec,' data follows'

! read local elements in this slice and output global DX elements
    open(unit=23,file='materials_file',status='old',action='read')
    do ispec=1,nspec
! beware: elements may not be listed in increasing order, they can appear in any order
      read(23,*) ispec_read,imat_read
      imat(ispec_read) = imat_read
    enddo
    close(23)

  do ispec=1,nspec
    val_color = imat(ispec) ! use material property read to color the elements
    write(11,*) val_color
  enddo

! define OpenDX field
    write(11,*) 'attribute "dep" string "connections"'
    write(11,*) 'object "irregular positions irregular connections" class field'
    write(11,*) 'component "positions" value 1'
    write(11,*) 'component "connections" value 2'
    write(11,*) 'component "data" value 3'
    write(11,*) 'end'

  close(11)

  end program convert_CUBIT_SPECFEM_to_DX

