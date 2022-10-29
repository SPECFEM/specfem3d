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

  program convert_mesh_files_ASCII_binary

! Dimitri Komatitsch, CNRS Marseille, France, March 2016

! convert mesh files from ASCII to binary to speed up future I/Os

  implicit none

! this code is for HEX8 only for now
  integer, parameter :: NGNOD = 8

  integer :: nspec,npoin
  integer :: ispec,ipoin
  integer :: ipoin_read,ispec_loop
  integer :: i1,i2,i3,i4,i5,i6,i7,i8

  double precision :: xread,yread,zread

  double precision, dimension(:), allocatable, target :: x,y,z

  integer, dimension(:), allocatable :: imaterial

  integer, dimension(:,:), allocatable :: ibool

  print *,'Converting mesh files from ASCII to binary...'
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

  print *,'Total number of points in the mesh read = ',npoin

! ************* read mesh elements *************

! open SPECFEM3D_Cartesian topology file to read the mesh elements
  open(unit=23,file='mesh_file',status='old',action='read')
  read(23,*) nspec
  allocate(ibool(NGNOD,nspec))

! loop on the whole mesh
  do ispec_loop = 1,nspec

    read(23,*) ispec,i1,i2,i3,i4,i5,i6,i7,i8

! store the ibool() array read
    ibool(1,ispec) = i1
    ibool(2,ispec) = i2
    ibool(3,ispec) = i3
    ibool(4,ispec) = i4
    ibool(5,ispec) = i5
    ibool(6,ispec) = i6
    ibool(7,ispec) = i7
    ibool(8,ispec) = i8

  enddo

  close(23)

  print *,'Total number of elements in the mesh read = ',nspec

! read the materials file
  allocate(imaterial(nspec))
  open(unit=23,file='materials_file',status='old',action='read')
! loop on the whole mesh
  do ispec_loop = 1,nspec
    read(23,*) ispec,i1
! store the imaterial() array read
    imaterial(ispec) = i1
  enddo
  close(23)

! write the new points in binary format
  open(unit=23,file='nodes_coords_file.bin',form='unformatted',status='unknown',action='write')
  write(23) npoin
  write(23) x
  write(23) y
  write(23) z
  close(23)

! write the new mesh elements in binary format
  open(unit=23,file='mesh_file.bin',form='unformatted',status='unknown',action='write')
  write(23) nspec
  write(23) ibool
  close(23)

! write the new material properties in binary format
  open(unit=23,file='materials_file.bin',form='unformatted',status='unknown',action='write')
  write(23) imaterial
  close(23)

  end program convert_mesh_files_ASCII_binary

