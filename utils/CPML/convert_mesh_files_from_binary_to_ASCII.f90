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

  program convert_mesh_files_binary_ASCII

! Dimitri Komatitsch, CNRS Marseille, France, March 2016

! convert mesh files back from binary to ASCII

  implicit none

! this code is for HEX8 only for now
  integer, parameter :: NGNOD = 8

  integer :: nspec,npoin
  integer :: ispec,ipoin

  double precision, dimension(:), allocatable, target :: x,y,z

  integer, dimension(:), allocatable :: imaterial

  integer, dimension(:,:), allocatable :: ibool

  print *,'Converting mesh files back from binary to ASCII...'
  print *

! read the new points in binary format
  open(unit=23,file='nodes_coords_file.bin',form='unformatted',status='old',action='read')
  read(23) npoin
  allocate(x(npoin))
  allocate(y(npoin))
  allocate(z(npoin))
  read(23) x
  read(23) y
  read(23) z
  close(23)

! read the new mesh elements in binary format
  open(unit=23,file='mesh_file.bin',form='unformatted',status='old',action='read')
  read(23) nspec
  allocate(ibool(NGNOD,nspec))
  read(23) ibool
  close(23)

! read the new material properties in binary format
  open(unit=23,file='materials_file.bin',form='unformatted',status='old',action='read')
  allocate(imaterial(nspec))
  read(23) imaterial
  close(23)

  print *,'Total number of elements in the mesh read = ',nspec
  print *,'Total number of points in the mesh read = ',npoin

! ************* write mesh points and elements *************

! open SPECFEM3D_Cartesian mesh file to write the points
  open(unit=23,file='nodes_coords_file',status='unknown',action='write')
  write(23,*) npoin
  do ipoin = 1,npoin
    write(23,*) ipoin,sngl(x(ipoin)),sngl(y(ipoin)),sngl(z(ipoin))
  enddo
  close(23)

! open SPECFEM3D_Cartesian topology file to write the mesh elements
  open(unit=23,file='mesh_file',status='unknown',action='write')
  write(23,*) nspec
! loop on the whole mesh
  do ispec = 1,nspec
    write(23,"(i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9)") ispec, &
                      ibool(1,ispec),ibool(2,ispec),ibool(3,ispec),ibool(4,ispec), &
                      ibool(5,ispec),ibool(6,ispec),ibool(7,ispec),ibool(8,ispec)
  enddo
  close(23)

! write the materials file
  open(unit=23,file='materials_file',status='unknown',action='write')
! loop on the whole mesh
  do ispec = 1,nspec
    write(23,*) ispec,imaterial(ispec)
  enddo
  close(23)

  end program convert_mesh_files_binary_ASCII

