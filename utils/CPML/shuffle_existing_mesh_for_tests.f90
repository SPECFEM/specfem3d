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

  program shuffle_existing_mesh_for_tests

! Dimitri Komatitsch, CNRS Marseille, France, March 2016

! rotate some of the elements of an existing mesh in order to fully test the CPML extrusion program I wrote
! to make sure all possible element orientations are encountered and thus tested

  implicit none

! this code is for HEX8 only

  integer :: nspec,ispec,ispec_loop,i1,i2,i3,i4,i5,i6,i7,i8

! open SPECFEM3D_Cartesian topology file to read the mesh elements
  open(unit=23,file='mesh_file',status='old',action='read')
  open(unit=24,file='mesh_file_shuffled',status='unknown',action='write')

  read(23,*) nspec
  write(24,*) nspec

! loop on the whole mesh
  do ispec_loop = 1,nspec

    read(23,*) ispec,i1,i2,i3,i4,i5,i6,i7,i8

! implement shuffling of the list of points, one possibility for each of the six faces of the reference cube
    if (mod(ispec_loop,6) == 0) then
      write(24,"(i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9)") ispec,i1,i2,i3,i4,i5,i6,i7,i8
    else if (mod(ispec_loop,6) == 1) then
      write(24,"(i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9)") ispec,i5,i1,i4,i8,i6,i2,i3,i7
    else if (mod(ispec_loop,6) == 2) then
      write(24,"(i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9)") ispec,i8,i7,i6,i5,i4,i3,i2,i1
    else if (mod(ispec_loop,6) == 3) then
      write(24,"(i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9)") ispec,i3,i2,i6,i7,i4,i1,i5,i8
    else if (mod(ispec_loop,6) == 4) then
      write(24,"(i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9)") ispec,i4,i3,i7,i8,i1,i2,i6,i5
    else if (mod(ispec_loop,6) == 5) then
      write(24,"(i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9)") ispec,i1,i5,i6,i2,i4,i8,i7,i3
    else
      stop 'incorrect value of the shuffling index'
    endif

  enddo

  close(23)
  close(24)

  end program shuffle_existing_mesh_for_tests

