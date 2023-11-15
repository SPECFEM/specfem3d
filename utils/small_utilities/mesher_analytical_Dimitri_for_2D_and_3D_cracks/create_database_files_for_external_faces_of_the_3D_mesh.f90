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

  program create_database_files_for_external_faces_of_the_3D_mesh

! Dimitri Komatitsch, CNRS Marseille, France, April 2013 and February 2017 and September 2018

  implicit none

  integer :: nspec,npoin
  integer :: ispec,ipoin
  integer :: ipoin_read,ispec_loop,iformat,NGNOD,ia
  integer :: i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,i21,i22,i23,i24,i25,i26,i27
  integer :: count_faces_found

  double precision, dimension(:), allocatable :: x,y,z

  double precision :: xread,yread,zread,xmin,xmax,ymin,ymax,zmin,zmax,limit,size_of_model

  logical :: already_found_a_face

  integer, dimension(:,:), allocatable :: ibool

! to make sure coordinate roundoff problems do not occur, use a tolerance of 0.5%
  double precision, parameter :: SMALL_PERCENTAGE_TOLERANCE = 1.005d0

  double precision, parameter :: SMALL_RELATIVE_VALUE = 0.5d-3

! read the mesh files in ASCII format (that is the standard case)
  iformat = 1

! the mesh contains HEX8 elements
  NGNOD = 8

! open SPECFEM3D_Cartesian mesh file to read the points
  if (iformat == 1) then
    open(unit=23,file='DATA_3D/nodes_coords_file',status='old',action='read')
    read(23,*) npoin
  else
    open(unit=23,file='DATA_3D/nodes_coords_file.bin',form='unformatted',status='old',action='read')
    read(23) npoin
  endif
  allocate(x(npoin))
  allocate(y(npoin))
  allocate(z(npoin))
  if (iformat == 1) then
    do ipoin = 1,npoin
      read(23,*) ipoin_read,xread,yread,zread
      x(ipoin_read) = xread
      y(ipoin_read) = yread
      z(ipoin_read) = zread
    enddo
  else
    read(23) x
    read(23) y
    read(23) z
  endif
  close(23)

! compute the min and max values of each coordinate
   xmin = minval(x)
   xmax = maxval(x)

   ymin = minval(y)
   ymax = maxval(y)

   zmin = minval(z)
   zmax = maxval(z)

  print *,'Xmin and Xmax of the mesh read = ',xmin,xmax
  print *,'Ymin and Ymax of the mesh read = ',ymin,ymax
  print *,'Zmin and Zmax of the mesh read = ',zmin,zmax
  print *

! ************* read mesh elements *************

! open SPECFEM3D_Cartesian topology file to read the mesh elements
  if (iformat == 1) then
    open(unit=23,file='DATA_3D/mesh_file',status='old',action='read')
    read(23,*) nspec
  else
    open(unit=23,file='DATA_3D/mesh_file.bin',form='unformatted',status='old',action='read')
    read(23) nspec
  endif

  allocate(ibool(NGNOD,nspec))

! loop on the whole mesh
  if (iformat == 1) then
    do ispec_loop = 1,nspec
      read(23,*) ispec,(ibool(ia,ispec), ia = 1,NGNOD)
    enddo
  else
    read(23) ibool
  endif

  close(23)

  print *,'Total number of elements in the mesh read = ',nspec
  print *

! ************* generate "absorbing_surface_file_xmin" *************

! first count the number of faces that are along that edge

  count_faces_found = 0

! Xmin
  size_of_model = xmax - xmin
  limit = xmin + SMALL_RELATIVE_VALUE*size_of_model

! loop on the whole mesh
  do ispec = 1,nspec

    i1 = ibool(1,ispec)
    i2 = ibool(2,ispec)
    i3 = ibool(3,ispec)
    i4 = ibool(4,ispec)
    i5 = ibool(5,ispec)
    i6 = ibool(6,ispec)
    i7 = ibool(7,ispec)
    i8 = ibool(8,ispec)

    if (.true.) then

      already_found_a_face = .false.

! test face 1 (bottom)
      if (x(i1) < limit .and. x(i2) < limit .and. x(i3) < limit .and. x(i4) < limit) then
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 2 (top)
      if (x(i5) < limit .and. x(i6) < limit .and. x(i7) < limit .and. x(i8) < limit) then
        if (already_found_a_face) stop 'error: element with two faces on the same edge found!'
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 3 (left)
      if (x(i1) < limit .and. x(i4) < limit .and. x(i8) < limit .and. x(i5) < limit) then
        if (already_found_a_face) stop 'error: element with two faces on the same edge found!'
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 4 (right)
      if (x(i2) < limit .and. x(i3) < limit .and. x(i7) < limit .and. x(i6) < limit) then
        if (already_found_a_face) stop 'error: element with two faces on the same edge found!'
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 5 (front)
      if (x(i1) < limit .and. x(i2) < limit .and. x(i6) < limit .and. x(i5) < limit) then
        if (already_found_a_face) stop 'error: element with two faces on the same edge found!'
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 6 (back)
      if (x(i4) < limit .and. x(i3) < limit .and. x(i7) < limit .and. x(i8) < limit) then
        if (already_found_a_face) stop 'error: element with two faces on the same edge found!'
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

    endif

  enddo

  print *,'found ',count_faces_found,' full faces on face Xmin'

!-----------------------------

! open SPECFEM3D_Cartesian topology file to read the mesh elements
  open(unit=24,file='DATA_3D/absorbing_surface_file_xmin',status='unknown',action='write')

! write the total number of face elements
  write(24,*) count_faces_found

! loop on the whole mesh
  do ispec = 1,nspec

    i1 = ibool(1,ispec)
    i2 = ibool(2,ispec)
    i3 = ibool(3,ispec)
    i4 = ibool(4,ispec)
    i5 = ibool(5,ispec)
    i6 = ibool(6,ispec)
    i7 = ibool(7,ispec)
    i8 = ibool(8,ispec)

    if (NGNOD == 27) then
       i9 = ibool(9,ispec)
      i10 = ibool(10,ispec)
      i11 = ibool(11,ispec)
      i12 = ibool(12,ispec)
      i13 = ibool(13,ispec)
      i14 = ibool(14,ispec)
      i15 = ibool(15,ispec)
      i16 = ibool(16,ispec)
      i17 = ibool(17,ispec)
      i18 = ibool(18,ispec)
      i19 = ibool(19,ispec)
      i20 = ibool(20,ispec)
      i21 = ibool(21,ispec)
      i22 = ibool(22,ispec)
      i23 = ibool(23,ispec)
      i24 = ibool(24,ispec)
      i25 = ibool(25,ispec)
      i26 = ibool(26,ispec)
      i27 = ibool(27,ispec)
    endif

    if (.true.) then

! for the six faces below it is important to make sure we write the four points
! in an order for which the normal to the face points outwards

! test face 1 (bottom)
      if (x(i1) < limit .and. x(i2) < limit .and. x(i3) < limit .and. x(i4) < limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i4,i3,i2,i1
        else
          write(24,"(10(i9,1x))") ispec,i4,i3,i2,i1,i11,i10,i9,i12,i21
        endif
      endif

! test face 2 (top)
      if (x(i5) < limit .and. x(i6) < limit .and. x(i7) < limit .and. x(i8) < limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i5,i6,i7,i8
        else
          write(24,"(10(i9,1x))") ispec,i5,i6,i7,i8,i17,i18,i19,i20,i26
        endif
      endif

! test face 3 (left)
      if (x(i1) < limit .and. x(i4) < limit .and. x(i8) < limit .and. x(i5) < limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i1,i5,i8,i4
        else
          write(24,"(10(i9,1x))") ispec,i1,i5,i8,i4,i13,i20,i16,i12,i25
        endif
      endif

! test face 4 (right)
      if (x(i2) < limit .and. x(i3) < limit .and. x(i7) < limit .and. x(i6) < limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i2,i3,i7,i6
        else
          write(24,"(10(i9,1x))") ispec,i2,i3,i7,i6,i10,i15,i18,i14,i23
        endif
      endif

! test face 5 (front)
      if (x(i1) < limit .and. x(i2) < limit .and. x(i6) < limit .and. x(i5) < limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i1,i2,i6,i5
        else
          write(24,"(10(i9,1x))") ispec,i1,i2,i6,i5,i9,i14,i17,i13,i22
        endif
      endif

! test face 6 (back)
      if (x(i4) < limit .and. x(i3) < limit .and. x(i7) < limit .and. x(i8) < limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i3,i4,i8,i7
        else
          write(24,"(10(i9,1x))") ispec,i3,i4,i8,i7,i11,i16,i19,i15,i24
        endif
      endif

    endif

  enddo

  close(24)

  print *,'File "absorbing_surface_file_xmin" has been successfully created'
  print *

! ************* generate "absorbing_surface_file_xmax" *************

! first count the number of faces that are along that edge

  count_faces_found = 0

! Xmax
  size_of_model = xmax - xmin
  limit = xmax - SMALL_RELATIVE_VALUE*size_of_model

! loop on the whole mesh
  do ispec = 1,nspec

    i1 = ibool(1,ispec)
    i2 = ibool(2,ispec)
    i3 = ibool(3,ispec)
    i4 = ibool(4,ispec)
    i5 = ibool(5,ispec)
    i6 = ibool(6,ispec)
    i7 = ibool(7,ispec)
    i8 = ibool(8,ispec)

    if (.true.) then

      already_found_a_face = .false.

! test face 1 (bottom)
      if (x(i1) > limit .and. x(i2) > limit .and. x(i3) > limit .and. x(i4) > limit) then
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 2 (top)
      if (x(i5) > limit .and. x(i6) > limit .and. x(i7) > limit .and. x(i8) > limit) then
        if (already_found_a_face) stop 'error: element with two faces on the same edge found!'
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 3 (left)
      if (x(i1) > limit .and. x(i4) > limit .and. x(i8) > limit .and. x(i5) > limit) then
        if (already_found_a_face) stop 'error: element with two faces on the same edge found!'
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 4 (right)
      if (x(i2) > limit .and. x(i3) > limit .and. x(i7) > limit .and. x(i6) > limit) then
        if (already_found_a_face) stop 'error: element with two faces on the same edge found!'
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 5 (front)
      if (x(i1) > limit .and. x(i2) > limit .and. x(i6) > limit .and. x(i5) > limit) then
        if (already_found_a_face) stop 'error: element with two faces on the same edge found!'
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 6 (back)
      if (x(i4) > limit .and. x(i3) > limit .and. x(i7) > limit .and. x(i8) > limit) then
        if (already_found_a_face) stop 'error: element with two faces on the same edge found!'
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

    endif

  enddo

  print *,'found ',count_faces_found,' full faces on face Xmax'

!-----------------------------

  open(unit=24,file='DATA_3D/absorbing_surface_file_xmax',status='unknown',action='write')

! write the total number of face elements
  write(24,*) count_faces_found

! loop on the whole mesh
  do ispec = 1,nspec

    i1 = ibool(1,ispec)
    i2 = ibool(2,ispec)
    i3 = ibool(3,ispec)
    i4 = ibool(4,ispec)
    i5 = ibool(5,ispec)
    i6 = ibool(6,ispec)
    i7 = ibool(7,ispec)
    i8 = ibool(8,ispec)

    if (NGNOD == 27) then
       i9 = ibool(9,ispec)
      i10 = ibool(10,ispec)
      i11 = ibool(11,ispec)
      i12 = ibool(12,ispec)
      i13 = ibool(13,ispec)
      i14 = ibool(14,ispec)
      i15 = ibool(15,ispec)
      i16 = ibool(16,ispec)
      i17 = ibool(17,ispec)
      i18 = ibool(18,ispec)
      i19 = ibool(19,ispec)
      i20 = ibool(20,ispec)
      i21 = ibool(21,ispec)
      i22 = ibool(22,ispec)
      i23 = ibool(23,ispec)
      i24 = ibool(24,ispec)
      i25 = ibool(25,ispec)
      i26 = ibool(26,ispec)
      i27 = ibool(27,ispec)
    endif

    if (.true.) then

! for the six faces below it is important to make sure we write the four points
! in an order for which the normal to the face points outwards

! test face 1 (bottom)
      if (x(i1) > limit .and. x(i2) > limit .and. x(i3) > limit .and. x(i4) > limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i4,i3,i2,i1
        else
          write(24,"(10(i9,1x))") ispec,i4,i3,i2,i1,i11,i10,i9,i12,i21
        endif
      endif

! test face 2 (top)
      if (x(i5) > limit .and. x(i6) > limit .and. x(i7) > limit .and. x(i8) > limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i5,i6,i7,i8
        else
          write(24,"(10(i9,1x))") ispec,i5,i6,i7,i8,i17,i18,i19,i20,i26
        endif
      endif

! test face 3 (left)
      if (x(i1) > limit .and. x(i4) > limit .and. x(i8) > limit .and. x(i5) > limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i1,i5,i8,i4
        else
          write(24,"(10(i9,1x))") ispec,i1,i5,i8,i4,i13,i20,i16,i12,i25
        endif
      endif

! test face 4 (right)
      if (x(i2) > limit .and. x(i3) > limit .and. x(i7) > limit .and. x(i6) > limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i2,i3,i7,i6
        else
          write(24,"(10(i9,1x))") ispec,i2,i3,i7,i6,i10,i15,i18,i14,i23
        endif
      endif

! test face 5 (front)
      if (x(i1) > limit .and. x(i2) > limit .and. x(i6) > limit .and. x(i5) > limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i1,i2,i6,i5
        else
          write(24,"(10(i9,1x))") ispec,i1,i2,i6,i5,i9,i14,i17,i13,i22
        endif
      endif

! test face 6 (back)
      if (x(i4) > limit .and. x(i3) > limit .and. x(i7) > limit .and. x(i8) > limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i3,i4,i8,i7
        else
          write(24,"(10(i9,1x))") ispec,i3,i4,i8,i7,i11,i16,i19,i15,i24
        endif
      endif

    endif

  enddo

  close(24)

  print *,'File "absorbing_surface_file_xmax" has been successfully created'
  print *

! ************* generate "absorbing_surface_file_ymin" *************

! first count the number of faces that are along that edge

  count_faces_found = 0

! Ymin
  size_of_model = ymax - ymin
  limit = ymin + SMALL_RELATIVE_VALUE*size_of_model

! loop on the whole mesh
  do ispec = 1,nspec

    i1 = ibool(1,ispec)
    i2 = ibool(2,ispec)
    i3 = ibool(3,ispec)
    i4 = ibool(4,ispec)
    i5 = ibool(5,ispec)
    i6 = ibool(6,ispec)
    i7 = ibool(7,ispec)
    i8 = ibool(8,ispec)

    if (.true.) then

      already_found_a_face = .false.

! test face 1 (bottom)
      if (y(i1) < limit .and. y(i2) < limit .and. y(i3) < limit .and. y(i4) < limit) then
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 2 (top)
      if (y(i5) < limit .and. y(i6) < limit .and. y(i7) < limit .and. y(i8) < limit) then
        if (already_found_a_face) stop 'error: element with two faces on the same edge found!'
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 3 (left)
      if (y(i1) < limit .and. y(i4) < limit .and. y(i8) < limit .and. y(i5) < limit) then
        if (already_found_a_face) stop 'error: element with two faces on the same edge found!'
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 4 (right)
      if (y(i2) < limit .and. y(i3) < limit .and. y(i7) < limit .and. y(i6) < limit) then
        if (already_found_a_face) stop 'error: element with two faces on the same edge found!'
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 5 (front)
      if (y(i1) < limit .and. y(i2) < limit .and. y(i6) < limit .and. y(i5) < limit) then
        if (already_found_a_face) stop 'error: element with two faces on the same edge found!'
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 6 (back)
      if (y(i4) < limit .and. y(i3) < limit .and. y(i7) < limit .and. y(i8) < limit) then
        if (already_found_a_face) stop 'error: element with two faces on the same edge found!'
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

    endif

  enddo

  print *,'found ',count_faces_found,' full faces on face Ymin'

!-----------------------------

  open(unit=24,file='DATA_3D/absorbing_surface_file_ymin',status='unknown',action='write')

! write the total number of face elements
  write(24,*) count_faces_found

! loop on the whole mesh
  do ispec = 1,nspec

    i1 = ibool(1,ispec)
    i2 = ibool(2,ispec)
    i3 = ibool(3,ispec)
    i4 = ibool(4,ispec)
    i5 = ibool(5,ispec)
    i6 = ibool(6,ispec)
    i7 = ibool(7,ispec)
    i8 = ibool(8,ispec)

    if (NGNOD == 27) then
       i9 = ibool(9,ispec)
      i10 = ibool(10,ispec)
      i11 = ibool(11,ispec)
      i12 = ibool(12,ispec)
      i13 = ibool(13,ispec)
      i14 = ibool(14,ispec)
      i15 = ibool(15,ispec)
      i16 = ibool(16,ispec)
      i17 = ibool(17,ispec)
      i18 = ibool(18,ispec)
      i19 = ibool(19,ispec)
      i20 = ibool(20,ispec)
      i21 = ibool(21,ispec)
      i22 = ibool(22,ispec)
      i23 = ibool(23,ispec)
      i24 = ibool(24,ispec)
      i25 = ibool(25,ispec)
      i26 = ibool(26,ispec)
      i27 = ibool(27,ispec)
    endif

    if (.true.) then

! for the six faces below it is important to make sure we write the four points
! in an order for which the normal to the face points outwards

! test face 1 (bottom)
      if (y(i1) < limit .and. y(i2) < limit .and. y(i3) < limit .and. y(i4) < limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i4,i3,i2,i1
        else
          write(24,"(10(i9,1x))") ispec,i4,i3,i2,i1,i11,i10,i9,i12,i21
        endif
      endif

! test face 2 (top)
      if (y(i5) < limit .and. y(i6) < limit .and. y(i7) < limit .and. y(i8) < limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i5,i6,i7,i8
        else
          write(24,"(10(i9,1x))") ispec,i5,i6,i7,i8,i17,i18,i19,i20,i26
        endif
      endif

! test face 3 (left)
      if (y(i1) < limit .and. y(i4) < limit .and. y(i8) < limit .and. y(i5) < limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i1,i5,i8,i4
        else
          write(24,"(10(i9,1x))") ispec,i1,i5,i8,i4,i13,i20,i16,i12,i25
        endif
      endif

! test face 4 (right)
      if (y(i2) < limit .and. y(i3) < limit .and. y(i7) < limit .and. y(i6) < limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i2,i3,i7,i6
        else
          write(24,"(10(i9,1x))") ispec,i2,i3,i7,i6,i10,i15,i18,i14,i23
        endif
      endif

! test face 5 (front)
      if (y(i1) < limit .and. y(i2) < limit .and. y(i6) < limit .and. y(i5) < limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i1,i2,i6,i5
        else
          write(24,"(10(i9,1x))") ispec,i1,i2,i6,i5,i9,i14,i17,i13,i22
        endif
      endif

! test face 6 (back)
      if (y(i4) < limit .and. y(i3) < limit .and. y(i7) < limit .and. y(i8) < limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i3,i4,i8,i7
        else
          write(24,"(10(i9,1x))") ispec,i3,i4,i8,i7,i11,i16,i19,i15,i24
        endif
      endif

    endif

  enddo

  close(24)

  print *,'File "absorbing_surface_file_ymin" has been successfully created'
  print *

! ************* generate "absorbing_surface_file_ymax" *************

! first count the number of faces that are along that edge

  count_faces_found = 0

! Ymax
  size_of_model = ymax - ymin
  limit = ymax - SMALL_RELATIVE_VALUE*size_of_model

! loop on the whole mesh
  do ispec = 1,nspec

    i1 = ibool(1,ispec)
    i2 = ibool(2,ispec)
    i3 = ibool(3,ispec)
    i4 = ibool(4,ispec)
    i5 = ibool(5,ispec)
    i6 = ibool(6,ispec)
    i7 = ibool(7,ispec)
    i8 = ibool(8,ispec)

    if (.true.) then

      already_found_a_face = .false.

! test face 1 (bottom)
      if (y(i1) > limit .and. y(i2) > limit .and. y(i3) > limit .and. y(i4) > limit) then
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 2 (top)
      if (y(i5) > limit .and. y(i6) > limit .and. y(i7) > limit .and. y(i8) > limit) then
        if (already_found_a_face) stop 'error: element with two faces on the same edge found!'
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 3 (left)
      if (y(i1) > limit .and. y(i4) > limit .and. y(i8) > limit .and. y(i5) > limit) then
        if (already_found_a_face) stop 'error: element with two faces on the same edge found!'
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 4 (right)
      if (y(i2) > limit .and. y(i3) > limit .and. y(i7) > limit .and. y(i6) > limit) then
        if (already_found_a_face) stop 'error: element with two faces on the same edge found!'
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 5 (front)
      if (y(i1) > limit .and. y(i2) > limit .and. y(i6) > limit .and. y(i5) > limit) then
        if (already_found_a_face) stop 'error: element with two faces on the same edge found!'
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 6 (back)
      if (y(i4) > limit .and. y(i3) > limit .and. y(i7) > limit .and. y(i8) > limit) then
        if (already_found_a_face) stop 'error: element with two faces on the same edge found!'
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

    endif

  enddo

  print *,'found ',count_faces_found,' full faces on face Ymax'

!-----------------------------

  open(unit=24,file='DATA_3D/absorbing_surface_file_ymax',status='unknown',action='write')

! write the total number of face elements
  write(24,*) count_faces_found

! loop on the whole mesh
  do ispec = 1,nspec

    i1 = ibool(1,ispec)
    i2 = ibool(2,ispec)
    i3 = ibool(3,ispec)
    i4 = ibool(4,ispec)
    i5 = ibool(5,ispec)
    i6 = ibool(6,ispec)
    i7 = ibool(7,ispec)
    i8 = ibool(8,ispec)

    if (NGNOD == 27) then
       i9 = ibool(9,ispec)
      i10 = ibool(10,ispec)
      i11 = ibool(11,ispec)
      i12 = ibool(12,ispec)
      i13 = ibool(13,ispec)
      i14 = ibool(14,ispec)
      i15 = ibool(15,ispec)
      i16 = ibool(16,ispec)
      i17 = ibool(17,ispec)
      i18 = ibool(18,ispec)
      i19 = ibool(19,ispec)
      i20 = ibool(20,ispec)
      i21 = ibool(21,ispec)
      i22 = ibool(22,ispec)
      i23 = ibool(23,ispec)
      i24 = ibool(24,ispec)
      i25 = ibool(25,ispec)
      i26 = ibool(26,ispec)
      i27 = ibool(27,ispec)
    endif

    if (.true.) then

! for the six faces below it is important to make sure we write the four points
! in an order for which the normal to the face points outwards

! test face 1 (bottom)
      if (y(i1) > limit .and. y(i2) > limit .and. y(i3) > limit .and. y(i4) > limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i4,i3,i2,i1
        else
          write(24,"(10(i9,1x))") ispec,i4,i3,i2,i1,i11,i10,i9,i12,i21
        endif
      endif

! test face 2 (top)
      if (y(i5) > limit .and. y(i6) > limit .and. y(i7) > limit .and. y(i8) > limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i5,i6,i7,i8
        else
          write(24,"(10(i9,1x))") ispec,i5,i6,i7,i8,i17,i18,i19,i20,i26
        endif
      endif

! test face 3 (left)
      if (y(i1) > limit .and. y(i4) > limit .and. y(i8) > limit .and. y(i5) > limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i1,i5,i8,i4
        else
          write(24,"(10(i9,1x))") ispec,i1,i5,i8,i4,i13,i20,i16,i12,i25
        endif
      endif

! test face 4 (right)
      if (y(i2) > limit .and. y(i3) > limit .and. y(i7) > limit .and. y(i6) > limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i2,i3,i7,i6
        else
          write(24,"(10(i9,1x))") ispec,i2,i3,i7,i6,i10,i15,i18,i14,i23
        endif
      endif

! test face 5 (front)
      if (y(i1) > limit .and. y(i2) > limit .and. y(i6) > limit .and. y(i5) > limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i1,i2,i6,i5
        else
          write(24,"(10(i9,1x))") ispec,i1,i2,i6,i5,i9,i14,i17,i13,i22
        endif
      endif

! test face 6 (back)
      if (y(i4) > limit .and. y(i3) > limit .and. y(i7) > limit .and. y(i8) > limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i3,i4,i8,i7
        else
          write(24,"(10(i9,1x))") ispec,i3,i4,i8,i7,i11,i16,i19,i15,i24
        endif
      endif

    endif

  enddo

  close(24)

  print *,'File "absorbing_surface_file_ymax" has been successfully created'
  print *

! ************* generate "absorbing_surface_file_bottom" *************

! first count the number of faces that are along that edge

  count_faces_found = 0

! Zmin
  size_of_model = zmax - zmin
  limit = zmin + SMALL_RELATIVE_VALUE*size_of_model

! loop on the whole mesh
  do ispec = 1,nspec

    i1 = ibool(1,ispec)
    i2 = ibool(2,ispec)
    i3 = ibool(3,ispec)
    i4 = ibool(4,ispec)
    i5 = ibool(5,ispec)
    i6 = ibool(6,ispec)
    i7 = ibool(7,ispec)
    i8 = ibool(8,ispec)

    if (.true.) then

      already_found_a_face = .false.

! test face 1 (bottom)
      if (z(i1) < limit .and. z(i2) < limit .and. z(i3) < limit .and. z(i4) < limit) then
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 2 (top)
      if (z(i5) < limit .and. z(i6) < limit .and. z(i7) < limit .and. z(i8) < limit) then
        if (already_found_a_face) stop 'error: element with two faces on the same edge found!'
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 3 (left)
      if (z(i1) < limit .and. z(i4) < limit .and. z(i8) < limit .and. z(i5) < limit) then
        if (already_found_a_face) stop 'error: element with two faces on the same edge found!'
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 4 (right)
      if (z(i2) < limit .and. z(i3) < limit .and. z(i7) < limit .and. z(i6) < limit) then
        if (already_found_a_face) stop 'error: element with two faces on the same edge found!'
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 5 (front)
      if (z(i1) < limit .and. z(i2) < limit .and. z(i6) < limit .and. z(i5) < limit) then
        if (already_found_a_face) stop 'error: element with two faces on the same edge found!'
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 6 (back)
      if (z(i4) < limit .and. z(i3) < limit .and. z(i7) < limit .and. z(i8) < limit) then
        if (already_found_a_face) stop 'error: element with two faces on the same edge found!'
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

    endif

  enddo

  print *,'found ',count_faces_found,' full faces on face Zmin'

!-----------------------------

  open(unit=24,file='DATA_3D/absorbing_surface_file_bottom',status='unknown',action='write')

! write the total number of face elements
  write(24,*) count_faces_found

! loop on the whole mesh
  do ispec = 1,nspec

    i1 = ibool(1,ispec)
    i2 = ibool(2,ispec)
    i3 = ibool(3,ispec)
    i4 = ibool(4,ispec)
    i5 = ibool(5,ispec)
    i6 = ibool(6,ispec)
    i7 = ibool(7,ispec)
    i8 = ibool(8,ispec)

    if (NGNOD == 27) then
       i9 = ibool(9,ispec)
      i10 = ibool(10,ispec)
      i11 = ibool(11,ispec)
      i12 = ibool(12,ispec)
      i13 = ibool(13,ispec)
      i14 = ibool(14,ispec)
      i15 = ibool(15,ispec)
      i16 = ibool(16,ispec)
      i17 = ibool(17,ispec)
      i18 = ibool(18,ispec)
      i19 = ibool(19,ispec)
      i20 = ibool(20,ispec)
      i21 = ibool(21,ispec)
      i22 = ibool(22,ispec)
      i23 = ibool(23,ispec)
      i24 = ibool(24,ispec)
      i25 = ibool(25,ispec)
      i26 = ibool(26,ispec)
      i27 = ibool(27,ispec)
    endif

    if (.true.) then

! for the six faces below it is important to make sure we write the four points
! in an order for which the normal to the face points outwards

! test face 1 (bottom)
      if (z(i1) < limit .and. z(i2) < limit .and. z(i3) < limit .and. z(i4) < limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i4,i3,i2,i1
        else
          write(24,"(10(i9,1x))") ispec,i4,i3,i2,i1,i11,i10,i9,i12,i21
        endif
      endif

! test face 2 (top)
      if (z(i5) < limit .and. z(i6) < limit .and. z(i7) < limit .and. z(i8) < limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i5,i6,i7,i8
        else
          write(24,"(10(i9,1x))") ispec,i5,i6,i7,i8,i17,i18,i19,i20,i26
        endif
      endif

! test face 3 (left)
      if (z(i1) < limit .and. z(i4) < limit .and. z(i8) < limit .and. z(i5) < limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i1,i5,i8,i4
        else
          write(24,"(10(i9,1x))") ispec,i1,i5,i8,i4,i13,i20,i16,i12,i25
        endif
      endif

! test face 4 (right)
      if (z(i2) < limit .and. z(i3) < limit .and. z(i7) < limit .and. z(i6) < limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i2,i3,i7,i6
        else
          write(24,"(10(i9,1x))") ispec,i2,i3,i7,i6,i10,i15,i18,i14,i23
        endif
      endif

! test face 5 (front)
      if (z(i1) < limit .and. z(i2) < limit .and. z(i6) < limit .and. z(i5) < limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i1,i2,i6,i5
        else
          write(24,"(10(i9,1x))") ispec,i1,i2,i6,i5,i9,i14,i17,i13,i22
        endif
      endif

! test face 6 (back)
      if (z(i4) < limit .and. z(i3) < limit .and. z(i7) < limit .and. z(i8) < limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i3,i4,i8,i7
        else
          write(24,"(10(i9,1x))") ispec,i3,i4,i8,i7,i11,i16,i19,i15,i24
        endif
      endif

    endif

  enddo

  close(24)

  print *,'File "absorbing_surface_file_bottom" has been successfully created'
  print *

! ************* generate "free_or_absorbing_surface_file_zmax" *************

! first count the number of faces that are along that edge

  count_faces_found = 0

! Zmax
  size_of_model = zmax - zmin
  limit = zmax - SMALL_RELATIVE_VALUE*size_of_model

! loop on the whole mesh
  do ispec = 1,nspec

    i1 = ibool(1,ispec)
    i2 = ibool(2,ispec)
    i3 = ibool(3,ispec)
    i4 = ibool(4,ispec)
    i5 = ibool(5,ispec)
    i6 = ibool(6,ispec)
    i7 = ibool(7,ispec)
    i8 = ibool(8,ispec)

      already_found_a_face = .false.

! test face 1 (bottom)
      if (z(i1) > limit .and. z(i2) > limit .and. z(i3) > limit .and. z(i4) > limit) then
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 2 (top)
      if (z(i5) > limit .and. z(i6) > limit .and. z(i7) > limit .and. z(i8) > limit) then
        if (already_found_a_face) stop 'error: element with two faces on the same edge found!'
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 3 (left)
      if (z(i1) > limit .and. z(i4) > limit .and. z(i8) > limit .and. z(i5) > limit) then
        if (already_found_a_face) stop 'error: element with two faces on the same edge found!'
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 4 (right)
      if (z(i2) > limit .and. z(i3) > limit .and. z(i7) > limit .and. z(i6) > limit) then
        if (already_found_a_face) stop 'error: element with two faces on the same edge found!'
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 5 (front)
      if (z(i1) > limit .and. z(i2) > limit .and. z(i6) > limit .and. z(i5) > limit) then
        if (already_found_a_face) stop 'error: element with two faces on the same edge found!'
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 6 (back)
      if (z(i4) > limit .and. z(i3) > limit .and. z(i7) > limit .and. z(i8) > limit) then
        if (already_found_a_face) stop 'error: element with two faces on the same edge found!'
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

  enddo

  print *,'found ',count_faces_found,' full faces on face Zmax'

!-----------------------------

  open(unit=24,file='DATA_3D/free_or_absorbing_surface_file_zmax',status='unknown',action='write')

! write the total number of face elements
  write(24,*) count_faces_found

! loop on the whole mesh
  do ispec = 1,nspec

    i1 = ibool(1,ispec)
    i2 = ibool(2,ispec)
    i3 = ibool(3,ispec)
    i4 = ibool(4,ispec)
    i5 = ibool(5,ispec)
    i6 = ibool(6,ispec)
    i7 = ibool(7,ispec)
    i8 = ibool(8,ispec)

    if (NGNOD == 27) then
       i9 = ibool(9,ispec)
      i10 = ibool(10,ispec)
      i11 = ibool(11,ispec)
      i12 = ibool(12,ispec)
      i13 = ibool(13,ispec)
      i14 = ibool(14,ispec)
      i15 = ibool(15,ispec)
      i16 = ibool(16,ispec)
      i17 = ibool(17,ispec)
      i18 = ibool(18,ispec)
      i19 = ibool(19,ispec)
      i20 = ibool(20,ispec)
      i21 = ibool(21,ispec)
      i22 = ibool(22,ispec)
      i23 = ibool(23,ispec)
      i24 = ibool(24,ispec)
      i25 = ibool(25,ispec)
      i26 = ibool(26,ispec)
      i27 = ibool(27,ispec)
    endif

! for the six faces below it is important to make sure we write the four points
! in an order for which the normal to the face points outwards

! test face 1 (bottom)
      if (z(i1) > limit .and. z(i2) > limit .and. z(i3) > limit .and. z(i4) > limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i4,i3,i2,i1
        else
          write(24,"(10(i9,1x))") ispec,i4,i3,i2,i1,i11,i10,i9,i12,i21
        endif
      endif

! test face 2 (top)
      if (z(i5) > limit .and. z(i6) > limit .and. z(i7) > limit .and. z(i8) > limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i5,i6,i7,i8
        else
          write(24,"(10(i9,1x))") ispec,i5,i6,i7,i8,i17,i18,i19,i20,i26
        endif
      endif

! test face 3 (left)
      if (z(i1) > limit .and. z(i4) > limit .and. z(i8) > limit .and. z(i5) > limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i1,i5,i8,i4
        else
          write(24,"(10(i9,1x))") ispec,i1,i5,i8,i4,i13,i20,i16,i12,i25
        endif
      endif

! test face 4 (right)
      if (z(i2) > limit .and. z(i3) > limit .and. z(i7) > limit .and. z(i6) > limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i2,i3,i7,i6
        else
          write(24,"(10(i9,1x))") ispec,i2,i3,i7,i6,i10,i15,i18,i14,i23
        endif
      endif

! test face 5 (front)
      if (z(i1) > limit .and. z(i2) > limit .and. z(i6) > limit .and. z(i5) > limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i1,i2,i6,i5
        else
          write(24,"(10(i9,1x))") ispec,i1,i2,i6,i5,i9,i14,i17,i13,i22
        endif
      endif

! test face 6 (back)
      if (z(i4) > limit .and. z(i3) > limit .and. z(i7) > limit .and. z(i8) > limit) then
        if (NGNOD == 8) then
          write(24,*) ispec,i3,i4,i8,i7
        else
          write(24,"(10(i9,1x))") ispec,i3,i4,i8,i7,i11,i16,i19,i15,i24
        endif
      endif

  enddo

  close(24)

  print *,'File "free_or_absorbing_surface_file_zmax" has been successfully created'
  print *

  end program create_database_files_for_external_faces_of_the_3D_mesh

