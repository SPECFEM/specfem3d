!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
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

  program add_CPML_layers_to_a_given_mesh

! Dimitri Komatitsch, CNRS Marseille, France, February 2016

! add PML layers around an existing mesh (i.e. create new elements and new points)

  implicit none

! this is for HEX8; the code below also works for HEX27, it then just uses the first 8 points of an element
! to determine if it belongs to a CPML layer
  integer, parameter :: NGNOD = 8

! number of PML layers to add on each side of the mesh
  integer :: NUMBER_OF_PML_LAYERS_TO_ADD

! size of each PML element to add when they are added on the Xmin or Xmax faces
  double precision :: SIZE_OF_X_ELEMENT_TO_ADD

  integer :: nspec,npoin,npoin_new,nspec_new,count_elem_faces_to_extend_xmin,iextend
  integer :: ispec,ipoin
  integer :: ipoin_read,ispec_loop
  integer :: i1,i2,i3,i4,i5,i6,i7,i8,elem_counter,ibool_counter,ia
  integer :: p1,p2,p3,p4

  double precision, dimension(:), allocatable :: x,y,z,xnew,ynew,znew

  integer, dimension(:), allocatable :: imaterial,imaterial_new

  integer, dimension(:,:), allocatable :: ibool,ibool_new

  double precision :: xread,yread,zread,xmin,xmax,ymin,ymax,zmin,zmax,limit,xsize,ysize,zsize

! logical :: ALSO_ADD_ON_THE_TOP_SURFACE
  logical :: need_to_extend_this_element

  double precision, parameter :: SMALL_RELATIVE_VALUE = 1.d-5

  print *
  print *,'IMPORTANT: it is your responsibility to make sure that the input mesh'
  print *,'that this code will read in SPECFEM3D format from files "nodes_coords_file" and "mesh_file"'
  print *,'has flat outer edges aligned with the coordinate grid axes (X, Y and/or Z),'
  print *,'so that this code can external CPML elements to it.'
  print *,'This code does NOT check that (because it cannot, in any easy way).'
  print *,'The mesh does not need to be structured nor regular, any non-structured'
  print *,'mesh is fine as long as it has flat outer faces, parallel to the axes.'
  print *

!! DK DK to Vadim: il faudra que tu remettes ca mais je l'enleve pour l'instant car dans mon exemple
!! DK DK to Vadim: je ne fais que le cas Xmin pour l'instant
! print *,'1 = use a free surface at the top of the mesh (most classical option)'
! print *,'2 = use a CPML absorbing layer at the top of the mesh (less classical option)'
! print *,'3 = exit'
! read(*,*) iflag
! if(iflag /= 1 .and. iflag /= 2) stop 'exiting...'
! if(iflag == 1) then
!   ALSO_ADD_ON_THE_TOP_SURFACE = .false.
! else
!   ALSO_ADD_ON_THE_TOP_SURFACE = .true.
! endif
! print *

  print *,'enter the X size (in meters) of each CPML element to add on the Xmin and Xmax faces:'
  read(*,*) SIZE_OF_X_ELEMENT_TO_ADD
  if(SIZE_OF_X_ELEMENT_TO_ADD <= 0.d0) stop 'NUMBER_OF_PML_LAYERS_TO_ADD must be > 0'
  print *

  print *,'enter the number of PML layers to add on each side of the mesh (usually 3, can also be 4):'
  read(*,*) NUMBER_OF_PML_LAYERS_TO_ADD
  if(NUMBER_OF_PML_LAYERS_TO_ADD < 1) stop 'NUMBER_OF_PML_LAYERS_TO_ADD must be >= 1'
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

! compute the min and max values of each coordinate
  xmin = minval(x)
  xmax = maxval(x)

  ymin = minval(y)
  ymax = maxval(y)

  zmin = minval(z)
  zmax = maxval(z)

  xsize = xmax - xmin
  ysize = ymax - ymin
  zsize = zmax - zmin

  print *,'Xmin and Xmax of the mesh read = ',xmin,xmax
  print *,'Ymin and Ymax of the mesh read = ',ymin,ymax
  print *,'Zmin and Zmax of the mesh read = ',zmin,zmax
  print *

  print *,'Size of the mesh read along X = ',xsize
  print *,'Size of the mesh read along Y = ',ysize
  print *,'Size of the mesh read along Z = ',zsize
  print *

! ************* read mesh elements and generate CPML flags *************

! open SPECFEM3D_Cartesian topology file to read the mesh elements
  open(unit=23,file='mesh_file',status='old',action='read')
  read(23,*) nspec

  count_elem_faces_to_extend_xmin = 0

! loop on the whole mesh
  do ispec_loop = 1,nspec

    read(23,*) ispec,i1,i2,i3,i4,i5,i6,i7,i8

! detect elements belonging to the Xmin face
    limit = xmin + xsize * SMALL_RELATIVE_VALUE

! test face 1 (bottom)
    if(x(i1) < limit .and. x(i2) < limit .and. x(i3) < limit .and. x(i4) < limit) &
      count_elem_faces_to_extend_xmin = count_elem_faces_to_extend_xmin + 1

! test face 2 (top)
    if(x(i5) < limit .and. x(i6) < limit .and. x(i7) < limit .and. x(i8) < limit) &
      count_elem_faces_to_extend_xmin = count_elem_faces_to_extend_xmin + 1

! test face 3 (left)
    if(x(i1) < limit .and. x(i4) < limit .and. x(i5) < limit .and. x(i8) < limit) &
      count_elem_faces_to_extend_xmin = count_elem_faces_to_extend_xmin + 1

! test face 4 (right)
    if(x(i2) < limit .and. x(i3) < limit .and. x(i7) < limit .and. x(i6) < limit) &
      count_elem_faces_to_extend_xmin = count_elem_faces_to_extend_xmin + 1

! test face 5 (front)
    if(x(i1) < limit .and. x(i2) < limit .and. x(i6) < limit .and. x(i5) < limit) &
      count_elem_faces_to_extend_xmin = count_elem_faces_to_extend_xmin + 1

! test face 6 (back)
    if(x(i4) < limit .and. x(i3) < limit .and. x(i7) < limit .and. x(i8) < limit) &
      count_elem_faces_to_extend_xmin = count_elem_faces_to_extend_xmin + 1

  enddo

  close(23)

  print *,'Total number of elements in the mesh before extension = ',nspec
  print *,'Number of element faces to extend  = ',count_elem_faces_to_extend_xmin
! we will add NUMBER_OF_PML_LAYERS_TO_ADD to each of the element faces detected that need to be extended
  nspec_new = nspec + count_elem_faces_to_extend_xmin * NUMBER_OF_PML_LAYERS_TO_ADD
! and each of these elements will have NGNOD points
! (some of them shared with other elements, but we do not care because they will be removed automatically by xdecompose_mesh)
  npoin_new = npoin + count_elem_faces_to_extend_xmin * NUMBER_OF_PML_LAYERS_TO_ADD * NGNOD
  print *,'Total number of elements in the mesh after extension = ',nspec_new
  print *

! allocate arrays with the right size of the future extended mesh
  allocate(imaterial(nspec_new))
  allocate(ibool(NGNOD,nspec_new))

! clear the arrays
  imaterial(:) = 0
  ibool(:,:) = 0

! now that the ibool array has been allocated with the right future size, read and store the topology of the original mesh
! open SPECFEM3D_Cartesian topology file to read the mesh elements
  open(unit=23,file='mesh_file',status='old',action='read')
  read(23,*) nspec

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

! read the materials file
  open(unit=23,file='materials_file',status='old',action='read')
! loop on the whole mesh
  do ispec_loop = 1,nspec
    read(23,*) ispec,i1
! store the imaterial() array read
    imaterial(ispec) = i1
  enddo
  close(23)

! allocate a new set of points, with multiples
  allocate(xnew(npoin_new))
  allocate(ynew(npoin_new))
  allocate(znew(npoin_new))

! copy the original points into the new set
  xnew(1:npoin) = x(1:npoin)
  ynew(1:npoin) = y(1:npoin)
  znew(1:npoin) = z(1:npoin)

! allocate a new set of elements, i.e. a new ibool()
  allocate(ibool_new(NGNOD,nspec_new))
  allocate(imaterial_new(nspec_new))

! copy the original elements into the new set
  ibool_new(:,:) = ibool(:,:)
  imaterial_new(:) = imaterial(:)

! position after which to start to create the new elements
  elem_counter = nspec
  ibool_counter = npoin

! loop on the whole original mesh
  do ispec = 1,nspec

    i1 = ibool(1,ispec)
    i2 = ibool(2,ispec)
    i3 = ibool(3,ispec)
    i4 = ibool(4,ispec)
    i5 = ibool(5,ispec)
    i6 = ibool(6,ispec)
    i7 = ibool(7,ispec)
    i8 = ibool(8,ispec)

! detect elements belonging to the Xmin face
    limit = xmin + xsize * SMALL_RELATIVE_VALUE

! reset flag
    need_to_extend_this_element = .false.

! test face 1 (bottom)
    if(x(i1) < limit .and. x(i2) < limit .and. x(i3) < limit .and. x(i4) < limit) then
      need_to_extend_this_element = .true.
      p1 = i1
      p2 = i2
      p3 = i3
      p4 = i4
    endif

! test face 2 (top)
    if(x(i5) < limit .and. x(i6) < limit .and. x(i7) < limit .and. x(i8) < limit) then
      need_to_extend_this_element = .true.
      p1 = i5
      p2 = i6
      p3 = i7
      p4 = i8
    endif

! test face 3 (left)
    if(x(i1) < limit .and. x(i4) < limit .and. x(i5) < limit .and. x(i8) < limit) then
      need_to_extend_this_element = .true.
      p1 = i1
      p2 = i4
      p3 = i5
      p4 = i8
    endif

! test face 4 (right)
    if(x(i2) < limit .and. x(i3) < limit .and. x(i7) < limit .and. x(i6) < limit) then
      need_to_extend_this_element = .true.
      p1 = i2
      p2 = i3
      p3 = i7
      p4 = i6
    endif

! test face 5 (front)
    if(x(i1) < limit .and. x(i2) < limit .and. x(i6) < limit .and. x(i5) < limit) then
      need_to_extend_this_element = .true.
      p1 = i1
      p2 = i2
      p3 = i6
      p4 = i5
    endif

! test face 6 (back)
    if(x(i4) < limit .and. x(i3) < limit .and. x(i7) < limit .and. x(i8) < limit) then
      need_to_extend_this_element = .true.
      p1 = i4
      p2 = i3
      p3 = i7
      p4 = i8
    endif

    if(need_to_extend_this_element) then

! create the NUMBER_OF_PML_LAYERS_TO_ADD new elements, to the left of the face found because we are on Xmin

! very important remark: it is OK to create duplicates of the mesh points in the loop below (i.e. not to tell the code
! that many of these points created are in fact shared between adjacent elements) because "xdecompose_mesh" will remove
! them automatically later on, thus no need to remove them here; this makes this PML mesh extrusion code much simpler to write.

      do iextend = 1,NUMBER_OF_PML_LAYERS_TO_ADD

        ! create a new element
        elem_counter = elem_counter + 1

        ! use the same material property for the extended elements as for the element being extended
        imaterial_new(elem_counter) = imaterial(ispec)

        ! create a new point
        ibool_counter = ibool_counter + 1
        ibool_new(1,elem_counter) = ibool_counter
! we use iextend for the first 4 points because we are on Xmin, to avoid a negative Jacobian (mirror symmetry of the element)
! on Xmax here we would use (iextend-1) instead
        xnew(ibool_counter) = x(p1) - SIZE_OF_X_ELEMENT_TO_ADD*iextend
        ynew(ibool_counter) = y(p1)
        znew(ibool_counter) = z(p1)

        ! create a new point
        ibool_counter = ibool_counter + 1
        ibool_new(2,elem_counter) = ibool_counter
        xnew(ibool_counter) = x(p2) - SIZE_OF_X_ELEMENT_TO_ADD*iextend
        ynew(ibool_counter) = y(p2)
        znew(ibool_counter) = z(p2)

        ! create a new point
        ibool_counter = ibool_counter + 1
        ibool_new(3,elem_counter) = ibool_counter
        xnew(ibool_counter) = x(p3) - SIZE_OF_X_ELEMENT_TO_ADD*iextend
        ynew(ibool_counter) = y(p3)
        znew(ibool_counter) = z(p3)

        ! create a new point
        ibool_counter = ibool_counter + 1
        ibool_new(4,elem_counter) = ibool_counter
        xnew(ibool_counter) = x(p4) - SIZE_OF_X_ELEMENT_TO_ADD*iextend
        ynew(ibool_counter) = y(p4)
        znew(ibool_counter) = z(p4)

        ! create a new point
        ibool_counter = ibool_counter + 1
        ibool_new(5,elem_counter) = ibool_counter
! we use (iextend-1) for the last 4 points because we are on Xmin, to avoid a negative Jacobian (mirror symmetry of the element)
! on Xmax here we would use iextend instead
        xnew(ibool_counter) = x(p1) - SIZE_OF_X_ELEMENT_TO_ADD*(iextend-1)
        ynew(ibool_counter) = y(p1)
        znew(ibool_counter) = z(p1)

        ! create a new point
        ibool_counter = ibool_counter + 1
        ibool_new(6,elem_counter) = ibool_counter
        xnew(ibool_counter) = x(p2) - SIZE_OF_X_ELEMENT_TO_ADD*(iextend-1)
        ynew(ibool_counter) = y(p2)
        znew(ibool_counter) = z(p2)

        ! create a new point
        ibool_counter = ibool_counter + 1
        ibool_new(7,elem_counter) = ibool_counter
        xnew(ibool_counter) = x(p3) - SIZE_OF_X_ELEMENT_TO_ADD*(iextend-1)
        ynew(ibool_counter) = y(p3)
        znew(ibool_counter) = z(p3)

        ! create a new point
        ibool_counter = ibool_counter + 1
        ibool_new(8,elem_counter) = ibool_counter
        xnew(ibool_counter) = x(p4) - SIZE_OF_X_ELEMENT_TO_ADD*(iextend-1)
        ynew(ibool_counter) = y(p4)
        znew(ibool_counter) = z(p4)

      enddo
    endif

  enddo

! write the new points (overwrite the old file)
  open(unit=23,file='nodes_coords_file',status='old',action='write')
  write(23,*) npoin_new
  do ipoin = 1,npoin_new
    write(23,*) ipoin,sngl(xnew(ipoin)),sngl(ynew(ipoin)),sngl(znew(ipoin))
  enddo
  close(23)

! write the new mesh elements (overwrite the old file)
  open(unit=23,file='mesh_file',status='old',action='write')
  write(23,*) nspec_new
! loop on the whole mesh
  do ispec = 1,nspec_new
    write(23,"(i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9)") ispec,(ibool_new(ia,ispec), ia = 1,NGNOD)
  enddo
  close(23)

! write the new material properties (overwrite the old file)
  open(unit=23,file='materials_file',status='old',action='write')
! loop on the whole mesh
  do ispec = 1,nspec_new
    write(23,*) ispec,imaterial_new(ispec)
  enddo
  close(23)

  end program add_CPML_layers_to_a_given_mesh

