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

! Dimitri Komatitsch, CNRS Marseille, France, March 2016

! add PML layers around an existing mesh (i.e. create new elements and new points)

  implicit none

! this code is for HEX8 only for now
  integer, parameter :: NGNOD = 8

! we are in 3D
  integer, parameter :: NDIM = 3

! number of GLL points in each direction, to check for negative Jacobians
  integer, parameter :: NGLLX = 5,NGLLY = NGLLX,NGLLZ = NGLLX

! number of PML layers to add on each side of the mesh
  integer :: NUMBER_OF_PML_LAYERS_TO_ADD

! size of each PML element to add when they are added on the Xmin and Xmax faces, Ymin and Ymax faces, Zmin and/or Zmax faces
  double precision :: SIZE_OF_X_ELEMENT_TO_ADD,SIZE_OF_Y_ELEMENT_TO_ADD,SIZE_OF_Z_ELEMENT_TO_ADD
  double precision :: SIZE_OF_XMIN_ELEMENT_TO_ADD,SIZE_OF_YMIN_ELEMENT_TO_ADD,SIZE_OF_ZMIN_ELEMENT_TO_ADD
  double precision :: SIZE_OF_XMAX_ELEMENT_TO_ADD,SIZE_OF_YMAX_ELEMENT_TO_ADD,SIZE_OF_ZMAX_ELEMENT_TO_ADD

  integer :: nspec,npoin,npoin_new,nspec_new,count_elem_faces_to_extend,iextend
  integer :: factor_x,factor_y,factor_z
  integer :: ispec,ipoin,iloop_on_X_Y_Z_faces,iloop_on_min_face_then_max_face
  integer :: ipoin_read,ispec_loop
  integer :: i1,i2,i3,i4,i5,i6,i7,i8,elem_counter,ibool_counter,ia,iflag,iformat,icompute_size
  integer :: p1,p2,p3,p4

  double precision, dimension(:), allocatable, target :: x,y,z
  double precision, dimension(:), allocatable :: x_new,y_new,z_new
  double precision, dimension(:), pointer :: coord_to_use1,coord_to_use2,coord_to_use3

  integer, dimension(:), allocatable :: imaterial,imaterial_new

  integer, dimension(:,:), allocatable :: ibool,ibool_new

! Gauss-Lobatto-Legendre points of integration, to check for negative Jacobians
  double precision xigll(NGLLX)
  double precision yigll(NGLLY)
  double precision zigll(NGLLZ)

! 3D shape function derivatives, to check for negative Jacobians
  double precision dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ)

  double precision, dimension(NGNOD) :: xelm,yelm,zelm

  double precision :: xread,yread,zread,xmin,xmax,ymin,ymax,zmin,zmax,limit,xsize,ysize,zsize
  double precision :: value_min,value_max,value_size,sum_of_distances,mean_distance

  logical :: ALSO_ADD_ON_THE_TOP_SURFACE,need_to_extend_this_element,found_a_negative_Jacobian1,found_a_negative_Jacobian2

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

  print *,'1 = read the mesh files in ASCII format (that is the standard case)'
  print *,'2 = read the mesh files in binary format (that is much faster, if your mesh is available in that format)'
  print *,'  (if not, you can run xconvert_mesh_files_from_ASCII_to_binary)'
  print *,'3 = exit'
  read(*,*) iformat
  if (iformat /= 1 .and. iformat /= 2) stop 'exiting...'

  print *,'enter the number of PML layers to add on each side of the mesh (usually 3, can also be 4):'
  read(*,*) NUMBER_OF_PML_LAYERS_TO_ADD
  if (NUMBER_OF_PML_LAYERS_TO_ADD < 1) stop 'NUMBER_OF_PML_LAYERS_TO_ADD must be >= 1'
  print *

  print *,'1 = use a free surface at the top of the mesh i.e. do not add a CPML layer at the top (most classical option)'
  print *,'2 = add a CPML absorbing layer at the top of the mesh (less classical option)'
  print *,'3 = exit'
  read(*,*) iflag
  if (iflag /= 1 .and. iflag /= 2) stop 'exiting...'
  if (iflag == 1) then
    ALSO_ADD_ON_THE_TOP_SURFACE = .false.
  else
    ALSO_ADD_ON_THE_TOP_SURFACE = .true.
  endif
  print *

  print *,'1 = compute the size of the PML elements to add automatically using the average size of edge elements'
  print *,'2 = enter the size of the PML elements to add manually'
  print *,'3 = exit'
  read(*,*) icompute_size
  if (icompute_size /= 1 .and. icompute_size /= 2) stop 'exiting...'

  if (icompute_size == 2) then

  print *,'enter the X size (in meters) of each CPML element to add on the Xmin face:'
  read(*,*) SIZE_OF_XMIN_ELEMENT_TO_ADD
  if (SIZE_OF_XMIN_ELEMENT_TO_ADD <= 0.d0) stop 'SIZE_OF_XMIN_ELEMENT_TO_ADD must be > 0'
  print *

  print *,'enter the X size (in meters) of each CPML element to add on the Xmax face:'
  read(*,*) SIZE_OF_XMAX_ELEMENT_TO_ADD
  if (SIZE_OF_XMAX_ELEMENT_TO_ADD <= 0.d0) stop 'SIZE_OF_XMAX_ELEMENT_TO_ADD must be > 0'
  print *

  print *,'enter the Y size (in meters) of each CPML element to add on the Ymin faces:'
  read(*,*) SIZE_OF_YMIN_ELEMENT_TO_ADD
  if (SIZE_OF_YMIN_ELEMENT_TO_ADD <= 0.d0) stop 'SIZE_OF_YMIN_ELEMENT_TO_ADD must be > 0'
  print *

  print *,'enter the Y size (in meters) of each CPML element to add on the Ymax faces:'
  read(*,*) SIZE_OF_YMAX_ELEMENT_TO_ADD
  if (SIZE_OF_YMAX_ELEMENT_TO_ADD <= 0.d0) stop 'SIZE_OF_YMAX_ELEMENT_TO_ADD must be > 0'
  print *

  print *,'enter the Z size (in meters) of each CPML element to add on the Zmin faces:'
  read(*,*) SIZE_OF_ZMIN_ELEMENT_TO_ADD
  if (SIZE_OF_ZMIN_ELEMENT_TO_ADD <= 0.d0) stop 'SIZE_OF_ZMIN_ELEMENT_TO_ADD must be > 0'
  print *

  if (ALSO_ADD_ON_THE_TOP_SURFACE) then
    print *,'enter the Z size (in meters) of each CPML element to add on the Zmax faces:'
    read(*,*) SIZE_OF_ZMAX_ELEMENT_TO_ADD
    if (SIZE_OF_ZMAX_ELEMENT_TO_ADD <= 0.d0) stop 'SIZE_OF_ZMAX_ELEMENT_TO_ADD must be > 0'
    print *
  endif

  endif

! hardwire GLL point location values to avoid having to link with a long library to compute them
  xigll(:) = (/ -1.d0 , -0.654653670707977d0 , 0.d0 , 0.654653670707977d0 , 1.d0 /)

! compute the derivatives of the 3D shape functions for a 8-node element
  call get_shape3D(dershape3D,xigll,yigll,zigll,NGNOD,NGLLX,NGLLY,NGLLZ,NDIM)

! open SPECFEM3D_Cartesian mesh file to read the points
  if (iformat == 1) then
    open(unit=23,file='nodes_coords_file',status='old',action='read')
    read(23,*) npoin
  else
    open(unit=23,file='nodes_coords_file.bin',form='unformatted',status='old',action='read')
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

! ************* read mesh elements *************

! open SPECFEM3D_Cartesian topology file to read the mesh elements
  if (iformat == 1) then
    open(unit=23,file='mesh_file',status='old',action='read')
    read(23,*) nspec
  else
    open(unit=23,file='mesh_file.bin',form='unformatted',status='old',action='read')
    read(23) nspec
  endif

  allocate(ibool(NGNOD,nspec))

! loop on the whole mesh
  if (iformat == 1) then
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
  else
    read(23) ibool
  endif

  close(23)

  print *,'Total number of elements in the mesh read = ',nspec

! read the materials file
  allocate(imaterial(nspec))
  if (iformat == 1) then
    open(unit=23,file='materials_file',status='old',action='read')
  else
    open(unit=23,file='materials_file.bin',form='unformatted',status='old',action='read')
  endif
! loop on the whole mesh
  if (iformat == 1) then
    do ispec_loop = 1,nspec
      read(23,*) ispec,i1
! store the imaterial() array read
      imaterial(ispec) = i1
    enddo
  else
    read(23) imaterial
  endif
  close(23)

! we need to extend/extrude the existing mesh by adding CPML elements
! along the X faces, then along the Y faces, then along the Z faces.

! loop on the three sets of faces to first add CPML elements along X, then along Y, then along Z
  do iloop_on_X_Y_Z_faces = 1,NDIM

    do iloop_on_min_face_then_max_face = 1,2  ! 1 is min face and 2 is max face (Xmin then Xmax, Ymin then Ymax, or Zmin then Zmax)

! do not add a CPML layer on the top surface if the user asked not to
      if (iloop_on_X_Y_Z_faces == 3 .and. iloop_on_min_face_then_max_face == 2 .and. .not. ALSO_ADD_ON_THE_TOP_SURFACE) cycle

    print *
    print *,'********************************************************************'
    if (iloop_on_X_Y_Z_faces == 1) then
      print *,'adding CPML elements along one of the two X faces of the existing mesh'
    else if (iloop_on_X_Y_Z_faces == 2) then
      print *,'adding CPML elements along one of the two Y faces of the existing mesh'
    else if (iloop_on_X_Y_Z_faces == 3) then
      print *,'adding CPML elements along one of the two Z faces of the existing mesh'
    else
      stop 'wrong index in loop on faces'
    endif
    print *,'********************************************************************'
    print *

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

  count_elem_faces_to_extend = 0

    if (iloop_on_X_Y_Z_faces == 1) then ! Xmin or Xmax face
      value_min = xmin
      value_max = xmax
      value_size = xsize
      coord_to_use1 => x  ! make coordinate array to use point to array x()
      coord_to_use2 => y
      coord_to_use3 => z
    else if (iloop_on_X_Y_Z_faces == 2) then ! Ymin or Ymax face
      value_min = ymin
      value_max = ymax
      value_size = ysize
      coord_to_use1 => y  ! make coordinate array to use point to array y()
      coord_to_use2 => x
      coord_to_use3 => z
    else if (iloop_on_X_Y_Z_faces == 3) then ! Zmin or Zmax face
      value_min = zmin
      value_max = zmax
      value_size = zsize
      coord_to_use1 => z  ! make coordinate array to use point to array z()
      coord_to_use2 => x
      coord_to_use3 => y
    else
      stop 'wrong index in loop on faces'
    endif

  if (minval(ibool) /= 1) stop 'error in minval(ibool)'

  sum_of_distances = 0.d0

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

      if (iloop_on_min_face_then_max_face == 1) then ! min face

! detect elements belonging to the min face
    limit = value_min + value_size * SMALL_RELATIVE_VALUE

! test face 1 (bottom)
    if (coord_to_use1(i1) < limit .and. coord_to_use1(i2) < limit .and. &
       coord_to_use1(i3) < limit .and. coord_to_use1(i4) < limit) then
      count_elem_faces_to_extend = count_elem_faces_to_extend + 1
      sum_of_distances = sum_of_distances + &
          sqrt((coord_to_use3(i1) - coord_to_use3(i2))**2 + (coord_to_use2(i1) - coord_to_use2(i2))**2)
    endif

! test face 2 (top)
    if (coord_to_use1(i5) < limit .and. coord_to_use1(i6) < limit .and. &
       coord_to_use1(i7) < limit .and. coord_to_use1(i8) < limit) then
      count_elem_faces_to_extend = count_elem_faces_to_extend + 1
      sum_of_distances = sum_of_distances + &
          sqrt((coord_to_use3(i5) - coord_to_use3(i6))**2 + (coord_to_use2(i5) - coord_to_use2(i6))**2)
    endif

! test face 3 (left)
    if (coord_to_use1(i1) < limit .and. coord_to_use1(i4) < limit .and. &
       coord_to_use1(i8) < limit .and. coord_to_use1(i5) < limit) then
      count_elem_faces_to_extend = count_elem_faces_to_extend + 1
      sum_of_distances = sum_of_distances + &
          sqrt((coord_to_use3(i1) - coord_to_use3(i4))**2 + (coord_to_use2(i1) - coord_to_use2(i4))**2)
    endif

! test face 4 (right)
    if (coord_to_use1(i2) < limit .and. coord_to_use1(i3) < limit .and. &
       coord_to_use1(i7) < limit .and. coord_to_use1(i6) < limit) then
      count_elem_faces_to_extend = count_elem_faces_to_extend + 1
      sum_of_distances = sum_of_distances + &
          sqrt((coord_to_use3(i2) - coord_to_use3(i3))**2 + (coord_to_use2(i2) - coord_to_use2(i3))**2)
    endif

! test face 5 (front)
    if (coord_to_use1(i1) < limit .and. coord_to_use1(i2) < limit .and. &
       coord_to_use1(i6) < limit .and. coord_to_use1(i5) < limit) then
      count_elem_faces_to_extend = count_elem_faces_to_extend + 1
      sum_of_distances = sum_of_distances + &
          sqrt((coord_to_use3(i1) - coord_to_use3(i2))**2 + (coord_to_use2(i1) - coord_to_use2(i2))**2)
    endif

! test face 6 (back)
    if (coord_to_use1(i4) < limit .and. coord_to_use1(i3) < limit .and. &
       coord_to_use1(i7) < limit .and. coord_to_use1(i8) < limit) then
      count_elem_faces_to_extend = count_elem_faces_to_extend + 1
      sum_of_distances = sum_of_distances + &
          sqrt((coord_to_use3(i4) - coord_to_use3(i3))**2 + (coord_to_use2(i4) - coord_to_use2(i3))**2)
    endif

      else ! max face

! detect elements belonging to the max face
    limit = value_max - value_size * SMALL_RELATIVE_VALUE

! test face 1 (bottom)
    if (coord_to_use1(i1) > limit .and. coord_to_use1(i2) > limit .and. &
       coord_to_use1(i3) > limit .and. coord_to_use1(i4) > limit) then
      count_elem_faces_to_extend = count_elem_faces_to_extend + 1
      sum_of_distances = sum_of_distances + &
          sqrt((coord_to_use3(i1) - coord_to_use3(i2))**2 + (coord_to_use2(i1) - coord_to_use2(i2))**2)
    endif

! test face 2 (top)
    if (coord_to_use1(i5) > limit .and. coord_to_use1(i6) > limit .and. &
       coord_to_use1(i7) > limit .and. coord_to_use1(i8) > limit) then
      count_elem_faces_to_extend = count_elem_faces_to_extend + 1
      sum_of_distances = sum_of_distances + &
          sqrt((coord_to_use3(i5) - coord_to_use3(i6))**2 + (coord_to_use2(i5) - coord_to_use2(i6))**2)
    endif

! test face 3 (left)
    if (coord_to_use1(i1) > limit .and. coord_to_use1(i4) > limit .and. &
       coord_to_use1(i8) > limit .and. coord_to_use1(i5) > limit) then
      count_elem_faces_to_extend = count_elem_faces_to_extend + 1
      sum_of_distances = sum_of_distances + &
          sqrt((coord_to_use3(i1) - coord_to_use3(i4))**2 + (coord_to_use2(i1) - coord_to_use2(i4))**2)
    endif

! test face 4 (right)
    if (coord_to_use1(i2) > limit .and. coord_to_use1(i3) > limit .and. &
       coord_to_use1(i7) > limit .and. coord_to_use1(i6) > limit) then
      count_elem_faces_to_extend = count_elem_faces_to_extend + 1
      sum_of_distances = sum_of_distances + &
          sqrt((coord_to_use3(i2) - coord_to_use3(i3))**2 + (coord_to_use2(i2) - coord_to_use2(i3))**2)
    endif

! test face 5 (front)
    if (coord_to_use1(i1) > limit .and. coord_to_use1(i2) > limit .and. &
       coord_to_use1(i6) > limit .and. coord_to_use1(i5) > limit) then
      count_elem_faces_to_extend = count_elem_faces_to_extend + 1
      sum_of_distances = sum_of_distances + &
          sqrt((coord_to_use3(i1) - coord_to_use3(i2))**2 + (coord_to_use2(i1) - coord_to_use2(i2))**2)
    endif

! test face 6 (back)
    if (coord_to_use1(i4) > limit .and. coord_to_use1(i3) > limit .and. &
       coord_to_use1(i7) > limit .and. coord_to_use1(i8) > limit) then
      count_elem_faces_to_extend = count_elem_faces_to_extend + 1
      sum_of_distances = sum_of_distances + &
          sqrt((coord_to_use3(i4) - coord_to_use3(i3))**2 + (coord_to_use2(i4) - coord_to_use2(i3))**2)
    endif

      endif

  enddo

  print *,'Total number of elements in the mesh before extension = ',nspec
  print *,'Number of element faces to extend  = ',count_elem_faces_to_extend
  if (count_elem_faces_to_extend == 0) stop 'error: number of element faces to extend detected is zero!'
! we will add NUMBER_OF_PML_LAYERS_TO_ADD to each of the element faces detected that need to be extended
  nspec_new = nspec + count_elem_faces_to_extend * NUMBER_OF_PML_LAYERS_TO_ADD
! and each of these elements will have NGNOD points
! (some of them shared with other elements, but we do not care because they will be removed automatically by xdecompose_mesh)
  npoin_new = npoin + count_elem_faces_to_extend * NUMBER_OF_PML_LAYERS_TO_ADD * NGNOD
  print *,'Total number of elements in the mesh after extension = ',nspec_new
  if (icompute_size == 1) then
    mean_distance = sum_of_distances / dble(count_elem_faces_to_extend)
    print *,'Computed mean size of the elements to extend = ',mean_distance
  endif
  print *

! allocate a new set of elements, i.e. a new ibool()
! allocate arrays with the right size of the future extended mesh
  allocate(imaterial_new(nspec_new))
  allocate(ibool_new(NGNOD,nspec_new))

! clear the arrays
  imaterial_new(:) = 0
  ibool_new(:,:) = 0

! copy the original arrays into the first part i.e. the non-extended part of the new ones
  ibool_new(:,1:nspec) = ibool(:,1:nspec)
  imaterial_new(1:nspec) = imaterial(1:nspec)

  if (minval(ibool) /= 1) stop 'error in minval(ibool)'

! allocate a new set of points, with multiples
  allocate(x_new(npoin_new))
  allocate(y_new(npoin_new))
  allocate(z_new(npoin_new))

! copy the original points into the new set
  x_new(1:npoin) = x(1:npoin)
  y_new(1:npoin) = y(1:npoin)
  z_new(1:npoin) = z(1:npoin)

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

! reset flag
    need_to_extend_this_element = .false.

      if (iloop_on_min_face_then_max_face == 1) then ! min face

! detect elements belonging to the min face
    limit = value_min + value_size * SMALL_RELATIVE_VALUE

! test face 1 (bottom)
    if (coord_to_use1(i1) < limit .and. coord_to_use1(i2) < limit .and. &
       coord_to_use1(i3) < limit .and. coord_to_use1(i4) < limit) then
      need_to_extend_this_element = .true.
      p1 = i1
      p2 = i2
      p3 = i3
      p4 = i4
    endif

! test face 2 (top)
    if (coord_to_use1(i5) < limit .and. coord_to_use1(i6) < limit .and. &
       coord_to_use1(i7) < limit .and. coord_to_use1(i8) < limit) then
      need_to_extend_this_element = .true.
      p1 = i5
      p2 = i6
      p3 = i7
      p4 = i8
    endif

! test face 3 (left)
    if (coord_to_use1(i1) < limit .and. coord_to_use1(i4) < limit .and. &
       coord_to_use1(i5) < limit .and. coord_to_use1(i8) < limit) then
      need_to_extend_this_element = .true.
      p1 = i1
      p2 = i4
      p3 = i8
      p4 = i5
    endif

! test face 4 (right)
    if (coord_to_use1(i2) < limit .and. coord_to_use1(i3) < limit .and. &
       coord_to_use1(i7) < limit .and. coord_to_use1(i6) < limit) then
      need_to_extend_this_element = .true.
      p1 = i2
      p2 = i3
      p3 = i7
      p4 = i6
    endif

! test face 5 (front)
    if (coord_to_use1(i1) < limit .and. coord_to_use1(i2) < limit .and. &
       coord_to_use1(i6) < limit .and. coord_to_use1(i5) < limit) then
      need_to_extend_this_element = .true.
      p1 = i1
      p2 = i2
      p3 = i6
      p4 = i5
    endif

! test face 6 (back)
    if (coord_to_use1(i4) < limit .and. coord_to_use1(i3) < limit .and. &
       coord_to_use1(i7) < limit .and. coord_to_use1(i8) < limit) then
      need_to_extend_this_element = .true.
      p1 = i4
      p2 = i3
      p3 = i7
      p4 = i8
    endif

      else ! max face

! detect elements belonging to the max face
    limit = value_max - value_size * SMALL_RELATIVE_VALUE

! test face 1 (bottom)
    if (coord_to_use1(i1) > limit .and. coord_to_use1(i2) > limit .and. &
       coord_to_use1(i3) > limit .and. coord_to_use1(i4) > limit) then
      need_to_extend_this_element = .true.
      p1 = i1
      p2 = i2
      p3 = i3
      p4 = i4
    endif

! test face 2 (top)
    if (coord_to_use1(i5) > limit .and. coord_to_use1(i6) > limit .and. &
       coord_to_use1(i7) > limit .and. coord_to_use1(i8) > limit) then
      need_to_extend_this_element = .true.
      p1 = i5
      p2 = i6
      p3 = i7
      p4 = i8
    endif

! test face 3 (left)
    if (coord_to_use1(i1) > limit .and. coord_to_use1(i4) > limit .and. &
       coord_to_use1(i5) > limit .and. coord_to_use1(i8) > limit) then
      need_to_extend_this_element = .true.
      p1 = i1
      p2 = i4
      p3 = i8
      p4 = i5
    endif

! test face 4 (right)
    if (coord_to_use1(i2) > limit .and. coord_to_use1(i3) > limit .and. &
       coord_to_use1(i7) > limit .and. coord_to_use1(i6) > limit) then
      need_to_extend_this_element = .true.
      p1 = i2
      p2 = i3
      p3 = i7
      p4 = i6
    endif

! test face 5 (front)
    if (coord_to_use1(i1) > limit .and. coord_to_use1(i2) > limit .and. &
       coord_to_use1(i6) > limit .and. coord_to_use1(i5) > limit) then
      need_to_extend_this_element = .true.
      p1 = i1
      p2 = i2
      p3 = i6
      p4 = i5
    endif

! test face 6 (back)
    if (coord_to_use1(i4) > limit .and. coord_to_use1(i3) > limit .and. &
       coord_to_use1(i7) > limit .and. coord_to_use1(i8) > limit) then
      need_to_extend_this_element = .true.
      p1 = i4
      p2 = i3
      p3 = i7
      p4 = i8
    endif

      endif

    if (need_to_extend_this_element) then

! create the NUMBER_OF_PML_LAYERS_TO_ADD new elements

! very important remark: it is OK to create duplicates of the mesh points in the loop below (i.e. not to tell the code
! that many of these points created are in fact shared between adjacent elements) because "xdecompose_mesh" will remove
! them automatically later on, thus no need to remove them here; this makes this PML mesh extrusion code much simpler to write.

    factor_x = 0
    factor_y = 0
    factor_z = 0

    SIZE_OF_X_ELEMENT_TO_ADD = 0.d0
    SIZE_OF_Y_ELEMENT_TO_ADD = 0.d0
    SIZE_OF_Z_ELEMENT_TO_ADD = 0.d0

    if (iloop_on_X_Y_Z_faces == 1) then  ! Xmin or Xmax
      if (iloop_on_min_face_then_max_face == 1) then ! min face
        factor_x = -1
        if (icompute_size == 1) SIZE_OF_XMIN_ELEMENT_TO_ADD = mean_distance
        SIZE_OF_X_ELEMENT_TO_ADD = SIZE_OF_XMIN_ELEMENT_TO_ADD
      else ! max face
        factor_x = +1
        if (icompute_size == 1) SIZE_OF_XMAX_ELEMENT_TO_ADD = mean_distance
        SIZE_OF_X_ELEMENT_TO_ADD = SIZE_OF_XMAX_ELEMENT_TO_ADD
      endif
    else if (iloop_on_X_Y_Z_faces == 2) then
      if (iloop_on_min_face_then_max_face == 1) then ! min face
        factor_y = -1
        if (icompute_size == 1) SIZE_OF_YMIN_ELEMENT_TO_ADD = mean_distance
        SIZE_OF_Y_ELEMENT_TO_ADD = SIZE_OF_YMIN_ELEMENT_TO_ADD
      else ! max face
        factor_y = +1
        if (icompute_size == 1) SIZE_OF_YMAX_ELEMENT_TO_ADD = mean_distance
        SIZE_OF_Y_ELEMENT_TO_ADD = SIZE_OF_YMAX_ELEMENT_TO_ADD
      endif
    else if (iloop_on_X_Y_Z_faces == 3) then
      if (iloop_on_min_face_then_max_face == 1) then ! min face
        factor_z = -1
        if (icompute_size == 1) SIZE_OF_ZMIN_ELEMENT_TO_ADD = mean_distance
        SIZE_OF_Z_ELEMENT_TO_ADD = SIZE_OF_ZMIN_ELEMENT_TO_ADD
      else ! max face
        factor_z = +1
        if (icompute_size == 1) SIZE_OF_ZMAX_ELEMENT_TO_ADD = mean_distance
        SIZE_OF_Z_ELEMENT_TO_ADD = SIZE_OF_ZMAX_ELEMENT_TO_ADD
      endif
    else
      stop 'wrong index in loop on faces'
    endif

      do iextend = 1,NUMBER_OF_PML_LAYERS_TO_ADD

        ! create a new element
        elem_counter = elem_counter + 1

        ! use the same material property for the extended elements as for the element being extended
        imaterial_new(elem_counter) = imaterial(ispec)

        ! create a new point
        ibool_counter = ibool_counter + 1
        ibool_new(1,elem_counter) = ibool_counter
        x_new(ibool_counter) = x(p1) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*iextend
        y_new(ibool_counter) = y(p1) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*iextend
        z_new(ibool_counter) = z(p1) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*iextend

        ! create a new point
        ibool_counter = ibool_counter + 1
        ibool_new(2,elem_counter) = ibool_counter
        x_new(ibool_counter) = x(p2) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*iextend
        y_new(ibool_counter) = y(p2) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*iextend
        z_new(ibool_counter) = z(p2) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*iextend

        ! create a new point
        ibool_counter = ibool_counter + 1
        ibool_new(3,elem_counter) = ibool_counter
        x_new(ibool_counter) = x(p3) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*iextend
        y_new(ibool_counter) = y(p3) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*iextend
        z_new(ibool_counter) = z(p3) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*iextend

        ! create a new point
        ibool_counter = ibool_counter + 1
        ibool_new(4,elem_counter) = ibool_counter
        x_new(ibool_counter) = x(p4) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*iextend
        y_new(ibool_counter) = y(p4) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*iextend
        z_new(ibool_counter) = z(p4) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*iextend

        ! create a new point
        ibool_counter = ibool_counter + 1
        ibool_new(5,elem_counter) = ibool_counter
        x_new(ibool_counter) = x(p1) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-1)
        y_new(ibool_counter) = y(p1) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-1)
        z_new(ibool_counter) = z(p1) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-1)

        ! create a new point
        ibool_counter = ibool_counter + 1
        ibool_new(6,elem_counter) = ibool_counter
        x_new(ibool_counter) = x(p2) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-1)
        y_new(ibool_counter) = y(p2) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-1)
        z_new(ibool_counter) = z(p2) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-1)

        ! create a new point
        ibool_counter = ibool_counter + 1
        ibool_new(7,elem_counter) = ibool_counter
        x_new(ibool_counter) = x(p3) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-1)
        y_new(ibool_counter) = y(p3) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-1)
        z_new(ibool_counter) = z(p3) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-1)

        ! create a new point
        ibool_counter = ibool_counter + 1
        ibool_new(8,elem_counter) = ibool_counter
        x_new(ibool_counter) = x(p4) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-1)
        y_new(ibool_counter) = y(p4) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-1)
        z_new(ibool_counter) = z(p4) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-1)

! now we need to test if the element created is flipped i.e. it has a negative Jacobian,
! and if so we will use the mirrored version of that element, which will then have a positive Jacobian

! check the element for a negative Jacobian
      xelm(1) = x_new(ibool_new(1,elem_counter))
      xelm(2) = x_new(ibool_new(2,elem_counter))
      xelm(3) = x_new(ibool_new(3,elem_counter))
      xelm(4) = x_new(ibool_new(4,elem_counter))
      xelm(5) = x_new(ibool_new(5,elem_counter))
      xelm(6) = x_new(ibool_new(6,elem_counter))
      xelm(7) = x_new(ibool_new(7,elem_counter))
      xelm(8) = x_new(ibool_new(8,elem_counter))

      yelm(1) = y_new(ibool_new(1,elem_counter))
      yelm(2) = y_new(ibool_new(2,elem_counter))
      yelm(3) = y_new(ibool_new(3,elem_counter))
      yelm(4) = y_new(ibool_new(4,elem_counter))
      yelm(5) = y_new(ibool_new(5,elem_counter))
      yelm(6) = y_new(ibool_new(6,elem_counter))
      yelm(7) = y_new(ibool_new(7,elem_counter))
      yelm(8) = y_new(ibool_new(8,elem_counter))

      zelm(1) = z_new(ibool_new(1,elem_counter))
      zelm(2) = z_new(ibool_new(2,elem_counter))
      zelm(3) = z_new(ibool_new(3,elem_counter))
      zelm(4) = z_new(ibool_new(4,elem_counter))
      zelm(5) = z_new(ibool_new(5,elem_counter))
      zelm(6) = z_new(ibool_new(6,elem_counter))
      zelm(7) = z_new(ibool_new(7,elem_counter))
      zelm(8) = z_new(ibool_new(8,elem_counter))

      call calc_jacobian(xelm,yelm,zelm,dershape3D,found_a_negative_jacobian1,NDIM,NGNOD,NGLLX,NGLLY,NGLLZ)

! check the mirrored (i.e. flipped/swapped) element for a negative Jacobian
! either this one or the non-mirrored one above should be OK, and thus we will select it
      xelm(1) = x_new(ibool_new(5,elem_counter))
      xelm(2) = x_new(ibool_new(6,elem_counter))
      xelm(3) = x_new(ibool_new(7,elem_counter))
      xelm(4) = x_new(ibool_new(8,elem_counter))
      xelm(5) = x_new(ibool_new(1,elem_counter))
      xelm(6) = x_new(ibool_new(2,elem_counter))
      xelm(7) = x_new(ibool_new(3,elem_counter))
      xelm(8) = x_new(ibool_new(4,elem_counter))

      yelm(1) = y_new(ibool_new(5,elem_counter))
      yelm(2) = y_new(ibool_new(6,elem_counter))
      yelm(3) = y_new(ibool_new(7,elem_counter))
      yelm(4) = y_new(ibool_new(8,elem_counter))
      yelm(5) = y_new(ibool_new(1,elem_counter))
      yelm(6) = y_new(ibool_new(2,elem_counter))
      yelm(7) = y_new(ibool_new(3,elem_counter))
      yelm(8) = y_new(ibool_new(4,elem_counter))

      zelm(1) = z_new(ibool_new(5,elem_counter))
      zelm(2) = z_new(ibool_new(6,elem_counter))
      zelm(3) = z_new(ibool_new(7,elem_counter))
      zelm(4) = z_new(ibool_new(8,elem_counter))
      zelm(5) = z_new(ibool_new(1,elem_counter))
      zelm(6) = z_new(ibool_new(2,elem_counter))
      zelm(7) = z_new(ibool_new(3,elem_counter))
      zelm(8) = z_new(ibool_new(4,elem_counter))

      call calc_jacobian(xelm,yelm,zelm,dershape3D,found_a_negative_jacobian2,NDIM,NGNOD,NGLLX,NGLLY,NGLLZ)

! this should never happen, it is just a safety test
      if (found_a_negative_jacobian1 .and. found_a_negative_jacobian2) &
        stop 'error: found a negative Jacobian that could not be fixed'

! this should also never happen, it is just a second safety test
      if (.not. found_a_negative_jacobian1 .and. .not. found_a_negative_jacobian2) &
        stop 'strange error: both the element created and its mirrored version have a positive Jacobian!'

! if we have found that the original element has a negative Jacobian and its mirrored element is fine,
! swap the points so that we use that mirrored element in the final mesh saved to disk instead of the original one
      if (found_a_negative_jacobian1) then
        i1 = ibool_new(5,elem_counter)
        i2 = ibool_new(6,elem_counter)
        i3 = ibool_new(7,elem_counter)
        i4 = ibool_new(8,elem_counter)
        i5 = ibool_new(1,elem_counter)
        i6 = ibool_new(2,elem_counter)
        i7 = ibool_new(3,elem_counter)
        i8 = ibool_new(4,elem_counter)

        ibool_new(1,elem_counter) = i1
        ibool_new(2,elem_counter) = i2
        ibool_new(3,elem_counter) = i3
        ibool_new(4,elem_counter) = i4
        ibool_new(5,elem_counter) = i5
        ibool_new(6,elem_counter) = i6
        ibool_new(7,elem_counter) = i7
        ibool_new(8,elem_counter) = i8
      endif

      enddo
    endif

  enddo

  if (minval(ibool_new) /= 1) stop 'error in minval(ibool_new)'

! deallocate the original arrays
  deallocate(x,y,z)
  deallocate(ibool)
  deallocate(imaterial)

! reallocate them with the new size
  allocate(x(npoin_new))
  allocate(y(npoin_new))
  allocate(z(npoin_new))
  allocate(imaterial(nspec_new))
  allocate(ibool(NGNOD,nspec_new))

! make the new ones become the old ones, to prepare for the next iteration of the two nested loops we are in,
! i.e. to make sure the next loop will extend the mesh from the new arrays rather than from the old ones
  x(:) = x_new(:)
  y(:) = y_new(:)
  z(:) = z_new(:)
  imaterial(:) = imaterial_new(:)
  ibool(:,:) = ibool_new(:,:)

! the new number of elements and points becomes the old one, for the same reason
  nspec = nspec_new
  npoin = npoin_new

! deallocate the new ones, to make sure they can be allocated again in the next iteration of the nested loops we are in
  deallocate(x_new,y_new,z_new)
  deallocate(ibool_new)
  deallocate(imaterial_new)

  if (minval(ibool) /= 1) stop 'error in minval(ibool)'

    enddo ! of iloop_on_min_face_then_max_face loop on Xmin then Xmax, or Ymin then Ymax, or Zmin then Zmax

! end of loop on the three sets of faces to first add CPML elements along X, then along Y, then along Z
  enddo

  if (iformat == 1) then ! write the output in ASCII format

! write the new points (overwrite the old file)
  open(unit=23,file='nodes_coords_file',status='old',action='write')
  write(23,*) npoin
  do ipoin = 1,npoin
    write(23,*) ipoin,sngl(x(ipoin)),sngl(y(ipoin)),sngl(z(ipoin))
  enddo
  close(23)

! write the new mesh elements (overwrite the old file)
  open(unit=23,file='mesh_file',status='old',action='write')
  write(23,*) nspec
! loop on the whole mesh
  do ispec = 1,nspec
    write(23,"(i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9)") ispec,(ibool(ia,ispec), ia = 1,NGNOD)
  enddo
  close(23)

! write the new material properties (overwrite the old file)
  open(unit=23,file='materials_file',status='old',action='write')
! loop on the whole mesh
  do ispec = 1,nspec
    write(23,*) ispec,imaterial(ispec)
  enddo
  close(23)

  else ! write the output in binary format

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

  endif

! output information for the next code (xconvert_external_layers_of_a_given_mesh_to_CPML_layers)
  print *
  print *,'Here are the values to use as input in the next code, xconvert_external_layers_of_a_given_mesh_to_CPML_layers:'
  print *,'  (also saved in file values_to_use_for_convert_external_layers_of_a_given_mesh_to_CPML_layers.txt)'
  print *
  print *,'THICKNESS_OF_XMIN_PML = ',sngl(SIZE_OF_XMIN_ELEMENT_TO_ADD * NUMBER_OF_PML_LAYERS_TO_ADD)
  print *,'THICKNESS_OF_XMAX_PML = ',sngl(SIZE_OF_XMAX_ELEMENT_TO_ADD * NUMBER_OF_PML_LAYERS_TO_ADD)
  print *,'THICKNESS_OF_YMIN_PML = ',sngl(SIZE_OF_YMIN_ELEMENT_TO_ADD * NUMBER_OF_PML_LAYERS_TO_ADD)
  print *,'THICKNESS_OF_YMAX_PML = ',sngl(SIZE_OF_YMAX_ELEMENT_TO_ADD * NUMBER_OF_PML_LAYERS_TO_ADD)
  print *,'THICKNESS_OF_ZMIN_PML = ',sngl(SIZE_OF_ZMIN_ELEMENT_TO_ADD * NUMBER_OF_PML_LAYERS_TO_ADD)
  if (ALSO_ADD_ON_THE_TOP_SURFACE) &
      print *,'THICKNESS_OF_ZMAX_PML = ',sngl(SIZE_OF_ZMAX_ELEMENT_TO_ADD * NUMBER_OF_PML_LAYERS_TO_ADD)
  print *

! save the thickness values to a text file
  open(unit=23,file='values_to_use_for_convert_external_layers_of_a_given_mesh_to_CPML_layers.txt',status='unknown',action='write')
  write(23,*) SIZE_OF_XMIN_ELEMENT_TO_ADD * NUMBER_OF_PML_LAYERS_TO_ADD
  write(23,*) SIZE_OF_XMAX_ELEMENT_TO_ADD * NUMBER_OF_PML_LAYERS_TO_ADD
  write(23,*) SIZE_OF_YMIN_ELEMENT_TO_ADD * NUMBER_OF_PML_LAYERS_TO_ADD
  write(23,*) SIZE_OF_YMAX_ELEMENT_TO_ADD * NUMBER_OF_PML_LAYERS_TO_ADD
  write(23,*) SIZE_OF_ZMIN_ELEMENT_TO_ADD * NUMBER_OF_PML_LAYERS_TO_ADD
  if (ALSO_ADD_ON_THE_TOP_SURFACE) then
    write(23,*) SIZE_OF_ZMAX_ELEMENT_TO_ADD * NUMBER_OF_PML_LAYERS_TO_ADD
  else
! convention (negative value) to say that this Zmax absorbing edge is turned off
    write(23,*) -1
  endif
  close(23)

  end program add_CPML_layers_to_a_given_mesh

!
!=====================================================================
!

  subroutine calc_jacobian(xelm,yelm,zelm,dershape3D,found_a_negative_jacobian,NDIM,NGNOD,NGLLX,NGLLY,NGLLZ)

  implicit none

  integer :: NDIM,NGNOD,NGLLX,NGLLY,NGLLZ

  logical :: found_a_negative_jacobian

  double precision, dimension(NGNOD) :: xelm,yelm,zelm
  double precision dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ)

  integer i,j,k,ia
  double precision xxi,xeta,xgamma,yxi,yeta,ygamma,zxi,zeta,zgamma
  double precision jacobian

  double precision, parameter :: ZERO = 0.d0

  found_a_negative_jacobian = .false.

! do k=1,NGLLZ
!   do j=1,NGLLY
!     do i=1,NGLLX
! for this CPML mesh extrusion routine it is sufficient to test the 8 corners of each element to reduce the cost
! because we just want to detect if the element is flipped or not, and if so flip it back
! do k=1,NGLLZ,NGLLZ-1
!   do j=1,NGLLY,NGLLY-1
!     do i=1,NGLLX,NGLLX-1
! it is even sufficient to test a single corner
  do k=1,1
    do j=1,1
      do i=1,1

      xxi = ZERO
      xeta = ZERO
      xgamma = ZERO
      yxi = ZERO
      yeta = ZERO
      ygamma = ZERO
      zxi = ZERO
      zeta = ZERO
      zgamma = ZERO

      do ia=1,NGNOD
        xxi = xxi + dershape3D(1,ia,i,j,k)*xelm(ia)
        xeta = xeta + dershape3D(2,ia,i,j,k)*xelm(ia)
        xgamma = xgamma + dershape3D(3,ia,i,j,k)*xelm(ia)

        yxi = yxi + dershape3D(1,ia,i,j,k)*yelm(ia)
        yeta = yeta + dershape3D(2,ia,i,j,k)*yelm(ia)
        ygamma = ygamma + dershape3D(3,ia,i,j,k)*yelm(ia)

        zxi = zxi + dershape3D(1,ia,i,j,k)*zelm(ia)
        zeta = zeta + dershape3D(2,ia,i,j,k)*zelm(ia)
        zgamma = zgamma + dershape3D(3,ia,i,j,k)*zelm(ia)
      enddo

      jacobian = xxi*(yeta*zgamma-ygamma*zeta) - xeta*(yxi*zgamma-ygamma*zxi) + xgamma*(yxi*zeta-yeta*zxi)

! check that the Jacobian transform is invertible, i.e. that the Jacobian never becomes negative or null
      if (jacobian <= ZERO) found_a_negative_jacobian = .true.

      enddo
    enddo
  enddo

  end subroutine calc_jacobian

!
!=====================================================================
!

! 3D shape functions for 8-node element

  subroutine get_shape3D(dershape3D,xigll,yigll,zigll,NGNOD,NGLLX,NGLLY,NGLLZ,NDIM)

  implicit none

  integer NGNOD,NGLLX,NGLLY,NGLLZ,NDIM

! Gauss-Lobatto-Legendre points of integration
  double precision xigll(NGLLX)
  double precision yigll(NGLLY)
  double precision zigll(NGLLZ)

! 3D shape functions and their derivatives
  double precision dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ)

  integer i,j,k

! location of the nodes of the 3D hexahedra elements
  double precision xi,eta,gamma
  double precision ra1,ra2,rb1,rb2,rc1,rc2

  double precision, parameter :: ONE = 1.d0
  double precision, parameter :: ONE_EIGHTH = 0.125d0

! check that the parameter file is correct
  if (NGNOD /= 8) stop 'volume elements should have 8 control nodes'

! ***
! *** create 3D shape functions and jacobian
! ***

  do i=1,NGLLX
    do j=1,NGLLY
      do k=1,NGLLZ

        xi = xigll(i)
        eta = yigll(j)
        gamma = zigll(k)

          ra1 = ONE + xi
          ra2 = ONE - xi

          rb1 = ONE + eta
          rb2 = ONE - eta

          rc1 = ONE + gamma
          rc2 = ONE - gamma

          dershape3D(1,1,i,j,k) = - ONE_EIGHTH*rb2*rc2
          dershape3D(1,2,i,j,k) = ONE_EIGHTH*rb2*rc2
          dershape3D(1,3,i,j,k) = ONE_EIGHTH*rb1*rc2
          dershape3D(1,4,i,j,k) = - ONE_EIGHTH*rb1*rc2
          dershape3D(1,5,i,j,k) = - ONE_EIGHTH*rb2*rc1
          dershape3D(1,6,i,j,k) = ONE_EIGHTH*rb2*rc1
          dershape3D(1,7,i,j,k) = ONE_EIGHTH*rb1*rc1
          dershape3D(1,8,i,j,k) = - ONE_EIGHTH*rb1*rc1

          dershape3D(2,1,i,j,k) = - ONE_EIGHTH*ra2*rc2
          dershape3D(2,2,i,j,k) = - ONE_EIGHTH*ra1*rc2
          dershape3D(2,3,i,j,k) = ONE_EIGHTH*ra1*rc2
          dershape3D(2,4,i,j,k) = ONE_EIGHTH*ra2*rc2
          dershape3D(2,5,i,j,k) = - ONE_EIGHTH*ra2*rc1
          dershape3D(2,6,i,j,k) = - ONE_EIGHTH*ra1*rc1
          dershape3D(2,7,i,j,k) = ONE_EIGHTH*ra1*rc1
          dershape3D(2,8,i,j,k) = ONE_EIGHTH*ra2*rc1

          dershape3D(3,1,i,j,k) = - ONE_EIGHTH*ra2*rb2
          dershape3D(3,2,i,j,k) = - ONE_EIGHTH*ra1*rb2
          dershape3D(3,3,i,j,k) = - ONE_EIGHTH*ra1*rb1
          dershape3D(3,4,i,j,k) = - ONE_EIGHTH*ra2*rb1
          dershape3D(3,5,i,j,k) = ONE_EIGHTH*ra2*rb2
          dershape3D(3,6,i,j,k) = ONE_EIGHTH*ra1*rb2
          dershape3D(3,7,i,j,k) = ONE_EIGHTH*ra1*rb1
          dershape3D(3,8,i,j,k) = ONE_EIGHTH*ra2*rb1

      enddo
    enddo
  enddo

  end subroutine get_shape3D

