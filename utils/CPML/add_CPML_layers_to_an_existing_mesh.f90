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
! The full text of the license is available in file "LICENSE".
!
!========================================================================

  program add_CPML_layers_to_a_given_mesh

! Dimitri Komatitsch, CNRS Marseille, France, March 2016 and February 2017

! add PML layers around an existing mesh (i.e. create new elements and new points)

  implicit none

! we are in 3D
  integer, parameter :: NDIM = 3

! number of GLL points in each direction, to check for negative Jacobians
  integer, parameter :: NGLLX = 5,NGLLY = NGLLX,NGLLZ = NGLLX

! number of PML and non-PML layers to add on each side of the mesh
  integer :: NUMBER_OF_PML_LAYERS_TO_ADD,NUMBER_OF_TRANSITION_LAYERS_TO_ADD,NUMBER_OF_LAYERS_TO_ADD_IN_THIS_STEP

! size of each PML element to add when they are added on the Xmin and Xmax faces, Ymin and Ymax faces, Zmin and/or Zmax faces
  double precision :: SIZE_OF_X_ELEMENT_TO_ADD,SIZE_OF_Y_ELEMENT_TO_ADD,SIZE_OF_Z_ELEMENT_TO_ADD
  double precision :: SIZE_OF_XMIN_ELEMENT_TO_ADD,SIZE_OF_YMIN_ELEMENT_TO_ADD,SIZE_OF_ZMIN_ELEMENT_TO_ADD
  double precision :: SIZE_OF_XMAX_ELEMENT_TO_ADD,SIZE_OF_YMAX_ELEMENT_TO_ADD,SIZE_OF_ZMAX_ELEMENT_TO_ADD

  integer :: nspec,npoin,npoin_new_max,npoin_new_real,npointot,nspec_new,count_elem_faces_to_extend,iextend
  integer :: factor_x,factor_y,factor_z
  integer :: ispec,ipoin,iloop_on_X_Y_Z_faces,iloop_on_min_face_then_max_face
  integer :: ipoin_read,ispec_loop
  integer :: i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,i21,i22,i23,i24,i25,i26,i27
  integer :: p1,p2,p3,p4,p5,p6,p7,p8,p9
  integer :: elem_counter,ia,iflag,iformat,icompute_size,istep,itype_of_hex,NGNOD,icheck_for_negative_Jacobians

  double precision, dimension(:), allocatable, target :: x,y,z
  double precision, dimension(:), allocatable :: x_new,y_new,z_new,xp,yp,zp
  double precision, dimension(:,:), allocatable :: x_copy,y_copy,z_copy
  double precision, dimension(:), pointer :: coord_to_use1,coord_to_use2,coord_to_use3

  integer, dimension(:), allocatable :: imaterial,imaterial_new,locval

  logical, dimension(:), allocatable :: ifseg

  integer :: ieoff,ilocnum,iglobnum

  integer, dimension(:,:), allocatable :: ibool,ibool_new

  integer :: count_def_mat,count_undef_mat,ier,idummy,num_mat,aniso_flag,idomain_id,imaterial_values_to_start_from

  real(kind=4) :: vp,vs,rho,qkappa,qmu

  character(len=256) :: line

! Gauss-Lobatto-Legendre points of integration, to check for negative Jacobians
  double precision xigll(NGLLX)
  double precision yigll(NGLLY)
  double precision zigll(NGLLZ)

! 3D shape function derivatives, to check for negative Jacobians
  double precision, dimension(:,:,:,:,:), allocatable :: dershape3D

  double precision, dimension(:), allocatable :: xelm,yelm,zelm

  double precision :: jacobian

  double precision, dimension(:), allocatable :: x_tmp1,y_tmp1,z_tmp1,x_tmp2,y_tmp2,z_tmp2

  double precision :: xmin,xmax,ymin,ymax,zmin,zmax,limit,xsize,ysize,zsize
  double precision :: value_min,value_max,value_size,sum_of_distances,mean_distance,very_small_distance

  logical :: ADD_ON_THE_XMIN_SURFACE,ADD_ON_THE_XMAX_SURFACE
  logical :: ADD_ON_THE_YMIN_SURFACE,ADD_ON_THE_YMAX_SURFACE
  logical :: ADD_ON_THE_ZMIN_SURFACE,ADD_ON_THE_ZMAX_SURFACE

  logical :: need_to_extend_this_element,found_a_negative_Jacobian1,found_a_negative_Jacobian2

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
  print *

  print *,'1 = the mesh contains HEX8 elements'
  print *,'2 = the mesh contains HEX27 elements'
  print *,'3 = exit'
  read(*,*) itype_of_hex
  if (itype_of_hex /= 1 .and. itype_of_hex /= 2) stop 'exiting...'
  if (itype_of_hex == 1) then
    NGNOD = 8
  else
    NGNOD = 27
  endif
  print *

  allocate(x_tmp1(NGNOD))
  allocate(y_tmp1(NGNOD))
  allocate(z_tmp1(NGNOD))
  allocate(x_tmp2(NGNOD))
  allocate(y_tmp2(NGNOD))
  allocate(z_tmp2(NGNOD))

  print *,'enter the number of PML layers to add on each side of the mesh (usually 3, can also be 4):'
  read(*,*) NUMBER_OF_PML_LAYERS_TO_ADD
  if (NUMBER_OF_PML_LAYERS_TO_ADD < 1) stop 'NUMBER_OF_PML_LAYERS_TO_ADD must be >= 1'
  print *

  print *,'enter the number of non-PML transition layers to add between each side of the mesh and the new PMLs &
                &that will be created (usually 0, 1 or 2) (if you have no idea, just enter 0):'
  read(*,*) NUMBER_OF_TRANSITION_LAYERS_TO_ADD
  if (NUMBER_OF_TRANSITION_LAYERS_TO_ADD < 0) stop 'NUMBER_OF_TRANSITION_LAYERS_TO_ADD must be >= 0'
  print *

  ADD_ON_THE_XMIN_SURFACE = .true.
  ADD_ON_THE_XMAX_SURFACE = .true.
  ADD_ON_THE_YMIN_SURFACE = .true.
  ADD_ON_THE_YMAX_SURFACE = .true.
  ADD_ON_THE_ZMIN_SURFACE = .true.
  ADD_ON_THE_ZMAX_SURFACE = .true.

  print *,'1 = use a free surface at the top of the mesh i.e. do not add a CPML layer at the top (most classical option)'
  print *,'2 = add a CPML absorbing layer at the top of the mesh (less classical option)'
  print *,'3 = exit'
  read(*,*) iflag
  if (iflag /= 1 .and. iflag /= 2) stop 'exiting...'
  if (iflag == 1) then
    ADD_ON_THE_ZMAX_SURFACE = .false.
  else
    ADD_ON_THE_ZMAX_SURFACE = .true.
  endif
  print *

  print *,'1 = compute the size of the PML elements to add automatically using the average size of edge elements'
  print *,'2 = enter the size of the PML elements to add manually'
  print *,'3 = exit'
  read(*,*) icompute_size
  if (icompute_size /= 1 .and. icompute_size /= 2) stop 'exiting...'

  if (icompute_size == 2) then

  print *,'enter the X size (in meters) of each CPML element to add on the Xmin face (enter -1 to turn that PML layer off):'
  read(*,*) SIZE_OF_XMIN_ELEMENT_TO_ADD
  if (SIZE_OF_XMIN_ELEMENT_TO_ADD <= 0.d0) then
    SIZE_OF_XMIN_ELEMENT_TO_ADD = 0.d0
    ADD_ON_THE_XMIN_SURFACE = .false.
  endif
  print *

  print *,'enter the X size (in meters) of each CPML element to add on the Xmax face (enter -1 to turn that PML layer off):'
  read(*,*) SIZE_OF_XMAX_ELEMENT_TO_ADD
  if (SIZE_OF_XMAX_ELEMENT_TO_ADD <= 0.d0) then
    SIZE_OF_XMAX_ELEMENT_TO_ADD = 0.d0
    ADD_ON_THE_XMAX_SURFACE = .false.
  endif
  print *

  print *,'enter the Y size (in meters) of each CPML element to add on the Ymin faces (enter -1 to turn that PML layer off):'
  read(*,*) SIZE_OF_YMIN_ELEMENT_TO_ADD
  if (SIZE_OF_YMIN_ELEMENT_TO_ADD <= 0.d0) then
    SIZE_OF_YMIN_ELEMENT_TO_ADD = 0.d0
    ADD_ON_THE_YMIN_SURFACE = .false.
  endif
  print *

  print *,'enter the Y size (in meters) of each CPML element to add on the Ymax faces (enter -1 to turn that PML layer off):'
  read(*,*) SIZE_OF_YMAX_ELEMENT_TO_ADD
  if (SIZE_OF_YMAX_ELEMENT_TO_ADD <= 0.d0) then
    SIZE_OF_YMAX_ELEMENT_TO_ADD = 0.d0
    ADD_ON_THE_YMAX_SURFACE = .false.
  endif
  print *

  print *,'enter the Z size (in meters) of each CPML element to add on the Zmin faces (enter -1 to turn that PML layer off):'
  read(*,*) SIZE_OF_ZMIN_ELEMENT_TO_ADD
  if (SIZE_OF_ZMIN_ELEMENT_TO_ADD <= 0.d0) then
    SIZE_OF_ZMIN_ELEMENT_TO_ADD = 0.d0
    ADD_ON_THE_ZMIN_SURFACE = .false.
  endif
  print *

  if (ADD_ON_THE_ZMAX_SURFACE) then
    print *,'enter the Z size (in meters) of each CPML element to add on the Zmax faces (enter -1 to turn that PML layer off):'
    read(*,*) SIZE_OF_ZMAX_ELEMENT_TO_ADD
    if (SIZE_OF_ZMAX_ELEMENT_TO_ADD <= 0.d0) then
      SIZE_OF_ZMAX_ELEMENT_TO_ADD = 0.d0
      ADD_ON_THE_ZMAX_SURFACE = .false.
    endif
    print *
  endif

  endif

! check that we need to add at least one PML, otherwise this code is useless
  if (.not. ADD_ON_THE_XMIN_SURFACE .and. .not. ADD_ON_THE_XMAX_SURFACE &
      .and. .not. ADD_ON_THE_YMIN_SURFACE .and. .not. ADD_ON_THE_YMAX_SURFACE &
      .and. .not. ADD_ON_THE_ZMIN_SURFACE .and. .not. ADD_ON_THE_ZMAX_SURFACE) &
    stop 'Error: the purpose of this code is to add at least one PML, but you have added none'

  print *,'1 = do not check the input mesh for negative Jacobians (normally safe because it should contain none)'
  print *,'2 = check the input mesh for negative Jacobians (does not hurt, usually takes less than a minute)'
  print *,'3 = exit'
  read(*,*) icheck_for_negative_Jacobians
  if (icheck_for_negative_Jacobians /= 1 .and. icheck_for_negative_Jacobians /= 2) stop 'exiting...'
  print *

! hardwire GLL point location values to avoid having to link with a long library to compute them
  xigll(:) = (/ -1.d0 , -0.654653670707977d0 , 0.d0 , 0.654653670707977d0 , 1.d0 /)
  yigll(:) = xigll(:)
  zigll(:) = xigll(:)

  allocate(dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ))

! compute the derivatives of the 3D shape functions for a 8-node or 27-node element
  call get_shape3D(dershape3D,xigll,yigll,zigll,NGNOD,NGLLX,NGLLY,NGLLZ,NDIM)

  allocate(xelm(NGNOD))
  allocate(yelm(NGNOD))
  allocate(zelm(NGNOD))

! The format of nummaterial_velocity_file must be:

! #(1)material_domain_id #(2)material_id  #(3)rho  #(4)vp   #(5)vs   #(6)Q_kappa   #(7)Q_mu  #(8)anisotropy_flag
!
! where
!     material_domain_id : 1=acoustic / 2=elastic / 3=poroelastic
!     material_id        : POSITIVE integer identifier corresponding to the identifier of material block
!     rho                : density
!     vp                 : P-velocity
!     vs                 : S-velocity
!     Q_kappa            : 9999 = no Q_kappa attenuation
!     Q_mu               : 9999 = no Q_mu attenuation
!     anisotropy_flag    : 0=no anisotropy/ 1,2,.. check with implementation in aniso_model.f90
!
! example:
!  2   1 2300 2800 1500 9999.0 9999.0 0

! Note that when poroelastic material, this file is a dummy except for material_domain_id & material_id,
! and that poroelastic materials are actually read from nummaterial_poroelastic_file, because CUBIT
! cannot support more than 10 attributes

!or

! #(1)material_domain_id #(2)material_id  tomography elastic  #(3)tomography_filename #(4)positive_unique_number
!
! where
!     material_domain_id : 1=acoustic / 2=elastic / 3=poroelastic
!     material_id        : NEGATIVE integer identifier corresponding to the identifier of material block
!     tomography_filename: filename of the tomography file
!     positive_unique_number: a positive unique identifier
!
! example:
!  2  -1 tomography elastic tomo.xyz 1

! read the file a first time to count the total number of materials (both defined and tomographic)
  count_def_mat = 0
  count_undef_mat = 0
  open(unit=98, file='nummaterial_velocity_file', status='old', action='read',form='formatted',iostat=ier)
  if (ier /= 0) stop 'Error opening nummaterial_velocity_file'

  ! counts materials (defined/undefined)
  ier = 0
  do while (ier == 0)
    ! note: format #material_domain_id #material_id #...
    read(98,'(A)',iostat=ier) line
    if (ier /= 0) exit

    ! skip empty/comment lines
    if (len_trim(line) == 0) cycle
    if (line(1:1) == '#' .or. line(1:1) == '!') cycle

    read(line,*,iostat=ier) idummy,num_mat
    if (ier /= 0) exit

    ! checks non-zero material id
    if (num_mat == 0) stop 'Error in nummaterial_velocity_file: material id 0 found. Material ids must be non-zero.'
    if (num_mat > 0) then
      ! positive materials_id: velocity values will be defined
      count_def_mat = count_def_mat + 1
    else
      ! negative materials_id: undefined material properties yet
      count_undef_mat = count_undef_mat + 1
    endif
  enddo
  close(98)

! read the file a second time to create the materials for the extruded PML layers (both defined and tomographic)
  open(unit=98, file='nummaterial_velocity_file', status='old', action='read',form='formatted',iostat=ier)
  if (ier /= 0) stop 'Error opening nummaterial_velocity_file'
  open(unit=99, file='nummaterial_velocity_file_new', status='unknown', action='write',form='formatted',iostat=ier)
  if (ier /= 0) stop 'Error opening nummaterial_velocity_file_new'

  ! counts materials (defined/undefined)
  ier = 0
  do while (ier == 0)
    ! note: format #material_domain_id #material_id #...
    read(98,'(A)',iostat=ier) line
    if (ier /= 0) exit

    ! skip empty/comment lines
    if (len_trim(line) == 0) cycle
    if (line(1:1) == '#' .or. line(1:1) == '!') cycle

    read(line,*,iostat=ier) idummy,num_mat
    if (ier /= 0) exit

    ! checks non-zero material id
    if (num_mat > 0) then
      ! positive materials_id: velocity values will be defined
      ! reads in defined material properties
      read(line,*) idomain_id,num_mat,rho,vp,vs,qkappa,qmu,aniso_flag

      ! sanity check: Q factor cannot be equal to zero, thus convert to 9999 to indicate no attenuation
      ! if users have used 0 to indicate that instead
      if (qkappa <= 0.000001) qkappa = 9999.
      if (qmu <= 0.000001) qmu = 9999.

      ! checks material_id bounds
      if (num_mat < 1 .or. num_mat > count_def_mat) &
        stop "Error in nummaterial_velocity_file: material id invalid for defined materials."

      ! copy the original material to the new file
      write(99,200) idomain_id,num_mat,rho,vp,vs,qkappa,qmu,aniso_flag

      ! create new material for the PML, with attenuation off
      write(99,200) idomain_id,num_mat + count_def_mat,rho,vp,vs,9999.,9999.,aniso_flag

      ! create new material for the transition layer that has PML off (if it exists)
      ! in this transition layer we purposely put less attenuation in order to smooth out the transition with PML,
      ! since currently PMLs have no attenuation
      if (NUMBER_OF_TRANSITION_LAYERS_TO_ADD > 0) &
        write(99,200) idomain_id,num_mat + 2*count_def_mat,rho,vp,vs,min(3.*qkappa,9999.),min(3.*qmu,9999.),aniso_flag

    else
      ! negative materials_id: undefined material properties yet
      ! just copy the line read to the new file
      write(99,*) line(1:len_trim(line))
    endif
  enddo
  close(98)
  close(99)

 200 format(i6,1x,i6,1x,f14.6,1x,f14.6,1x,f14.6,1x,f14.6,1x,f14.6,1x,i6)

! replace the old file with the new one
  call system('mv -v nummaterial_velocity_file_new nummaterial_velocity_file')

! perform two steps: first we create the transition layer (if any), then we create the PML layers
! we do this in two separate steps in order to be able to easily create a transition layer that uses
! intermediate Q values, otherwise if doing it in a single step there are difficult problems in the PML mesh corners.
! The only drawback of using two steps is that the mesh is read from disk and then written to disk twice instead of once.

  do istep = 1,2

    if (istep == 1) then
      ! do nothing in the first step if the number of layers to create is equal to zero
      if (NUMBER_OF_TRANSITION_LAYERS_TO_ADD == 0) cycle
      NUMBER_OF_LAYERS_TO_ADD_IN_THIS_STEP = NUMBER_OF_TRANSITION_LAYERS_TO_ADD
      print *
      print *,'about to create the transition layer, will first need to read the input mesh files...'
      print *
    else
      NUMBER_OF_LAYERS_TO_ADD_IN_THIS_STEP = NUMBER_OF_PML_LAYERS_TO_ADD
      print *
      print *,'about to create the PML layer, will first need to read the input mesh files...'
      print *
    endif

    if (NGNOD == 8) then
      print *,'reading HEX8 input mesh file...'
      print *
    else
      print *,'reading HEX27 input mesh file...'
      print *
    endif

! ************* read mesh points *************

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
        read(23,*) ipoin_read,x(ipoin_read),y(ipoin_read),z(ipoin_read)
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
        read(23,*) ispec,(ibool(ia,ispec), ia = 1,NGNOD)
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
        read(23,*) ispec,imaterial(ispec)
      enddo
    else
      read(23) imaterial
    endif
    close(23)

! check the mesh read to make sure it contains no negative Jacobians
    if (icheck_for_negative_Jacobians == 2) then
      do ispec = 1,nspec
!   check the element for a negative Jacobian
        do ia = 1,NGNOD
          xelm(ia) = x(ibool(ia,ispec))
          yelm(ia) = y(ibool(ia,ispec))
          zelm(ia) = z(ibool(ia,ispec))
        enddo

        call calc_jacobian(xelm,yelm,zelm,dershape3D,found_a_negative_jacobian1,NDIM,NGNOD,NGLLX,NGLLY,NGLLZ,jacobian)
        if (found_a_negative_jacobian1) then
          print *,'detected an element with negative Jacobian in the input mesh: element ',ispec, &
                     ' in which the Jacobian is ',jacobian
          stop 'error: the mesh read contains a negative Jacobian!'
        endif
      enddo
      print *,'input mesh successfully checked for negative Jacobians, it contains none'
      print *
    endif

! we need to extend/extrude the existing mesh by adding CPML elements
! along the X faces, then along the Y faces, then along the Z faces.

! loop on the three sets of faces to first add CPML elements along X, then along Y, then along Z
    do iloop_on_X_Y_Z_faces = 1,NDIM

! 1 is min face and 2 is max face (Xmin then Xmax, Ymin then Ymax, or Zmin then Zmax)
      do iloop_on_min_face_then_max_face = 1,2


! do not add a CPML layer on a given surface if the user asked not to
        if (iloop_on_X_Y_Z_faces == 1 .and. iloop_on_min_face_then_max_face == 1 .and. .not. ADD_ON_THE_XMIN_SURFACE) cycle
        if (iloop_on_X_Y_Z_faces == 1 .and. iloop_on_min_face_then_max_face == 2 .and. .not. ADD_ON_THE_XMAX_SURFACE) cycle
        if (iloop_on_X_Y_Z_faces == 2 .and. iloop_on_min_face_then_max_face == 1 .and. .not. ADD_ON_THE_YMIN_SURFACE) cycle
        if (iloop_on_X_Y_Z_faces == 2 .and. iloop_on_min_face_then_max_face == 2 .and. .not. ADD_ON_THE_YMAX_SURFACE) cycle
        if (iloop_on_X_Y_Z_faces == 3 .and. iloop_on_min_face_then_max_face == 1 .and. .not. ADD_ON_THE_ZMIN_SURFACE) cycle
        if (iloop_on_X_Y_Z_faces == 3 .and. iloop_on_min_face_then_max_face == 2 .and. .not. ADD_ON_THE_ZMAX_SURFACE) cycle

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

! we can use the 8 corners of the element only for the test here, even if the element is HEX27
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

! we will add NUMBER_OF_LAYERS_TO_ADD_IN_THIS_STEP to each of the element faces detected that need to be extended
        nspec_new = nspec + count_elem_faces_to_extend * NUMBER_OF_LAYERS_TO_ADD_IN_THIS_STEP
        print *,'Total number of elements in the mesh after extension = ',nspec_new

! and each of these elements will have NGNOD points
! (some of them shared with other elements, but we will remove the multiples below, thus here it is a maximum
        npoin_new_max = npoin + count_elem_faces_to_extend * NUMBER_OF_LAYERS_TO_ADD_IN_THIS_STEP * NGNOD

        mean_distance = sum_of_distances / dble(count_elem_faces_to_extend)
        very_small_distance = mean_distance / 10000.d0
        if (icompute_size == 1) print *,'Computed mean size of the elements to extend = ',mean_distance
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
        allocate(x_new(npoin_new_max))
        allocate(y_new(npoin_new_max))
        allocate(z_new(npoin_new_max))

! copy the original points into the new set
        x_new(1:npoin) = x(1:npoin)
        y_new(1:npoin) = y(1:npoin)
        z_new(1:npoin) = z(1:npoin)

! position after which to start to create the new elements
        elem_counter = nspec
        npoin_new_real = npoin

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

              if (NGNOD == 27) then
                p5 = i9
                p6 = i10
                p7 = i11
                p8 = i12
                p9 = i21
              endif
            endif

! test face 2 (top)
            if (coord_to_use1(i5) < limit .and. coord_to_use1(i6) < limit .and. &
               coord_to_use1(i7) < limit .and. coord_to_use1(i8) < limit) then
              need_to_extend_this_element = .true.
              p1 = i5
              p2 = i6
              p3 = i7
              p4 = i8

              if (NGNOD == 27) then
                p5 = i17
                p6 = i18
                p7 = i19
                p8 = i20
                p9 = i26
              endif
            endif

! test face 3 (left)
            if (coord_to_use1(i1) < limit .and. coord_to_use1(i4) < limit .and. &
               coord_to_use1(i5) < limit .and. coord_to_use1(i8) < limit) then
              need_to_extend_this_element = .true.
              p1 = i1
              p2 = i4
              p3 = i8
              p4 = i5

              if (NGNOD == 27) then
                p5 = i12
                p6 = i16
                p7 = i20
                p8 = i13
                p9 = i25
              endif
            endif

! test face 4 (right)
            if (coord_to_use1(i2) < limit .and. coord_to_use1(i3) < limit .and. &
               coord_to_use1(i7) < limit .and. coord_to_use1(i6) < limit) then
              need_to_extend_this_element = .true.
              p1 = i2
              p2 = i3
              p3 = i7
              p4 = i6

              if (NGNOD == 27) then
                p5 = i10
                p6 = i15
                p7 = i18
                p8 = i14
                p9 = i23
              endif
            endif

! test face 5 (front)
            if (coord_to_use1(i1) < limit .and. coord_to_use1(i2) < limit .and. &
               coord_to_use1(i6) < limit .and. coord_to_use1(i5) < limit) then
              need_to_extend_this_element = .true.
              p1 = i1
              p2 = i2
              p3 = i6
              p4 = i5

              if (NGNOD == 27) then
                p5 = i9
                p6 = i14
                p7 = i17
                p8 = i13
                p9 = i22
              endif
            endif

! test face 6 (back)
            if (coord_to_use1(i4) < limit .and. coord_to_use1(i3) < limit .and. &
               coord_to_use1(i7) < limit .and. coord_to_use1(i8) < limit) then
              need_to_extend_this_element = .true.
              p1 = i4
              p2 = i3
              p3 = i7
              p4 = i8

              if (NGNOD == 27) then
                p5 = i11
                p6 = i15
                p7 = i19
                p8 = i16
                p9 = i24
              endif
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

              if (NGNOD == 27) then
                p5 = i9
                p6 = i10
                p7 = i11
                p8 = i12
                p9 = i21
              endif
            endif

! test face 2 (top)
            if (coord_to_use1(i5) > limit .and. coord_to_use1(i6) > limit .and. &
               coord_to_use1(i7) > limit .and. coord_to_use1(i8) > limit) then
              need_to_extend_this_element = .true.
              p1 = i5
              p2 = i6
              p3 = i7
              p4 = i8

              if (NGNOD == 27) then
                p5 = i17
                p6 = i18
                p7 = i19
                p8 = i20
                p9 = i26
              endif
            endif

! test face 3 (left)
            if (coord_to_use1(i1) > limit .and. coord_to_use1(i4) > limit .and. &
               coord_to_use1(i5) > limit .and. coord_to_use1(i8) > limit) then
              need_to_extend_this_element = .true.
              p1 = i1
              p2 = i4
              p3 = i8
              p4 = i5

              if (NGNOD == 27) then
                p5 = i12
                p6 = i16
                p7 = i20
                p8 = i13
                p9 = i25
              endif
            endif

! test face 4 (right)
            if (coord_to_use1(i2) > limit .and. coord_to_use1(i3) > limit .and. &
               coord_to_use1(i7) > limit .and. coord_to_use1(i6) > limit) then
              need_to_extend_this_element = .true.
              p1 = i2
              p2 = i3
              p3 = i7
              p4 = i6

              if (NGNOD == 27) then
                p5 = i10
                p6 = i15
                p7 = i18
                p8 = i14
                p9 = i23
              endif
            endif

! test face 5 (front)
            if (coord_to_use1(i1) > limit .and. coord_to_use1(i2) > limit .and. &
               coord_to_use1(i6) > limit .and. coord_to_use1(i5) > limit) then
              need_to_extend_this_element = .true.
              p1 = i1
              p2 = i2
              p3 = i6
              p4 = i5

              if (NGNOD == 27) then
                p5 = i9
                p6 = i14
                p7 = i17
                p8 = i13
                p9 = i22
              endif
            endif

! test face 6 (back)
            if (coord_to_use1(i4) > limit .and. coord_to_use1(i3) > limit .and. &
               coord_to_use1(i7) > limit .and. coord_to_use1(i8) > limit) then
              need_to_extend_this_element = .true.
              p1 = i4
              p2 = i3
              p3 = i7
              p4 = i8

              if (NGNOD == 27) then
                p5 = i11
                p6 = i15
                p7 = i19
                p8 = i16
                p9 = i24
              endif
            endif

          endif

          if (need_to_extend_this_element) then

! create the NUMBER_OF_LAYERS_TO_ADD_IN_THIS_STEP new elements

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

            do iextend = 1,NUMBER_OF_LAYERS_TO_ADD_IN_THIS_STEP

              ! create a new element
              elem_counter = elem_counter + 1

              ! create the material property for the extended elements
              imaterial_values_to_start_from = imaterial(ispec)

              ! if external tomographic model, do not change the value
              if (imaterial_values_to_start_from < 0) then
                imaterial_new(elem_counter) = imaterial_values_to_start_from

              ! else change the value to turn off attenuation in the PMLs
              else

                ! map back this starting value to the interval [1:count_def_mat]
                do while (imaterial_values_to_start_from > count_def_mat)
                  imaterial_values_to_start_from = imaterial_values_to_start_from - count_def_mat
                enddo

                if (istep == 1) then
                  ! use new material for the transition layer that has PML off (if it exists)
                  imaterial_new(elem_counter) = imaterial_values_to_start_from + 2*count_def_mat
                else
                  ! use new material for the PML, with attenuation off
                  imaterial_new(elem_counter) = imaterial_values_to_start_from + count_def_mat
                endif

              endif

! first create the normal element

              ! bottom face of HEX8
              x_tmp1(1) = x(p1) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*iextend
              y_tmp1(1) = y(p1) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*iextend
              z_tmp1(1) = z(p1) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*iextend

              x_tmp1(2) = x(p2) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*iextend
              y_tmp1(2) = y(p2) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*iextend
              z_tmp1(2) = z(p2) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*iextend

              x_tmp1(3) = x(p3) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*iextend
              y_tmp1(3) = y(p3) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*iextend
              z_tmp1(3) = z(p3) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*iextend

              x_tmp1(4) = x(p4) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*iextend
              y_tmp1(4) = y(p4) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*iextend
              z_tmp1(4) = z(p4) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*iextend

              ! top face of HEX8
              x_tmp1(5) = x(p1) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-1)
              y_tmp1(5) = y(p1) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-1)
              z_tmp1(5) = z(p1) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-1)

              x_tmp1(6) = x(p2) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-1)
              y_tmp1(6) = y(p2) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-1)
              z_tmp1(6) = z(p2) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-1)

              x_tmp1(7) = x(p3) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-1)
              y_tmp1(7) = y(p3) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-1)
              z_tmp1(7) = z(p3) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-1)

              x_tmp1(8) = x(p4) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-1)
              y_tmp1(8) = y(p4) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-1)
              z_tmp1(8) = z(p4) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-1)

              if (NGNOD == 27) then

                ! remaining points of bottom face of HEX27
                x_tmp1(9) = x(p5) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*iextend
                y_tmp1(9) = y(p5) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*iextend
                z_tmp1(9) = z(p5) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*iextend

                x_tmp1(10) = x(p6) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*iextend
                y_tmp1(10) = y(p6) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*iextend
                z_tmp1(10) = z(p6) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*iextend

                x_tmp1(11) = x(p7) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*iextend
                y_tmp1(11) = y(p7) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*iextend
                z_tmp1(11) = z(p7) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*iextend

                x_tmp1(12) = x(p8) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*iextend
                y_tmp1(12) = y(p8) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*iextend
                z_tmp1(12) = z(p8) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*iextend

                x_tmp1(21) = x(p9) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*iextend
                y_tmp1(21) = y(p9) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*iextend
                z_tmp1(21) = z(p9) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*iextend

                ! remaining points of top face of HEX27
                x_tmp1(17) = x(p5) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-1)
                y_tmp1(17) = y(p5) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-1)
                z_tmp1(17) = z(p5) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-1)

                x_tmp1(18) = x(p6) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-1)
                y_tmp1(18) = y(p6) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-1)
                z_tmp1(18) = z(p6) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-1)

                x_tmp1(19) = x(p7) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-1)
                y_tmp1(19) = y(p7) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-1)
                z_tmp1(19) = z(p7) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-1)

                x_tmp1(20) = x(p8) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-1)
                y_tmp1(20) = y(p8) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-1)
                z_tmp1(20) = z(p8) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-1)

                x_tmp1(26) = x(p9) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-1)
                y_tmp1(26) = y(p9) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-1)
                z_tmp1(26) = z(p9) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-1)

                ! remaining points of middle cutplane (middle "face") of HEX27
                x_tmp1(13) = x(p1) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-0.5d0)
                y_tmp1(13) = y(p1) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-0.5d0)
                z_tmp1(13) = z(p1) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-0.5d0)

                x_tmp1(14) = x(p2) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-0.5d0)
                y_tmp1(14) = y(p2) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-0.5d0)
                z_tmp1(14) = z(p2) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-0.5d0)

                x_tmp1(15) = x(p3) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-0.5d0)
                y_tmp1(15) = y(p3) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-0.5d0)
                z_tmp1(15) = z(p3) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-0.5d0)

                x_tmp1(16) = x(p4) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-0.5d0)
                y_tmp1(16) = y(p4) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-0.5d0)
                z_tmp1(16) = z(p4) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-0.5d0)

                x_tmp1(22) = x(p5) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-0.5d0)
                y_tmp1(22) = y(p5) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-0.5d0)
                z_tmp1(22) = z(p5) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-0.5d0)

                x_tmp1(23) = x(p6) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-0.5d0)
                y_tmp1(23) = y(p6) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-0.5d0)
                z_tmp1(23) = z(p6) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-0.5d0)

                x_tmp1(24) = x(p7) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-0.5d0)
                y_tmp1(24) = y(p7) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-0.5d0)
                z_tmp1(24) = z(p7) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-0.5d0)

                x_tmp1(25) = x(p8) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-0.5d0)
                y_tmp1(25) = y(p8) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-0.5d0)
                z_tmp1(25) = z(p8) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-0.5d0)

                x_tmp1(27) = x(p9) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-0.5d0)
                y_tmp1(27) = y(p9) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-0.5d0)
                z_tmp1(27) = z(p9) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-0.5d0)

              endif

! then create the mirrored element

              ! bottom face of HEX8
              x_tmp2(1) = x(p1) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-1)
              y_tmp2(1) = y(p1) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-1)
              z_tmp2(1) = z(p1) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-1)

              x_tmp2(2) = x(p2) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-1)
              y_tmp2(2) = y(p2) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-1)
              z_tmp2(2) = z(p2) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-1)

              x_tmp2(3) = x(p3) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-1)
              y_tmp2(3) = y(p3) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-1)
              z_tmp2(3) = z(p3) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-1)

              x_tmp2(4) = x(p4) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-1)
              y_tmp2(4) = y(p4) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-1)
              z_tmp2(4) = z(p4) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-1)

              ! top face of HEX8
              x_tmp2(5) = x(p1) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*iextend
              y_tmp2(5) = y(p1) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*iextend
              z_tmp2(5) = z(p1) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*iextend

              x_tmp2(6) = x(p2) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*iextend
              y_tmp2(6) = y(p2) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*iextend
              z_tmp2(6) = z(p2) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*iextend

              x_tmp2(7) = x(p3) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*iextend
              y_tmp2(7) = y(p3) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*iextend
              z_tmp2(7) = z(p3) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*iextend

              x_tmp2(8) = x(p4) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*iextend
              y_tmp2(8) = y(p4) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*iextend
              z_tmp2(8) = z(p4) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*iextend

              if (NGNOD == 27) then

                ! remaining points of bottom face of HEX27
                x_tmp2(9) = x(p5) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-1)
                y_tmp2(9) = y(p5) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-1)
                z_tmp2(9) = z(p5) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-1)

                x_tmp2(10) = x(p6) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-1)
                y_tmp2(10) = y(p6) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-1)
                z_tmp2(10) = z(p6) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-1)

                x_tmp2(11) = x(p7) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-1)
                y_tmp2(11) = y(p7) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-1)
                z_tmp2(11) = z(p7) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-1)

                x_tmp2(12) = x(p8) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-1)
                y_tmp2(12) = y(p8) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-1)
                z_tmp2(12) = z(p8) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-1)

                x_tmp2(21) = x(p9) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-1)
                y_tmp2(21) = y(p9) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-1)
                z_tmp2(21) = z(p9) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-1)

                ! remaining points of top face of HEX27
                x_tmp2(17) = x(p5) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*iextend
                y_tmp2(17) = y(p5) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*iextend
                z_tmp2(17) = z(p5) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*iextend

                x_tmp2(18) = x(p6) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*iextend
                y_tmp2(18) = y(p6) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*iextend
                z_tmp2(18) = z(p6) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*iextend

                x_tmp2(19) = x(p7) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*iextend
                y_tmp2(19) = y(p7) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*iextend
                z_tmp2(19) = z(p7) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*iextend

                x_tmp2(20) = x(p8) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*iextend
                y_tmp2(20) = y(p8) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*iextend
                z_tmp2(20) = z(p8) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*iextend

                x_tmp2(26) = x(p9) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*iextend
                y_tmp2(26) = y(p9) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*iextend
                z_tmp2(26) = z(p9) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*iextend

                ! remaining points of middle cutplane (middle "face") of HEX27
                x_tmp2(13) = x(p1) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-0.5d0)
                y_tmp2(13) = y(p1) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-0.5d0)
                z_tmp2(13) = z(p1) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-0.5d0)

                x_tmp2(14) = x(p2) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-0.5d0)
                y_tmp2(14) = y(p2) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-0.5d0)
                z_tmp2(14) = z(p2) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-0.5d0)

                x_tmp2(15) = x(p3) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-0.5d0)
                y_tmp2(15) = y(p3) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-0.5d0)
                z_tmp2(15) = z(p3) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-0.5d0)

                x_tmp2(16) = x(p4) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-0.5d0)
                y_tmp2(16) = y(p4) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-0.5d0)
                z_tmp2(16) = z(p4) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-0.5d0)

                x_tmp2(22) = x(p5) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-0.5d0)
                y_tmp2(22) = y(p5) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-0.5d0)
                z_tmp2(22) = z(p5) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-0.5d0)

                x_tmp2(23) = x(p6) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-0.5d0)
                y_tmp2(23) = y(p6) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-0.5d0)
                z_tmp2(23) = z(p6) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-0.5d0)

                x_tmp2(24) = x(p7) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-0.5d0)
                y_tmp2(24) = y(p7) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-0.5d0)
                z_tmp2(24) = z(p7) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-0.5d0)

                x_tmp2(25) = x(p8) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-0.5d0)
                y_tmp2(25) = y(p8) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-0.5d0)
                z_tmp2(25) = z(p8) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-0.5d0)

                x_tmp2(27) = x(p9) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-0.5d0)
                y_tmp2(27) = y(p9) + factor_y*SIZE_OF_Y_ELEMENT_TO_ADD*(iextend-0.5d0)
                z_tmp2(27) = z(p9) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-0.5d0)

              endif

! now we need to test if the element created is flipped i.e. it has a negative Jacobian,
! and if so we will use the mirrored version of that element, which will then have a positive Jacobian

! check the element for a negative Jacobian
              do ia = 1,NGNOD
                xelm(ia) = x_tmp1(ia)
                yelm(ia) = y_tmp1(ia)
                zelm(ia) = z_tmp1(ia)
              enddo
              call calc_jacobian(xelm,yelm,zelm,dershape3D,found_a_negative_jacobian1,NDIM,NGNOD,NGLLX,NGLLY,NGLLZ,jacobian)

! check the mirrored (i.e. flipped/swapped) element for a negative Jacobian
! either this one or the non-mirrored one above should be OK, and thus we will select it
              do ia = 1,NGNOD
                xelm(ia) = x_tmp2(ia)
                yelm(ia) = y_tmp2(ia)
                zelm(ia) = z_tmp2(ia)
              enddo
              call calc_jacobian(xelm,yelm,zelm,dershape3D,found_a_negative_jacobian2,NDIM,NGNOD,NGLLX,NGLLY,NGLLZ,jacobian)

              ! this should never happen, it is just a safety test
              if (found_a_negative_jacobian1 .and. found_a_negative_jacobian2) &
                stop 'error: found a negative Jacobian that could not be fixed'

              ! this should also never happen, it is just a second safety test
              if (.not. found_a_negative_jacobian1 .and. .not. found_a_negative_jacobian2) &
                stop 'strange error: both the element created and its mirrored version have a positive Jacobian!'

! if we have found that the original element has a negative Jacobian and its mirrored element is fine,
! swap the points so that we use that mirrored element in the final mesh saved to disk instead of the original one
! i.e. implement a mirror symmetry here
              if (found_a_negative_jacobian1) then
                do ia = 1,NGNOD
                  npoin_new_real = npoin_new_real + 1
                  ibool_new(ia,elem_counter) = npoin_new_real
                  x_new(npoin_new_real) = x_tmp2(ia)
                  y_new(npoin_new_real) = y_tmp2(ia)
                  z_new(npoin_new_real) = z_tmp2(ia)
                enddo
              else
                do ia = 1,NGNOD
                  npoin_new_real = npoin_new_real + 1
                  ibool_new(ia,elem_counter) = npoin_new_real
                  x_new(npoin_new_real) = x_tmp1(ia)
                  y_new(npoin_new_real) = y_tmp1(ia)
                  z_new(npoin_new_real) = z_tmp1(ia)
                enddo
              endif

            enddo
          endif

        enddo

        if (minval(ibool_new) /= 1) stop 'error in minval(ibool_new)'
        if (maxval(ibool_new) > npoin_new_max) stop 'error in maxval(ibool_new)'

! deallocate the original arrays
        deallocate(x,y,z)
        deallocate(ibool)
        deallocate(imaterial)

! reallocate them with the new size
        allocate(x(npoin_new_real))
        allocate(y(npoin_new_real))
        allocate(z(npoin_new_real))
        allocate(imaterial(nspec_new))
        allocate(ibool(NGNOD,nspec_new))

! make the new ones become the old ones, to prepare for the next iteration of the two nested loops we are in,
! i.e. to make sure the next loop will extend the mesh from the new arrays rather than from the old ones
        x(:) = x_new(1:npoin_new_real)
        y(:) = y_new(1:npoin_new_real)
        z(:) = z_new(1:npoin_new_real)
        imaterial(:) = imaterial_new(:)
        ibool(:,:) = ibool_new(:,:)

! the new number of elements and points becomes the old one, for the same reason
        nspec = nspec_new
        npoin = npoin_new_real

! deallocate the new ones, to make sure they can be allocated again in the next iteration of the nested loops we are in
        deallocate(x_new,y_new,z_new)
        deallocate(ibool_new)
        deallocate(imaterial_new)

        if (minval(ibool) /= 1) stop 'error in minval(ibool)'
        if (maxval(ibool) > npoin) stop 'error in maxval(ibool)'

      enddo ! of iloop_on_min_face_then_max_face loop on Xmin then Xmax, or Ymin then Ymax, or Zmin then Zmax

! end of loop on the three sets of faces to first add CPML elements along X, then along Y, then along Z
    enddo

! we must remove all point multiples here using a fast sorting routine

    print *,'beginning of multiple point removal based on sorting...'

    npointot = nspec*NGNOD

    allocate(locval(npointot))
    allocate(ifseg(npointot))

    allocate(xp(npointot))
    allocate(yp(npointot))
    allocate(zp(npointot))

    allocate(x_copy(NGNOD,nspec))
    allocate(y_copy(NGNOD,nspec))
    allocate(z_copy(NGNOD,nspec))

    do ispec=1,nspec
      ieoff = NGNOD * (ispec-1)
      ilocnum = 0
      do ia = 1,NGNOD
        ilocnum = ilocnum + 1
        xp(ilocnum+ieoff) = x(ibool(ia,ispec))
        yp(ilocnum+ieoff) = y(ibool(ia,ispec))
        zp(ilocnum+ieoff) = z(ibool(ia,ispec))

! create a copy, since the sorting below will destroy the arrays
        x_copy(ia,ispec) = x(ibool(ia,ispec))
        y_copy(ia,ispec) = y(ibool(ia,ispec))
        z_copy(ia,ispec) = z(ibool(ia,ispec))
      enddo
    enddo

! gets ibool indexing from local (NGNOD points) to global points
    call get_global(NDIM,npointot,xp,yp,zp,ibool,locval,ifseg,npoin,very_small_distance)

    print *,'done with multiple point removal based on sorting'
    print *,'found a total of ',npoin,' unique mesh points'

! unique global point locations
    do ispec = 1, nspec
      do ia = 1,NGNOD
        iglobnum = ibool(ia,ispec)
        x(iglobnum) = x_copy(ia,ispec)
        y(iglobnum) = y_copy(ia,ispec)
        z(iglobnum) = z_copy(ia,ispec)
      enddo
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
      if (NGNOD == 8) then
        do ispec = 1,nspec
          write(23,"(9(i9,1x))") ispec,(ibool(ia,ispec), ia = 1,NGNOD)
        enddo
      else
        do ispec = 1,nspec
          write(23,"(28(i9,1x))") ispec,(ibool(ia,ispec), ia = 1,NGNOD)
        enddo
      endif
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

    deallocate(x)
    deallocate(y)
    deallocate(z)
    deallocate(imaterial)
    deallocate(ibool)

    deallocate(locval)
    deallocate(ifseg)

    deallocate(xp)
    deallocate(yp)
    deallocate(zp)

    deallocate(x_copy)
    deallocate(y_copy)
    deallocate(z_copy)

  enddo ! of do istep = 1,2

! output information for the next code (xconvert_external_layers_of_a_given_mesh_to_CPML_layers)
  print *
  print *,'Here are the values to use as input in the next code, xconvert_external_layers_of_a_given_mesh_to_CPML_layers:'
  print *,'  (also saved in file values_to_use_for_convert_external_layers_of_a_given_mesh_to_CPML_layers.txt)'
  print *
  if (ADD_ON_THE_XMIN_SURFACE) print *,'THICKNESS_OF_XMIN_PML = ',sngl(SIZE_OF_XMIN_ELEMENT_TO_ADD * NUMBER_OF_PML_LAYERS_TO_ADD)
  if (ADD_ON_THE_XMAX_SURFACE) print *,'THICKNESS_OF_XMAX_PML = ',sngl(SIZE_OF_XMAX_ELEMENT_TO_ADD * NUMBER_OF_PML_LAYERS_TO_ADD)
  if (ADD_ON_THE_YMIN_SURFACE) print *,'THICKNESS_OF_YMIN_PML = ',sngl(SIZE_OF_YMIN_ELEMENT_TO_ADD * NUMBER_OF_PML_LAYERS_TO_ADD)
  if (ADD_ON_THE_YMAX_SURFACE) print *,'THICKNESS_OF_YMAX_PML = ',sngl(SIZE_OF_YMAX_ELEMENT_TO_ADD * NUMBER_OF_PML_LAYERS_TO_ADD)
  if (ADD_ON_THE_ZMIN_SURFACE) print *,'THICKNESS_OF_ZMIN_PML = ',sngl(SIZE_OF_ZMIN_ELEMENT_TO_ADD * NUMBER_OF_PML_LAYERS_TO_ADD)
  if (ADD_ON_THE_ZMAX_SURFACE) print *,'THICKNESS_OF_ZMAX_PML = ',sngl(SIZE_OF_ZMAX_ELEMENT_TO_ADD * NUMBER_OF_PML_LAYERS_TO_ADD)
  print *

! save the thickness values to a text file
  open(unit=23,file='values_to_use_for_convert_external_layers_of_a_given_mesh_to_CPML_layers.txt', &
       status='unknown',action='write')

  if (ADD_ON_THE_XMIN_SURFACE) then
    write(23,*) SIZE_OF_XMIN_ELEMENT_TO_ADD * NUMBER_OF_PML_LAYERS_TO_ADD
  else
! convention (negative value) to say that this absorbing edge is turned off
    write(23,*) -1
  endif

  if (ADD_ON_THE_XMAX_SURFACE) then
    write(23,*) SIZE_OF_XMAX_ELEMENT_TO_ADD * NUMBER_OF_PML_LAYERS_TO_ADD
  else
! convention (negative value) to say that this absorbing edge is turned off
    write(23,*) -1
  endif

  if (ADD_ON_THE_YMIN_SURFACE) then
    write(23,*) SIZE_OF_YMIN_ELEMENT_TO_ADD * NUMBER_OF_PML_LAYERS_TO_ADD
  else
! convention (negative value) to say that this absorbing edge is turned off
    write(23,*) -1
  endif

  if (ADD_ON_THE_YMAX_SURFACE) then
    write(23,*) SIZE_OF_YMAX_ELEMENT_TO_ADD * NUMBER_OF_PML_LAYERS_TO_ADD
  else
! convention (negative value) to say that this absorbing edge is turned off
    write(23,*) -1
  endif

  if (ADD_ON_THE_ZMIN_SURFACE) then
    write(23,*) SIZE_OF_ZMIN_ELEMENT_TO_ADD * NUMBER_OF_PML_LAYERS_TO_ADD
  else
! convention (negative value) to say that this absorbing edge is turned off
    write(23,*) -1
  endif

  if (ADD_ON_THE_ZMAX_SURFACE) then
    write(23,*) SIZE_OF_ZMAX_ELEMENT_TO_ADD * NUMBER_OF_PML_LAYERS_TO_ADD
  else
! convention (negative value) to say that this absorbing edge is turned off
    write(23,*) -1
  endif

  close(23)

  end program add_CPML_layers_to_a_given_mesh

!
!=====================================================================
!

! 3D shape functions for 8-node or 27-node element

  subroutine get_shape3D(dershape3D,xigll,yigll,zigll,NGNOD,NGLLX,NGLLY,NGLLZ,NDIM)

  implicit none

! a few useful constants
  double precision, parameter :: ZERO = 0.d0, ONE = 1.d0, TWO = 2.d0, HALF = 0.5d0

! very small value
  double precision, parameter :: TINYVAL = 1.d-9

  integer NGNOD,NGLLX,NGLLY,NGLLZ,NDIM

! Gauss-Lobatto-Legendre points of integration
  double precision xigll(NGLLX)
  double precision yigll(NGLLY)
  double precision zigll(NGLLZ)

! 3D shape function derivatives
  double precision dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ)

  integer i,j,k,ia

! location of the nodes of the 3D hexahedra elements
  double precision xi,eta,gamma
  double precision ra1,ra2,rb1,rb2,rc1,rc2

! for checking the 3D shape functions
  double precision sumdershapexi,sumdershapeeta,sumdershapegamma

  double precision, parameter :: ONE_EIGHTH = 0.125d0

! check that the parameter file is correct
  if (NGNOD /= 8 .and. NGNOD /= 27) stop 'volume elements should have 8 or 27 control nodes'

! ***
! *** create 3D shape functions and jacobian
! ***

  do i=1,NGLLX
    do j=1,NGLLY
      do k=1,NGLLZ

        xi = xigll(i)
        eta = yigll(j)
        gamma = zigll(k)

        !--- case of a 3D 8-node element (Dhatt-Touzot p. 115)
        if (NGNOD == 8) then

          ra1 = one + xi
          ra2 = one - xi

          rb1 = one + eta
          rb2 = one - eta

          rc1 = one + gamma
          rc2 = one - gamma

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

        else

          ! note: put further initialization for NGNOD == 27 into subroutine
          !       to avoid compilation errors in case NGNOD == 8
          call get_shape3D_27(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ,dershape3D,xi,eta,gamma,i,j,k)

        endif

      enddo
    enddo
  enddo

!--- check the shape functions and their derivatives

  do i=1,NGLLX
    do j=1,NGLLY
      do k=1,NGLLZ

        sumdershapexi = ZERO
        sumdershapeeta = ZERO
        sumdershapegamma = ZERO

        do ia=1,NGNOD
          sumdershapexi = sumdershapexi + dershape3D(1,ia,i,j,k)
          sumdershapeeta = sumdershapeeta + dershape3D(2,ia,i,j,k)
          sumdershapegamma = sumdershapegamma + dershape3D(3,ia,i,j,k)
        enddo

        ! sum of shape functions should be one
        ! sum of derivative of shape functions should be zero
        if (abs(sumdershapexi) > TINYVAL) stop 'error in xi derivative of 3D shape functions'
        if (abs(sumdershapeeta) > TINYVAL) stop 'error in eta derivative of 3D shape functions'
        if (abs(sumdershapegamma) > TINYVAL) stop 'error in gamma derivative of 3D shape functions'

      enddo
    enddo
  enddo

  end subroutine get_shape3D

!
!-------------------------------------------------------------------------------------------------
!

!--- case of a 3D 27-node element

  subroutine get_shape3D_27(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ,dershape3D,xi,eta,gamma,i,j,k)

  implicit none

! a few useful constants
  double precision, parameter :: ZERO = 0.d0, ONE = 1.d0, TWO = 2.d0, HALF = 0.5d0

  integer :: NDIM,NGNOD,NGLLX,NGLLY,NGLLZ,i,j,k

! 3D shape function derivatives
  double precision dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ)

! location of the nodes of the 3D hexahedra elements
  double precision xi,eta,gamma
  double precision l1xi,l2xi,l3xi,l1eta,l2eta,l3eta,l1gamma,l2gamma,l3gamma
  double precision l1pxi,l2pxi,l3pxi,l1peta,l2peta,l3peta,l1pgamma,l2pgamma,l3pgamma

  l1xi=HALF*xi*(xi-ONE)
  l2xi=ONE-xi**2
  l3xi=HALF*xi*(xi+ONE)

  l1pxi=xi-HALF
  l2pxi=-TWO*xi
  l3pxi=xi+HALF

  l1eta=HALF*eta*(eta-ONE)
  l2eta=ONE-eta**2
  l3eta=HALF*eta*(eta+ONE)

  l1peta=eta-HALF
  l2peta=-TWO*eta
  l3peta=eta+HALF

  l1gamma=HALF*gamma*(gamma-ONE)
  l2gamma=ONE-gamma**2
  l3gamma=HALF*gamma*(gamma+ONE)

  l1pgamma=gamma-HALF
  l2pgamma=-TWO*gamma
  l3pgamma=gamma+HALF

  ! corner nodes

  dershape3D(1,1,i,j,k)=l1pxi*l1eta*l1gamma
  dershape3D(1,2,i,j,k)=l3pxi*l1eta*l1gamma
  dershape3D(1,3,i,j,k)=l3pxi*l3eta*l1gamma
  dershape3D(1,4,i,j,k)=l1pxi*l3eta*l1gamma
  dershape3D(1,5,i,j,k)=l1pxi*l1eta*l3gamma
  dershape3D(1,6,i,j,k)=l3pxi*l1eta*l3gamma
  dershape3D(1,7,i,j,k)=l3pxi*l3eta*l3gamma
  dershape3D(1,8,i,j,k)=l1pxi*l3eta*l3gamma

  dershape3D(2,1,i,j,k)=l1xi*l1peta*l1gamma
  dershape3D(2,2,i,j,k)=l3xi*l1peta*l1gamma
  dershape3D(2,3,i,j,k)=l3xi*l3peta*l1gamma
  dershape3D(2,4,i,j,k)=l1xi*l3peta*l1gamma
  dershape3D(2,5,i,j,k)=l1xi*l1peta*l3gamma
  dershape3D(2,6,i,j,k)=l3xi*l1peta*l3gamma
  dershape3D(2,7,i,j,k)=l3xi*l3peta*l3gamma
  dershape3D(2,8,i,j,k)=l1xi*l3peta*l3gamma

  dershape3D(3,1,i,j,k)=l1xi*l1eta*l1pgamma
  dershape3D(3,2,i,j,k)=l3xi*l1eta*l1pgamma
  dershape3D(3,3,i,j,k)=l3xi*l3eta*l1pgamma
  dershape3D(3,4,i,j,k)=l1xi*l3eta*l1pgamma
  dershape3D(3,5,i,j,k)=l1xi*l1eta*l3pgamma
  dershape3D(3,6,i,j,k)=l3xi*l1eta*l3pgamma
  dershape3D(3,7,i,j,k)=l3xi*l3eta*l3pgamma
  dershape3D(3,8,i,j,k)=l1xi*l3eta*l3pgamma

  ! midside nodes

  dershape3D(1,9,i,j,k)=l2pxi*l1eta*l1gamma
  dershape3D(1,10,i,j,k)=l3pxi*l2eta*l1gamma
  dershape3D(1,11,i,j,k)=l2pxi*l3eta*l1gamma
  dershape3D(1,12,i,j,k)=l1pxi*l2eta*l1gamma
  dershape3D(1,13,i,j,k)=l1pxi*l1eta*l2gamma
  dershape3D(1,14,i,j,k)=l3pxi*l1eta*l2gamma
  dershape3D(1,15,i,j,k)=l3pxi*l3eta*l2gamma
  dershape3D(1,16,i,j,k)=l1pxi*l3eta*l2gamma
  dershape3D(1,17,i,j,k)=l2pxi*l1eta*l3gamma
  dershape3D(1,18,i,j,k)=l3pxi*l2eta*l3gamma
  dershape3D(1,19,i,j,k)=l2pxi*l3eta*l3gamma
  dershape3D(1,20,i,j,k)=l1pxi*l2eta*l3gamma

  dershape3D(2,9,i,j,k)=l2xi*l1peta*l1gamma
  dershape3D(2,10,i,j,k)=l3xi*l2peta*l1gamma
  dershape3D(2,11,i,j,k)=l2xi*l3peta*l1gamma
  dershape3D(2,12,i,j,k)=l1xi*l2peta*l1gamma
  dershape3D(2,13,i,j,k)=l1xi*l1peta*l2gamma
  dershape3D(2,14,i,j,k)=l3xi*l1peta*l2gamma
  dershape3D(2,15,i,j,k)=l3xi*l3peta*l2gamma
  dershape3D(2,16,i,j,k)=l1xi*l3peta*l2gamma
  dershape3D(2,17,i,j,k)=l2xi*l1peta*l3gamma
  dershape3D(2,18,i,j,k)=l3xi*l2peta*l3gamma
  dershape3D(2,19,i,j,k)=l2xi*l3peta*l3gamma
  dershape3D(2,20,i,j,k)=l1xi*l2peta*l3gamma

  dershape3D(3,9,i,j,k)=l2xi*l1eta*l1pgamma
  dershape3D(3,10,i,j,k)=l3xi*l2eta*l1pgamma
  dershape3D(3,11,i,j,k)=l2xi*l3eta*l1pgamma
  dershape3D(3,12,i,j,k)=l1xi*l2eta*l1pgamma
  dershape3D(3,13,i,j,k)=l1xi*l1eta*l2pgamma
  dershape3D(3,14,i,j,k)=l3xi*l1eta*l2pgamma
  dershape3D(3,15,i,j,k)=l3xi*l3eta*l2pgamma
  dershape3D(3,16,i,j,k)=l1xi*l3eta*l2pgamma
  dershape3D(3,17,i,j,k)=l2xi*l1eta*l3pgamma
  dershape3D(3,18,i,j,k)=l3xi*l2eta*l3pgamma
  dershape3D(3,19,i,j,k)=l2xi*l3eta*l3pgamma
  dershape3D(3,20,i,j,k)=l1xi*l2eta*l3pgamma

  ! side center nodes

  dershape3D(1,21,i,j,k)=l2pxi*l2eta*l1gamma
  dershape3D(1,22,i,j,k)=l2pxi*l1eta*l2gamma
  dershape3D(1,23,i,j,k)=l3pxi*l2eta*l2gamma
  dershape3D(1,24,i,j,k)=l2pxi*l3eta*l2gamma
  dershape3D(1,25,i,j,k)=l1pxi*l2eta*l2gamma
  dershape3D(1,26,i,j,k)=l2pxi*l2eta*l3gamma

  dershape3D(2,21,i,j,k)=l2xi*l2peta*l1gamma
  dershape3D(2,22,i,j,k)=l2xi*l1peta*l2gamma
  dershape3D(2,23,i,j,k)=l3xi*l2peta*l2gamma
  dershape3D(2,24,i,j,k)=l2xi*l3peta*l2gamma
  dershape3D(2,25,i,j,k)=l1xi*l2peta*l2gamma
  dershape3D(2,26,i,j,k)=l2xi*l2peta*l3gamma

  dershape3D(3,21,i,j,k)=l2xi*l2eta*l1pgamma
  dershape3D(3,22,i,j,k)=l2xi*l1eta*l2pgamma
  dershape3D(3,23,i,j,k)=l3xi*l2eta*l2pgamma
  dershape3D(3,24,i,j,k)=l2xi*l3eta*l2pgamma
  dershape3D(3,25,i,j,k)=l1xi*l2eta*l2pgamma
  dershape3D(3,26,i,j,k)=l2xi*l2eta*l3pgamma

  ! center node

  dershape3D(1,27,i,j,k)=l2pxi*l2eta*l2gamma
  dershape3D(2,27,i,j,k)=l2xi*l2peta*l2gamma
  dershape3D(3,27,i,j,k)=l2xi*l2eta*l2pgamma

  end subroutine get_shape3D_27

!
!=====================================================================
!

  subroutine calc_jacobian(xelm,yelm,zelm,dershape3D,found_a_negative_jacobian,NDIM,NGNOD,NGLLX,NGLLY,NGLLZ,jacobian)

  implicit none

  integer :: NDIM,NGNOD,NGLLX,NGLLY,NGLLZ

  logical :: found_a_negative_jacobian

  double precision, dimension(NGNOD) :: xelm,yelm,zelm
  double precision dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ)

  integer :: i,j,k,ia
  double precision :: xxi,xeta,xgamma,yxi,yeta,ygamma,zxi,zeta,zgamma
  double precision :: jacobian

  double precision, parameter :: ZERO = 0.d0

  found_a_negative_jacobian = .false.

! do k=1,NGLLZ
!   do j=1,NGLLY
!     do i=1,NGLLX
! for this CPML mesh extrusion routine it is sufficient to test the 8 corners of each element to reduce the cost
! because we just want to detect if the element is flipped or not, and if so flip it back
  do k=1,NGLLZ,NGLLZ-1
    do j=1,NGLLY,NGLLY-1
      do i=1,NGLLX,NGLLX-1

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
      if (jacobian <= ZERO) then
        found_a_negative_jacobian = .true.
        return
      endif

      enddo
    enddo
  enddo

  end subroutine calc_jacobian

! ------------------------------------------------------------------

! subroutines to sort indexing arrays based on geometrical coordinates instead of based on topology (because that is much faster)

  subroutine get_global(NDIM,npointot,xp,yp,zp,iglob,locval,ifseg,nglob,SMALLVALTOL)

! this routine MUST be in double precision to avoid sensitivity
! to roundoff errors in the coordinates of the points

! non-structured global numbering software provided by Paul F. Fischer

! leave sorting subroutines in same source file to allow for inlining

  implicit none

  integer NDIM
  integer npointot
  integer nglob
  integer iglob(npointot),locval(npointot)
  logical ifseg(npointot)
  double precision xp(npointot),yp(npointot),zp(npointot)

! geometry tolerance parameter to calculate number of independent grid points
! small value for double precision and to avoid sensitivity to roundoff
  double precision SMALLVALTOL

  integer :: ier

  integer, dimension(:), allocatable :: ninseg,idummy

! dynamically allocate arrays
  allocate(ninseg(npointot),stat=ier)
  if (ier /= 0) stop 'error allocating array ninseg'
  allocate(idummy(npointot),stat=ier)
  if (ier /= 0) stop 'error allocating array idummy'

  call sort_array_coordinates(NDIM,npointot,xp,yp,zp,idummy,iglob,locval,ifseg, &
                              nglob,ninseg,SMALLVALTOL)

! deallocate arrays
  deallocate(ninseg)
  deallocate(idummy)

  end subroutine get_global

! ------------------------------------------------------------------

! subroutines to sort indexing arrays based on geometrical coordinates instead of based on topology (because that is much faster)

  subroutine sort_array_coordinates(NDIM,npointot,x,y,z,ibool,iglob,locval,ifseg,nglob,ninseg,xtol)

! this routine MUST be in double precision to avoid sensitivity
! to roundoff errors in the coordinates of the points
!
! returns: sorted indexing array (ibool),  reordering array (iglob) & number of global points (nglob)

  implicit none

  integer, intent(in) :: NDIM,npointot
  double precision, dimension(npointot), intent(inout) :: x, y, z
  integer, dimension(npointot), intent(inout) :: ibool

  integer, dimension(npointot), intent(out) :: iglob, locval, ninseg
  logical, dimension(npointot), intent(out) :: ifseg
  integer, intent(out) :: nglob
  double precision, intent(in) :: xtol

  ! local parameters
  integer :: i, j
  integer :: nseg, ioff, iseg, ig

  ! establish initial pointers
  do i = 1,npointot
    locval(i) = i
  enddo

  ifseg(:) = .false.

  nseg = 1
  ifseg(1) = .true.
  ninseg(1) = npointot

  do j = 1,NDIM

    ! sort within each segment
    ioff = 1
    do iseg = 1,nseg
      if (j == 1) then
        ! sort on X
        call heap_sort_multi(ninseg(iseg), x(ioff), y(ioff), z(ioff), ibool(ioff), locval(ioff))
      else if (j == 2) then
        ! then sort on Y for a sublist of given constant X
        call heap_sort_multi(ninseg(iseg), y(ioff), x(ioff), z(ioff), ibool(ioff), locval(ioff))
      else
        ! then sort on Z for a sublist of given constant X and Y
        call heap_sort_multi(ninseg(iseg), z(ioff), x(ioff), y(ioff), ibool(ioff), locval(ioff))
      endif
      ioff = ioff + ninseg(iseg)
    enddo

    ! check for jumps in current coordinate
    ! define a tolerance, normalized radius is 1., so let's use a small value
    if (j == 1) then
      do i = 2,npointot
        if (dabs(x(i) - x(i-1)) > xtol) ifseg(i) = .true.
      enddo
    else if (j == 2) then
      do i = 2,npointot
        if (dabs(y(i) - y(i-1)) > xtol) ifseg(i) = .true.
      enddo
    else
      do i = 2,npointot
        if (dabs(z(i) - z(i-1)) > xtol) ifseg(i) = .true.
      enddo
    endif

    ! count up number of different segments
    nseg = 0
    do i = 1,npointot
      if (ifseg(i)) then
        nseg = nseg + 1
        ninseg(nseg) = 1
      else
        ninseg(nseg) = ninseg(nseg) + 1
      endif
    enddo

  enddo

  ! assign global node numbers (now sorted lexicographically)
  ig = 0
  do i = 1,npointot
    ! eliminate the multiples by using a single (new) point number for all the points that have the same X Y Z after sorting
    if (ifseg(i)) ig = ig + 1
    iglob(locval(i)) = ig
  enddo

  nglob = ig

  end subroutine sort_array_coordinates

!
!--------------------
!


! sorting routine left here for inlining

! -------------------- library for sorting routine ------------------

! sorting routines put here in same file to allow for inlining

! this directive avoids triggering a random bug in Intel ifort v13 (in the compiler, not in SPECFEM),
! fixed in later versions of Intel ifort, which also ignore this directive because it was discontinued
!$DIR NOOPTIMIZE
  subroutine heap_sort_multi(N, dx, dy, dz, ia, ib)

  implicit none
  integer, intent(in) :: N
  double precision, dimension(N), intent(inout) :: dx
  double precision, dimension(N), intent(inout) :: dy
  double precision, dimension(N), intent(inout) :: dz
  integer, dimension(N), intent(inout) :: ia
  integer, dimension(N), intent(inout) :: ib

  integer :: i

  ! checks if anything to do
  if (N < 2) return

  ! builds heap
  do i = N/2, 1, -1
    call heap_sort_siftdown(i, n)
  enddo

  ! sorts array
  do i = N, 2, -1
    ! swaps last and first entry in this section
    call dswap(dx, 1, i)
    call dswap(dy, 1, i)
    call dswap(dz, 1, i)
    call iswap(ia, 1, i)
    call iswap(ib, 1, i)
    call heap_sort_siftdown(1, i - 1)
  enddo

  contains

! this directive avoids triggering a random bug in Intel ifort v13 (in the compiler, not in SPECFEM),
! fixed in later versions of Intel ifort, which also ignore this directive because it was discontinued
!$DIR NOOPTIMIZE
    subroutine dswap(A, i, j)

    double precision, dimension(:), intent(inout) :: A
    integer, intent(in) :: i
    integer, intent(in) :: j

    double precision :: tmp

    tmp = A(i)
    A(i) = A(j)
    A(j) = tmp

    end subroutine

    subroutine iswap(A, i, j)

    integer, dimension(:), intent(inout) :: A
    integer, intent(in) :: i
    integer, intent(in) :: j

    integer :: tmp

    tmp = A(i)
    A(i) = A(j)
    A(j) = tmp

    end subroutine

! this directive avoids triggering a random bug in Intel ifort v13 (in the compiler, not in SPECFEM),
! fixed in later versions of Intel ifort, which also ignore this directive because it was discontinued
!$DIR NOOPTIMIZE
    subroutine heap_sort_siftdown(start, bottom)

    integer, intent(in) :: start
    integer, intent(in) :: bottom

    integer :: i, j
    double precision :: xtmp, ytmp, ztmp
    integer :: atmp, btmp

    i = start
    xtmp = dx(i)
    ytmp = dy(i)
    ztmp = dz(i)
    atmp = ia(i)
    btmp = ib(i)

    j = 2 * i
    do while (j <= bottom)
      ! chooses larger value first in this section
      if (j < bottom) then
        if (dx(j) <= dx(j+1)) j = j + 1
      endif

      ! checks if section already smaller than initial value
      if (dx(j) < xtmp) exit

      dx(i) = dx(j)
      dy(i) = dy(j)
      dz(i) = dz(j)
      ia(i) = ia(j)
      ib(i) = ib(j)
      i = j
      j = 2 * i
    enddo

    dx(i) = xtmp
    dy(i) = ytmp
    dz(i) = ztmp
    ia(i) = atmp
    ib(i) = btmp

    end subroutine heap_sort_siftdown

  end subroutine heap_sort_multi

