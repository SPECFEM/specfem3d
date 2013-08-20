!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
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

  subroutine get_coupling_surfaces(myrank, &
                        nspec,ibool,NPROC, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        num_interfaces_ext_mesh,max_interface_size_ext_mesh, &
                        my_neighbours_ext_mesh)

! determines coupling surface for acoustic-elastic domains
! based on ispec_is_acoustic, ispec_is_elastic and ispec_is_poroelastic arrays

  use create_regions_mesh_ext_par
  implicit none

! number of spectral elements in each block
  integer :: myrank,nspec,NPROC

! arrays with the mesh
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

! MPI communication
  integer :: num_interfaces_ext_mesh,max_interface_size_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: my_neighbours_ext_mesh
  integer, dimension(NGLLX*NGLLX*max_interface_size_ext_mesh,num_interfaces_ext_mesh) :: &
            ibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh

  ! local parameters
  integer, dimension(:), allocatable :: elastic_flag,acoustic_flag,poroelastic_flag
  integer, dimension(:,:), allocatable :: ibool_interfaces_ext_mesh_dummy
  integer :: max_nibool_interfaces_ext_mesh
  integer :: count_elastic,count_acoustic,count_poroelastic
  integer :: ispec,i,j,k,iglob,ier,inum

  ! initializes number of coupling faces
  num_coupling_ac_el_faces = 0
  num_coupling_ac_po_faces = 0
  num_coupling_el_po_faces = 0

  ! sets flags for acoustic / elastic / poroelastic on global points
  allocate(acoustic_flag(nglob_dummy),stat=ier)
  if( ier /= 0 ) stop 'error allocating array acoustic_flag'
  allocate(elastic_flag(nglob_dummy),stat=ier)
  if( ier /= 0 ) stop 'error allocating array elastic_flag'
  allocate(poroelastic_flag(nglob_dummy),stat=ier)
  if( ier /= 0 ) stop 'error allocating array poroelastic_flag'

  acoustic_flag(:) = 0
  elastic_flag(:) = 0
  poroelastic_flag(:) = 0

  count_acoustic = 0
  count_elastic = 0
  count_poroelastic = 0

  do ispec = 1, nspec
    ! counts elements
    if( ispec_is_acoustic(ispec) ) count_acoustic = count_acoustic + 1
    if( ispec_is_elastic(ispec) ) count_elastic = count_elastic + 1
    if( ispec_is_poroelastic(ispec) ) count_poroelastic = count_poroelastic + 1

    ! sets flags on global points
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          ! global index
          iglob = ibool(i,j,k,ispec)
          ! sets acoustic flag
          if( ispec_is_acoustic(ispec) ) acoustic_flag(iglob) =  myrank+1
          ! sets elastic flag
          if( ispec_is_elastic(ispec) ) elastic_flag(iglob) =  myrank+1
          ! sets poroelastic flag
          if( ispec_is_poroelastic(ispec) ) poroelastic_flag(iglob) =  myrank+1
        enddo
      enddo
    enddo
  enddo
  call sum_all_i(count_acoustic,inum)
  if( myrank == 0 ) then
    write(IMAIN,*) '     total acoustic elements   :',inum
  endif
  call sum_all_i(count_elastic,inum)
  if( myrank == 0 ) then
    write(IMAIN,*) '     total elastic elements    :',inum
  endif
  call sum_all_i(count_poroelastic,inum)
  if( myrank == 0 ) then
    write(IMAIN,*) '     total poroelastic elements:',inum
  endif

  ! collects contributions from different MPI partitions
  ! sets up MPI communications
  max_nibool_interfaces_ext_mesh = maxval( nibool_interfaces_ext_mesh(:) )
  allocate(ibool_interfaces_ext_mesh_dummy(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
  if( ier /= 0 ) stop 'error allocating array ibool_interfaces_ext_mesh_dummy'
  do i = 1, num_interfaces_ext_mesh
     ibool_interfaces_ext_mesh_dummy(:,i) = ibool_interfaces_ext_mesh(1:max_nibool_interfaces_ext_mesh,i)
  enddo
  ! sums acoustic flags
  if( ACOUSTIC_SIMULATION ) then
    call assemble_MPI_scalar_i_blocking(NPROC,nglob_dummy,acoustic_flag, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh_dummy,&
                        my_neighbours_ext_mesh)
  endif

  ! sums elastic flags
  if( ELASTIC_SIMULATION ) then
    call assemble_MPI_scalar_i_blocking(NPROC,nglob_dummy,elastic_flag, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh_dummy,&
                        my_neighbours_ext_mesh)
  endif

  ! sums poroelastic flags
  if( POROELASTIC_SIMULATION ) then
    call assemble_MPI_scalar_i_blocking(NPROC,nglob_dummy,poroelastic_flag, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh_dummy,&
                        my_neighbours_ext_mesh)
  endif

  ! determines common faces between different domains

  ! acoustic - elastic domain coupling
  if( ACOUSTIC_SIMULATION .and. ELASTIC_SIMULATION ) then
    call get_coupling_surfaces_ac_el(myrank,nspec,ibool,elastic_flag)
  endif

  ! acoustic - poroelastic domain coupling
  if( ACOUSTIC_SIMULATION .and. POROELASTIC_SIMULATION ) then
    call get_coupling_surfaces_ac_poro(myrank,nspec,ibool,acoustic_flag)
  endif

  ! elastic - poroelastic domain coupling
  if( ELASTIC_SIMULATION .and. POROELASTIC_SIMULATION ) then
    call get_coupling_surfaces_el_poro(myrank,nspec,ibool,elastic_flag)
  endif

  ! frees temporary arrays
  deallocate(acoustic_flag,elastic_flag,poroelastic_flag)

  end subroutine get_coupling_surfaces

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_coupling_surfaces_ac_el(myrank,nspec,ibool,elastic_flag)

! determines coupling surface for acoustic-elastic domains

  use generate_databases_par, only: NGNOD2D
  use create_regions_mesh_ext_par
  implicit none

! number of spectral elements in each block
  integer :: myrank,nspec

! arrays with the mesh
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

  integer,dimension(nglob_dummy) :: elastic_flag

! local parameters
  ! (assumes NGLLX=NGLLY=NGLLZ)
  real(kind=CUSTOM_REAL),dimension(NGNOD2D_FOUR_CORNERS) :: xcoord,ycoord,zcoord
  real(kind=CUSTOM_REAL) :: jacobian2Dw_face(NGLLX,NGLLY)
  real(kind=CUSTOM_REAL) :: normal_face(NDIM,NGLLX,NGLLY)
  real(kind=CUSTOM_REAL),dimension(:,:,:),allocatable :: tmp_normal
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: tmp_jacobian2Dw
  integer :: ijk_face(3,NGLLX,NGLLY)
  integer,dimension(:,:,:),allocatable :: tmp_ijk
  integer,dimension(:),allocatable :: tmp_ispec
  logical, dimension(:), allocatable :: mask_ibool

  integer,dimension(NGNOD2D_FOUR_CORNERS) :: iglob_corners_ref
  integer :: ispec,i,j,k,igll,ier,iglob,iglob_midpoint
  integer :: inum,iface_ref

  ! corners indices of reference cube faces
  integer,dimension(3,4),parameter :: iface1_corner_ijk = &
             reshape( (/ 1,1,1, 1,NGLLY,1, 1,NGLLY,NGLLZ, 1,1,NGLLZ /),(/3,4/))   ! xmin
  integer,dimension(3,4),parameter :: iface2_corner_ijk = &
             reshape( (/ NGLLX,1,1, NGLLX,NGLLY,1, NGLLX,NGLLY,NGLLZ, NGLLX,1,NGLLZ  /),(/3,4/))   ! xmax
  integer,dimension(3,4),parameter :: iface3_corner_ijk = &
             reshape( (/ 1,1,1, 1,1,NGLLZ, NGLLX,1,NGLLZ, NGLLX,1,1  /),(/3,4/))   ! ymin
  integer,dimension(3,4),parameter :: iface4_corner_ijk = &
             reshape( (/ 1,NGLLY,1, NGLLX,NGLLY,1, NGLLX,NGLLY,NGLLZ, 1,NGLLY,NGLLZ /),(/3,4/))   ! ymax
  integer,dimension(3,4),parameter :: iface5_corner_ijk = &
             reshape( (/ 1,1,1, 1,NGLLY,1, NGLLX,NGLLY,1, NGLLX,1,1 /),(/3,4/))  ! bottom
  integer,dimension(3,4),parameter :: iface6_corner_ijk = &
             reshape( (/ 1,1,NGLLZ, NGLLX,1,NGLLZ, NGLLX,NGLLY,NGLLZ, 1,NGLLY,NGLLZ  /),(/3,4/))   ! top
  integer,dimension(3,4,6),parameter :: iface_all_corner_ijk = &
             reshape( (/ iface1_corner_ijk,iface2_corner_ijk, &
                 iface3_corner_ijk,iface4_corner_ijk, &
                 iface5_corner_ijk,iface6_corner_ijk /),(/3,4,6/))   ! all faces
  ! midpoint indices for each face (xmin,xmax,ymin,ymax,zmin,zmax)
  integer,dimension(3,6),parameter :: iface_all_midpointijk = &
             reshape( (/ 1,2,2, NGLLX,2,2, 2,1,2, 2,NGLLY,2, 2,2,1, 2,2,NGLLZ  /),(/3,6/))

  ! allocates temporary arrays
  allocate(tmp_normal(NDIM,NGLLSQUARE,nspec*6),stat=ier)
  if( ier /= 0 ) stop 'error allocating array tmp_normal'
  allocate(tmp_jacobian2Dw(NGLLSQUARE,nspec*6),stat=ier)
  if( ier /= 0 ) stop 'error allocating array tmp_jacobian2Dw'
  allocate(tmp_ijk(3,NGLLSQUARE,nspec*6),stat=ier)
  if( ier /= 0 ) stop 'error allocating array tmp_ijk'
  allocate(tmp_ispec(nspec*6),stat=ier)
  if( ier /= 0 ) stop 'error allocating array tmp_ispec'
  tmp_ispec(:) = 0
  tmp_ijk(:,:,:) = 0
  tmp_normal(:,:,:) = 0.0
  tmp_jacobian2Dw(:,:) = 0.0

  allocate(mask_ibool(nglob_dummy),stat=ier)
  if( ier /= 0 ) stop 'error allocating array mask_ibool'
  mask_ibool(:) = .false.

  ! loops over all element faces and
  ! counts number of coupling faces between acoustic and elastic elements
  inum = 0

  ! coupling surfaces: takes point of view from acoustic elements, i.e. if element is acoustic
  !                               and has an elastic neighbor, then this acoustic element is added to
  !                               the coupling elements and its corresponding coupling surface
  !                               ( no matter in which MPI partition it is)
  ! note: we use acoustic elements as reference elements because we will need
  !          density from acoustic element when coupling pressure in case of gravity
  do ispec=1,nspec

    if( ispec_is_acoustic(ispec) ) then

      ! loops over each face
      do iface_ref= 1, 6

        ! takes indices of corners of reference face
        call get_element_corners(ispec,iface_ref,xcoord,ycoord,zcoord,iglob_corners_ref, &
                                ibool,nspec,nglob_dummy,xstore_dummy,ystore_dummy,zstore_dummy, &
                                iface_all_corner_ijk)

        ! checks if face is has an elastic side
        ! note: this check also returns true if layer is only 1 element thick and face points are in touch with elastic elements,
        !       but faces are actually between acoustic elements;
        !       we therefore check again the midpoints and reset them once counted
        if( elastic_flag( iglob_corners_ref(1) ) >= 1 .and. &
           elastic_flag( iglob_corners_ref(2) ) >= 1 .and. &
           elastic_flag( iglob_corners_ref(3) ) >= 1 .and. &
           elastic_flag( iglob_corners_ref(4) ) >= 1) then

          ! reference midpoint on face (used to avoid redundant face counting)
          i = iface_all_midpointijk(1,iface_ref)
          j = iface_all_midpointijk(2,iface_ref)
          k = iface_all_midpointijk(3,iface_ref)
          iglob_midpoint = ibool(i,j,k,ispec)

          ! gets face GLL points i,j,k indices from element face
          call get_element_face_gll_indices(iface_ref,ijk_face,NGLLX,NGLLY)

          ! checks if points on this face are masked already
          if( .not. mask_ibool(iglob_midpoint) .and. &
              (elastic_flag(iglob_midpoint) >= 1) ) then

            ! gets face GLL 2Djacobian, weighted from element face
            call get_jacobian_boundary_face(myrank,nspec, &
                    xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob_dummy, &
                    dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
                    wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                    ispec,iface_ref,jacobian2Dw_face,normal_face,NGLLX,NGLLY,NGNOD2D)

            ! normal convention: points away from acoustic, reference element
            !                                switch normal direction if necessary
            do j=1,NGLLY
              do i=1,NGLLX
                ! directs normals such that they point outwards of element
                call get_element_face_normal(ispec,iface_ref,xcoord,ycoord,zcoord, &
                                            ibool,nspec,nglob_dummy, &
                                            xstore_dummy,ystore_dummy,zstore_dummy, &
                                            normal_face(:,i,j) )
                ! makes sure that it always points away from acoustic element,
                ! otherwise switch direction
                ! note: this should not happen, since we only loop over acoustic elements
                if( ispec_is_elastic(ispec) ) stop 'error acoustic-elastic coupling surface'
              enddo
            enddo

            ! stores informations about this face
            inum = inum + 1
            tmp_ispec(inum) = ispec
            igll = 0
            do j=1,NGLLY
              do i=1,NGLLX
                ! adds all gll points on this face
                igll = igll + 1

                ! do we need to store local i,j,k,ispec info? or only global indices iglob?
                tmp_ijk(:,igll,inum) = ijk_face(:,i,j)

                ! stores weighted jacobian and normals
                tmp_jacobian2Dw(igll,inum) = jacobian2Dw_face(i,j)
                tmp_normal(:,igll,inum) = normal_face(:,i,j)

                ! masks global points ( to avoid redundant counting of faces)
                iglob = ibool(ijk_face(1,i,j),ijk_face(2,i,j),ijk_face(3,i,j),ispec)
                mask_ibool(iglob) = .true.
              enddo
            enddo
          else
            ! assumes to be already collected by lower rank process, masks face points
            do j=1,NGLLY
              do i=1,NGLLX
                iglob = ibool(ijk_face(1,i,j),ijk_face(2,i,j),ijk_face(3,i,j),ispec)
                mask_ibool(iglob) = .true.
              enddo
            enddo
          endif ! mask_ibool
        endif ! elastic_flag
      enddo ! iface_ref
    endif ! ispec_is_acoustic
  enddo ! ispec

! stores completed coupling face informations
!
! note: no need to store material parameters on these coupling points
!          for acoustic-elastic interface
  num_coupling_ac_el_faces = inum
  allocate(coupling_ac_el_normal(NDIM,NGLLSQUARE,num_coupling_ac_el_faces),stat=ier)
  if( ier /= 0 ) stop 'error allocating array coupling_ac_el_normal'
  allocate(coupling_ac_el_jacobian2Dw(NGLLSQUARE,num_coupling_ac_el_faces),stat=ier)
  if( ier /= 0 ) stop 'error allocating array coupling_ac_el_jacobian2Dw'
  allocate(coupling_ac_el_ijk(3,NGLLSQUARE,num_coupling_ac_el_faces),stat=ier)
  if( ier /= 0 ) stop 'error allocating array coupling_ac_el_ijk'
  allocate(coupling_ac_el_ispec(num_coupling_ac_el_faces),stat=ier)
  if( ier /= 0 ) stop 'error allocating array coupling_ac_el_ispec'
  do inum = 1,num_coupling_ac_el_faces
    coupling_ac_el_normal(:,:,inum) = tmp_normal(:,:,inum)
    coupling_ac_el_jacobian2Dw(:,inum) = tmp_jacobian2Dw(:,inum)
    coupling_ac_el_ijk(:,:,inum) = tmp_ijk(:,:,inum)
    coupling_ac_el_ispec(inum) = tmp_ispec(inum)
  enddo

! user output
  call sum_all_i(num_coupling_ac_el_faces,inum)
  if( myrank == 0 ) then
    write(IMAIN,*) '     acoustic-elastic coupling    : total number of faces = ',inum
  endif

  deallocate(mask_ibool)

  end subroutine get_coupling_surfaces_ac_el


!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_coupling_surfaces_ac_poro(myrank,nspec,ibool,acoustic_flag)

! determines coupling surface for acoustic-poroelastic domains

  use generate_databases_par, only: NGNOD2D
  use create_regions_mesh_ext_par
  implicit none

! number of spectral elements in each block
  integer :: myrank,nspec

! arrays with the mesh
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

  integer,dimension(nglob_dummy) :: acoustic_flag

! local parameters
  ! (assumes NGLLX=NGLLY=NGLLZ)
  real(kind=CUSTOM_REAL),dimension(NGNOD2D_FOUR_CORNERS) :: xcoord,ycoord,zcoord
  real(kind=CUSTOM_REAL) :: jacobian2Dw_face(NGLLX,NGLLY)
  real(kind=CUSTOM_REAL) :: normal_face(NDIM,NGLLX,NGLLY)
  real(kind=CUSTOM_REAL),dimension(:,:,:),allocatable :: tmp_normal
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: tmp_jacobian2Dw
  integer :: ijk_face(3,NGLLX,NGLLY)
  integer,dimension(:,:,:),allocatable :: tmp_ijk
  integer,dimension(:),allocatable :: tmp_ispec

  integer,dimension(NGNOD2D_FOUR_CORNERS) :: iglob_corners_ref
  integer :: ispec,i,j,igll,ier
  integer :: inum,iface_ref

  ! corners indices of reference cube faces
  integer,dimension(3,4),parameter :: iface1_corner_ijk = &
             reshape( (/ 1,1,1, 1,NGLLY,1, 1,NGLLY,NGLLZ, 1,1,NGLLZ /),(/3,4/))   ! xmin
  integer,dimension(3,4),parameter :: iface2_corner_ijk = &
             reshape( (/ NGLLX,1,1, NGLLX,NGLLY,1, NGLLX,NGLLY,NGLLZ, NGLLX,1,NGLLZ  /),(/3,4/))   ! xmax
  integer,dimension(3,4),parameter :: iface3_corner_ijk = &
             reshape( (/ 1,1,1, 1,1,NGLLZ, NGLLX,1,NGLLZ, NGLLX,1,1  /),(/3,4/))   ! ymin
  integer,dimension(3,4),parameter :: iface4_corner_ijk = &
             reshape( (/ 1,NGLLY,1, NGLLX,NGLLY,1, NGLLX,NGLLY,NGLLZ, 1,NGLLY,NGLLZ /),(/3,4/))   ! ymax
  integer,dimension(3,4),parameter :: iface5_corner_ijk = &
             reshape( (/ 1,1,1, 1,NGLLY,1, NGLLX,NGLLY,1, NGLLX,1,1 /),(/3,4/))  ! bottom
  integer,dimension(3,4),parameter :: iface6_corner_ijk = &
             reshape( (/ 1,1,NGLLZ, NGLLX,1,NGLLZ, NGLLX,NGLLY,NGLLZ, 1,NGLLY,NGLLZ  /),(/3,4/))   ! top
  integer,dimension(3,4,6),parameter :: iface_all_corner_ijk = &
             reshape( (/ iface1_corner_ijk,iface2_corner_ijk, &
                 iface3_corner_ijk,iface4_corner_ijk, &
                 iface5_corner_ijk,iface6_corner_ijk /),(/3,4,6/))   ! all faces
  ! midpoint indices for each face (xmin,xmax,ymin,ymax,zmin,zmax)
  integer,dimension(3,6),parameter :: iface_all_midpointijk = &
             reshape( (/ 1,2,2, NGLLX,2,2, 2,1,2, 2,NGLLY,2, 2,2,1, 2,2,NGLLZ  /),(/3,6/))   ! top

  ! allocates temporary arrays
  allocate(tmp_normal(NDIM,NGLLSQUARE,nspec*6),stat=ier)
  if( ier /= 0 ) stop 'error allocating array tmp_normal'
  allocate(tmp_jacobian2Dw(NGLLSQUARE,nspec*6),stat=ier)
  if( ier /= 0 ) stop 'error allocating array tmp_jacobian2Dw'
  allocate(tmp_ijk(3,NGLLSQUARE,nspec*6),stat=ier)
  if( ier /= 0 ) stop 'error allocating array tmp_ijk'
  allocate(tmp_ispec(nspec*6),stat=ier)
  if( ier /= 0 ) stop 'error allocating array tmp_ispec'
  tmp_ispec(:) = 0
  tmp_ijk(:,:,:) = 0
  tmp_normal(:,:,:) = 0.0
  tmp_jacobian2Dw(:,:) = 0.0

  ! loops over all element faces and
  ! counts number of coupling faces between acoustic and poroelastic elements
  inum = 0
  do ispec=1,nspec

    if(ispec_is_poroelastic(ispec)) then

      ! loops over each face
      do iface_ref= 1, 6

        ! takes indices of corners of reference face
        call get_element_corners(ispec,iface_ref,xcoord,ycoord,zcoord,iglob_corners_ref, &
                                ibool,nspec,nglob_dummy,xstore_dummy,ystore_dummy,zstore_dummy, &
                                iface_all_corner_ijk)

        ! checks if face has acoustic side
        if( acoustic_flag( iglob_corners_ref(1) ) >= 1 .and. &
           acoustic_flag( iglob_corners_ref(2) ) >= 1 .and. &
           acoustic_flag( iglob_corners_ref(3) ) >= 1 .and. &
           acoustic_flag( iglob_corners_ref(4) ) >= 1) then

          ! gets face GLL points i,j,k indices from poroelastic element face
          call get_element_face_gll_indices(iface_ref,ijk_face,NGLLX,NGLLY)

          ! gets face GLL 2Djacobian, weighted from element face
          call get_jacobian_boundary_face(myrank,nspec, &
                    xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob_dummy, &
                    dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
                    wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                    ispec,iface_ref,jacobian2Dw_face,normal_face,NGLLX,NGLLY,NGNOD2D)

          ! normal convention: points away from acoustic, reference element
          !                                switch normal direction if necessary
          do j=1,NGLLY
            do i=1,NGLLX
              ! directs normals such that they point outwards of element
              call get_element_face_normal(ispec,iface_ref,xcoord,ycoord,zcoord, &
                                          ibool,nspec,nglob_dummy, &
                                          xstore_dummy,ystore_dummy,zstore_dummy, &
                                          normal_face(:,i,j) )
              ! reverse the sign, we know we are in a poroelastic element
              normal_face(:,i,j) = - normal_face(:,i,j)
            enddo
          enddo

          ! stores informations about this face
          inum = inum + 1
          tmp_ispec(inum) = ispec
          igll = 0
          do j=1,NGLLY
            do i=1,NGLLX
              ! adds all gll points on this face
              igll = igll + 1

              ! we need to store local i,j,k,ispec info
              tmp_ijk(:,igll,inum) = ijk_face(:,i,j)

              ! stores weighted jacobian and normals
              tmp_jacobian2Dw(igll,inum) = jacobian2Dw_face(i,j)
              tmp_normal(:,igll,inum) = normal_face(:,i,j)
            enddo
          enddo
        endif ! acoustic_flag
      enddo ! iface_ref
    endif ! ispec_is_poroelastic
  enddo ! ispec

! stores completed coupling face informations
!
! note: for this coupling we need to have access to porous properties. The construction is such
! that i,j,k, & face correspond to poroelastic interface. Note that the normal
! is pointing outward the acoustic element
  num_coupling_ac_po_faces = inum
  allocate(coupling_ac_po_normal(NDIM,NGLLSQUARE,num_coupling_ac_po_faces),stat=ier)
  if( ier /= 0 ) stop 'error allocating array coupling_ac_po_normal'
  allocate(coupling_ac_po_jacobian2Dw(NGLLSQUARE,num_coupling_ac_po_faces),stat=ier)
  if( ier /= 0 ) stop 'error allocating array coupling_ac_po_jacobian2Dw'
  allocate(coupling_ac_po_ijk(3,NGLLSQUARE,num_coupling_ac_po_faces),stat=ier)
  if( ier /= 0 ) stop 'error allocating array coupling_ac_po_ijk'
  allocate(coupling_ac_po_ispec(num_coupling_ac_po_faces),stat=ier)
  if( ier /= 0 ) stop 'error allocating array coupling_ac_po_ispec'
  do inum = 1,num_coupling_ac_po_faces
    coupling_ac_po_normal(:,:,inum) = tmp_normal(:,:,inum)
    coupling_ac_po_jacobian2Dw(:,inum) = tmp_jacobian2Dw(:,inum)
    coupling_ac_po_ijk(:,:,inum) = tmp_ijk(:,:,inum)
    coupling_ac_po_ispec(inum) = tmp_ispec(inum)
  enddo

! user output
  call sum_all_i(num_coupling_ac_po_faces,inum)
  if( myrank == 0 ) then
    write(IMAIN,*) '     acoustic-poroelastic coupling: total number of faces = ',inum
  endif

  end subroutine get_coupling_surfaces_ac_poro

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_coupling_surfaces_el_poro(myrank,nspec,ibool,elastic_flag)

! determines coupling surface for elastic-poroelastic domains

  use generate_databases_par, only: NGNOD2D
  use create_regions_mesh_ext_par
  implicit none

! number of spectral elements in each block
  integer :: myrank,nspec

! arrays with the mesh
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

  integer,dimension(nglob_dummy) :: elastic_flag

! local parameters
  ! (assumes NGLLX=NGLLY=NGLLZ)
  real(kind=CUSTOM_REAL),dimension(NGNOD2D_FOUR_CORNERS) :: xcoord,ycoord,zcoord
  real(kind=CUSTOM_REAL) :: jacobian2Dw_face(NGLLX,NGLLY)
  real(kind=CUSTOM_REAL) :: normal_face(NDIM,NGLLX,NGLLY)
  real(kind=CUSTOM_REAL),dimension(:,:,:),allocatable :: tmp_normal
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: tmp_jacobian2Dw
  integer :: ijk_face_po(3,NGLLX,NGLLY), ijk_face_el(3,NGLLX,NGLLY)
  integer,dimension(:,:,:),allocatable :: tmp_ijk,tmp_ijk_el
  integer,dimension(:),allocatable :: tmp_ispec,tmp_ispec_el

  integer,dimension(NGNOD2D_FOUR_CORNERS) :: iglob_corners_ref,iglob_corners_ref_el
  integer :: ispec,i,j,k,igll,ier
  integer :: ispec_el,ispec_ref_el
  integer :: inum,iface_ref,iface_ref_el,iface_el,icorner

  ! corners indices of reference cube faces
  integer,dimension(3,4),parameter :: iface1_corner_ijk = &
             reshape( (/ 1,1,1, 1,NGLLY,1, 1,NGLLY,NGLLZ, 1,1,NGLLZ /),(/3,4/))   ! xmin
  integer,dimension(3,4),parameter :: iface2_corner_ijk = &
             reshape( (/ NGLLX,1,1, NGLLX,NGLLY,1, NGLLX,NGLLY,NGLLZ, NGLLX,1,NGLLZ  /),(/3,4/))   ! xmax
  integer,dimension(3,4),parameter :: iface3_corner_ijk = &
             reshape( (/ 1,1,1, 1,1,NGLLZ, NGLLX,1,NGLLZ, NGLLX,1,1  /),(/3,4/))   ! ymin
  integer,dimension(3,4),parameter :: iface4_corner_ijk = &
             reshape( (/ 1,NGLLY,1, NGLLX,NGLLY,1, NGLLX,NGLLY,NGLLZ, 1,NGLLY,NGLLZ /),(/3,4/))   ! ymax
  integer,dimension(3,4),parameter :: iface5_corner_ijk = &
             reshape( (/ 1,1,1, 1,NGLLY,1, NGLLX,NGLLY,1, NGLLX,1,1 /),(/3,4/))  ! bottom
  integer,dimension(3,4),parameter :: iface6_corner_ijk = &
             reshape( (/ 1,1,NGLLZ, NGLLX,1,NGLLZ, NGLLX,NGLLY,NGLLZ, 1,NGLLY,NGLLZ  /),(/3,4/))   ! top
  integer,dimension(3,4,6),parameter :: iface_all_corner_ijk = &
             reshape( (/ iface1_corner_ijk,iface2_corner_ijk, &
                 iface3_corner_ijk,iface4_corner_ijk, &
                 iface5_corner_ijk,iface6_corner_ijk /),(/3,4,6/))   ! all faces
  ! midpoint indices for each face (xmin,xmax,ymin,ymax,zmin,zmax)
  integer,dimension(3,6),parameter :: iface_all_midpointijk = &
             reshape( (/ 1,2,2, NGLLX,2,2, 2,1,2, 2,NGLLY,2, 2,2,1, 2,2,NGLLZ  /),(/3,6/))   ! top

! allocates temporary arrays
  allocate(tmp_normal(NDIM,NGLLSQUARE,nspec*6),stat=ier)
  if( ier /= 0 ) stop 'error allocating array tmp_normal'
  allocate(tmp_jacobian2Dw(NGLLSQUARE,nspec*6),stat=ier)
  if( ier /= 0 ) stop 'error allocating array tmp_jacobian2Dw'
  allocate(tmp_ijk(3,NGLLSQUARE,nspec*6),stat=ier)
  if( ier /= 0 ) stop 'error allocating array tmp_ijk'
  allocate(tmp_ijk_el(3,NGLLSQUARE,nspec*6),stat=ier)
  if( ier /= 0 ) stop 'error allocating array tmp_ijk_el'
  allocate(tmp_ispec(nspec*6),stat=ier)
  if( ier /= 0 ) stop 'error allocating array tmp_ispec'
  allocate(tmp_ispec_el(nspec*6),stat=ier)
  if( ier /= 0 ) stop 'error allocating array tmp_ispec_el'
  tmp_ispec(:) = 0
  tmp_ispec_el(:) = 0
  tmp_ijk(:,:,:) = 0
  tmp_ijk_el(:,:,:) = 0
  tmp_normal(:,:,:) = 0.0
  tmp_jacobian2Dw(:,:) = 0.0

  ! loops over all element faces and
  ! counts number of coupling faces between elastic and poroelastic elements
  inum = 0
  do ispec=1,nspec

    if(ispec_is_poroelastic(ispec)) then

      ! loops over each face
      do iface_ref= 1, 6

        ! takes indices of corners of reference face
        call get_element_corners(ispec,iface_ref,xcoord,ycoord,zcoord,iglob_corners_ref, &
                                ibool,nspec,nglob_dummy,xstore_dummy,ystore_dummy,zstore_dummy, &
                                iface_all_corner_ijk)

        ! checks if face has elastic side
        if( elastic_flag( iglob_corners_ref(1) ) >= 1 .and. &
           elastic_flag( iglob_corners_ref(2) ) >= 1 .and. &
           elastic_flag( iglob_corners_ref(3) ) >= 1 .and. &
           elastic_flag( iglob_corners_ref(4) ) >= 1) then

          ! need to find elastic element for coupling
          !
          ! note: this assumes that both, elastic and poroelastic element, are in the same
          !          partition; check with decomposition that this is valid for this mesh partitioning
          do ispec_el=1,nspec
            if(ispec_is_elastic(ispec_el))then
              do iface_el=6,1,-1
                ! takes indices of corners of reference face
                do icorner = 1,NGNOD2D_FOUR_CORNERS
                  i = iface_all_corner_ijk(1,icorner,iface_el)
                  j = iface_all_corner_ijk(2,icorner,iface_el)
                  k = iface_all_corner_ijk(3,icorner,iface_el)
                  ! global reference indices
                  iglob_corners_ref_el(icorner) = ibool(i,j,k,ispec_el)
                enddo

                if ( (iglob_corners_ref(1) == iglob_corners_ref_el(3)) .and. &
                  (iglob_corners_ref(3) == iglob_corners_ref_el(1)) ) then

                  iface_ref_el = iface_el ![Christina Morency]: for some reason this shows a wrong orientation
                                          ! but the calculation is ok.
                  ispec_ref_el = ispec_el

                  ! gets face GLL points i,j,k indices from poroelastic element face
                  call get_element_face_gll_indices(iface_ref,ijk_face_po,NGLLX,NGLLY)
                  ! gets face GLL points i,j,k indices from elastic element face
                  call get_element_face_gll_indices(iface_ref_el,ijk_face_el,NGLLX,NGLLY)

                  ! gets face GLL 2Djacobian, weighted from element face
                  call get_jacobian_boundary_face(myrank,nspec, &
                          xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob_dummy, &
                          dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
                          wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                          ispec,iface_ref,jacobian2Dw_face,normal_face,NGLLX,NGLLY,NGNOD2D)

                  ! normal convention: points away from poroelastic, reference element
                  do j=1,NGLLY
                    do i=1,NGLLX
                      ! directs normals such that they point outwards of poroelastic element
                      call get_element_face_normal(ispec,iface_ref,xcoord,ycoord,zcoord, &
                                                  ibool,nspec,nglob_dummy, &
                                                  xstore_dummy,ystore_dummy,zstore_dummy, &
                                                  normal_face(:,i,j) )
                    enddo
                  enddo

                  ! stores informations about this face
                  inum = inum + 1
                  tmp_ispec(inum) = ispec
                  tmp_ispec_el(inum) = ispec_ref_el
                  igll = 0
                  do j=1,NGLLY
                    do i=1,NGLLX
                      ! adds all gll points on this face
                      igll = igll + 1

                      ! we need to store local i,j,k,ispec info
                      tmp_ijk(:,igll,inum) = ijk_face_po(:,i,j)
                      tmp_ijk_el(:,igll,inum) = ijk_face_el(:,NGLLY-j+1,NGLLX-i+1)

                      ! stores weighted jacobian and normals
                      tmp_jacobian2Dw(igll,inum) = jacobian2Dw_face(i,j)
                      tmp_normal(:,igll,inum) = normal_face(:,i,j)
                    enddo
                  enddo
                endif ! if

              enddo ! do iface_ref_el=1,6
            endif ! if(ispec_is_elastic(ispec_el))then
          enddo ! do ispec_el=1,nspec
        endif ! elastic_flag
      enddo ! iface_ref
    endif ! ispec_is_poroelastic
  enddo ! ispec

! stores completed coupling face informations
!
! note: for this coupling we need to have access to porous properties. The construction is such
! that i,j,k, & face correspond to poroelastic interface. Note that the normal
! is pointing outward the poroelastic element
  num_coupling_el_po_faces = inum
  allocate(coupling_el_po_normal(NDIM,NGLLSQUARE,num_coupling_el_po_faces),stat=ier)
  if( ier /= 0 ) stop 'error allocating array coupling_el_po_normal'
  allocate(coupling_el_po_jacobian2Dw(NGLLSQUARE,num_coupling_el_po_faces),stat=ier)
  if( ier /= 0 ) stop 'error allocating array coupling_el_po_jacobian2Dw'
  allocate(coupling_el_po_ijk(3,NGLLSQUARE,num_coupling_el_po_faces),stat=ier)
  if( ier /= 0 ) stop 'error allocating array coupling_el_po_ijk'
  allocate(coupling_po_el_ijk(3,NGLLSQUARE,num_coupling_el_po_faces),stat=ier)
  if( ier /= 0 ) stop 'error allocating array coupling_po_el_ijk'
  allocate(coupling_el_po_ispec(num_coupling_el_po_faces),stat=ier)
  if( ier /= 0 ) stop 'error allocating array coupling_el_po_ispec'
  allocate(coupling_po_el_ispec(num_coupling_el_po_faces),stat=ier)
  if( ier /= 0 ) stop 'error allocating array coupling_po_el_ispec'
  do inum = 1,num_coupling_el_po_faces
    coupling_el_po_normal(:,:,inum) = tmp_normal(:,:,inum)
    coupling_el_po_jacobian2Dw(:,inum) = tmp_jacobian2Dw(:,inum)
    coupling_el_po_ijk(:,:,inum) = tmp_ijk(:,:,inum)
    coupling_po_el_ijk(:,:,inum) = tmp_ijk_el(:,:,inum)
    coupling_el_po_ispec(inum) = tmp_ispec(inum)
    coupling_po_el_ispec(inum) = tmp_ispec_el(inum)
  enddo

! user output
  call sum_all_i(num_coupling_el_po_faces,inum)
  if( myrank == 0 ) then
    write(IMAIN,*) '     elastic-poroelastic coupling : total number of faces = ',inum
  endif

  end subroutine get_coupling_surfaces_el_poro

