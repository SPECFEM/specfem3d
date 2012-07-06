!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and University of Pau / CNRS / INRIA
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            November 2010
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

  subroutine save_moho_arrays( myrank,nglob,nspec, &
                        nspec2D_moho_ext,ibelm_moho,nodes_ibelm_moho, &
                        nodes_coords_ext_mesh,nnodes_ext_mesh,ibool )

  use create_regions_mesh_ext_par
  implicit none

  integer :: nspec2D_moho_ext
  integer, dimension(nspec2D_moho_ext) :: ibelm_moho
  integer, dimension(4,nspec2D_moho_ext) :: nodes_ibelm_moho

  integer :: myrank,nglob,nspec

  ! data from the external mesh
  integer :: nnodes_ext_mesh
  double precision, dimension(NDIM,nnodes_ext_mesh) :: nodes_coords_ext_mesh

  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

! local parameters
  ! Moho mesh
  real(CUSTOM_REAL), dimension(:,:,:),allocatable :: normal_moho_top
  real(CUSTOM_REAL), dimension(:,:,:),allocatable :: normal_moho_bot
  integer,dimension(:,:,:),allocatable :: ijk_moho_top, ijk_moho_bot
  integer,dimension(:),allocatable :: ibelm_moho_top, ibelm_moho_bot
  integer :: NSPEC2D_MOHO
  logical, dimension(:),allocatable :: is_moho_top, is_moho_bot

  real(kind=CUSTOM_REAL),dimension(NGNOD2D) :: xcoord,ycoord,zcoord
  real(kind=CUSTOM_REAL) :: jacobian2Dw_face(NGLLX,NGLLY)
  real(kind=CUSTOM_REAL) :: normal_face(NDIM,NGLLX,NGLLY)
  real(kind=CUSTOM_REAL),dimension(NDIM):: normal
  integer :: ijk_face(3,NGLLX,NGLLY)

  real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: iglob_normals
  integer,dimension(:),allocatable:: iglob_is_surface

  integer :: imoho_bot,imoho_top
  integer :: ispec2D,ispec,icorner,iface,i,j,k,igll,iglob
  integer :: iglob_midpoint,idirect,counter
  integer :: imoho_top_all,imoho_bot_all,imoho_all

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

  ! temporary arrays for passing information
  allocate(iglob_is_surface(nglob))
  allocate(iglob_normals(NDIM,nglob))
  iglob_is_surface = 0
  iglob_normals = 0._CUSTOM_REAL

  ! loops over given moho surface elements
  do ispec2D=1, nspec2D_moho_ext

    ! gets element id
    ispec = ibelm_moho(ispec2D)

    ! looks for i,j,k indices of GLL points on boundary face
    ! determines element face by given CUBIT corners
    ! (note: uses point locations rather than point indices to find the element face,
    !            because the indices refer no more to the newly indexed ibool array )
    do icorner=1,NGNOD2D
      xcoord(icorner) = nodes_coords_ext_mesh(1,nodes_ibelm_moho(icorner,ispec2D))
      ycoord(icorner) = nodes_coords_ext_mesh(2,nodes_ibelm_moho(icorner,ispec2D))
      zcoord(icorner) = nodes_coords_ext_mesh(3,nodes_ibelm_moho(icorner,ispec2D))
    enddo

    ! sets face id of reference element associated with this face
    call get_element_face_id(ispec,xcoord,ycoord,zcoord, &
                            ibool,nspec,nglob, &
                            xstore_dummy,ystore_dummy,zstore_dummy, &
                            iface)

    ! ijk indices of GLL points for face id
    call get_element_face_gll_indices(iface,ijk_face,NGLLX,NGLLZ)

    ! weighted jacobian and normal
    call get_jacobian_boundary_face(myrank,nspec, &
              xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob,&
              dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
              wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,&
              ispec,iface,jacobian2Dw_face,normal_face,NGLLX,NGLLZ)

    ! normal convention: points away from element
    ! switch normal direction if necessary
    do j=1,NGLLY
      do i=1,NGLLX
          call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                      ibool,nspec,nglob, &
                                      xstore_dummy,ystore_dummy,zstore_dummy, &
                                      normal_face(:,i,j) )
      enddo
    enddo

    ! stores information on global points on moho surface
    igll = 0
    do j=1,NGLLY
      do i=1,NGLLX
        iglob = ibool(ijk_face(1,i,j),ijk_face(2,i,j),ijk_face(3,i,j),ispec)
        ! sets flag
        iglob_is_surface(iglob) = ispec2D
        ! sets normals
        iglob_normals(:,iglob) = normal_face(:,i,j)
      enddo
    enddo
  enddo

  ! stores moho elements
  NSPEC2D_MOHO = nspec2D_moho_ext

  allocate(ibelm_moho_bot(NSPEC2D_MOHO))
  allocate(ibelm_moho_top(NSPEC2D_MOHO))
  allocate(normal_moho_top(NDIM,NGLLSQUARE,NSPEC2D_MOHO))
  allocate(normal_moho_bot(NDIM,NGLLSQUARE,NSPEC2D_MOHO))
  allocate(ijk_moho_bot(3,NGLLSQUARE,NSPEC2D_MOHO))
  allocate(ijk_moho_top(3,NGLLSQUARE,NSPEC2D_MOHO))
  ibelm_moho_bot = 0
  ibelm_moho_top = 0

  ! element flags
  allocate(is_moho_top(nspec))
  allocate(is_moho_bot(nspec))
  is_moho_top = .false.
  is_moho_bot = .false.

  ! finds spectral elements with moho surface
  imoho_top = 0
  imoho_bot = 0
  do ispec=1,nspec

    ! loops over each face
    do iface = 1,6
      ! checks if corners of face on surface
      counter = 0
      do icorner = 1,NGNOD2D
        i = iface_all_corner_ijk(1,icorner,iface)
        j = iface_all_corner_ijk(2,icorner,iface)
        k = iface_all_corner_ijk(3,icorner,iface)
        iglob = ibool(i,j,k,ispec)

        ! checks if point on surface
        if( iglob_is_surface(iglob) > 0 ) then
          counter = counter+1

          ! reference corner coordinates
          xcoord(icorner) = xstore_dummy(iglob)
          ycoord(icorner) = ystore_dummy(iglob)
          zcoord(icorner) = zstore_dummy(iglob)
        endif
      enddo

      ! stores moho informations
      if( counter == NGNOD2D ) then

        ! gets face GLL points i,j,k indices from element face
        call get_element_face_gll_indices(iface,ijk_face,NGLLX,NGLLY)

        ! re-computes face infos
        ! weighted jacobian and normal
        call get_jacobian_boundary_face(myrank,nspec, &
              xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob,&
              dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
              wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,&
              ispec,iface,jacobian2Dw_face,normal_face,NGLLX,NGLLZ)

        ! normal convention: points away from element
        ! switch normal direction if necessary
        do j=1,NGLLZ
          do i=1,NGLLX
            call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                      ibool,nspec,nglob, &
                                      xstore_dummy,ystore_dummy,zstore_dummy, &
                                      normal_face(:,i,j) )
          enddo
        enddo

        ! takes normal stored temporary on a face midpoint
        i = iface_all_midpointijk(1,iface)
        j = iface_all_midpointijk(2,iface)
        k = iface_all_midpointijk(3,iface)
        iglob_midpoint = ibool(i,j,k,ispec)
        normal(:) = iglob_normals(:,iglob_midpoint)

        ! determines whether normal points into element or not (top/bottom distinction)
        call get_element_face_normal_idirect(ispec,iface,xcoord,ycoord,zcoord, &
                              ibool,nspec,nglob, &
                              xstore_dummy,ystore_dummy,zstore_dummy, &
                              normal,idirect )

        ! takes moho surface element id given by id on midpoint
        ispec2D = iglob_is_surface(iglob_midpoint)

        ! sets face infos for bottom (normal points away from element)
        if( idirect == 1 ) then

          ! checks validity
          if( is_moho_bot( ispec) .eqv. .true. ) then
            print*,'error: moho surface geometry bottom'
            print*,'  does not allow for mulitple element faces in kernel computation'
            call exit_mpi(myrank,'error moho bottom elements')
          endif

          imoho_bot = imoho_bot + 1
          is_moho_bot(ispec) = .true.
          ibelm_moho_bot(ispec2D) = ispec

          ! stores on surface gll points (assuming NGLLX = NGLLY = NGLLZ)
          igll = 0
          do j=1,NGLLZ
            do i=1,NGLLX
              igll = igll+1
              ijk_moho_bot(:,igll,ispec2D) = ijk_face(:,i,j)
              normal_moho_bot(:,igll,ispec2D) = normal_face(:,i,j)
            enddo
          enddo

        ! sets face infos for top element
        else if( idirect == 2 ) then

          ! checks validity
          if( is_moho_top( ispec) .eqv. .true. ) then
            print*,'error: moho surface geometry top'
            print*,'  does not allow for mulitple element faces kernel computation'
            call exit_mpi(myrank,'error moho top elements')
          endif

          imoho_top = imoho_top + 1
          is_moho_top(ispec) = .true.
          ibelm_moho_top(ispec2D) = ispec

          ! gll points
          igll = 0
          do j=1,NGLLZ
            do i=1,NGLLX
              igll = igll+1
              ijk_moho_top(:,igll,ispec) = ijk_face(:,i,j)
              ! note: top elements have normal pointing into element
              normal_moho_top(:,igll,ispec) = - normal_face(:,i,j)
            enddo
          enddo
        endif

      endif ! counter

    enddo ! iface

    ! checks validity of top/bottom distinction
    if( is_moho_top(ispec) .and. is_moho_bot(ispec) ) then
      print*,'error: moho surface elements confusing'
      print*,'  element:',ispec,'has top and bottom surface'
      call exit_mpi(myrank,'error moho surface element')
    endif

  enddo ! ispec2D

  ! note: surface e.g. could be at the free-surface and have no top elements etc...
  ! user output
  call sum_all_i( imoho_top, imoho_top_all )
  call sum_all_i( imoho_bot, imoho_bot_all )
  call sum_all_i( NSPEC2D_MOHO, imoho_all )
  if( myrank == 0 ) then
    write(IMAIN,*) '********'
    write(IMAIN,*) 'Moho surface:'
    write(IMAIN,*) '    total surface elements: ',imoho_all
    write(IMAIN,*) '    top elements   :',imoho_top_all
    write(IMAIN,*) '    bottom elements:',imoho_bot_all
    write(IMAIN,*) '********'
  endif

  ! saves moho files: total number of elements, corner points, all points
  open(unit=27,file=prname(1:len_trim(prname))//'ibelm_moho.bin',status='unknown',form='unformatted')
  write(27) NSPEC2D_MOHO
  write(27) ibelm_moho_top
  write(27) ibelm_moho_bot
  write(27) ijk_moho_top
  write(27) ijk_moho_bot
  close(27)
  open(unit=27,file=prname(1:len_trim(prname))//'normal_moho.bin',status='unknown',form='unformatted')
  write(27) normal_moho_top
  write(27) normal_moho_bot
  close(27)
  open(unit=27,file=prname(1:len_trim(prname))//'is_moho.bin',status='unknown',form='unformatted')
  write(27) is_moho_top
  write(27) is_moho_bot
  close(27)

end subroutine save_moho_arrays
