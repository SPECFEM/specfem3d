!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
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

  subroutine get_absorbing_boundary(nspec,ibool, &
                                    nodes_coords_ext_mesh,nnodes_ext_mesh, &
                                    ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                                    nodes_ibelm_xmin,nodes_ibelm_xmax,nodes_ibelm_ymin,nodes_ibelm_ymax, &
                                    nodes_ibelm_bottom,nodes_ibelm_top, &
                                    nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                                    nspec2D_bottom,nspec2D_top)

! determines absorbing boundaries/free-surface, 2D jacobians, face normals for Stacey conditions

  use constants, only: myrank,NGLLX,NGLLY,NGLLZ,NDIM,NGNOD2D_FOUR_CORNERS,IMAIN

  use generate_databases_par, only: STACEY_INSTEAD_OF_FREE_SURFACE, PML_INSTEAD_OF_FREE_SURFACE, NGNOD2D, &
    STACEY_ABSORBING_CONDITIONS,PML_CONDITIONS,BOTTOM_FREE_SURFACE

  use create_regions_mesh_ext_par

  ! injection technique
  use constants, only: INJECTION_TECHNIQUE_IS_DSM
  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE,MESH_A_CHUNK_OF_THE_EARTH,INJECTION_TECHNIQUE_TYPE

  implicit none

  ! number of spectral elements in each block
  integer,intent(in) :: nspec

  ! arrays with the mesh
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool

  ! data from the external mesh
  integer,intent(in) :: nnodes_ext_mesh
  double precision, dimension(NDIM,nnodes_ext_mesh),intent(in) :: nodes_coords_ext_mesh

  ! absorbing boundaries (as defined in CUBIT)
  integer,intent(in)  :: nspec2D_xmin, nspec2D_xmax, nspec2D_ymin, nspec2D_ymax, NSPEC2D_BOTTOM, NSPEC2D_TOP
  ! element indices containing a boundary
  integer, dimension(nspec2D_xmin),intent(in) :: ibelm_xmin
  integer, dimension(nspec2D_xmax),intent(in) :: ibelm_xmax
  integer, dimension(nspec2D_ymin),intent(in) :: ibelm_ymin
  integer, dimension(nspec2D_ymax),intent(in) :: ibelm_ymax
  integer, dimension(NSPEC2D_BOTTOM),intent(in) :: ibelm_bottom
  integer, dimension(NSPEC2D_TOP),intent(in) :: ibelm_top

  ! corner node indices of boundary faces coming from CUBIT
  integer, dimension(NGNOD2D,nspec2D_xmin),intent(in) :: nodes_ibelm_xmin
  integer, dimension(NGNOD2D,nspec2D_xmax),intent(in) :: nodes_ibelm_xmax
  integer, dimension(NGNOD2D,nspec2D_ymin),intent(in) :: nodes_ibelm_ymin
  integer, dimension(NGNOD2D,nspec2D_ymax),intent(in) :: nodes_ibelm_ymax
  integer, dimension(NGNOD2D,NSPEC2D_BOTTOM),intent(in) :: nodes_ibelm_bottom
  integer, dimension(NGNOD2D,NSPEC2D_TOP),intent(in) :: nodes_ibelm_top

  ! local parameters
  ! (assumes NGLLX=NGLLY=NGLLZ)
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: jacobian2Dw_face
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY) :: normal_face
  real(kind=CUSTOM_REAL), dimension(NDIM) :: lnormal

  integer :: ijk_face(3,NGLLX,NGLLY)

  ! face corner locations
  real(kind=CUSTOM_REAL),dimension(NGNOD2D_FOUR_CORNERS) :: xcoord,ycoord,zcoord
  integer :: ispec,ispec2D,icorner,itop,iabsval,iface,igll,i,j,igllfree,ifree

  !! CD CD
  !! additional local parameters For coupling with DSM
  integer :: ier
  logical, dimension(:,:),allocatable :: iboun

  ! corner locations for faces
  real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: xcoord_iboun,ycoord_iboun,zcoord_iboun
  character(len=27) namefile

  ! sets flag in array iboun for elements with an absorbing boundary faces
  if (COUPLE_WITH_INJECTION_TECHNIQUE .or. MESH_A_CHUNK_OF_THE_EARTH) then
    ! allocate temporary flag array
    allocate(iboun(6,nspec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 900')
    allocate(xcoord_iboun(NGNOD2D,6,nspec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 901')
    allocate(ycoord_iboun(NGNOD2D,6,nspec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 902')
    allocate(zcoord_iboun(NGNOD2D,6,nspec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 903')
    if (ier /= 0) stop 'not enough memory to allocate arrays'
    iboun(:,:) = .false.
    xcoord_iboun(:,:,:) = 0.0_CUSTOM_REAL; ycoord_iboun(:,:,:) = 0.0_CUSTOM_REAL; zcoord_iboun(:,:,:) = 0.0_CUSTOM_REAL

    if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_DSM) then
       write(namefile,'(a17,i6.6,a4)') 'xmin_gll_for_dsm_',myrank,'.txt'
       open(123,file=namefile)
       write(123,*) nspec2D_xmin
    endif
  endif

  ! abs face counter
  iabsval = 0

  ! free surface face counter
  ifree = 0

  ! xmin
  if (myrank == 0) then
    write(IMAIN,*) '     boundary xmin   :',nspec2D_xmin
    call flush_IMAIN()
  endif

  ijk_face(:,:,:) = 0
  normal_face(:,:,:) = 0.0_CUSTOM_REAL
  jacobian2Dw_face(:,:) = 0.0_CUSTOM_REAL

  do ispec2D = 1, nspec2D_xmin
    ! sets element
    ispec = ibelm_xmin(ispec2D)

    ! looks for i,j,k indices of GLL points on boundary face
    ! determines element face by given CUBIT corners
    do icorner = 1,NGNOD2D_FOUR_CORNERS
      xcoord(icorner) = nodes_coords_ext_mesh(1,nodes_ibelm_xmin(icorner,ispec2D))
      ycoord(icorner) = nodes_coords_ext_mesh(2,nodes_ibelm_xmin(icorner,ispec2D))
      zcoord(icorner) = nodes_coords_ext_mesh(3,nodes_ibelm_xmin(icorner,ispec2D))
    enddo

    ! sets face id of reference element associated with this face
    call get_element_face_id(ispec,xcoord,ycoord,zcoord, &
                             ibool,nspec,nglob_unique, &
                             xstore_unique,ystore_unique,zstore_unique,iface)

    if (COUPLE_WITH_INJECTION_TECHNIQUE .or. MESH_A_CHUNK_OF_THE_EARTH) then
      iboun(iface,ispec) = .true.
    endif

    ! ijk indices of GLL points on face
    call get_element_face_gll_indices(iface,ijk_face,NGLLX,NGLLZ)

    ! weighted jacobian and normal
    call get_jacobian_boundary_face(nspec, &
                                    xstore_unique,ystore_unique,zstore_unique,ibool,nglob_unique, &
                                    dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
                                    wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                                    ispec,iface,jacobian2Dw_face,normal_face,NGLLX,NGLLZ,NGNOD2D)

    ! normal convention: points away from element
    ! switch normal direction if necessary
    do j = 1,NGLLZ
      do i = 1,NGLLX
        lnormal(:) = normal_face(:,i,j)
        call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                     ibool,nspec,nglob_unique, &
                                     xstore_unique,ystore_unique,zstore_unique, &
                                     lnormal )
        normal_face(:,i,j) = lnormal(:)

        ! writes out xmin boundary point locations
        if ( COUPLE_WITH_INJECTION_TECHNIQUE .and. INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_DSM ) then
          write(123,'(i10,3f20.10)') ispec, xstore_unique(ibool(i,j,1,ispec)), ystore_unique(ibool(i,j,1,ispec)), &
                                            zstore_unique(ibool(i,j,1,ispec))
        endif
      enddo
    enddo

    ! sets face infos
    iabsval = iabsval + 1

    ! checks counter
    if (iabsval > num_abs_boundary_faces) then
      print *,'Error xmin: rank',myrank,'has invalid iabsval counter ',iabsval,'out of ',num_abs_boundary_faces
      stop 'Invalid iabsval counter'
    endif

    abs_boundary_ispec(iabsval) = ispec

    ! GLL points -- assuming NGLLX = NGLLY = NGLLZ
    igll = 0
    do j = 1,NGLLZ
      do i = 1,NGLLX
        igll = igll+1
        abs_boundary_ijk(:,igll,iabsval) = ijk_face(:,i,j)
        abs_boundary_jacobian2Dw(igll,iabsval) = jacobian2Dw_face(i,j)
        abs_boundary_normal(:,igll,iabsval) = normal_face(:,i,j)
      enddo
    enddo
  enddo ! nspec2D_xmin

  if ( COUPLE_WITH_INJECTION_TECHNIQUE .and. (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_DSM) ) then
    close(123)
  endif

  ! xmax
  if (myrank == 0) then
    write(IMAIN,*) '     boundary xmax   :',nspec2D_xmax
    call flush_IMAIN()
  endif

  ijk_face(:,:,:) = 0
  normal_face(:,:,:) = 0.0_CUSTOM_REAL
  jacobian2Dw_face(:,:) = 0.0_CUSTOM_REAL

  do ispec2D = 1, nspec2D_xmax
    ! sets element
    ispec = ibelm_xmax(ispec2D)

    ! looks for i,j,k indices of GLL points on boundary face
    ! determines element face by given CUBIT corners
    do icorner = 1,NGNOD2D_FOUR_CORNERS
      xcoord(icorner) = nodes_coords_ext_mesh(1,nodes_ibelm_xmax(icorner,ispec2D))
      ycoord(icorner) = nodes_coords_ext_mesh(2,nodes_ibelm_xmax(icorner,ispec2D))
      zcoord(icorner) = nodes_coords_ext_mesh(3,nodes_ibelm_xmax(icorner,ispec2D))
    enddo

    ! sets face id of reference element associated with this face
    call get_element_face_id(ispec,xcoord,ycoord,zcoord, &
                             ibool,nspec,nglob_unique, &
                             xstore_unique,ystore_unique,zstore_unique,iface)

    if (COUPLE_WITH_INJECTION_TECHNIQUE .or. MESH_A_CHUNK_OF_THE_EARTH) then
      iboun(iface,ispec) = .true.
    endif

    ! ijk indices of GLL points on face
    call get_element_face_gll_indices(iface,ijk_face,NGLLX,NGLLZ)

    ! weighted jacobian and normal
    call get_jacobian_boundary_face(nspec, &
                                    xstore_unique,ystore_unique,zstore_unique,ibool,nglob_unique, &
                                    dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
                                    wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                                    ispec,iface,jacobian2Dw_face,normal_face,NGLLX,NGLLZ,NGNOD2D)

    ! normal convention: points away from element
    ! switch normal direction if necessary
    do j = 1,NGLLZ
      do i = 1,NGLLX
        lnormal(:) = normal_face(:,i,j)
        call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                     ibool,nspec,nglob_unique, &
                                     xstore_unique,ystore_unique,zstore_unique, &
                                     lnormal )
        normal_face(:,i,j) = lnormal(:)
      enddo
    enddo

    ! sets face infos
    iabsval = iabsval + 1

    ! checks counter
    if (iabsval > num_abs_boundary_faces) then
      print *,'Error xmax: rank',myrank,'has invalid iabsval counter ',iabsval,'out of ',num_abs_boundary_faces
      stop 'Invalid iabsval counter'
    endif

    abs_boundary_ispec(iabsval) = ispec

    ! GLL points -- assuming NGLLX = NGLLY = NGLLZ
    igll = 0
    do j = 1,NGLLZ
      do i = 1,NGLLX
        igll = igll+1
        abs_boundary_ijk(:,igll,iabsval) = ijk_face(:,i,j)
        abs_boundary_jacobian2Dw(igll,iabsval) = jacobian2Dw_face(i,j)
        abs_boundary_normal(:,igll,iabsval) = normal_face(:,i,j)
      enddo
    enddo
  enddo

  ! ymin
  if (myrank == 0) then
    write(IMAIN,*) '     boundary ymin   :',nspec2D_ymin
    call flush_IMAIN()
  endif

  ijk_face(:,:,:) = 0
  normal_face(:,:,:) = 0.0_CUSTOM_REAL
  jacobian2Dw_face(:,:) = 0.0_CUSTOM_REAL

  do ispec2D = 1, nspec2D_ymin
    ! sets element
    ispec = ibelm_ymin(ispec2D)

    ! looks for i,j,k indices of GLL points on boundary face
    ! determines element face by given CUBIT corners
    do icorner = 1,NGNOD2D_FOUR_CORNERS
      xcoord(icorner) = nodes_coords_ext_mesh(1,nodes_ibelm_ymin(icorner,ispec2D))
      ycoord(icorner) = nodes_coords_ext_mesh(2,nodes_ibelm_ymin(icorner,ispec2D))
      zcoord(icorner) = nodes_coords_ext_mesh(3,nodes_ibelm_ymin(icorner,ispec2D))
    enddo

    ! sets face id of reference element associated with this face
    call get_element_face_id(ispec,xcoord,ycoord,zcoord, &
                             ibool,nspec,nglob_unique, &
                             xstore_unique,ystore_unique,zstore_unique,iface)

    if (COUPLE_WITH_INJECTION_TECHNIQUE .or. MESH_A_CHUNK_OF_THE_EARTH) then
      iboun(iface,ispec) = .true.
    endif

    ! ijk indices of GLL points on face
    call get_element_face_gll_indices(iface,ijk_face,NGLLY,NGLLZ)

    ! weighted jacobian and normal
    call get_jacobian_boundary_face(nspec, &
                                    xstore_unique,ystore_unique,zstore_unique,ibool,nglob_unique, &
                                    dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
                                    wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                                    ispec,iface,jacobian2Dw_face,normal_face,NGLLY,NGLLZ,NGNOD2D)

    ! normal convention: points away from element
    ! switch normal direction if necessary
    do j = 1,NGLLZ
      do i = 1,NGLLY
        lnormal(:) = normal_face(:,i,j)
        call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                     ibool,nspec,nglob_unique, &
                                     xstore_unique,ystore_unique,zstore_unique, &
                                     lnormal )
        normal_face(:,i,j) = lnormal(:)
      enddo
    enddo

    ! sets face infos
    iabsval = iabsval + 1

    ! checks counter
    if (iabsval > num_abs_boundary_faces) then
      print *,'Error ymin: rank',myrank,'has invalid iabsval counter ',iabsval,'out of ',num_abs_boundary_faces
      stop 'Invalid iabsval counter'
    endif

    abs_boundary_ispec(iabsval) = ispec

    ! GLL points -- assuming NGLLX = NGLLY = NGLLZ
    igll = 0
    do j = 1,NGLLZ
      do i = 1,NGLLY
        igll = igll+1
        abs_boundary_ijk(:,igll,iabsval) = ijk_face(:,i,j)
        abs_boundary_jacobian2Dw(igll,iabsval) = jacobian2Dw_face(i,j)
        abs_boundary_normal(:,igll,iabsval) = normal_face(:,i,j)
      enddo
    enddo
  enddo

  ! ymax
  if (myrank == 0) then
    write(IMAIN,*) '     boundary ymax   :',nspec2D_ymax
    call flush_IMAIN()
  endif

  ijk_face(:,:,:) = 0
  normal_face(:,:,:) = 0.0_CUSTOM_REAL
  jacobian2Dw_face(:,:) = 0.0_CUSTOM_REAL

  do ispec2D = 1, nspec2D_ymax
    ! sets element
    ispec = ibelm_ymax(ispec2D)

    ! looks for i,j,k indices of GLL points on boundary face
    ! determines element face by given CUBIT corners
    do icorner = 1,NGNOD2D_FOUR_CORNERS
      xcoord(icorner) = nodes_coords_ext_mesh(1,nodes_ibelm_ymax(icorner,ispec2D))
      ycoord(icorner) = nodes_coords_ext_mesh(2,nodes_ibelm_ymax(icorner,ispec2D))
      zcoord(icorner) = nodes_coords_ext_mesh(3,nodes_ibelm_ymax(icorner,ispec2D))
    enddo

    ! sets face id of reference element associated with this face
    call get_element_face_id(ispec,xcoord,ycoord,zcoord, &
                             ibool,nspec,nglob_unique, &
                             xstore_unique,ystore_unique,zstore_unique,iface)

    if (COUPLE_WITH_INJECTION_TECHNIQUE .or. MESH_A_CHUNK_OF_THE_EARTH) then
      iboun(iface,ispec) = .true.
    endif

    ! ijk indices of GLL points on face
    call get_element_face_gll_indices(iface,ijk_face,NGLLY,NGLLZ)

    ! weighted jacobian and normal
    call get_jacobian_boundary_face(nspec, &
                                    xstore_unique,ystore_unique,zstore_unique,ibool,nglob_unique, &
                                    dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
                                    wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                                    ispec,iface,jacobian2Dw_face,normal_face,NGLLY,NGLLZ,NGNOD2D)

    ! normal convention: points away from element
    ! switch normal direction if necessary
    do j = 1,NGLLZ
      do i = 1,NGLLY
        lnormal(:) = normal_face(:,i,j)
        call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                     ibool,nspec,nglob_unique, &
                                     xstore_unique,ystore_unique,zstore_unique, &
                                     lnormal )
        normal_face(:,i,j) = lnormal(:)
      enddo
    enddo

    ! sets face infos
    iabsval = iabsval + 1

    ! checks counter
    if (iabsval > num_abs_boundary_faces) then
      print *,'Error ymax: rank',myrank,'has invalid iabsval counter ',iabsval,'out of ',num_abs_boundary_faces
      stop 'Invalid iabsval counter'
    endif

    abs_boundary_ispec(iabsval) = ispec

    ! GLL points -- assuming NGLLX = NGLLY = NGLLZ
    igll = 0
    do j = 1,NGLLY
      do i = 1,NGLLX
        igll = igll+1
        abs_boundary_ijk(:,igll,iabsval) = ijk_face(:,i,j)
        abs_boundary_jacobian2Dw(igll,iabsval) = jacobian2Dw_face(i,j)
        abs_boundary_normal(:,igll,iabsval) = normal_face(:,i,j)
      enddo
    enddo
  enddo

  ! bottom
  if (myrank == 0) then
    write(IMAIN,*) '     boundary bottom :',NSPEC2D_BOTTOM
    call flush_IMAIN()
  endif

  ijk_face(:,:,:) = 0
  normal_face(:,:,:) = 0.0_CUSTOM_REAL
  jacobian2Dw_face(:,:) = 0.0_CUSTOM_REAL

  do ispec2D = 1, NSPEC2D_BOTTOM
    ! sets element
    ispec = ibelm_bottom(ispec2D)

    ! looks for i,j,k indices of GLL points on boundary face
    ! determines element face by given CUBIT corners
    do icorner = 1,NGNOD2D_FOUR_CORNERS
      xcoord(icorner) = nodes_coords_ext_mesh(1,nodes_ibelm_bottom(icorner,ispec2D))
      ycoord(icorner) = nodes_coords_ext_mesh(2,nodes_ibelm_bottom(icorner,ispec2D))
      zcoord(icorner) = nodes_coords_ext_mesh(3,nodes_ibelm_bottom(icorner,ispec2D))
    enddo

    ! sets face id of reference element associated with this face
    call get_element_face_id(ispec,xcoord,ycoord,zcoord, &
                             ibool,nspec,nglob_unique, &
                             xstore_unique,ystore_unique,zstore_unique,iface)

    if (COUPLE_WITH_INJECTION_TECHNIQUE .or. MESH_A_CHUNK_OF_THE_EARTH) then
      iboun(iface,ispec) = .true.
    endif

    ! ijk indices of GLL points on face
    call get_element_face_gll_indices(iface,ijk_face,NGLLX,NGLLY)

    ! weighted jacobian and normal
    call get_jacobian_boundary_face(nspec, &
                                    xstore_unique,ystore_unique,zstore_unique,ibool,nglob_unique, &
                                    dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
                                    wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                                    ispec,iface,jacobian2Dw_face,normal_face,NGLLX,NGLLY,NGNOD2D)

    ! normal convention: points away from element
    ! switch normal direction if necessary
    do j = 1,NGLLY
      do i = 1,NGLLX
        lnormal(:) = normal_face(:,i,j)
        call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                     ibool,nspec,nglob_unique, &
                                     xstore_unique,ystore_unique,zstore_unique, &
                                     lnormal )
        normal_face(:,i,j) = lnormal(:)
      enddo
    enddo

    ! use bottom free surface instead of absorbing Stacey condition
    if (BOTTOM_FREE_SURFACE) then
       ! stores free surface
       ! sets face infos
       ifree = ifree + 1

       ! checks counter
       if (ifree > num_free_surface_faces) then
         print *,'Error bottom: rank',myrank,'has invalid ifree counter ',ifree,'out of ',num_free_surface_faces
         stop 'Invalid ifree counter'
       endif

       free_surface_ispec(ifree) = ispec

       ! GLL points -- assuming NGLLX = NGLLY = NGLLZ
       igllfree = 0
       do j = 1,NGLLY
          do i = 1,NGLLX
             igllfree = igllfree+1
             free_surface_ijk(:,igllfree,ifree) = ijk_face(:,i,j)
             free_surface_jacobian2Dw(igllfree,ifree) = jacobian2Dw_face(i,j)
             free_surface_normal(:,igllfree,ifree) = normal_face(:,i,j)
          enddo
       enddo
    else
       ! sets face infos
       iabsval = iabsval + 1

       ! checks counter
       if (iabsval > num_abs_boundary_faces) then
         print *,'Error bottom: rank',myrank,'has invalid iabsval counter ',iabsval,'out of ',num_abs_boundary_faces
         stop 'Invalid iabsval counter'
       endif

       abs_boundary_ispec(iabsval) = ispec

       ! GLL points -- assuming NGLLX = NGLLY = NGLLZ
       igll = 0
       do j = 1,NGLLY
          do i = 1,NGLLX
             igll = igll+1
             abs_boundary_ijk(:,igll,iabsval) = ijk_face(:,i,j)
             abs_boundary_jacobian2Dw(igll,iabsval) = jacobian2Dw_face(i,j)
             abs_boundary_normal(:,igll,iabsval) = normal_face(:,i,j)
          enddo
       enddo
    endif
  enddo   ! NSPEC2D_BOTTOM

  ! top
  if (myrank == 0) then
    write(IMAIN,*) '     boundary top    :',NSPEC2D_TOP
    call flush_IMAIN()
  endif

  ijk_face(:,:,:) = 0
  normal_face(:,:,:) = 0.0_CUSTOM_REAL
  jacobian2Dw_face(:,:) = 0.0_CUSTOM_REAL

  do ispec2D = 1, NSPEC2D_TOP
    ! sets element
    ispec = ibelm_top(ispec2D)

    ! looks for i,j,k indices of GLL points on boundary face
    ! determines element face by given CUBIT corners
    do icorner = 1,NGNOD2D_FOUR_CORNERS
      xcoord(icorner) = nodes_coords_ext_mesh(1,nodes_ibelm_top(icorner,ispec2D))
      ycoord(icorner) = nodes_coords_ext_mesh(2,nodes_ibelm_top(icorner,ispec2D))
      zcoord(icorner) = nodes_coords_ext_mesh(3,nodes_ibelm_top(icorner,ispec2D))
    enddo

    ! sets face id of reference element associated with this face
    call get_element_face_id(ispec,xcoord,ycoord,zcoord, &
                             ibool,nspec,nglob_unique, &
                             xstore_unique,ystore_unique,zstore_unique,iface)

    if (COUPLE_WITH_INJECTION_TECHNIQUE .or. MESH_A_CHUNK_OF_THE_EARTH) then
      iboun(iface,ispec) = .true.
    endif

    ! ijk indices of GLL points on face
    call get_element_face_gll_indices(iface,ijk_face,NGLLX,NGLLY)

    ! weighted jacobian and normal
    call get_jacobian_boundary_face(nspec, &
                                    xstore_unique,ystore_unique,zstore_unique,ibool,nglob_unique, &
                                    dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
                                    wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                                    ispec,iface,jacobian2Dw_face,normal_face,NGLLX,NGLLY,NGNOD2D)

    ! normal convention: points away from element
    ! switch normal direction if necessary
    do j = 1,NGLLY
      do i = 1,NGLLX
        lnormal(:) = normal_face(:,i,j)
        call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                     ibool,nspec,nglob_unique, &
                                     xstore_unique,ystore_unique,zstore_unique, &
                                     lnormal )
        normal_face(:,i,j) = lnormal(:)
      enddo
    enddo

    ! stores surface infos
    if (STACEY_ABSORBING_CONDITIONS) then
      if (.not. STACEY_INSTEAD_OF_FREE_SURFACE) then
        ! stores free surface
        ! sets face infos
        ifree = ifree + 1

        ! checks counter
        if (ifree > num_free_surface_faces) then
          print *,'Error top: rank',myrank,'has invalid ifree counter ',ifree,'out of ',num_free_surface_faces
          stop 'Invalid ifree counter'
        endif

        free_surface_ispec(ifree) = ispec

        ! GLL points -- assuming NGLLX = NGLLY = NGLLZ
        igllfree = 0
        do j = 1,NGLLY
          do i = 1,NGLLX
            igllfree = igllfree+1
            free_surface_ijk(:,igllfree,ifree) = ijk_face(:,i,j)
            free_surface_jacobian2Dw(igllfree,ifree) = jacobian2Dw_face(i,j)
            free_surface_normal(:,igllfree,ifree) = normal_face(:,i,j)
          enddo
        enddo
      else
        if (.not. BOTTOM_FREE_SURFACE) then
          ! stores free surface and adds it also to absorbing boundaries
          ! sets face infos
          ifree = ifree + 1

          ! checks counter
          if (ifree > num_free_surface_faces) then
            print *,'Error top: rank',myrank,'has invalid ifree counter ',ifree,'out of ',num_free_surface_faces
            stop 'Invalid ifree counter'
          endif

          free_surface_ispec(ifree) = ispec

          ! GLL points -- assuming NGLLX = NGLLY = NGLLZ
          igllfree = 0
          do j = 1,NGLLY
            do i = 1,NGLLX
              igllfree = igllfree+1
              free_surface_ijk(:,igllfree,ifree) = ijk_face(:,i,j)
              free_surface_jacobian2Dw(igllfree,ifree) = jacobian2Dw_face(i,j)
              free_surface_normal(:,igllfree,ifree) = normal_face(:,i,j)
            enddo
          enddo
        endif

        ! adds face infos to absorbing boundary surface
        iabsval = iabsval + 1

        ! checks counter
        if (iabsval > num_abs_boundary_faces) then
          print *,'Error top: rank',myrank,'has invalid iabsval counter ',iabsval,'out of ',num_abs_boundary_faces
          stop 'Invalid iabsval counter'
        endif

        abs_boundary_ispec(iabsval) = ispec

        ! GLL points -- assuming NGLLX = NGLLY = NGLLZ
        igll = 0
        do j = 1,NGLLY
          do i = 1,NGLLX
            igll = igll+1
            abs_boundary_ijk(:,igll,iabsval) = ijk_face(:,i,j)
            abs_boundary_jacobian2Dw(igll,iabsval) = jacobian2Dw_face(i,j)
            abs_boundary_normal(:,igll,iabsval) = normal_face(:,i,j)
          enddo
        enddo
      endif

    else if (PML_CONDITIONS) then
      if (.not. PML_INSTEAD_OF_FREE_SURFACE) then
        ! stores free surface
        ! sets face infos
        ifree = ifree + 1

        ! checks counter
        if (ifree > num_free_surface_faces) then
          print *,'Error top: rank',myrank,'has invalid ifree counter ',ifree,'out of ',num_free_surface_faces
          stop 'Invalid ifree counter'
        endif

        free_surface_ispec(ifree) = ispec

        ! GLL points -- assuming NGLLX = NGLLY = NGLLZ
        igllfree = 0
        do j = 1,NGLLY
          do i = 1,NGLLX
            igllfree = igllfree+1
            free_surface_ijk(:,igllfree,ifree) = ijk_face(:,i,j)
            free_surface_jacobian2Dw(igllfree,ifree) = jacobian2Dw_face(i,j)
            free_surface_normal(:,igllfree,ifree) = normal_face(:,i,j)
          enddo
        enddo
      else
        ! stores free surface and adds it also to absorbing boundaries
        ! sets face infos
        ifree = ifree + 1

        ! checks counter
        if (ifree > num_free_surface_faces) then
          print *,'Error top: rank',myrank,'has invalid ifree counter ',ifree,'out of ',num_free_surface_faces
          stop 'Invalid ifree counter'
        endif

        free_surface_ispec(ifree) = ispec

        ! GLL points -- assuming NGLLX = NGLLY = NGLLZ
        igllfree = 0
        do j = 1,NGLLY
          do i = 1,NGLLX
            igllfree = igllfree+1
            free_surface_ijk(:,igllfree,ifree) = ijk_face(:,i,j)
            free_surface_jacobian2Dw(igllfree,ifree) = jacobian2Dw_face(i,j)
            free_surface_normal(:,igllfree,ifree) = normal_face(:,i,j)
          enddo
        enddo

        ! adds face infos to absorbing boundary surface
        iabsval = iabsval + 1

        ! checks counter
        if (iabsval > num_abs_boundary_faces) then
          print *,'Error top: rank',myrank,'has invalid iabsval counter ',iabsval,'out of ',num_abs_boundary_faces
          stop 'Invalid iabsval counter'
        endif

        abs_boundary_ispec(iabsval) = ispec

        ! GLL points -- assuming NGLLX = NGLLY = NGLLZ
        igll = 0
        do j = 1,NGLLY
          do i = 1,NGLLX
            igll = igll+1
            abs_boundary_ijk(:,igll,iabsval) = ijk_face(:,i,j)
            abs_boundary_jacobian2Dw(igll,iabsval) = jacobian2Dw_face(i,j)
            abs_boundary_normal(:,igll,iabsval) = normal_face(:,i,j)
          enddo
        enddo
      endif

    else
      ! stores free surface
      ! sets face infos
      ifree = ifree + 1
      free_surface_ispec(ifree) = ispec

      ! GLL points -- assuming NGLLX = NGLLY = NGLLZ
      igllfree = 0
      do j = 1,NGLLY
        do i = 1,NGLLX
          igllfree = igllfree+1
          free_surface_ijk(:,igllfree,ifree) = ijk_face(:,i,j)
          free_surface_jacobian2Dw(igllfree,ifree) = jacobian2Dw_face(i,j)
          free_surface_normal(:,igllfree,ifree) = normal_face(:,i,j)
        enddo
      enddo
    endif
  enddo   ! NSPEC2D_TOP

  ! checks counters
  if (ifree /= num_free_surface_faces) then
    print *,'error number of free surface faces:',ifree,num_free_surface_faces
    stop 'error number of free surface faces'
  endif

  if (iabsval /= num_abs_boundary_faces) then
    print *,'error number of absorbing faces:',iabsval,num_abs_boundary_faces
    stop 'error number of absorbing faces'
  endif

  call sum_all_i(num_free_surface_faces,itop)
  call sum_all_i(num_abs_boundary_faces,iabsval)
  if (myrank == 0) then
    write(IMAIN,*) '     absorbing boundary:'
    write(IMAIN,*) '     total number of free faces = ',itop
    write(IMAIN,*) '     total number of faces      = ',iabsval
    if ((PML_CONDITIONS .and. PML_INSTEAD_OF_FREE_SURFACE) .or. &
       (STACEY_ABSORBING_CONDITIONS .and. STACEY_INSTEAD_OF_FREE_SURFACE)) then
       write(IMAIN,*) '     absorbing boundary includes free surface (i.e., top surface converted from free to absorbing)'
    endif

    ! when users set PML_CONDITIONS and PML_INSTEAD_OF_FREE_SURFACE to be .true. they should also
    ! provide a non-empty free_or_absorbing_surface_file_zmax file, since we need it to determine ibelm_top(),
    ! which is the outer boundary of top CPML or Stacey layer.
    if (((PML_CONDITIONS .and. PML_INSTEAD_OF_FREE_SURFACE) .or. &
         (STACEY_ABSORBING_CONDITIONS .and. STACEY_INSTEAD_OF_FREE_SURFACE)) .and. itop == 0) then
       print *,'the free_or_absorbing_surface_file_zmax contains no absorbing element, but Zmax absorption is turned on'
       stop 'error: number of Zmax absorbing elements cannot be zero in free_or_absorbing_surface_file_zmax in such a case'
    endif
    call flush_IMAIN()
  endif

  end subroutine get_absorbing_boundary

