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

  subroutine get_absorbing_boundary(myrank,nspec,ibool, &
                            nodes_coords_ext_mesh,nnodes_ext_mesh, &
                            ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                            nodes_ibelm_xmin,nodes_ibelm_xmax,nodes_ibelm_ymin,nodes_ibelm_ymax, &
                            nodes_ibelm_bottom,nodes_ibelm_top, &
                            nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                            nspec2D_bottom,nspec2D_top)

! determines absorbing boundaries/free-surface, 2D jacobians, face normals for Stacey conditions

  use generate_databases_par, only: STACEY_INSTEAD_OF_FREE_SURFACE, PML_INSTEAD_OF_FREE_SURFACE, NGNOD2D, &
                                      STACEY_ABSORBING_CONDITIONS,PML_CONDITIONS
  use create_regions_mesh_ext_par

  implicit none

! number of spectral elements in each block
  integer :: myrank,nspec

! arrays with the mesh
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

! data from the external mesh
  integer :: nnodes_ext_mesh
  double precision, dimension(NDIM,nnodes_ext_mesh) :: nodes_coords_ext_mesh

! absorbing boundaries (as defined in CUBIT)
  integer  :: nspec2D_xmin, nspec2D_xmax, nspec2D_ymin, nspec2D_ymax, NSPEC2D_BOTTOM, NSPEC2D_TOP
  ! element indices containing a boundary
  integer, dimension(nspec2D_xmin)  :: ibelm_xmin
  integer, dimension(nspec2D_xmax)  :: ibelm_xmax
  integer, dimension(nspec2D_ymin)  :: ibelm_ymin
  integer, dimension(nspec2D_ymax)  :: ibelm_ymax
  integer, dimension(NSPEC2D_BOTTOM)  :: ibelm_bottom
  integer, dimension(NSPEC2D_TOP)  :: ibelm_top

  ! corner node indices of boundary faces coming from CUBIT
  integer, dimension(NGNOD2D,nspec2D_xmin)  :: nodes_ibelm_xmin
  integer, dimension(NGNOD2D,nspec2D_xmax)  :: nodes_ibelm_xmax
  integer, dimension(NGNOD2D,nspec2D_ymin)  :: nodes_ibelm_ymin
  integer, dimension(NGNOD2D,nspec2D_ymax)  :: nodes_ibelm_ymax
  integer, dimension(NGNOD2D,NSPEC2D_BOTTOM)  :: nodes_ibelm_bottom
  integer, dimension(NGNOD2D,NSPEC2D_TOP)  :: nodes_ibelm_top

! local parameters
  ! (assumes NGLLX=NGLLY=NGLLZ)
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: jacobian2Dw_face
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY) :: normal_face
  real(kind=CUSTOM_REAL), dimension(NDIM) :: lnormal

  integer:: ijk_face(3,NGLLX,NGLLY)

  ! face corner locations
  real(kind=CUSTOM_REAL),dimension(NGNOD2D_FOUR_CORNERS) :: xcoord,ycoord,zcoord
  integer  :: ispec,ispec2D,icorner,itop,iabs,iface,igll,i,j,igllfree,ifree

  ! abs face counter
  iabs = 0

  ! xmin
  ijk_face(:,:,:) = 0
  normal_face(:,:,:) = 0.0_CUSTOM_REAL
  jacobian2Dw_face(:,:) = 0.0_CUSTOM_REAL
  do ispec2D = 1, nspec2D_xmin
    ! sets element
    ispec = ibelm_xmin(ispec2D)

    ! looks for i,j,k indices of GLL points on boundary face
    ! determines element face by given CUBIT corners
    do icorner=1,NGNOD2D_FOUR_CORNERS
      xcoord(icorner) = nodes_coords_ext_mesh(1,nodes_ibelm_xmin(icorner,ispec2D))
      ycoord(icorner) = nodes_coords_ext_mesh(2,nodes_ibelm_xmin(icorner,ispec2D))
      zcoord(icorner) = nodes_coords_ext_mesh(3,nodes_ibelm_xmin(icorner,ispec2D))
    enddo

    ! sets face id of reference element associated with this face
    call get_element_face_id(ispec,xcoord,ycoord,zcoord, &
                            ibool,nspec,nglob_dummy, &
                            xstore_dummy,ystore_dummy,zstore_dummy, &
                            iface)

    ! ijk indices of GLL points for face id
    call get_element_face_gll_indices(iface,ijk_face,NGLLX,NGLLZ)

    ! weighted jacobian and normal
    call get_jacobian_boundary_face(myrank,nspec, &
              xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob_dummy,&
              dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
              wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,&
              ispec,iface,jacobian2Dw_face,normal_face,NGLLX,NGLLZ,NGNOD2D)

    ! normal convention: points away from element
    ! switch normal direction if necessary
    do j=1,NGLLZ
      do i=1,NGLLX
        lnormal(:) = normal_face(:,i,j)
        call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                      ibool,nspec,nglob_dummy, &
                                      xstore_dummy,ystore_dummy,zstore_dummy, &
                                      lnormal )
        normal_face(:,i,j) = lnormal(:)
      enddo
    enddo

    ! sets face infos
    iabs = iabs + 1
    abs_boundary_ispec(iabs) = ispec

    ! gll points -- assuming NGLLX = NGLLY = NGLLZ
    igll = 0
    do j=1,NGLLZ
      do i=1,NGLLX
        igll = igll+1
        abs_boundary_ijk(:,igll,iabs) = ijk_face(:,i,j)
        abs_boundary_jacobian2Dw(igll,iabs) = jacobian2Dw_face(i,j)
        abs_boundary_normal(:,igll,iabs) = normal_face(:,i,j)
      enddo
    enddo

  enddo ! nspec2D_xmin

  ! xmax
  ijk_face(:,:,:) = 0
  normal_face(:,:,:) = 0.0_CUSTOM_REAL
  jacobian2Dw_face(:,:) = 0.0_CUSTOM_REAL
  do ispec2D = 1, nspec2D_xmax
    ! sets element
    ispec = ibelm_xmax(ispec2D)

    ! looks for i,j,k indices of GLL points on boundary face
    ! determines element face by given CUBIT corners
    do icorner=1,NGNOD2D_FOUR_CORNERS
      xcoord(icorner) = nodes_coords_ext_mesh(1,nodes_ibelm_xmax(icorner,ispec2D))
      ycoord(icorner) = nodes_coords_ext_mesh(2,nodes_ibelm_xmax(icorner,ispec2D))
      zcoord(icorner) = nodes_coords_ext_mesh(3,nodes_ibelm_xmax(icorner,ispec2D))
    enddo

    ! sets face id of reference element associated with this face
    call get_element_face_id(ispec,xcoord,ycoord,zcoord,&
                              ibool,nspec,nglob_dummy, &
                              xstore_dummy,ystore_dummy,zstore_dummy, &
                              iface )

    ! ijk indices of GLL points on face
    call get_element_face_gll_indices(iface,ijk_face,NGLLX,NGLLZ)

    ! weighted jacobian and normal
    call get_jacobian_boundary_face(myrank,nspec, &
              xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob_dummy,&
              dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
              wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,&
              ispec,iface,jacobian2Dw_face,normal_face,NGLLX,NGLLZ,NGNOD2D)

    ! normal convention: points away from element
    ! switch normal direction if necessary
    do j=1,NGLLZ
      do i=1,NGLLX
        lnormal(:) = normal_face(:,i,j)
        call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                      ibool,nspec,nglob_dummy, &
                                      xstore_dummy,ystore_dummy,zstore_dummy, &
                                      lnormal )
        normal_face(:,i,j) = lnormal(:)
      enddo
    enddo

    ! sets face infos
    iabs = iabs + 1
    abs_boundary_ispec(iabs) = ispec

    ! gll points -- assuming NGLLX = NGLLY = NGLLZ
    igll = 0
    do j=1,NGLLZ
      do i=1,NGLLX
        igll = igll+1
        abs_boundary_ijk(:,igll,iabs) = ijk_face(:,i,j)
        abs_boundary_jacobian2Dw(igll,iabs) = jacobian2Dw_face(i,j)
        abs_boundary_normal(:,igll,iabs) = normal_face(:,i,j)
      enddo
    enddo

  enddo

  ! ymin
  ijk_face(:,:,:) = 0
  normal_face(:,:,:) = 0.0_CUSTOM_REAL
  jacobian2Dw_face(:,:) = 0.0_CUSTOM_REAL
  do ispec2D = 1, nspec2D_ymin
    ! sets element
    ispec = ibelm_ymin(ispec2D)

    ! looks for i,j,k indices of GLL points on boundary face
    ! determines element face by given CUBIT corners
    do icorner=1,NGNOD2D_FOUR_CORNERS
      xcoord(icorner) = nodes_coords_ext_mesh(1,nodes_ibelm_ymin(icorner,ispec2D))
      ycoord(icorner) = nodes_coords_ext_mesh(2,nodes_ibelm_ymin(icorner,ispec2D))
      zcoord(icorner) = nodes_coords_ext_mesh(3,nodes_ibelm_ymin(icorner,ispec2D))
    enddo

    ! sets face id of reference element associated with this face
    call get_element_face_id(ispec,xcoord,ycoord,zcoord,&
                              ibool,nspec,nglob_dummy, &
                              xstore_dummy,ystore_dummy,zstore_dummy, &
                              iface )

    ! ijk indices of GLL points on face
    call get_element_face_gll_indices(iface,ijk_face,NGLLY,NGLLZ)

    ! weighted jacobian and normal
    call get_jacobian_boundary_face(myrank,nspec, &
              xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob_dummy,&
              dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
              wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,&
              ispec,iface,jacobian2Dw_face,normal_face,NGLLY,NGLLZ,NGNOD2D)

    ! normal convention: points away from element
    ! switch normal direction if necessary
    do j=1,NGLLZ
      do i=1,NGLLY
        lnormal(:) = normal_face(:,i,j)
        call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                      ibool,nspec,nglob_dummy, &
                                      xstore_dummy,ystore_dummy,zstore_dummy, &
                                      lnormal )
        normal_face(:,i,j) = lnormal(:)
      enddo
    enddo

    ! sets face infos
    iabs = iabs + 1
    abs_boundary_ispec(iabs) = ispec

    ! gll points -- assuming NGLLX = NGLLY = NGLLZ
    igll = 0
    do j=1,NGLLZ
      do i=1,NGLLY
        igll = igll+1
        abs_boundary_ijk(:,igll,iabs) = ijk_face(:,i,j)
        abs_boundary_jacobian2Dw(igll,iabs) = jacobian2Dw_face(i,j)
        abs_boundary_normal(:,igll,iabs) = normal_face(:,i,j)
      enddo
    enddo

  enddo

  ! ymax
  ijk_face(:,:,:) = 0
  normal_face(:,:,:) = 0.0_CUSTOM_REAL
  jacobian2Dw_face(:,:) = 0.0_CUSTOM_REAL
  do ispec2D = 1, nspec2D_ymax
    ! sets element
    ispec = ibelm_ymax(ispec2D)

    ! looks for i,j,k indices of GLL points on boundary face
    ! determines element face by given CUBIT corners
    do icorner=1,NGNOD2D_FOUR_CORNERS
      xcoord(icorner) = nodes_coords_ext_mesh(1,nodes_ibelm_ymax(icorner,ispec2D))
      ycoord(icorner) = nodes_coords_ext_mesh(2,nodes_ibelm_ymax(icorner,ispec2D))
      zcoord(icorner) = nodes_coords_ext_mesh(3,nodes_ibelm_ymax(icorner,ispec2D))
    enddo

    ! sets face id of reference element associated with this face
    call get_element_face_id(ispec,xcoord,ycoord,zcoord,&
                              ibool,nspec,nglob_dummy, &
                              xstore_dummy,ystore_dummy,zstore_dummy, &
                              iface )

    ! ijk indices of GLL points on face
    call get_element_face_gll_indices(iface,ijk_face,NGLLY,NGLLZ)

    ! weighted jacobian and normal
    call get_jacobian_boundary_face(myrank,nspec, &
              xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob_dummy, &
              dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
              wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
              ispec,iface,jacobian2Dw_face,normal_face,NGLLY,NGLLZ,NGNOD2D)

    ! normal convention: points away from element
    ! switch normal direction if necessary
    do j=1,NGLLZ
      do i=1,NGLLY
        lnormal(:) = normal_face(:,i,j)
        call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                      ibool,nspec,nglob_dummy, &
                                      xstore_dummy,ystore_dummy,zstore_dummy, &
                                      lnormal )
        normal_face(:,i,j) = lnormal(:)
      enddo
    enddo

    ! sets face infos
    iabs = iabs + 1
    abs_boundary_ispec(iabs) = ispec

    ! gll points -- assuming NGLLX = NGLLY = NGLLZ
    igll = 0
    do j=1,NGLLY
      do i=1,NGLLX
        igll = igll+1
        abs_boundary_ijk(:,igll,iabs) = ijk_face(:,i,j)
        abs_boundary_jacobian2Dw(igll,iabs) = jacobian2Dw_face(i,j)
        abs_boundary_normal(:,igll,iabs) = normal_face(:,i,j)
      enddo
    enddo

  enddo

  ! bottom
  ijk_face(:,:,:) = 0
  normal_face(:,:,:) = 0.0_CUSTOM_REAL
  jacobian2Dw_face(:,:) = 0.0_CUSTOM_REAL
  do ispec2D = 1, NSPEC2D_BOTTOM
    ! sets element
    ispec = ibelm_bottom(ispec2D)

    ! looks for i,j,k indices of GLL points on boundary face
    ! determines element face by given CUBIT corners
    do icorner=1,NGNOD2D_FOUR_CORNERS
      xcoord(icorner) = nodes_coords_ext_mesh(1,nodes_ibelm_bottom(icorner,ispec2D))
      ycoord(icorner) = nodes_coords_ext_mesh(2,nodes_ibelm_bottom(icorner,ispec2D))
      zcoord(icorner) = nodes_coords_ext_mesh(3,nodes_ibelm_bottom(icorner,ispec2D))
    enddo

    ! sets face id of reference element associated with this face
    call get_element_face_id(ispec,xcoord,ycoord,zcoord,&
                              ibool,nspec,nglob_dummy, &
                              xstore_dummy,ystore_dummy,zstore_dummy, &
                              iface )

    ! ijk indices of GLL points on face
    call get_element_face_gll_indices(iface,ijk_face,NGLLX,NGLLY)

    ! weighted jacobian and normal
    call get_jacobian_boundary_face(myrank,nspec, &
              xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob_dummy, &
              dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
              wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
              ispec,iface,jacobian2Dw_face,normal_face,NGLLX,NGLLY,NGNOD2D)

    ! normal convention: points away from element
    ! switch normal direction if necessary
    do j=1,NGLLY
      do i=1,NGLLX
        lnormal(:) = normal_face(:,i,j)
        call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                      ibool,nspec,nglob_dummy, &
                                      xstore_dummy,ystore_dummy,zstore_dummy, &
                                      lnormal )
        normal_face(:,i,j) = lnormal(:)
      enddo
    enddo

    ! sets face infos
    iabs = iabs + 1
    abs_boundary_ispec(iabs) = ispec

    ! gll points -- assuming NGLLX = NGLLY = NGLLZ
    igll = 0
    do j=1,NGLLY
      do i=1,NGLLX
        igll = igll+1
        abs_boundary_ijk(:,igll,iabs) = ijk_face(:,i,j)
        abs_boundary_jacobian2Dw(igll,iabs) = jacobian2Dw_face(i,j)
        abs_boundary_normal(:,igll,iabs) = normal_face(:,i,j)
      enddo
    enddo

  enddo

  ! top
  ijk_face(:,:,:) = 0
  normal_face(:,:,:) = 0.0_CUSTOM_REAL
  jacobian2Dw_face(:,:) = 0.0_CUSTOM_REAL
  ! free surface face counter
  ifree = 0

  do ispec2D = 1, NSPEC2D_TOP
    ! sets element
    ispec = ibelm_top(ispec2D)

    ! looks for i,j,k indices of GLL points on boundary face
    ! determines element face by given CUBIT corners
    do icorner=1,NGNOD2D_FOUR_CORNERS
      xcoord(icorner) = nodes_coords_ext_mesh(1,nodes_ibelm_top(icorner,ispec2D))
      ycoord(icorner) = nodes_coords_ext_mesh(2,nodes_ibelm_top(icorner,ispec2D))
      zcoord(icorner) = nodes_coords_ext_mesh(3,nodes_ibelm_top(icorner,ispec2D))
    enddo

    ! sets face id of reference element associated with this face
    call get_element_face_id(ispec,xcoord,ycoord,zcoord,&
                              ibool,nspec,nglob_dummy, &
                              xstore_dummy,ystore_dummy,zstore_dummy, &
                              iface )

    ! ijk indices of GLL points on face
    call get_element_face_gll_indices(iface,ijk_face,NGLLX,NGLLY)

    ! weighted jacobian and normal
    call get_jacobian_boundary_face(myrank,nspec, &
              xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob_dummy, &
              dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
              wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
              ispec,iface,jacobian2Dw_face,normal_face,NGLLX,NGLLY,NGNOD2D)

    ! normal convention: points away from element
    ! switch normal direction if necessary
    do j=1,NGLLY
      do i=1,NGLLX
        lnormal(:) = normal_face(:,i,j)
        call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                      ibool,nspec,nglob_dummy, &
                                      xstore_dummy,ystore_dummy,zstore_dummy, &
                                      lnormal )
        normal_face(:,i,j) = lnormal(:)
      enddo
    enddo

    ! stores surface infos
    if(STACEY_ABSORBING_CONDITIONS)then
       if( .not. STACEY_INSTEAD_OF_FREE_SURFACE ) then
         ! stores free surface
         ! sets face infos
         ifree = ifree + 1
         free_surface_ispec(ifree) = ispec

         ! gll points -- assuming NGLLX = NGLLY = NGLLZ
         igllfree = 0
         do j=1,NGLLY
           do i=1,NGLLX
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
         free_surface_ispec(ifree) = ispec

         ! gll points -- assuming NGLLX = NGLLY = NGLLZ
         igllfree = 0
         do j=1,NGLLY
           do i=1,NGLLX
             igllfree = igllfree+1
             free_surface_ijk(:,igllfree,ifree) = ijk_face(:,i,j)
             free_surface_jacobian2Dw(igllfree,ifree) = jacobian2Dw_face(i,j)
             free_surface_normal(:,igllfree,ifree) = normal_face(:,i,j)
           enddo
         enddo

         ! adds face infos to absorbing boundary surface
         iabs = iabs + 1
         abs_boundary_ispec(iabs) = ispec

         ! gll points -- assuming NGLLX = NGLLY = NGLLZ
         igll = 0
         do j=1,NGLLY
           do i=1,NGLLX
             igll = igll+1
             abs_boundary_ijk(:,igll,iabs) = ijk_face(:,i,j)
             abs_boundary_jacobian2Dw(igll,iabs) = jacobian2Dw_face(i,j)
             abs_boundary_normal(:,igll,iabs) = normal_face(:,i,j)
           enddo
         enddo
       endif

    else if(PML_CONDITIONS)then
       if( .not. PML_INSTEAD_OF_FREE_SURFACE ) then
         ! stores free surface
         ! sets face infos
         ifree = ifree + 1
         free_surface_ispec(ifree) = ispec

         ! gll points -- assuming NGLLX = NGLLY = NGLLZ
         igllfree = 0
         do j=1,NGLLY
           do i=1,NGLLX
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
         free_surface_ispec(ifree) = ispec

         ! gll points -- assuming NGLLX = NGLLY = NGLLZ
         igllfree = 0
         do j=1,NGLLY
           do i=1,NGLLX
             igllfree = igllfree+1
             free_surface_ijk(:,igllfree,ifree) = ijk_face(:,i,j)
             free_surface_jacobian2Dw(igllfree,ifree) = jacobian2Dw_face(i,j)
             free_surface_normal(:,igllfree,ifree) = normal_face(:,i,j)
           enddo
         enddo

         ! adds face infos to absorbing boundary surface
         iabs = iabs + 1
         abs_boundary_ispec(iabs) = ispec

         ! gll points -- assuming NGLLX = NGLLY = NGLLZ
         igll = 0
         do j=1,NGLLY
           do i=1,NGLLX
             igll = igll+1
             abs_boundary_ijk(:,igll,iabs) = ijk_face(:,i,j)
             abs_boundary_jacobian2Dw(igll,iabs) = jacobian2Dw_face(i,j)
             abs_boundary_normal(:,igll,iabs) = normal_face(:,i,j)
           enddo
         enddo
       endif

    else
      ! stores free surface
      ! sets face infos
      ifree = ifree + 1
      free_surface_ispec(ifree) = ispec

      ! gll points -- assuming NGLLX = NGLLY = NGLLZ
      igllfree = 0
      do j=1,NGLLY
       do i=1,NGLLX
         igllfree = igllfree+1
         free_surface_ijk(:,igllfree,ifree) = ijk_face(:,i,j)
         free_surface_jacobian2Dw(igllfree,ifree) = jacobian2Dw_face(i,j)
         free_surface_normal(:,igllfree,ifree) = normal_face(:,i,j)
       enddo
      enddo
    endif

  enddo

  ! checks counters
  if( ifree /= num_free_surface_faces ) then
    print*,'error number of free surface faces:',ifree,num_free_surface_faces
    stop 'error number of free surface faces'
  endif

  if( iabs /= num_abs_boundary_faces ) then
    print*,'error number of absorbing faces:',iabs,num_abs_boundary_faces
    stop 'error number of absorbing faces'
  endif

  call sum_all_i(num_free_surface_faces,itop)
  call sum_all_i(num_abs_boundary_faces,iabs)
  if( myrank == 0 ) then
    write(IMAIN,*) '     absorbing boundary:'
    write(IMAIN,*) '     total number of free faces = ',itop
    write(IMAIN,*) '     total number of faces = ',iabs
    if((PML_CONDITIONS .and. PML_INSTEAD_OF_FREE_SURFACE) .or. &
       (STACEY_ABSORBING_CONDITIONS .and. STACEY_INSTEAD_OF_FREE_SURFACE)) then
       write(IMAIN,*) '     absorbing boundary includes free surface (i.e., top surface converted from free to absorbing)'
    endif
! when users set PML_CONDITIONS and PML_INSTEAD_OF_FREE_SURFACE to be .true. they should also
! provide a non-empty free_or_absorbing_surface_file_zmax file, since we need it to determine ibelm_top(),
! which is the outer boundary of top CPML or Stacey layer.
    if( ((PML_CONDITIONS .and. PML_INSTEAD_OF_FREE_SURFACE) .or. &
         (STACEY_ABSORBING_CONDITIONS .and. STACEY_INSTEAD_OF_FREE_SURFACE)) .and. itop == 0 ) then
       print *,'the free_or_absorbing_surface_file_zmax contains no absorbing element, but Zmax absorption is turned on'
       stop 'error: number of Zmax absorbing elements cannot be zero in free_or_absorbing_surface_file_zmax in such a case'
    endif
  endif

  end subroutine get_absorbing_boundary

