!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  1 . 4
!               ---------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
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

  subroutine get_absorbing_boundary(myrank,nspec,nglob,ibool, &
                            nodes_coords_ext_mesh,nnodes_ext_mesh, &
                            ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                            nodes_ibelm_xmin,nodes_ibelm_xmax,nodes_ibelm_ymin,nodes_ibelm_ymax, &
                            nodes_ibelm_bottom,nodes_ibelm_top, &
                            nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                            nspec2D_bottom,nspec2D_top)

! determines absorbing boundaries/free-surface, 2D jacobians, face normals for Stacey conditions

  use create_regions_mesh_ext_par 
  implicit none

! number of spectral elements in each block
  integer :: myrank,nspec,nglob

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
  integer, dimension(4,nspec2D_xmin)  :: nodes_ibelm_xmin  
  integer, dimension(4,nspec2D_xmax)  :: nodes_ibelm_xmax
  integer, dimension(4,nspec2D_ymin)  :: nodes_ibelm_ymin
  integer, dimension(4,nspec2D_ymax)  :: nodes_ibelm_ymax
  integer, dimension(4,NSPEC2D_BOTTOM)  :: nodes_ibelm_bottom
  integer, dimension(4,NSPEC2D_TOP)  :: nodes_ibelm_top
  
! local parameters
  logical, dimension(:,:),allocatable :: iboun   ! pll 

  ! (assumes NGLLX=NGLLY=NGLLZ)
  real(kind=CUSTOM_REAL) :: jacobian2Dw_face(NGLLX,NGLLY)
  real(kind=CUSTOM_REAL) :: normal_face(NDIM,NGLLX,NGLLY)
  integer:: ijk_face(3,NGLLX,NGLLY)
  
  ! corner locations for faces
  real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: xcoord_iboun,ycoord_iboun,zcoord_iboun
  
  ! face corner locations
  real(kind=CUSTOM_REAL),dimension(NGNOD2D) :: xcoord,ycoord,zcoord    
  integer  :: ispec,ispec2D,icorner,ier,iabs,iface,igll,i,j,igllfree,ifree
  
! allocate temporary flag array
  allocate(iboun(6,nspec), &
          xcoord_iboun(NGNOD2D,6,nspec), &
          ycoord_iboun(NGNOD2D,6,nspec), &
          zcoord_iboun(NGNOD2D,6,nspec),stat=ier)
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')
  
! sets flag in array iboun for elements with an absorbing boundary faces
  iboun(:,:) = .false. 

! abs face counter  
  iabs = 0
  
  ! xmin   
  do ispec2D = 1, nspec2D_xmin 
    ! sets element 
    ispec = ibelm_xmin(ispec2D)
     
    !if(myrank == 0 ) print*,'xmin:',ispec2D,ispec
    
    ! looks for i,j,k indices of GLL points on boundary face
    ! determines element face by given CUBIT corners
    do icorner=1,NGNOD2D
      xcoord(icorner) = nodes_coords_ext_mesh(1,nodes_ibelm_xmin(icorner,ispec2D))
      ycoord(icorner) = nodes_coords_ext_mesh(2,nodes_ibelm_xmin(icorner,ispec2D))
      zcoord(icorner) = nodes_coords_ext_mesh(3,nodes_ibelm_xmin(icorner,ispec2D))
      !print*,'corner look:',icorner,xcoord(icorner),ycoord(icorner),zcoord(icorner)
    enddo
    
    ! sets face id of reference element associated with this face
    call get_element_face_id(ispec,xcoord,ycoord,zcoord, &
                            ibool,nspec,nglob, &
                            xstore_dummy,ystore_dummy,zstore_dummy, &
                            iface)

    iboun(iface,ispec) = .true. 

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
    do j=1,NGLLZ
      do i=1,NGLLX
          call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                      ibool,nspec,nglob, &
                                      xstore_dummy,ystore_dummy,zstore_dummy, &
                                      normal_face(:,i,j) )
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
  do ispec2D = 1, nspec2D_xmax 
    ! sets element
    ispec = ibelm_xmax(ispec2D)
     
    ! looks for i,j,k indices of GLL points on boundary face
    ! determines element face by given CUBIT corners
    do icorner=1,NGNOD2D
      xcoord(icorner) = nodes_coords_ext_mesh(1,nodes_ibelm_xmax(icorner,ispec2D))
      ycoord(icorner) = nodes_coords_ext_mesh(2,nodes_ibelm_xmax(icorner,ispec2D))
      zcoord(icorner) = nodes_coords_ext_mesh(3,nodes_ibelm_xmax(icorner,ispec2D))
    enddo
    
    ! sets face id of reference element associated with this face
    call get_element_face_id(ispec,xcoord,ycoord,zcoord,&
                              ibool,nspec,nglob, &
                              xstore_dummy,ystore_dummy,zstore_dummy, &
                              iface )   
    iboun(iface,ispec) = .true. 
                              
    ! ijk indices of GLL points on face
    call get_element_face_gll_indices(iface,ijk_face,NGLLX,NGLLZ)
    
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
  do ispec2D = 1, nspec2D_ymin 
    ! sets element 
    ispec = ibelm_ymin(ispec2D)
     
    ! looks for i,j,k indices of GLL points on boundary face
    ! determines element face by given CUBIT corners
    do icorner=1,NGNOD2D
      xcoord(icorner) = nodes_coords_ext_mesh(1,nodes_ibelm_ymin(icorner,ispec2D))
      ycoord(icorner) = nodes_coords_ext_mesh(2,nodes_ibelm_ymin(icorner,ispec2D))
      zcoord(icorner) = nodes_coords_ext_mesh(3,nodes_ibelm_ymin(icorner,ispec2D))
    enddo
    
    ! sets face id of reference element associated with this face
    call get_element_face_id(ispec,xcoord,ycoord,zcoord,&
                              ibool,nspec,nglob, &
                              xstore_dummy,ystore_dummy,zstore_dummy, &
                              iface )   
    iboun(iface,ispec) = .true. 
                              
    ! ijk indices of GLL points on face
    call get_element_face_gll_indices(iface,ijk_face,NGLLY,NGLLZ)

    ! weighted jacobian and normal                          
    call get_jacobian_boundary_face(myrank,nspec, & 
              xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob,&
              dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
              wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,&                                          
              ispec,iface,jacobian2Dw_face,normal_face,NGLLY,NGLLZ)                              

    ! normal convention: points away from element
    ! switch normal direction if necessary
    do j=1,NGLLZ
      do i=1,NGLLY
          call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                      ibool,nspec,nglob, &
                                      xstore_dummy,ystore_dummy,zstore_dummy, &
                                      normal_face(:,i,j) )
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
  do ispec2D = 1, nspec2D_ymax 
    ! sets element 
    ispec = ibelm_ymax(ispec2D)
     
    ! looks for i,j,k indices of GLL points on boundary face
    ! determines element face by given CUBIT corners
    do icorner=1,NGNOD2D
      xcoord(icorner) = nodes_coords_ext_mesh(1,nodes_ibelm_ymax(icorner,ispec2D))
      ycoord(icorner) = nodes_coords_ext_mesh(2,nodes_ibelm_ymax(icorner,ispec2D))
      zcoord(icorner) = nodes_coords_ext_mesh(3,nodes_ibelm_ymax(icorner,ispec2D))
    enddo
    
    ! sets face id of reference element associated with this face
    call get_element_face_id(ispec,xcoord,ycoord,zcoord,&
                              ibool,nspec,nglob, &
                              xstore_dummy,ystore_dummy,zstore_dummy, &
                              iface )   
    iboun(iface,ispec) = .true. 
                              
    ! ijk indices of GLL points on face
    call get_element_face_gll_indices(iface,ijk_face,NGLLY,NGLLZ)                              

    ! weighted jacobian and normal                          
    call get_jacobian_boundary_face(myrank,nspec, &
              xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob, &
              dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
              wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
              ispec,iface,jacobian2Dw_face,normal_face,NGLLY,NGLLZ) 

    ! normal convention: points away from element
    ! switch normal direction if necessary
    do j=1,NGLLZ
      do i=1,NGLLY
          call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                      ibool,nspec,nglob, &
                                      xstore_dummy,ystore_dummy,zstore_dummy, &
                                      normal_face(:,i,j) )
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
  do ispec2D = 1, NSPEC2D_BOTTOM
    ! sets element 
    ispec = ibelm_bottom(ispec2D)
     
    ! looks for i,j,k indices of GLL points on boundary face
    ! determines element face by given CUBIT corners
    do icorner=1,NGNOD2D
      xcoord(icorner) = nodes_coords_ext_mesh(1,nodes_ibelm_bottom(icorner,ispec2D))
      ycoord(icorner) = nodes_coords_ext_mesh(2,nodes_ibelm_bottom(icorner,ispec2D))
      zcoord(icorner) = nodes_coords_ext_mesh(3,nodes_ibelm_bottom(icorner,ispec2D))
    enddo
    
    ! sets face id of reference element associated with this face
    call get_element_face_id(ispec,xcoord,ycoord,zcoord,&
                              ibool,nspec,nglob, &
                              xstore_dummy,ystore_dummy,zstore_dummy, &
                              iface )   
    iboun(iface,ispec) = .true. 
                              
    ! ijk indices of GLL points on face
    call get_element_face_gll_indices(iface,ijk_face,NGLLX,NGLLY)
    
    ! weighted jacobian and normal                          
    call get_jacobian_boundary_face(myrank,nspec, &
              xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob, &
              dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
              wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
              ispec,iface,jacobian2Dw_face,normal_face,NGLLX,NGLLY) 

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
  ! free surface face counter
  ifree = 0
  do ispec2D = 1, NSPEC2D_TOP
    ! sets element 
    ispec = ibelm_top(ispec2D)
     
    ! looks for i,j,k indices of GLL points on boundary face
    ! determines element face by given CUBIT corners
    do icorner=1,NGNOD2D
      xcoord(icorner) = nodes_coords_ext_mesh(1,nodes_ibelm_top(icorner,ispec2D))
      ycoord(icorner) = nodes_coords_ext_mesh(2,nodes_ibelm_top(icorner,ispec2D))
      zcoord(icorner) = nodes_coords_ext_mesh(3,nodes_ibelm_top(icorner,ispec2D))
    enddo
    
    ! sets face id of reference element associated with this face
    call get_element_face_id(ispec,xcoord,ycoord,zcoord,&
                              ibool,nspec,nglob, &
                              xstore_dummy,ystore_dummy,zstore_dummy, &
                              iface )
    iboun(iface,ispec) = .true. 
                              
    ! ijk indices of GLL points on face
    call get_element_face_gll_indices(iface,ijk_face,NGLLX,NGLLY)

    ! weighted jacobian and normal                          
    call get_jacobian_boundary_face(myrank,nspec, &
              xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob, &
              dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
              wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
              ispec,iface,jacobian2Dw_face,normal_face,NGLLX,NGLLY) 

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

    ! stores surface infos
    if( .not. ABSORB_FREE_SURFACE ) then
      ! store for free surface
      !jacobian2D_top(:,:,ispec2D) = jacobian2Dw_face(:,:)
      !normal_top(:,:,:,ispec2D) = normal_face(:,:,:)  

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
      
      ! resets free surface 
      ifree = 1
      free_surface_ispec(:) = 0
      free_surface_ijk(:,:,:) = 0
      free_surface_jacobian2Dw(:,:) = 0.0
      free_surface_normal(:,:,:) = 0.0
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

  call sum_all_i(num_abs_boundary_faces,iabs)
  if( myrank == 0 ) then
    write(IMAIN,*) '     absorbing boundary:'
    write(IMAIN,*) '     total number of faces = ',iabs
    if( ABSORB_FREE_SURFACE ) then
    write(IMAIN,*) '     absorbing boundary includes free surface'
    endif
  endif

  end subroutine get_absorbing_boundary

