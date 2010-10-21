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

  subroutine get_coupling_surfaces(myrank, &
                        nspec,nglob,ibool,NPROC, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        num_interfaces_ext_mesh,max_interface_size_ext_mesh, &
                        my_neighbours_ext_mesh)
                            
! determines coupling surface for acoustic-elastic domains

  use create_regions_mesh_ext_par 
  implicit none

! number of spectral elements in each block
  integer :: myrank,nspec,nglob,NPROC

! arrays with the mesh
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

! MPI communication
  integer :: num_interfaces_ext_mesh,max_interface_size_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: my_neighbours_ext_mesh
  integer, dimension(NGLLX*NGLLX*max_interface_size_ext_mesh,num_interfaces_ext_mesh) :: &
            ibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh

! local parameters
  ! (assumes NGLLX=NGLLY=NGLLZ)
  real(kind=CUSTOM_REAL),dimension(NGNOD2D) :: xcoord,ycoord,zcoord    
  real(kind=CUSTOM_REAL) :: jacobian2Dw_face(NGLLX,NGLLY)
  real(kind=CUSTOM_REAL) :: normal_face(NDIM,NGLLX,NGLLY)
  real(kind=CUSTOM_REAL),dimension(:,:,:),allocatable :: tmp_normal
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: tmp_jacobian2Dw
  integer :: ijk_face(3,NGLLX,NGLLY)
  integer,dimension(:,:,:),allocatable :: tmp_ijk
  integer,dimension(:),allocatable :: tmp_ispec

  integer,dimension(NGNOD2D) :: iglob_corners_ref !,iglob_corners
  integer :: ispec,i,j,k,igll,ier,iglob
  integer :: inum,iface_ref,icorner,iglob_midpoint ! iface,ispec_neighbor
  integer :: count_elastic,count_acoustic
  
  ! mpi interface communication
  integer, dimension(:), allocatable :: elastic_flag,acoustic_flag,test_flag
  integer, dimension(:,:), allocatable :: ibool_interfaces_ext_mesh_dummy
  integer :: max_nibool_interfaces_ext_mesh
  logical, dimension(:), allocatable :: mask_ibool
  
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

  
  ! test vtk output
  !integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: gll_data
  !character(len=256):: prname_file
  
! allocates temporary arrays  
  allocate(tmp_normal(NDIM,NGLLSQUARE,nspec*6))
  allocate(tmp_jacobian2Dw(NGLLSQUARE,nspec*6))  
  allocate(tmp_ijk(3,NGLLSQUARE,nspec*6))
  allocate(tmp_ispec(nspec*6))
  tmp_ispec(:) = 0
  tmp_ijk(:,:,:) = 0
  tmp_normal(:,:,:) = 0.0
  tmp_jacobian2Dw(:,:) = 0.0
  
  ! sets flags for acoustic / elastic on global points
  allocate(elastic_flag(nglob),stat=ier)
  allocate(acoustic_flag(nglob),stat=ier)  
  allocate(test_flag(nglob),stat=ier)  
  allocate(mask_ibool(nglob),stat=ier)
  if( ier /= 0 ) stop 'error allocate flag array'  
  elastic_flag(:) = 0
  acoustic_flag(:) = 0
  test_flag(:) = 0
  count_elastic = 0
  count_acoustic = 0
  do ispec = 1, nspec
    ! counts elements
    if( ispec_is_elastic(ispec) ) count_elastic = count_elastic + 1
    if( ispec_is_acoustic(ispec) ) count_acoustic = count_acoustic + 1
    
    ! sets flags on global points
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          ! global index
          iglob = ibool(i,j,k,ispec)         
          ! sets elastic flag
          if( ispec_is_elastic(ispec) ) elastic_flag(iglob) =  myrank+1
          ! sets acoustic flag
          if( ispec_is_acoustic(ispec) ) acoustic_flag(iglob) =  myrank+1
          ! sets test flag
          test_flag(iglob) = myrank+1
        enddo
      enddo
    enddo
  enddo
  call sum_all_i(count_acoustic,inum)
  if( myrank == 0 ) then
    write(IMAIN,*) '     total acoustic elements:',inum
  endif   
  call sum_all_i(count_elastic,inum)
  if( myrank == 0 ) then
    write(IMAIN,*) '     total elastic elements :',inum
  endif   

  ! collects contributions from different MPI partitions
  ! sets up MPI communications
  max_nibool_interfaces_ext_mesh = maxval( nibool_interfaces_ext_mesh(:) )
  allocate(ibool_interfaces_ext_mesh_dummy(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
  if( ier /= 0 ) stop 'error allocating array'  
  do i = 1, num_interfaces_ext_mesh
     ibool_interfaces_ext_mesh_dummy(:,i) = ibool_interfaces_ext_mesh(1:max_nibool_interfaces_ext_mesh,i)
  enddo  
  ! sums elastic flags
  call assemble_MPI_scalar_i_ext_mesh(NPROC,nglob,elastic_flag, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh_dummy,&
                        my_neighbours_ext_mesh)
  ! sums acoustic flags
  call assemble_MPI_scalar_i_ext_mesh(NPROC,nglob,acoustic_flag, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh_dummy,&
                        my_neighbours_ext_mesh)

  ! sums test flags
  call assemble_MPI_scalar_i_ext_mesh(NPROC,nglob,test_flag, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh_dummy,&
                        my_neighbours_ext_mesh)

  ! loops over all element faces and 
  ! counts number of coupling faces between acoustic and elastic elements
  mask_ibool(:) = .false.
  inum = 0    
  do ispec=1,nspec

    ! loops over each face
    do iface_ref= 1, 6      

      ! takes indices of corners of reference face
      do icorner = 1,NGNOD2D
        i = iface_all_corner_ijk(1,icorner,iface_ref)
        j = iface_all_corner_ijk(2,icorner,iface_ref)
        k = iface_all_corner_ijk(3,icorner,iface_ref)
        ! global reference indices
        iglob_corners_ref(icorner) = ibool(i,j,k,ispec)

        ! reference corner coordinates
        xcoord(icorner) = xstore_dummy(iglob_corners_ref(icorner))
        ycoord(icorner) = ystore_dummy(iglob_corners_ref(icorner))
        zcoord(icorner) = zstore_dummy(iglob_corners_ref(icorner))                  
      enddo
      
      ! checks if face has acoustic side
      if( acoustic_flag( iglob_corners_ref(1) ) >= 1 .and. &
         acoustic_flag( iglob_corners_ref(2) ) >= 1 .and. &
         acoustic_flag( iglob_corners_ref(3) ) >= 1 .and. &
         acoustic_flag( iglob_corners_ref(4) ) >= 1) then        
        ! checks if face is has an elastic side 
        if( elastic_flag( iglob_corners_ref(1) ) >= 1 .and. &
           elastic_flag( iglob_corners_ref(2) ) >= 1 .and. &
           elastic_flag( iglob_corners_ref(3) ) >= 1 .and. &
           elastic_flag( iglob_corners_ref(4) ) >= 1) then

          ! reference midpoint on face (used to avoid redundant face counting)
          i = iface_all_midpointijk(1,iface_ref)
          j = iface_all_midpointijk(2,iface_ref)
          k = iface_all_midpointijk(3,iface_ref)      
          iglob_midpoint = ibool(i,j,k,ispec)

          ! checks if points on this face are masked already
          if( .not. mask_ibool(iglob_midpoint) ) then

            ! gets face GLL points i,j,k indices from element face
            call get_element_face_gll_indices(iface_ref,ijk_face,NGLLX,NGLLY)
            
            ! takes each element face only once, if it lies on an MPI interface
            ! note: this is not exactly load balanced
            !          lowest rank process collects as many faces as possible, second lowest as so forth
            if( (test_flag(iglob_midpoint) == myrank+1) .or. &
               (test_flag(iglob_midpoint) > 2*(myrank+1)) ) then
            
              ! gets face GLL 2Djacobian, weighted from element face
              call get_jacobian_boundary_face(myrank,nspec, &
                        xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob, &
                        dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                        ispec,iface_ref,jacobian2Dw_face,normal_face,NGLLX,NGLLY)

              ! normal convention: points away from acoustic, reference element
              !                                switch normal direction if necessary
              do j=1,NGLLY
                do i=1,NGLLX
                    ! directs normals such that they point outwards of element
                    call get_element_face_normal(ispec,iface_ref,xcoord,ycoord,zcoord, &
                                                ibool,nspec,nglob, &
                                                xstore_dummy,ystore_dummy,zstore_dummy, &
                                                normal_face(:,i,j) )
                    ! makes sure that it always points away from acoustic element, 
                    ! otherwise switch direction
                    if( ispec_is_elastic(ispec) ) normal_face(:,i,j) = - normal_face(:,i,j)
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
            endif ! test_flag
          endif ! mask_ibool          
        endif ! elastic_flag
      endif ! acoustic_flag
    enddo ! iface_ref
  enddo ! ispec
    
! stores completed coupling face informations  
! 
! note: no need to store material parameters on these coupling points 
!          for acoustic-elastic interface
  num_coupling_ac_el_faces = inum
  allocate(coupling_ac_el_normal(NDIM,NGLLSQUARE,num_coupling_ac_el_faces))
  allocate(coupling_ac_el_jacobian2Dw(NGLLSQUARE,num_coupling_ac_el_faces))
  allocate(coupling_ac_el_ijk(3,NGLLSQUARE,num_coupling_ac_el_faces))
  allocate(coupling_ac_el_ispec(num_coupling_ac_el_faces))
  do inum = 1,num_coupling_ac_el_faces
    coupling_ac_el_normal(:,:,inum) = tmp_normal(:,:,inum)
    coupling_ac_el_jacobian2Dw(:,inum) = tmp_jacobian2Dw(:,inum)
    coupling_ac_el_ijk(:,:,inum) = tmp_ijk(:,:,inum)
    coupling_ac_el_ispec(inum) = tmp_ispec(inum)    
  enddo

! user output
  call sum_all_i(num_coupling_ac_el_faces,inum)
  if( myrank == 0 ) then
    write(IMAIN,*) '     acoustic-elastic coupling:'
    write(IMAIN,*) '     total number of faces = ',inum
  endif  

  end subroutine get_coupling_surfaces

