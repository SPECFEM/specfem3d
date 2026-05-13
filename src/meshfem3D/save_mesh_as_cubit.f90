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


  subroutine save_mesh_files_as_cubit(nspec,nglob, &
                                      nodes_coords, ispec_material_id, &
                                      nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                                      ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top)

! subroutine to save meshes in case of a single MPI process

  use constants, only: NDIM,IMAIN,myrank,IIN_DB
  use constants_meshfem, only: NGLLX_M,NGLLY_M,NGLLZ_M
  use shared_parameters, only: NGNOD

  use meshfem_par, only: ibool, &
    NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
    NPROC_XI,NPROC_ETA, &
    NMATERIALS,material_properties, &
    nspec_CPML,CPML_to_spec,CPML_regions

  !! setting up wavefield discontinuity interface
  use shared_parameters, only: IS_WAVEFIELD_DISCONTINUITY

  implicit none

  integer, parameter :: IIN_database = IIN_DB

  ! number of spectral elements in each block
  integer, intent(in):: nspec

  ! number of vertices in each block
  integer, intent(in) :: nglob

  ! MPI Cartesian topology
  ! E for East (= XI_MIN), W for West (= XI_MAX), S for South (= ETA_MIN), N for North (= ETA_MAX)
  !integer, parameter :: W=1,E=2,S=3,N=4,NW=5,NE=6,SE=7,SW=8

  ! arrays with the mesh
  double precision, intent(in) :: nodes_coords(nglob,NDIM)

  integer, intent(in) :: ispec_material_id(nspec)

  ! boundary parameters locator
  integer, intent(in) :: nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax

  integer, intent(in) :: ibelm_xmin(NSPEC2DMAX_XMIN_XMAX),ibelm_xmax(NSPEC2DMAX_XMIN_XMAX)
  integer, intent(in) :: ibelm_ymin(NSPEC2DMAX_YMIN_YMAX),ibelm_ymax(NSPEC2DMAX_YMIN_YMAX)
  integer, intent(in) :: ibelm_bottom(NSPEC2D_BOTTOM)
  integer, intent(in) :: ibelm_top(NSPEC2D_TOP)

  ! local parameters
  integer :: i,ispec,iglob,ier
  integer :: mat_id,domain_id
  double precision  :: z_bottom
  integer :: nspec_CPML_total
  integer :: nnodes_mesh,inode
  integer, dimension(:),allocatable :: iglob_to_nodeid
  logical, dimension(:), allocatable :: mask_iglob

  ! safety check
  ! only for single process at the moment
  if (.not. (NPROC_XI == 1 .and. NPROC_ETA == 1)) then
    print *,'Error: SAVE_MESH_AS_CUBIT output requires NPROC_XI == NPROC_ETA == 1'
    print *,'       using NPROC_XI = ',NPROC_XI,' and NPROC_ETA = ',NPROC_ETA
    print *,'Please update your Mesh_Par_file and re-run the mesher...'
    stop 'Invalid NPROC_XI and/or NPROC_ETA for SAVE_MESH_AS_CUBIT output'
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  saving mesh as Cubit mesh in directory: MESH/'
    if (NGNOD /= 8) then
      write(IMAIN,*) '    (Using these Cubit files will require using NGNOD = 8 in Par_file for decomposer)'
    endif
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  z_bottom = 0.d0 ! will shift coordinates in z-direction

  ! nodes needed by ibool for outputting anchor nodes
  allocate(mask_iglob(nglob),iglob_to_nodeid(nglob),stat=ier)
  if (ier /= 0) stop 'Error allocating mask_iglob'
  mask_iglob(:) = .false.
  iglob_to_nodeid(:) = 0

  ! only output corner points as mesh, no internal points, otherwise decomposer will complain about unused nodes
  do ispec = 1,nspec
    ! corner points
    mask_iglob(ibool(1,1,1,ispec)) = .true.
    mask_iglob(ibool(NGLLZ_M,1,1,ispec)) = .true.
    mask_iglob(ibool(NGLLZ_M,NGLLZ_M,1,ispec)) = .true.
    mask_iglob(ibool(1,NGLLZ_M,1,ispec)) = .true.
    mask_iglob(ibool(1,1,NGLLZ_M,ispec)) = .true.
    mask_iglob(ibool(NGLLZ_M,1,NGLLZ_M,ispec)) = .true.
    mask_iglob(ibool(NGLLZ_M,NGLLZ_M,NGLLZ_M,ispec)) = .true.
    mask_iglob(ibool(1,NGLLZ_M,NGLLZ_M,ispec)) = .true.
  enddo

  ! corner points only
  nnodes_mesh = count(mask_iglob(:))
  if (nnodes_mesh < 1 .or. nnodes_mesh > nglob) then
    print *,'Error: nnodes_mesh ',nnodes_mesh,' is invalid for nglob = ',nglob
    stop 'Error nnodes_mesh invalid'
  endif

  ! maps iglob to new unique node IDs
  inode = 0
  do iglob = 1,nglob
    if (mask_iglob(iglob)) then
      inode = inode + 1
      iglob_to_nodeid(iglob) = inode
    endif
  enddo
  if (inode /= nnodes_mesh) then
    print *,'Error: inode count ',inode,' for nnodes_mesh ',nnodes_mesh,' is invalid; nglob = ',nglob
    stop 'Error inode count invalid'
  endif

  ! safety check just to re-visit this routine in case NGNOD changes
  if (NGNOD /= 8 .and. NGNOD /= 27) then
    stop 'Error invalid NGNOD for routine save_mesh_files_as_cubit()'
  endif

  ! outputs mesh as files in MESH/ directory
  ! (used for CUBIT models stored in a specfem-readable way; will need to run xdecompose_mesh with these files to continue)

  open(IIN_database, file = 'MESH/nummaterial_velocity_file',status='unknown',action='write',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening file: ','MESH/nummaterial_velocity_file'
    print *,'Please check if directory MESH/ exists for saving mesh files as cubit...'
    stop 'Error opening mesh file'
  endif

  do i = 1,NMATERIALS
    domain_id = int(material_properties(i,7))
    mat_id =  int(material_properties(i,8))
    if (domain_id > 0) then
      ! format: #domain_id #material_id #rho #vp #vs #Qkappa #Qmu #anisotropy_flag
      write(IIN_database,'(2i6,5f15.5,i6)') domain_id,mat_id,material_properties(i,1:5),0
    else
      write(*,*) 'STOP: undefined mat not yet implemented'
      stop
    endif
  enddo
  close(IIN_database)

  open(IIN_database,file='MESH/materials_file')
  do ispec = 1, nspec
    ! format: #ispec #material_id
    write(IIN_database,*) ispec,ispec_material_id(ispec)
  enddo

  open(IIN_database,file='MESH/nodes_coords_file')
  write(IIN_database,*) nnodes_mesh
  do iglob = 1,nglob
    ! point coordinates
    if (mask_iglob(iglob)) then
      ! format: #iglob #x #y #z(corrected by z_bottom)
      write(IIN_database,'(i14,3x,3(f20.5,1x))') iglob_to_nodeid(iglob), &
                                                 nodes_coords(iglob,1), &
                                                 nodes_coords(iglob,2), &
                                                 nodes_coords(iglob,3)-z_bottom
    endif
  enddo
  close(IIN_database)

  open(IIN_database,file='MESH/mesh_file')
  write(IIN_database,*) nspec
  do ispec = 1,nspec
    ! corner points
    ! format: #ispec #iglob1 #iglob2 #iglob3 #iglob4 #iglob5 #iglob6 #iglob7 #iglob8
    write(IIN_database,'(9i15)')  ispec, &
                                  iglob_to_nodeid(ibool(1,1,1,ispec)), &
                                  iglob_to_nodeid(ibool(NGLLZ_M,1,1,ispec)), &
                                  iglob_to_nodeid(ibool(NGLLZ_M,NGLLZ_M,1,ispec)), &
                                  iglob_to_nodeid(ibool(1,NGLLZ_M,1,ispec)), &
                                  iglob_to_nodeid(ibool(1,1,NGLLZ_M,ispec)), &
                                  iglob_to_nodeid(ibool(NGLLZ_M,1,NGLLZ_M,ispec)), &
                                  iglob_to_nodeid(ibool(NGLLZ_M,NGLLZ_M,NGLLZ_M,ispec)), &
                                  iglob_to_nodeid(ibool(1,NGLLZ_M,NGLLZ_M,ispec))
  enddo
  close(IIN_database)

  open(IIN_database,file='MESH/absorbing_surface_file_xmin')
  write(IIN_database,*) nspec2D_xmin
  do i = 1,nspec2D_xmin
    ! format: #boundary_face_id #iglob1 #iglob2 #iglob3 #iglob4
    write(IIN_database,'(5(i10,1x))') ibelm_xmin(i), &
                                      iglob_to_nodeid(ibool(1,1,1,ibelm_xmin(i))), &
                                      iglob_to_nodeid(ibool(1,NGLLY_M,1,ibelm_xmin(i))), &
                                      iglob_to_nodeid(ibool(1,1,NGLLZ_M,ibelm_xmin(i))), &
                                      iglob_to_nodeid(ibool(1,NGLLY_M,NGLLZ_M,ibelm_xmin(i)))
  enddo
  close(IIN_database)

  open(IIN_database,file='MESH/absorbing_surface_file_xmax')
  write(IIN_database,*) nspec2D_xmax
  do i = 1,nspec2D_xmax
    ! format: #boundary_face_id #iglob1 #iglob2 #iglob3 #iglob4
    write(IIN_database,'(5(i10,1x))') ibelm_xmax(i), &
                                      iglob_to_nodeid(ibool(NGLLX_M,1,1,ibelm_xmax(i))), &
                                      iglob_to_nodeid(ibool(NGLLX_M,NGLLY_M,1,ibelm_xmax(i))), &
                                      iglob_to_nodeid(ibool(NGLLX_M,1,NGLLZ_M,ibelm_xmax(i))), &
                                      iglob_to_nodeid(ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ibelm_xmax(i)))
  enddo
  close(IIN_database)

  open(IIN_database,file='MESH/absorbing_surface_file_ymin')
  write(IIN_database,*) nspec2D_ymin
  do i = 1,nspec2D_ymin
    ! format: #boundary_face_id #iglob1 #iglob2 #iglob3 #iglob4
    write(IIN_database,'(5(i10,1x))') ibelm_ymin(i), &
                                      iglob_to_nodeid(ibool(1,1,1,ibelm_ymin(i))), &
                                      iglob_to_nodeid(ibool(NGLLX_M,1,1,ibelm_ymin(i))), &
                                      iglob_to_nodeid(ibool(1,1,NGLLZ_M,ibelm_ymin(i))), &
                                      iglob_to_nodeid(ibool(NGLLX_M,1,NGLLZ_M,ibelm_ymin(i)))
  enddo
  close(IIN_database)

  open(IIN_database,file='MESH/absorbing_surface_file_ymax')
  write(IIN_database,*) nspec2D_ymax
  do i = 1,nspec2D_ymax
    ! format: #boundary_face_id #iglob1 #iglob2 #iglob3 #iglob4
    write(IIN_database,'(5(i10,1x))') ibelm_ymax(i), &
                                      iglob_to_nodeid(ibool(NGLLX_M,NGLLY_M,1,ibelm_ymax(i))), &
                                      iglob_to_nodeid(ibool(1,NGLLY_M,1,ibelm_ymax(i))), &
                                      iglob_to_nodeid(ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ibelm_ymax(i))), &
                                      iglob_to_nodeid(ibool(1,NGLLY_M,NGLLZ_M,ibelm_ymax(i)))
  enddo


  open(IIN_database,file='MESH/absorbing_surface_file_bottom')
  write(IIN_database,*) NSPEC2D_BOTTOM
  do i = 1,NSPEC2D_BOTTOM
    ! format: #boundary_face_id #iglob1 #iglob2 #iglob3 #iglob4
    write(IIN_database,'(5(i10,1x))') ibelm_bottom(i), &
                                      iglob_to_nodeid(ibool(1,1,1,ibelm_bottom(i))), &
                                      iglob_to_nodeid(ibool(NGLLX_M,1,1,ibelm_bottom(i))), &
                                      iglob_to_nodeid(ibool(NGLLX_M,NGLLY_M,1,ibelm_bottom(i))), &
                                      iglob_to_nodeid(ibool(1,NGLLY_M,1,ibelm_bottom(i)))
  enddo
  close(IIN_database)

  open(IIN_database,file='MESH/free_or_absorbing_surface_file_zmax')
  write(IIN_database,*) NSPEC2D_TOP
  do i = 1,NSPEC2D_TOP
    ! format: #boundary_face_id #iglob1 #iglob2 #iglob3 #iglob4
    write(IIN_database,'(5(i10,1x))') ibelm_top(i), &
                                      iglob_to_nodeid(ibool(1,1,NGLLZ_M,ibelm_top(i))), &
                                      iglob_to_nodeid(ibool(NGLLX_M,1,NGLLZ_M,ibelm_top(i))), &
                                      iglob_to_nodeid(ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ibelm_top(i))), &
                                      iglob_to_nodeid(ibool(1,NGLLY_M,NGLLZ_M,ibelm_top(i)))
  enddo
  close(IIN_database)

  ! CPML
  call sum_all_i(nspec_CPML,nspec_CPML_total)
  call synchronize_all()
  call bcast_all_singlei(nspec_CPML_total)
  call synchronize_all()

  if (nspec_CPML_total > 0) then
    open(IIN_database,file='MESH/absorbing_cpml_file')
    write(IIN_database,*) nspec_CPML
    do i = 1,nspec_CPML
      ! format: #cpml_ispec #cpml_region_id
      write(IIN_database,*) CPML_to_spec(i), CPML_regions(i)
    enddo
    close(IIN_database)
  endif

  ! setting up wavefield discontinuity interface
  if (IS_WAVEFIELD_DISCONTINUITY) then
    call write_wavefield_discontinuity_file()
  endif

  deallocate(mask_iglob,iglob_to_nodeid)

  end subroutine save_mesh_files_as_cubit



!---------------------------------------------------------------------------------------------------------------

  subroutine save_mesh_files_for_coupled_model(nspec, &
                                               nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                                               ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                                               xgrid,ygrid,zgrid)

! subroutine to save meshes in case of a single MPI process

  use constants, only: NGLLX, NGLLY, NGLLZ, NDIM, ZERO, IMAIN, myrank, &
    INJECTION_TECHNIQUE_IS_AXISEM
  use constants_meshfem, only: NGLLX_M,NGLLY_M,NGLLZ_M

  use shared_parameters, only: NGNOD,COUPLE_WITH_INJECTION_TECHNIQUE,INJECTION_TECHNIQUE_TYPE

  use meshfem_par, only: NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
    NPROC_XI,NPROC_ETA

  implicit none

  ! number of spectral elements in each block
  integer :: nspec

  ! MPI Cartesian topology
  ! E for East (= XI_MIN), W for West (= XI_MAX), S for South (= ETA_MIN), N for North (= ETA_MAX)
  !integer, parameter :: W=1,E=2,S=3,N=4,NW=5,NE=6,SE=7,SW=8

  !! VM VM add all GLL points for Axisem coupling
  double precision, dimension(NGLLX_M,NGLLY_M,NGLLZ_M,nspec) :: xgrid, ygrid, zgrid

  ! boundary parameters locator
  integer :: nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax

  integer :: ibelm_xmin(NSPEC2DMAX_XMIN_XMAX),ibelm_xmax(NSPEC2DMAX_XMIN_XMAX)
  integer :: ibelm_ymin(NSPEC2DMAX_YMIN_YMAX),ibelm_ymax(NSPEC2DMAX_YMIN_YMAX)
  integer :: ibelm_bottom(NSPEC2D_BOTTOM)
  integer :: ibelm_top(NSPEC2D_TOP)

  integer :: i,ispec

  double precision  :: z_bottom

  ! for axisem coupling case  ( only serial case for mesher use scotch after)
  integer, parameter :: nlayer = 12 !! (number of layer in the model iasp91, or ak135, or prem (one more layer than the model)

  ! GLL points
  double precision, parameter :: GAUSSALPHA = 0.d0, GAUSSBETA = 0.d0

  double precision   :: rotation_matrix(3,3)
  double precision   :: zlayer(nlayer), vpv(nlayer,4), vsv(nlayer,4), density(nlayer,4)
  integer            :: ilayer, updown(NGLLZ)

  !! GLL points
  double precision, dimension(NGLLX,NGLLY,NGLLZ) ::  longitud, latitud, radius
  double precision, dimension(NGLLX,NGLLY,NGLLZ) ::  xstore, ystore, zstore
   !! Element control points
  double precision, dimension(NGNOD) :: xelm, yelm, zelm
  double precision, dimension(NGLLZ) :: radius_Z ! for temporary copies

  !! 3D shape functions and their derivatives
  double precision :: shape3D(NGNOD,NGLLX,NGLLY,NGLLZ)
  double precision :: dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ)
  !! GLL points and weights of integration
  double precision :: xigll(NGLLX), yigll(NGLLY), zigll(NGLLZ)
  double precision :: wxgll(NGLLX), wygll(NGLLY), wzgll(NGLLZ)

  double precision  :: ANGULAR_WIDTH_ETA_RAD, ANGULAR_WIDTH_XI_RAD
  double precision  :: lat_center_chunk, lon_center_chunk, chunk_depth, chunk_azi
  double precision  :: radius_of_box_top

  integer :: ielm, j,k, imin,imax,jmin,jmax,kmin,kmax,icounter_pts,icounter_faces
  integer :: nel_lat, nel_lon, nel_depth
  logical :: buried_box

  character(len=10)  :: line
  character(len=250) :: model1D_file

  double precision,parameter  :: deg2rad = 3.141592653589793d0/180.d0

  ! safety check
  if (.not. COUPLE_WITH_INJECTION_TECHNIQUE) return

  ! adds files only in case of AxiSEM coupling
  if (INJECTION_TECHNIQUE_TYPE /= INJECTION_TECHNIQUE_IS_AXISEM) return

  ! only for single process at the moment
  if (.not. (NPROC_XI == 1 .and. NPROC_ETA == 1)) then
    print *,'Error: SAVE_MESH_AS_CUBIT output requires NPROC_XI == NPROC_ETA == 1'
    print *,'       using NPROC_XI = ',NPROC_XI,' and NPROC_ETA = ',NPROC_ETA
    print *,'Please update your Mesh_Par_file and re-run the mesher...'
    stop 'Invalid NPROC_XI and/or NPROC_ETA for SAVE_MESH_AS_CUBIT output'
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  AxiSEM coupling:'
    write(IMAIN,*) '    saving mesh files for coupled model in directory: MESH/'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

1000 format(3f30.10)

  z_bottom = 0.d0 ! will shift coordinates in z-direction

  !
  !--- set up coordinates of the Gauss-Lobatto-Legendre points
  !
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

  !
  !--- if number of points is odd, the middle abscissa is exactly zero
  !
  if (mod(NGLLX,2) /= 0) xigll((NGLLX - 1)/2 + 1) = ZERO
  if (mod(NGLLY,2) /= 0) yigll((NGLLY - 1)/2 + 1) = ZERO
  if (mod(NGLLZ,2) /= 0) zigll((NGLLZ - 1)/2 + 1) = ZERO

  !
  !--- get the 3-D shape functions
  !
  call get_shape3D(shape3D,dershape3D,xigll,yigll,zigll,NGNOD,NGLLX,NGLLY,NGLLZ)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    reading parameter for coupling file: MESH/ParFileMeshChunk'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  !! reading parameters for coupling
  open(90, file='MESH/ParFileMeshChunk',action='read')
  read(90,'(a)') line
  read(90,*) ANGULAR_WIDTH_XI_RAD, ANGULAR_WIDTH_ETA_RAD
  read(90,'(a)') line
  read(90,*) lon_center_chunk, lat_center_chunk, chunk_azi
  read(90,'(a)') line
  read(90,*) chunk_depth
  read(90,'(a)') line
  read(90,*) nel_lon,nel_lat, nel_depth
  read(90,'(a)') line
  read(90,'(a)') model1D_file
  read(90,'(a)') line
  read(90,*) buried_box
  if (buried_box) then
    read(90,'(a)') line
    read(90,*) radius_of_box_top
    radius_of_box_top =  radius_of_box_top * 1000.
  else
    radius_of_box_top = 6371000.
  endif
  model1D_file = 'MESH/'//trim(model1D_file)
  close(90)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    chunk : angular width xi/eta                = ',ANGULAR_WIDTH_XI_RAD,'/',ANGULAR_WIDTH_ETA_RAD
    write(IMAIN,*) '            center lon/lat/azi                  = ',lon_center_chunk,'/',lat_center_chunk,'/',chunk_azi
    write(IMAIN,*) '            depth                               = ',chunk_depth
    write(IMAIN,*) '            number of elements nlon/nlat/ndepth = ',nel_lon,'/',nel_lat,'/',nel_depth
    write(IMAIN,*) '            model 1D file                       = ',trim(model1D_file)
    write(IMAIN,*) '    box is buried : ',buried_box
    if (buried_box) then
      write(IMAIN,*) '                     radius of box top = ',radius_of_box_top
    endif
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! read 1D AxiSEM model
  call Read_dsm_model(model1D_file,vpv,vsv,density,zlayer,nlayer)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    creating input 1D model: ','MESH/model_1D.in'
    write(IMAIN,*) '                             number of layers = ',nlayer
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! modele 1D
  open(88,file='MESH/model_1D.in')
  write(88,*) nlayer,4
  do i = 1,nlayer
    write(88,*) zlayer(i)
    write(88,'(4f20.10)') vpv(i,:)
    write(88,'(4f20.10)') vsv(i,:)
    write(88,'(4f20.10)') density(i,:)
  enddo
  z_bottom = minval(zgrid(:,:,:,:))
  write(88,*)  radius_of_box_top + z_bottom                       !6371000.+z_bottom
  write(88,*)  lon_center_chunk,  lat_center_chunk,  chunk_azi
  close(88)

  ! compute rotation matrix
  call compute_rotation_matrix(rotation_matrix,lon_center_chunk,lat_center_chunk, chunk_azi)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    creating boundary files: ','MESH/list_ggl_boundary_spherical.txt'
    write(IMAIN,*) '                             ','MESH/list_ggl_boundary_Cartesian.txt'
    write(IMAIN,*) '                             ','MESH/flags_boundary.txt'
    write(IMAIN,*) '                             ','MESH/Nb_ielm_faces.txt'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! counters
  icounter_pts = 0
  icounter_faces = 0

  open(91, file = 'MESH/list_ggl_boundary_spherical.txt')
  open(92, file = 'MESH/list_ggl_boundary_Cartesian.txt')
  open(89, file = 'MESH/flags_boundary.txt')

  open(90, file = 'MESH/Nb_ielm_faces.txt')
  write(90,*)  nspec2D_xmin
  write(90,*)  nspec2D_xmax
  write(90,*)  nspec2D_ymin
  write(90,*)  nspec2D_ymax
  write(90,*)  nspec2D_bottom
  write(90,*)  nspec2D_top
  close(90)

  ! xmin
  do ielm = 1,nspec2D_xmin

    ispec = ibelm_xmin(ielm)

    write(89,*) ispec,ielm,1

    ! face counter
    icounter_faces = icounter_faces + 1

    xelm(1)=xgrid(1,1,1,ispec)
    xelm(2)=xgrid(NGLLX_M,1,1,ispec)
    xelm(3)=xgrid(NGLLX_M,NGLLY_M,1,ispec)
    xelm(4)=xgrid(1,NGLLY_M,1,ispec)
    xelm(5)=xgrid(1,1,NGLLZ_M,ispec)
    xelm(6)=xgrid(NGLLX_M,1,NGLLZ_M,ispec)
    xelm(7)=xgrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
    xelm(8)=xgrid(1,NGLLY_M,NGLLZ_M,ispec)

    yelm(1)=ygrid(1,1,1,ispec)
    yelm(2)=ygrid(NGLLX_M,1,1,ispec)
    yelm(3)=ygrid(NGLLX_M,NGLLY_M,1,ispec)
    yelm(4)=ygrid(1,NGLLY_M,1,ispec)
    yelm(5)=ygrid(1,1,NGLLZ_M,ispec)
    yelm(6)=ygrid(NGLLX_M,1,NGLLZ_M,ispec)
    yelm(7)=ygrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
    yelm(8)=ygrid(1,NGLLY_M,NGLLZ_M,ispec)

    zelm(1)=zgrid(1,1,1,ispec)
    zelm(2)=zgrid(NGLLX_M,1,1,ispec)
    zelm(3)=zgrid(NGLLX_M,NGLLY_M,1,ispec)
    zelm(4)=zgrid(1,NGLLY_M,1,ispec)
    zelm(5)=zgrid(1,1,NGLLZ_M,ispec)
    zelm(6)=zgrid(NGLLX_M,1,NGLLY_M,ispec)
    zelm(7)=zgrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
    zelm(8)=zgrid(1,NGLLY_M,NGLLZ_M,ispec)

    call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
    zstore(:,:,:) = zstore(:,:,:) + radius_of_box_top             ! 6371000.

    call Cartesian2spheric(xstore,ystore,zstore,rotation_matrix,longitud,latitud,radius,deg2rad)
    zstore(:,:,:) = zstore(:,:,:) - radius_of_box_top             ! 6371000.

    radius_Z(:) = radius(3,3,:) ! to avoid warning about temporary copies in routine call
    call find_layer_in_axisem_model(ilayer,updown,radius_Z,zlayer,nlayer)

    imin = 1
    imax = 1
    jmin = 1
    jmax = NGLLY
    kmin = 1
    kmax = NGLLZ

    do k = kmin,kmax
      do j = jmin,jmax
        do i = imin,imax
          write(92,'(3f25.10,i10,6i3)') xstore(i,j,k),ystore(i,j,k),zstore(i,j,k),ispec,i,j,k,1, &
                                        ilayer,updown(k)
          write(91,1000) radius(i,j,k), latitud(i,j,k), longitud(i,j,k)
          ! point counter
          icounter_pts = icounter_pts + 1
        enddo
      enddo
    enddo
  enddo

  ! xmax
  do ielm = 1,nspec2D_xmax

    ispec = ibelm_xmax(ielm)

    write(89,*) ispec,ielm,2

    ! face counter
    icounter_faces = icounter_faces + 1

    xelm(1)=xgrid(1,1,1,ispec)
    xelm(2)=xgrid(NGLLX_M,1,1,ispec)
    xelm(3)=xgrid(NGLLX_M,NGLLY_M,1,ispec)
    xelm(4)=xgrid(1,NGLLY_M,1,ispec)
    xelm(5)=xgrid(1,1,NGLLZ_M,ispec)
    xelm(6)=xgrid(NGLLX_M,1,NGLLZ_M,ispec)
    xelm(7)=xgrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
    xelm(8)=xgrid(1,NGLLY_M,NGLLZ_M,ispec)

    yelm(1)=ygrid(1,1,1,ispec)
    yelm(2)=ygrid(NGLLX_M,1,1,ispec)
    yelm(3)=ygrid(NGLLX_M,NGLLY_M,1,ispec)
    yelm(4)=ygrid(1,NGLLY_M,1,ispec)
    yelm(5)=ygrid(1,1,NGLLZ_M,ispec)
    yelm(6)=ygrid(NGLLX_M,1,NGLLZ_M,ispec)
    yelm(7)=ygrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
    yelm(8)=ygrid(1,NGLLY_M,NGLLZ_M,ispec)

    zelm(1)=zgrid(1,1,1,ispec)
    zelm(2)=zgrid(NGLLX_M,1,1,ispec)
    zelm(3)=zgrid(NGLLX_M,NGLLY_M,1,ispec)
    zelm(4)=zgrid(1,NGLLY_M,1,ispec)
    zelm(5)=zgrid(1,1,NGLLZ_M,ispec)
    zelm(6)=zgrid(NGLLX_M,1,NGLLZ_M,ispec)
    zelm(7)=zgrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
    zelm(8)=zgrid(1,NGLLY_M,NGLLZ_M,ispec)

    call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
    zstore(:,:,:) = zstore(:,:,:) + radius_of_box_top             ! 6371000.

    call Cartesian2spheric(xstore,ystore,zstore,rotation_matrix,longitud,latitud,radius,deg2rad)
    zstore(:,:,:) = zstore(:,:,:) - radius_of_box_top             ! 6371000.

    radius_Z(:) = radius(3,3,:) ! to avoid warning about temporary copies in routine call
    call find_layer_in_axisem_model(ilayer,updown,radius_Z,zlayer,nlayer)

    imin = NGLLX
    imax = NGLLX
    jmin = 1
    jmax = NGLLY
    kmin = 1
    kmax = NGLLZ

    do k = kmin,kmax
      do j = jmin,jmax
        do i = imin,imax
          write(92,'(3f25.10,i10,6i3)') xstore(i,j,k),ystore(i,j,k),zstore(i,j,k),ispec,i,j,k,2, &
                                        ilayer,updown(k)
          write(91,1000) radius(i,j,k), latitud(i,j,k), longitud(i,j,k)
          ! point counter
          icounter_pts = icounter_pts + 1
        enddo
      enddo
    enddo
  enddo

  ! ymin
  do ielm = 1,nspec2D_ymin

    ispec = ibelm_ymin(ielm)

    write(89,*) ispec,ielm,3

    ! face counter
    icounter_faces = icounter_faces + 1

    xelm(1)=xgrid(1,1,1,ispec)
    xelm(2)=xgrid(NGLLX_M,1,1,ispec)
    xelm(3)=xgrid(NGLLX_M,NGLLY_M,1,ispec)
    xelm(4)=xgrid(1,NGLLY_M,1,ispec)
    xelm(5)=xgrid(1,1,NGLLZ_M,ispec)
    xelm(6)=xgrid(NGLLX_M,1,NGLLZ_M,ispec)
    xelm(7)=xgrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
    xelm(8)=xgrid(1,NGLLY_M,NGLLZ_M,ispec)

    yelm(1)=ygrid(1,1,1,ispec)
    yelm(2)=ygrid(NGLLX_M,1,1,ispec)
    yelm(3)=ygrid(NGLLX_M,NGLLY_M,1,ispec)
    yelm(4)=ygrid(1,NGLLY_M,1,ispec)
    yelm(5)=ygrid(1,1,NGLLZ_M,ispec)
    yelm(6)=ygrid(NGLLX_M,1,NGLLZ_M,ispec)
    yelm(7)=ygrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
    yelm(8)=ygrid(1,NGLLY_M,NGLLZ_M,ispec)

    zelm(1)=zgrid(1,1,1,ispec)
    zelm(2)=zgrid(NGLLX_M,1,1,ispec)
    zelm(3)=zgrid(NGLLX_M,NGLLY_M,1,ispec)
    zelm(4)=zgrid(1,NGLLY_M,1,ispec)
    zelm(5)=zgrid(1,1,NGLLZ_M,ispec)
    zelm(6)=zgrid(NGLLX_M,1,NGLLZ_M,ispec)
    zelm(7)=zgrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
    zelm(8)=zgrid(1,NGLLY_M,NGLLZ_M,ispec)

    call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
    zstore(:,:,:) = zstore(:,:,:) + radius_of_box_top           ! 6371000.

    call Cartesian2spheric(xstore,ystore,zstore,rotation_matrix,longitud,latitud,radius,deg2rad)
    zstore(:,:,:) = zstore(:,:,:) - radius_of_box_top           ! 6371000.

    radius_Z(:) = radius(3,3,:) ! to avoid warning about temporary copies in routine call
    call find_layer_in_axisem_model(ilayer,updown,radius_Z,zlayer,nlayer)

    imin = 1
    imax = NGLLX
    jmin = 1
    jmax = 1
    kmin = 1
    kmax = NGLLZ

    do k = kmin,kmax
      do j = jmin,jmax
        do i = imin,imax
          write(92,'(3f25.10,i10,6i3)') xstore(i,j,k),ystore(i,j,k),zstore(i,j,k),ispec,i,j,k,3, &
                                        ilayer,updown(k)
          write(91,1000) radius(i,j,k), latitud(i,j,k), longitud(i,j,k)
          ! point counter
          icounter_pts = icounter_pts + 1
        enddo
      enddo
    enddo
  enddo

  ! ymax
  do ielm = 1,nspec2D_ymax

    ispec = ibelm_ymax(ielm)

    write(89,*) ispec,ielm,4

    ! face counter
    icounter_faces = icounter_faces + 1

    xelm(1)=xgrid(1,1,1,ispec)
    xelm(2)=xgrid(NGLLX_M,1,1,ispec)
    xelm(3)=xgrid(NGLLX_M,NGLLY_M,1,ispec)
    xelm(4)=xgrid(1,NGLLY_M,1,ispec)
    xelm(5)=xgrid(1,1,NGLLZ_M,ispec)
    xelm(6)=xgrid(NGLLX_M,1,NGLLZ_M,ispec)
    xelm(7)=xgrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
    xelm(8)=xgrid(1,NGLLY_M,NGLLZ_M,ispec)

    yelm(1)=ygrid(1,1,1,ispec)
    yelm(2)=ygrid(NGLLX_M,1,1,ispec)
    yelm(3)=ygrid(NGLLX_M,NGLLY_M,1,ispec)
    yelm(4)=ygrid(1,NGLLY_M,1,ispec)
    yelm(5)=ygrid(1,1,NGLLZ_M,ispec)
    yelm(6)=ygrid(NGLLX_M,1,NGLLZ_M,ispec)
    yelm(7)=ygrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
    yelm(8)=ygrid(1,NGLLY_M,NGLLZ_M,ispec)

    zelm(1)=zgrid(1,1,1,ispec)
    zelm(2)=zgrid(NGLLX_M,1,1,ispec)
    zelm(3)=zgrid(NGLLX_M,NGLLY_M,1,ispec)
    zelm(4)=zgrid(1,NGLLY_M,1,ispec)
    zelm(5)=zgrid(1,1,NGLLZ_M,ispec)
    zelm(6)=zgrid(NGLLX_M,1,NGLLZ_M,ispec)
    zelm(7)=zgrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
    zelm(8)=zgrid(1,NGLLY_M,NGLLZ_M,ispec)

    call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
    zstore(:,:,:) = zstore(:,:,:) + radius_of_box_top           ! 6371000.

    call Cartesian2spheric(xstore,ystore,zstore,rotation_matrix,longitud,latitud,radius,deg2rad)
    zstore(:,:,:) = zstore(:,:,:) - radius_of_box_top           ! 6371000.

    radius_Z(:) = radius(3,3,:) ! to avoid warning about temporary copies in routine call
    call find_layer_in_axisem_model(ilayer,updown,radius_Z,zlayer,nlayer)

    imin = 1
    imax = NGLLX
    jmin = NGLLY
    jmax = NGLLY
    kmin = 1
    kmax = NGLLZ

    do k = kmin,kmax
      do j = jmin,jmax
        do i = imin,imax
          write(92,'(3f25.10,i10,6i3)') xstore(i,j,k),ystore(i,j,k),zstore(i,j,k),ispec,i,j,k,4, &
                                        ilayer,updown(k)
          write(91,1000) radius(i,j,k), latitud(i,j,k), longitud(i,j,k)
          ! point counter
          icounter_pts = icounter_pts + 1
        enddo
      enddo
    enddo
  enddo

  ! bottom
  do ielm = 1,nspec2D_BOTTOM

    ispec = ibelm_bottom(ielm)

    write(89,*) ispec,ielm,5

    ! face counter
    icounter_faces = icounter_faces + 1

    xelm(1)=xgrid(1,1,1,ispec)
    xelm(2)=xgrid(NGLLX_M,1,1,ispec)
    xelm(3)=xgrid(NGLLX_M,NGLLY_M,1,ispec)
    xelm(4)=xgrid(1,NGLLY_M,1,ispec)
    xelm(5)=xgrid(1,1,NGLLZ_M,ispec)
    xelm(6)=xgrid(NGLLX_M,1,NGLLZ_M,ispec)
    xelm(7)=xgrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
    xelm(8)=xgrid(1,NGLLY_M,NGLLZ_M,ispec)

    yelm(1)=ygrid(1,1,1,ispec)
    yelm(2)=ygrid(NGLLX_M,1,1,ispec)
    yelm(3)=ygrid(NGLLX_M,NGLLY_M,1,ispec)
    yelm(4)=ygrid(1,NGLLY_M,1,ispec)
    yelm(5)=ygrid(1,1,NGLLZ_M,ispec)
    yelm(6)=ygrid(NGLLX_M,1,NGLLZ_M,ispec)
    yelm(7)=ygrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
    yelm(8)=ygrid(1,NGLLY_M,NGLLZ_M,ispec)

    zelm(1)=zgrid(1,1,1,ispec)
    zelm(2)=zgrid(NGLLX_M,1,1,ispec)
    zelm(3)=zgrid(NGLLX_M,NGLLY_M,1,ispec)
    zelm(4)=zgrid(1,NGLLY_M,1,ispec)
    zelm(5)=zgrid(1,1,NGLLZ_M,ispec)
    zelm(6)=zgrid(NGLLX_M,1,NGLLZ_M,ispec)
    zelm(7)=zgrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
    zelm(8)=zgrid(1,NGLLY_M,NGLLZ_M,ispec)

    call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
    zstore(:,:,:) = zstore(:,:,:) + radius_of_box_top           ! 6371000.

    call Cartesian2spheric(xstore,ystore,zstore,rotation_matrix,longitud,latitud,radius,deg2rad)
    zstore(:,:,:) = zstore(:,:,:) - radius_of_box_top           ! 6371000.

    radius_Z(:) = radius(3,3,:) ! to avoid warning about temporary copies in routine call
    call find_layer_in_axisem_model(ilayer,updown,radius_Z,zlayer,nlayer)

    imin = 1
    imax = NGLLX
    jmin = 1
    jmax = NGLLY
    kmin = 1
    kmax = 1

    do k = kmin,kmax
      do j = jmin,jmax
        do i = imin,imax
          write(92,'(3f25.10,i10,6i3)') xstore(i,j,k),ystore(i,j,k),zstore(i,j,k),ispec,i,j,k,5, &
                                        ilayer,updown(k)
          write(91,1000) radius(i,j,k), latitud(i,j,k), longitud(i,j,k)
          ! point counter
          icounter_pts = icounter_pts + 1
        enddo
      enddo
    enddo
  enddo

  if (buried_box) then
    ! top
    do ielm = 1,nspec2D_TOP

      ispec = ibelm_top(ielm)

      write(89,*) ispec,ielm,6

      ! face counter
      icounter_faces = icounter_faces + 1

      xelm(1)=xgrid(1,1,1,ispec)
      xelm(2)=xgrid(NGLLX_M,1,1,ispec)
      xelm(3)=xgrid(NGLLX_M,NGLLY_M,1,ispec)
      xelm(4)=xgrid(1,NGLLY_M,1,ispec)
      xelm(5)=xgrid(1,1,NGLLZ_M,ispec)
      xelm(6)=xgrid(NGLLX_M,1,NGLLZ_M,ispec)
      xelm(7)=xgrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
      xelm(8)=xgrid(1,NGLLY_M,NGLLZ_M,ispec)

      yelm(1)=ygrid(1,1,1,ispec)
      yelm(2)=ygrid(NGLLX_M,1,1,ispec)
      yelm(3)=ygrid(NGLLX_M,NGLLY_M,1,ispec)
      yelm(4)=ygrid(1,NGLLY_M,1,ispec)
      yelm(5)=ygrid(1,1,NGLLZ_M,ispec)
      yelm(6)=ygrid(NGLLX_M,1,NGLLZ_M,ispec)
      yelm(7)=ygrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
      yelm(8)=ygrid(1,NGLLY_M,NGLLZ_M,ispec)

      zelm(1)=zgrid(1,1,1,ispec)
      zelm(2)=zgrid(NGLLX_M,1,1,ispec)
      zelm(3)=zgrid(NGLLX_M,NGLLY_M,1,ispec)
      zelm(4)=zgrid(1,NGLLY_M,1,ispec)
      zelm(5)=zgrid(1,1,NGLLZ_M,ispec)
      zelm(6)=zgrid(NGLLX_M,1,NGLLZ_M,ispec)
      zelm(7)=zgrid(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
      zelm(8)=zgrid(1,NGLLY_M,NGLLZ_M,ispec)

      call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
      zstore(:,:,:) = zstore(:,:,:) + radius_of_box_top            ! 6371000.

      call Cartesian2spheric(xstore,ystore,zstore,rotation_matrix,longitud,latitud,radius,deg2rad)
      zstore(:,:,:) = zstore(:,:,:) - radius_of_box_top            ! 6371000.

      radius_Z(:) = radius(3,3,:) ! to avoid warning about temporary copies in routine call
      call find_layer_in_axisem_model(ilayer,updown,radius_Z,zlayer,nlayer)

      imin = 1
      imax = NGLLX
      jmin = 1
      jmax = NGLLY
      kmin = NGLLZ
      kmax = NGLLZ

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax
            write(92,'(3f25.10,i10,6i3)') xstore(i,j,k),ystore(i,j,k),zstore(i,j,k),ispec,i,j,k,6, &
                                          ilayer,updown(k)
            write(91,1000) radius(i,j,k), latitud(i,j,k), longitud(i,j,k)
            ! point counter
            icounter_pts = icounter_pts + 1
          enddo
        enddo
      enddo
    enddo
  endif

  close(89)
  close(91)
  close(92)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    total number of boundary faces = ',icounter_faces
    write(IMAIN,*) '    total number of boundary points = ',icounter_pts
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine save_mesh_files_for_coupled_model
