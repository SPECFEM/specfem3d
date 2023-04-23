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


! creation of arrays for movie and shakemap routines for external meshes

  subroutine setup_movie_meshes()

  use specfem_par
  use specfem_par_movie

  implicit none

  ! face corner indices
  integer,dimension(NGNOD2D_FOUR_CORNERS),parameter :: iorderi = (/ 1,NGLLX,NGLLX,1 /)
  integer,dimension(NGNOD2D_FOUR_CORNERS),parameter :: iorderj = (/ 1,1,NGLLY,NGLLY /)
  integer :: i,j,k,iglob,ispec,ispec2D,ia,ier
  integer :: iface,ipoin,imin,imax,jmin,jmax,kmin,kmax
  character(len=MAX_STRING_LEN) :: filename
  integer :: nfaces_m,npoin,npoin_elem
  real(kind=CUSTOM_REAL),dimension(1):: dummy

  !! DK DK create a copy of the iface number to avoid an (erroneous) warning when compiling on a Cray with crayftn
  integer :: iface_copy

  ! number of points per element surface
  if (USE_HIGHRES_FOR_MOVIES) then
    npoin_elem = NGLLX*NGLLY
  else
    npoin_elem = NGNOD2D_FOUR_CORNERS
  endif

  ! number of faces on surface
  if (nfaces_surface == 0) then
    ! dummy arrays
    nfaces_m = 1
  else
    nfaces_m = nfaces_surface
  endif

  ! total number of surface points
  npoin = npoin_elem * nfaces_m

  ! surface elements
  allocate(faces_surface_ispec(nfaces_m),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2098')
  allocate(faces_surface_ibool(npoin_elem,nfaces_m),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2099')
  if (ier /= 0) stop 'error allocating array faces_surface_ispec'
  faces_surface_ispec(:) = 0
  faces_surface_ibool(:,:) = 0

  ! point locations
  allocate(store_val_x(npoin),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2100')
  allocate(store_val_y(npoin),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2101')
  allocate(store_val_z(npoin),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2102')
  if (ier /= 0) stop 'error allocating location arrays for highres movie'

  ! surface movie data
  if (MOVIE_SURFACE) then
    allocate(store_val_ux(npoin),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2103')
    allocate(store_val_uy(npoin),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2104')
    allocate(store_val_uz(npoin),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2105')
    if (ier /= 0) stop 'error allocating arrays for highres movie'
  endif

  ! shakemap data
  if (CREATE_SHAKEMAP) then
    allocate(shakemap_ux(npoin),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2106')
    allocate(shakemap_uy(npoin),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2107')
    allocate(shakemap_uz(npoin),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2108')
    if (ier /= 0) stop 'error allocating arrays for highres shakemap'
    ! initializes shakemap values
    shakemap_ux(:) = 0.0_CUSTOM_REAL
    shakemap_uy(:) = 0.0_CUSTOM_REAL
    shakemap_uz(:) = 0.0_CUSTOM_REAL
  endif

  ! number of surface faces for all partitions together
  call sum_all_i(nfaces_surface,nfaces_surface_glob_ext_mesh)

  ! arrays used for collected/gathered fields
  if (myrank == 0) then
    ! all point locations
    allocate(store_val_x_all(npoin_elem*nfaces_surface_glob_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2109')
    allocate(store_val_y_all(npoin_elem*nfaces_surface_glob_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2110')
    allocate(store_val_z_all(npoin_elem*nfaces_surface_glob_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2111')
    if (ier /= 0) stop 'error allocating arrays for highres movie'

    ! surface movie
    if (MOVIE_SURFACE) then
      allocate(store_val_ux_all(npoin_elem*nfaces_surface_glob_ext_mesh),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2112')
      allocate(store_val_uy_all(npoin_elem*nfaces_surface_glob_ext_mesh),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2113')
      allocate(store_val_uz_all(npoin_elem*nfaces_surface_glob_ext_mesh),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2114')
      if (ier /= 0) stop 'error allocating arrays for highres movie'
    endif

    ! shakemap
    if (CREATE_SHAKEMAP) then
      allocate(shakemap_ux_all(npoin_elem*nfaces_surface_glob_ext_mesh),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2115')
      allocate(shakemap_uy_all(npoin_elem*nfaces_surface_glob_ext_mesh),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2116')
      allocate(shakemap_uz_all(npoin_elem*nfaces_surface_glob_ext_mesh),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2117')
      if (ier /= 0) stop 'error allocating arrays for highres movie'
      shakemap_ux_all(:) = 0.0_CUSTOM_REAL
      shakemap_uy_all(:) = 0.0_CUSTOM_REAL
      shakemap_uz_all(:) = 0.0_CUSTOM_REAL
    endif
  endif

  ! arrays for collecting movies and shakemaps
  allocate(nfaces_perproc_surface(NPROC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2118')
  allocate(faces_surface_offset(NPROC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2119')
  if (ier /= 0) stop 'error allocating array for movie faces'

  ! number of faces per slice
  call gather_all_singlei(nfaces_surface,nfaces_perproc_surface,NPROC)

  ! array offsets
  faces_surface_offset(1) = 0
  do i = 2, NPROC
    faces_surface_offset(i) = sum(nfaces_perproc_surface(1:i-1))
  enddo
  faces_surface_offset(:) = faces_surface_offset(:)*npoin_elem

  ! determines surface elements and indexes for movie/shakemap
  ! stores global indices of GLL points on the surface to array faces_surface_ibool
  select case (MOVIE_TYPE)
  case (1,3)
    ! only for top, free surface
    ! check
    if (num_free_surface_faces /= nfaces_surface) then
      print *,'Error: rank ',myrank,'has invalid num_free_surface_faces and nfaces_surface values for surface movie'
      print *,'  num_free_surface_faces = ',num_free_surface_faces,'nfaces_surface = ',nfaces_surface
      stop 'Error num_free_surface_faces and nfaces_surface for surface MOVIE_TYPE == 1'
    endif

    ! determines movie indexing for all faces on top, free surface
    do iface = 1,num_free_surface_faces
      ispec = free_surface_ispec(iface)

      ! copies element addressing for movie surface
      faces_surface_ispec(iface) = ispec

      ! high_resolution
      if (USE_HIGHRES_FOR_MOVIES) then
        do ipoin = 1, npoin_elem
          i = free_surface_ijk(1,ipoin,iface)
          j = free_surface_ijk(2,ipoin,iface)
          k = free_surface_ijk(3,ipoin,iface)
          iglob = ibool(i,j,k,ispec)
          ! sets indices
          faces_surface_ibool(ipoin,iface) = iglob
        enddo
      else
        imin = minval(free_surface_ijk(1,:,iface))
        imax = maxval(free_surface_ijk(1,:,iface))
        jmin = minval(free_surface_ijk(2,:,iface))
        jmax = maxval(free_surface_ijk(2,:,iface))
        kmin = minval(free_surface_ijk(3,:,iface))
        kmax = maxval(free_surface_ijk(3,:,iface))
        do ipoin = 1,npoin_elem
          ! corner points
          if (imin == imax) then
            iglob = ibool(imin,iorderi(ipoin),iorderj(ipoin),ispec)
          else if (jmin == jmax) then
            iglob = ibool(iorderi(ipoin),jmin,iorderj(ipoin),ispec)
          else
            iglob = ibool(iorderi(ipoin),iorderj(ipoin),kmin,ispec)
          endif
          ! sets indices
          faces_surface_ibool(ipoin,iface) = iglob
        enddo
      endif
    enddo

  case (2)
    ! all external, outer surfaces
    ! face counter
    iface = 0
    ! for all mesh elements touching an external, outer surface
    do ispec = 1, NSPEC_AB
      if (ispec_is_surface_external_mesh(ispec)) then
        ! determines indexing for all faces on a outer surface
        !! DK DK create a copy of the ispec number to avoid an (erroneous) warning when compiling on a Cray with crayftn
        iface_copy = iface
        call setup_movie_face_indices(ispec,iface_copy)
        iface = iface_copy !! DK DK this also only to avoid the (erroneous) warning when compiling on a Cray with crayftn
      endif
    enddo

    ! checks number of faces
    if (iface /= nfaces_surface) then
      print *,'Error: number of movie faces found ',iface,'should be',nfaces_surface
      call exit_mpi(myrank,'Error number of faces')
    endif

  case default
    ! not recognized type
    call exit_MPI(myrank,'Invalid MOVIE_TYPE value for surface movie and/or shakemap, must be 1, 2 or 3')

  end select

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'mesh surfaces:'
    write(IMAIN,*) '  high-resolution  = ',USE_HIGHRES_FOR_MOVIES
    write(IMAIN,*) '  nfaces per slice = ',nfaces_perproc_surface(:)
    write(IMAIN,*)
    write(IMAIN,*) '  total number of surface faces is ',nfaces_surface_glob_ext_mesh
    write(IMAIN,*)
    if (PLOT_CROSS_SECTIONS) then
      write(IMAIN,*) '  used on cross-sections'
      write(IMAIN,*)
    endif
    ! shakemaps
    if (CREATE_SHAKEMAP) then
      write(IMAIN,*) 'shakemap:'
      select case (MOVIE_TYPE)
      case (1)
        write(IMAIN,*) '  movie type 1: horizontal peak-ground values'
      case (2)
        write(IMAIN,*) '  movie type 2: maximum length of particle vector'
      case (3)
        write(IMAIN,*) '  movie type 3: geometric mean of horizontal peak-ground values'
      case default
        print *,'Error: invalid MOVIE_TYPE value, must be 1, 2 or 3 for CREATE_SHAKEMAP'
        stop 'Error invalid MOVIE_TYPE for CREATE_SHAKEMAP'
      end select
      write(IMAIN,*)
    endif
    ! surface movie
    if (MOVIE_SURFACE) then
      write(IMAIN,*) 'surface movies:'
      if (SAVE_DISPLACEMENT) then
        write(IMAIN,*) '  saving: particle displacements'
      else
        write(IMAIN,*) '  saving: particle velocities'
      endif
      write(IMAIN,*) '  number of steps between frames = ',NTSTEP_BETWEEN_FRAMES
      write(IMAIN,*)
    endif
    call flush_IMAIN()

    ! updates number of surface elements in an include file for the movies
    if (nfaces_surface_glob_ext_mesh > 0) then
      filename = OUTPUT_FILES(1:len_trim(OUTPUT_FILES)) // '/surface_from_mesher.h'
      open(unit=IOUT,file=trim(filename),status='unknown')
      write(IOUT,*) '!'
      write(IOUT,*) '! this is the parameter file for static compilation for movie creation'
      write(IOUT,*) '!'
      write(IOUT,*) '! number of elements containing surface faces '
      write(IOUT,*) '! ---------------'
      write(IOUT,*)
      write(IOUT,*) 'integer,parameter :: NSPEC_SURFACE_EXT_MESH = ',nfaces_surface_glob_ext_mesh
      write(IOUT,*)
      close(IOUT)
    endif
  endif

  ! for gathering movie data
  ! all GLL points on surface
  nfaces_perproc_surface(:) = nfaces_perproc_surface(:) * npoin_elem
  nfaces_surface_points = nfaces_surface * npoin_elem
  nfaces_surface_glob_points = nfaces_surface_glob_ext_mesh * npoin_elem

  ! sets surface point locations
  do ispec2D = 1,nfaces_surface
    do ipoin = 1,npoin_elem
      iglob = faces_surface_ibool(ipoin,ispec2D)
      ia = npoin_elem * (ispec2D - 1) + ipoin
      ! x,y,z point coordinates
      store_val_x(ia) = xstore(iglob)
      store_val_y(ia) = ystore(iglob)
      store_val_z(ia) = zstore(iglob)
    enddo
  enddo
  ! main process collects all info
  ! collects locations only once
  ! main collects all
  if (myrank == 0) then
    call gatherv_all_cr(store_val_x,nfaces_surface_points, &
       store_val_x_all,nfaces_perproc_surface,faces_surface_offset, &
       nfaces_surface_glob_points,NPROC)
    call gatherv_all_cr(store_val_y,nfaces_surface_points, &
       store_val_y_all,nfaces_perproc_surface,faces_surface_offset, &
       nfaces_surface_glob_points,NPROC)
    call gatherv_all_cr(store_val_z,nfaces_surface_points, &
       store_val_z_all,nfaces_perproc_surface,faces_surface_offset, &
       nfaces_surface_glob_points,NPROC)
  else
    ! secondarys just send
    call gatherv_all_cr(store_val_x,nfaces_surface_points, &
       dummy,nfaces_perproc_surface,faces_surface_offset, &
       1,NPROC)
    call gatherv_all_cr(store_val_y,nfaces_surface_points, &
       dummy,nfaces_perproc_surface,faces_surface_offset, &
       1,NPROC)
    call gatherv_all_cr(store_val_z,nfaces_surface_points, &
       dummy,nfaces_perproc_surface,faces_surface_offset, &
       1,NPROC)
  endif

  end subroutine setup_movie_meshes

!================================================================

  subroutine setup_movie_face_indices(ispec,iface)

! sets up arrays faces_surface_ibool,faces_surface_ispec

  use constants, only: NGLLX,NGLLY,NGLLZ
  use specfem_par, only: USE_HIGHRES_FOR_MOVIES,ibool,iglob_is_surface_external_mesh

  use specfem_par_movie, only: faces_surface_ibool,faces_surface_ispec

  implicit none

  integer, intent(in) :: ispec
  integer, intent(inout) :: iface

  ! local parameters
  integer :: iglob,i,j,k,ipoin

  ! checks if interior point on element surface is on an external mesh surface

  ! zmin face
  iglob = ibool(2,2,1,ispec)
  if (iglob_is_surface_external_mesh(iglob)) then
    iface = iface + 1
    faces_surface_ispec(iface) = ispec
    if (USE_HIGHRES_FOR_MOVIES) then
      ipoin = 0
      do j = NGLLY, 1, -1
        do i = 1, NGLLX
          ipoin = ipoin+1
          faces_surface_ibool(ipoin,iface) = ibool(i,j,1,ispec)
        enddo
      enddo
    else
      ! only corners
      faces_surface_ibool(1,iface) = ibool(1,1,1,ispec)
      faces_surface_ibool(2,iface) = ibool(1,NGLLY,1,ispec)
      faces_surface_ibool(3,iface) = ibool(NGLLX,NGLLY,1,ispec)
      faces_surface_ibool(4,iface) = ibool(NGLLX,1,1,ispec)
    endif
  endif

  ! zmax face
  iglob = ibool(2,2,NGLLZ,ispec)
  if (iglob_is_surface_external_mesh(iglob)) then
    iface = iface + 1
    faces_surface_ispec(iface) = ispec
    if (USE_HIGHRES_FOR_MOVIES) then
      ipoin = 0
      do j = 1, NGLLY
        do i = 1, NGLLX
          ipoin = ipoin+1
          faces_surface_ibool(ipoin,iface) = ibool(i,j,NGLLZ,ispec)
        enddo
      enddo
    else
      ! only corners
      faces_surface_ibool(1,iface) = ibool(1,1,NGLLZ,ispec)
      faces_surface_ibool(2,iface) = ibool(NGLLX,1,NGLLZ,ispec)
      faces_surface_ibool(3,iface) = ibool(NGLLX,NGLLY,NGLLZ,ispec)
      faces_surface_ibool(4,iface) = ibool(1,NGLLY,NGLLZ,ispec)
    endif
  endif

  ! ymin face
  iglob = ibool(2,1,2,ispec)
  if (iglob_is_surface_external_mesh(iglob)) then
    iface = iface + 1
    faces_surface_ispec(iface) = ispec
    if (USE_HIGHRES_FOR_MOVIES) then
      ipoin = 0
      do k = 1, NGLLZ
        do i = 1, NGLLX
          ipoin = ipoin+1
          faces_surface_ibool(ipoin,iface) = ibool(i,1,k,ispec)
        enddo
      enddo
    else
      ! only corners
      faces_surface_ibool(1,iface) = ibool(1,1,1,ispec)
      faces_surface_ibool(2,iface) = ibool(NGLLX,1,1,ispec)
      faces_surface_ibool(3,iface) = ibool(NGLLX,1,NGLLZ,ispec)
      faces_surface_ibool(4,iface) = ibool(1,1,NGLLZ,ispec)
    endif
  endif

  ! ymax face
  iglob = ibool(2,NGLLY,2,ispec)
  if (iglob_is_surface_external_mesh(iglob)) then
    iface = iface + 1
    faces_surface_ispec(iface) = ispec
    if (USE_HIGHRES_FOR_MOVIES) then
      ipoin = 0
      do k = 1, NGLLZ
        do i = NGLLX, 1, -1
          ipoin = ipoin+1
          faces_surface_ibool(ipoin,iface) = ibool(i,NGLLY,k,ispec)
        enddo
      enddo
    else
      ! only corners
      faces_surface_ibool(1,iface) = ibool(NGLLX,NGLLY,1,ispec)
      faces_surface_ibool(2,iface) = ibool(1,NGLLY,1,ispec)
      faces_surface_ibool(3,iface) = ibool(1,NGLLY,NGLLZ,ispec)
      faces_surface_ibool(4,iface) = ibool(NGLLX,NGLLY,NGLLZ,ispec)
    endif
  endif

  ! xmin face
  iglob = ibool(1,2,2,ispec)
  if (iglob_is_surface_external_mesh(iglob)) then
    iface = iface + 1
    faces_surface_ispec(iface) = ispec
    if (USE_HIGHRES_FOR_MOVIES) then
      ipoin = 0
      do k = 1, NGLLZ
        do j = NGLLY, 1, -1
          ipoin = ipoin+1
          faces_surface_ibool(ipoin,iface) = ibool(1,j,k,ispec)
        enddo
      enddo
    else
      ! only corners
      faces_surface_ibool(1,iface) = ibool(1,NGLLY,1,ispec)
      faces_surface_ibool(2,iface) = ibool(1,1,1,ispec)
      faces_surface_ibool(3,iface) = ibool(1,1,NGLLZ,ispec)
      faces_surface_ibool(4,iface) = ibool(1,NGLLY,NGLLZ,ispec)
    endif
  endif

  ! xmax face
  iglob = ibool(NGLLX,2,2,ispec)
  if (iglob_is_surface_external_mesh(iglob)) then
    iface = iface + 1
    faces_surface_ispec(iface) = ispec
    if (USE_HIGHRES_FOR_MOVIES) then
      ipoin = 0
      do k = 1, NGLLZ
        do j = 1, NGLLY
          ipoin = ipoin+1
          faces_surface_ibool(ipoin,iface) = ibool(NGLLX,j,k,ispec)
        enddo
      enddo
    else
      ! only corners
      faces_surface_ibool(1,iface) = ibool(NGLLX,1,1,ispec)
      faces_surface_ibool(2,iface) = ibool(NGLLX,NGLLY,1,ispec)
      faces_surface_ibool(3,iface) = ibool(NGLLX,NGLLY,NGLLZ,ispec)
      faces_surface_ibool(4,iface) = ibool(NGLLX,1,NGLLZ,ispec)
    endif
  endif

  end subroutine setup_movie_face_indices
