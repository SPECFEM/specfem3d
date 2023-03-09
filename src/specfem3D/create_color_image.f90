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

  module image_PNM_par

  use constants, only: CUSTOM_REAL

  ! ----------------------------------------------
  ! USER PARAMETER
  ! non linear display to enhance small amplitudes in color images
  real(kind=CUSTOM_REAL), parameter :: POWER_DISPLAY_COLOR = 0.30_CUSTOM_REAL

  ! amplitude threshold
  real(kind=CUSTOM_REAL),parameter :: image_cutsnaps  = 1.e-2

  ! create temporary image files in binary PNM P6 format (smaller)
  ! or ASCII PNM P3 format (easier to edit)
  logical, parameter :: BINARY_FILE = .true.

  !! defaults
  ! image data output:
  !   type = 1 : displ/velocity x-component
  !   type = 2 : displ/velocity y-component
  !   type = 3 : displ/velocity z-component
  !   type = 4 : displ/velocity norm
  !   type = 5 : only model vp
  !   type = 6 : only model vs
  !   type = 7 : only model rho
  integer :: PNM_IMAGE_TYPE = 4

  ! cross-section surface
  ! cross-section origin point
  real(kind=CUSTOM_REAL) :: PNM_section_xorg = 0.0
  real(kind=CUSTOM_REAL) :: PNM_section_yorg = 0.0
  real(kind=CUSTOM_REAL) :: PNM_section_zorg = 0.0

  ! cross-section surface normal
  real(kind=CUSTOM_REAL) :: PNM_section_nx = 0.0
  real(kind=CUSTOM_REAL) :: PNM_section_ny = 0.0
  real(kind=CUSTOM_REAL) :: PNM_section_nz = 1.0 ! for horizontal cross-section

  ! cross-section (in-plane) horizontal-direction
  real(kind=CUSTOM_REAL) :: PNM_section_hdirx = 1.0
  real(kind=CUSTOM_REAL) :: PNM_section_hdiry = 0.0
  real(kind=CUSTOM_REAL) :: PNM_section_hdirz = 0.0

  ! cross-section (in-plane) vertical-direction
  real(kind=CUSTOM_REAL) :: PNM_section_vdirx = 0.0
  real(kind=CUSTOM_REAL) :: PNM_section_vdiry = 1.0
  real(kind=CUSTOM_REAL) :: PNM_section_vdirz = 0.0

  ! background for wavefields
  integer :: BACKGROUND_TYPE = 1   ! use 1 == vp, 2 == vs, 3 == rho as gray background

  ! color palette
  integer :: PNM_COLOR_PALETTE = 1     ! use 1 == blue-red, 2 == rainbow(blue-cyan-green-yellow-orange-red)
                                       !     3==tomo, 4==terrain

  ! ----------------------------------------------

  ! image data
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable:: image_color_background_display
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable:: image_color_data

  integer,dimension(:,:),allocatable :: iglob_image_color
  integer,dimension(:,:),allocatable :: ispec_image_color

  ! pixel data
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: data_pixel_recv
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: data_pixel_send
  integer,dimension(:),allocatable :: num_pixel_loc
  integer,dimension(:),allocatable :: nb_pixel_per_proc
  integer,dimension(:,:),allocatable :: num_pixel_recv
  integer :: NX_IMAGE_color,NZ_IMAGE_color
  integer :: nb_pixel_loc

  end module image_PNM_par

!=============================================================

  subroutine write_PNM_initialize()

  use image_PNM_par

  use specfem_par, only: NGLOB_AB,NSPEC_AB,NPROC,ibool,xstore,ystore,zstore, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,my_neighbors_ext_mesh, &
                        ibool_interfaces_ext_mesh,myrank

  use constants, only: HUGEVAL,NGLLX,NGLLY,NGLLZ,IMAIN

  implicit none
  ! local parameters
  ! image sizes
  real(kind=CUSTOM_REAL):: xmin_color_image_loc,xmax_color_image_loc
  real(kind=CUSTOM_REAL):: xmin_color_image,xmax_color_image
  real(kind=CUSTOM_REAL):: zmin_color_image_loc,zmax_color_image_loc
  real(kind=CUSTOM_REAL):: zmin_color_image,zmax_color_image
  ! image pixels
  real(kind=CUSTOM_REAL):: size_pixel_horizontal,size_pixel_vertical
  real(kind=CUSTOM_REAL):: dist_pixel,dist_min_pixel
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: dist_pixel_image,dist_pixel_recv
  real(kind=CUSTOM_REAL):: pixel_midpoint_x,pixel_midpoint_z,x_loc,z_loc,xtmp,ztmp
  real(kind=CUSTOM_REAL):: ratio
  real(kind=CUSTOM_REAL):: distance_x1,distance_x2,distance_z1,distance_z2
  integer :: npgeo,npgeo_glob
  integer :: i,j,k,iproc,iglob,ispec,ier
  integer :: nx0,nz0
  ! data from mesh
  real(kind=CUSTOM_REAL),dimension(:),allocatable:: xcoord,zcoord
  integer,dimension(:),allocatable :: iglob_coord,ispec_coord
  logical,dimension(:),allocatable:: ispec_is_image_surface,iglob_is_image_surface
  integer :: num_iglob_image_surface
  integer :: countval,locval(1),irank
  integer :: zoom_factor = 2
  logical :: zoom
  integer, dimension(1,0:NPROC-1) :: tmp_pixel_per_proc

  ! checks image type
  if (PNM_IMAGE_TYPE > 7 .or. PNM_IMAGE_TYPE < 1) then
    call exit_mpi(myrank,'That type is not implemented for PNM images yet')
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '********'
    write(IMAIN,*) 'PNM image:'
    ! type = 1 : velocity V_x component
    if (PNM_IMAGE_TYPE == 1) write(IMAIN,*) '  type: velocity V_x component'
    ! type = 2 : velocity V_y component
    if (PNM_IMAGE_TYPE == 2) write(IMAIN,*) '  type: velocity V_y component'
    ! type = 3 : velocity V_z component
    if (PNM_IMAGE_TYPE == 3) write(IMAIN,*) '  type: velocity V_z component'
    ! type = 4 : velocity V norm
    if (PNM_IMAGE_TYPE == 4) write(IMAIN,*) '  type: velocity norm'
    ! only vp
    if (PNM_IMAGE_TYPE == 5) write(IMAIN,*) '  type: vp model'
    ! only vs
    if (PNM_IMAGE_TYPE == 6) write(IMAIN,*) '  type: vs model'
    ! only rho
    if (PNM_IMAGE_TYPE == 7) write(IMAIN,*) '  type: rho model'
  endif

  ! finds global points on image surface
  allocate(ispec_is_image_surface(NSPEC_AB),iglob_is_image_surface(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1766')
  if (ier /= 0) call exit_mpi(myrank,'error allocating image ispec and iglob')

  call detect_surface_PNM_image(NPROC,NGLOB_AB,NSPEC_AB,ibool, &
                                ispec_is_image_surface, &
                                iglob_is_image_surface, &
                                num_iglob_image_surface, &
                                num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                nibool_interfaces_ext_mesh,my_neighbors_ext_mesh, &
                                ibool_interfaces_ext_mesh, &
                                xstore,ystore,zstore,myrank)

  ! extracts points on surface
  allocate(xcoord(num_iglob_image_surface),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1767')
  allocate(zcoord(num_iglob_image_surface),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1768')
  allocate(iglob_coord(num_iglob_image_surface),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1769')
  allocate(ispec_coord(num_iglob_image_surface),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1770')
  if (ier /= 0) call exit_mpi(myrank,'error allocating xyz image coordinates')
  xcoord(:) = 0.0_CUSTOM_REAL; zcoord(:) = 0.0_CUSTOM_REAL
  iglob_coord(:) = 0; ispec_coord(:) = 0

  countval = 0
  do ispec=1,NSPEC_AB
    if (ispec_is_image_surface(ispec)) then
      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX
            iglob = ibool(i,j,k,ispec)
            if (iglob_is_image_surface(iglob)) then
              countval = countval + 1
              ! coordinates with respect to horizontal and vertical direction
              xcoord(countval)= xstore(iglob)*PNM_section_hdirx &
                              + ystore(iglob)*PNM_section_hdiry &
                              + zstore(iglob)*PNM_section_hdirz
              zcoord(countval)= xstore(iglob)*PNM_section_vdirx &
                              + ystore(iglob)*PNM_section_vdiry &
                              + zstore(iglob)*PNM_section_vdirz
              iglob_coord(countval) = iglob
              ispec_coord(countval) = ispec

              ! reset iglob flag
              iglob_is_image_surface(iglob) = .false.
            endif
          enddo
        enddo
      enddo
    endif
  enddo

  if (countval /= num_iglob_image_surface) call exit_mpi(myrank,'error image point number')

  ! horizontal size of the image
  xmin_color_image_loc = minval( xcoord(:) )
  xmax_color_image_loc = maxval( xcoord(:) )

  ! vertical size
  zmin_color_image_loc = minval( zcoord(:) )
  zmax_color_image_loc = maxval( zcoord(:) )

  ! global values
  xmin_color_image = xmin_color_image_loc
  xmax_color_image = xmax_color_image_loc
  zmin_color_image = zmin_color_image_loc
  zmax_color_image = zmax_color_image_loc

  ! global number of points on image slice
  npgeo = num_iglob_image_surface
  npgeo_glob = npgeo

  !MPI for all processes
  ! takes minimum of all process and stores it in xmin_color_image
  call min_all_all_cr(xmin_color_image_loc,xmin_color_image)
  call min_all_all_cr(zmin_color_image_loc,zmin_color_image)
  call max_all_all_cr(xmax_color_image_loc,xmax_color_image)
  call max_all_all_cr(zmax_color_image_loc,zmax_color_image)
  call sum_all_all_i(npgeo,npgeo_glob)

  ! compute number of pixels in the horizontal direction and pixels in the vertical
  ! direction based on ratio of sizes
  ratio = abs(xmax_color_image - xmin_color_image)/abs(zmax_color_image - zmin_color_image)
  NX_IMAGE_color = nint( sqrt( ratio * dble(npgeo_glob) ) )
  NZ_IMAGE_color = nint( dble(npgeo_glob) / NX_IMAGE_color )

  ! convert pixel sizes to even numbers because easier to reduce size,
  ! create MPEG movies in postprocessing
  NX_IMAGE_color = 2 * (NX_IMAGE_color / 2)
  NZ_IMAGE_color = 2 * (NZ_IMAGE_color / 2)

  ! check that image size is not too big
  if (NX_IMAGE_color > 4096 .or. NZ_IMAGE_color > 4096) then
    ! half of it
    NX_IMAGE_color = NX_IMAGE_color / 2
    NZ_IMAGE_color = NZ_IMAGE_color / 2
    ! even numbers
    NX_IMAGE_color = 2 * (NX_IMAGE_color / 2)
    NZ_IMAGE_color = 2 * (NZ_IMAGE_color / 2)
  endif

  ! ...and not too small
  zoom = .false.
  if (NX_IMAGE_color < 600 .or. NZ_IMAGE_color < 600) then
    ! increases it
    nx0 = NX_IMAGE_color
    nz0 = NZ_IMAGE_color
    do while (NX_IMAGE_color < 600 .or. NZ_IMAGE_color < 600)
      zoom_factor = zoom_factor + 2
      NX_IMAGE_color = nx0 * zoom_factor
      NZ_IMAGE_color = nz0 * zoom_factor
      zoom = .true.
    enddo
  endif

  ! create all the pixels
  if (NX_IMAGE_color /= 0) then
    size_pixel_horizontal = (xmax_color_image - xmin_color_image) / dble(NX_IMAGE_color)
  else
    size_pixel_horizontal = 0.0
  endif

  if (NZ_IMAGE_color /= 0) then
    size_pixel_vertical = (zmax_color_image - zmin_color_image) / dble(NZ_IMAGE_color)
  else
    size_pixel_vertical = 0.0
  endif

  if (myrank == 0) then
    write(IMAIN,*) '  image points: ',npgeo_glob
    write(IMAIN,*) '  image xmin/xmax: ',xmin_color_image,'/',xmax_color_image
    write(IMAIN,*) '  image zmin/zmax: ',zmin_color_image,'/',zmax_color_image
    write(IMAIN,*) '  pixel numbers: ',NX_IMAGE_color,' x ',NZ_IMAGE_color
    write(IMAIN,*) '  pixel sizes  : ',size_pixel_horizontal,' x ',size_pixel_vertical
  endif

  ! allocate an array for the grid point that corresponds to a given image data point
  allocate(iglob_image_color(NX_IMAGE_color,NZ_IMAGE_color),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1771')
  allocate(ispec_image_color(NX_IMAGE_color,NZ_IMAGE_color),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1772')
  if (ier /= 0) call exit_mpi(myrank,'error allocating iglob_image_color')
  allocate(dist_pixel_image(NX_IMAGE_color,NZ_IMAGE_color),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1773')
  if (ier /= 0) call exit_mpi(myrank,'error allocating dist pixel image')
  iglob_image_color(:,:) = -1
  ispec_image_color(:,:) = 0
  dist_pixel_image(:,:) = HUGEVAL

  if (zoom) then
    distance_x1 = zoom_factor*size_pixel_horizontal
    distance_x2 = (zoom_factor+1)*size_pixel_horizontal
    distance_z1 = zoom_factor*size_pixel_vertical
    distance_z2 = (zoom_factor+1)*size_pixel_vertical
  else
    distance_x1 = 0.0
    distance_x2 = 2.0*size_pixel_horizontal
    distance_z1 = 0.0
    distance_z2 = 2.0*size_pixel_vertical
  endif

  do j=1,NZ_IMAGE_color
    do i=1,NX_IMAGE_color
      ! calculates midpoint of pixel
      xtmp = xmin_color_image + (i-1)*size_pixel_horizontal
      ztmp = zmin_color_image + (j-1)*size_pixel_vertical
      pixel_midpoint_x =  xtmp + 0.5*size_pixel_horizontal
      pixel_midpoint_z =  ztmp + 0.5*size_pixel_vertical

      ! avoid points on image border rim
      if (pixel_midpoint_x < xmin_color_image_loc &
        .or. pixel_midpoint_x > xmax_color_image_loc) cycle
      if (pixel_midpoint_z < zmin_color_image_loc &
        .or. pixel_midpoint_z > zmax_color_image_loc) cycle

      ! looks for closest point to midpoint of pixel
      dist_min_pixel = HUGEVAL
      do iglob=1,num_iglob_image_surface
        ! point location with respect to image surface
        x_loc = xcoord(iglob)
        z_loc = zcoord(iglob)

        ! checks if inside pixel range for larger numbers of points, minimizing computation time
        if (x_loc < xtmp - distance_x1 .or. x_loc > xtmp + distance_x2) cycle
        if (z_loc < ztmp - distance_z1 .or. z_loc > ztmp + distance_z2) cycle

        ! stores closest iglob
        x_loc = pixel_midpoint_x - x_loc
        z_loc = pixel_midpoint_z - z_loc
        dist_pixel = x_loc*x_loc + z_loc*z_loc
        if (dist_pixel < dist_min_pixel) then
          dist_min_pixel = dist_pixel
          dist_pixel_image(i,j) = dist_min_pixel
          iglob_image_color(i,j) = iglob_coord(iglob)
          ispec_image_color(i,j) = ispec_coord(iglob)
        endif
      enddo
    enddo
  enddo
  deallocate(xcoord,zcoord,iglob_coord,ispec_coord)

  ! gather info from other processes as well
  allocate(dist_pixel_recv(NX_IMAGE_color,0:NPROC-1),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1774')
  if (ier /= 0) call exit_mpi(myrank,'error allocating dist pixel recv')
  dist_pixel_recv(:,:) = HUGEVAL
  nb_pixel_loc = 0
  do j=1,NZ_IMAGE_color
    ! compares with other processes
    call gather_all_all_cr(dist_pixel_image(:,j),dist_pixel_recv,NX_IMAGE_color,NPROC)

    ! selects entries
    do i=1,NX_IMAGE_color
      ! note: minimum location will be between 1 and NPROC
      locval = minloc(dist_pixel_recv(i,:))
      irank = locval(1) - 1
      ! store only own best points
      if (irank == myrank .and. dist_pixel_recv(i,irank) < HUGEVAL) then
        ! increases count
        nb_pixel_loc = nb_pixel_loc + 1
      else
        ! resets index
        iglob_image_color(i,j) = -1
        ispec_image_color(i,j) = 0
      endif
    enddo
  enddo
  deallocate(dist_pixel_recv,dist_pixel_image)

  ! creating and filling array num_pixel_loc with the positions of each colored
  ! pixel owned by the local process (useful for parallel jobs)
  if (nb_pixel_loc > 0) then
    allocate(num_pixel_loc(nb_pixel_loc),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1775')
    if (ier /= 0) stop 'error allocating array num_pixel_loc'
  endif
  nb_pixel_loc = 0
  do i = 1, NX_IMAGE_color
    do j = 1, NZ_IMAGE_color
      if (iglob_image_color(i,j) /= -1) then
        nb_pixel_loc = nb_pixel_loc + 1
        num_pixel_loc(nb_pixel_loc) = (j-1)*NX_IMAGE_color + i
      endif
    enddo
  enddo
  ! checks if array is allocated
  if (nb_pixel_loc > 0) then
    if (.not. allocated(num_pixel_loc)) call exit_MPI(myrank,'error num_pixel_loc allocation')
  endif

  ! filling array iglob_image_color, containing info on which process owns which pixels.
  allocate(nb_pixel_per_proc(0:NPROC-1),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1776')
  if (ier /= 0) stop 'error allocating array nb_pixel_per_proc'
  nb_pixel_per_proc(:) = 0

  call gather_all_singlei(nb_pixel_loc,tmp_pixel_per_proc,NPROC)
  nb_pixel_per_proc(:) = tmp_pixel_per_proc(1,:)

  ! allocates receiving array
  if (myrank == 0) then
    allocate( num_pixel_recv(maxval(nb_pixel_per_proc(:)),0:NPROC-1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1777')
    if (ier /= 0) stop 'error allocating array num_pixel_recv'
  endif
  ! fills iglob_image_color index array
  if (NPROC > 1) then
    if (myrank == 0) then
      do iproc = 1, NPROC-1
        if (nb_pixel_per_proc(iproc) > 0) then
          call recv_i(num_pixel_recv(:,iproc),nb_pixel_per_proc(iproc),iproc,42)
          ! stores proc number instead where iglob_image_color wouldn't be defined (=-1)
          do k = 1, nb_pixel_per_proc(iproc)
            j = ceiling(real(num_pixel_recv(k,iproc)) / real(NX_IMAGE_color))
            i = num_pixel_recv(k,iproc) - (j-1)*NX_IMAGE_color
            iglob_image_color(i,j) = iproc
          enddo
        endif
      enddo
    else
      if (nb_pixel_loc > 0) then
        call send_i(num_pixel_loc(:),nb_pixel_loc,0,42)
      endif
    endif
  endif

  ! allocate an array for image data
  allocate(image_color_data(NX_IMAGE_color,NZ_IMAGE_color),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1778')
  allocate(image_color_background_display(NX_IMAGE_color,NZ_IMAGE_color),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1779')
  if (ier /= 0) call exit_mpi(myrank,'error allocating image data')
  image_color_data(:,:) = 0._CUSTOM_REAL
  image_color_background_display(:,:) = 0._CUSTOM_REAL

  if (myrank == 0) then
    allocate( data_pixel_recv(maxval(nb_pixel_per_proc(:))),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1780')
    if (ier /= 0) stop 'error allocating array data_pixel_recv'
    data_pixel_recv(:) = 0._CUSTOM_REAL
  endif
  if (nb_pixel_loc > 0) then
    allocate(data_pixel_send(nb_pixel_loc),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1781')
    if (ier /= 0) call exit_mpi(myrank,'error allocating image send data')
    data_pixel_send(:) = 0._CUSTOM_REAL
  endif

  ! handles wavefield background data
  if (PNM_IMAGE_TYPE >= 1 .and. PNM_IMAGE_TYPE <= 4) then
    call write_PNM_background()
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '******** '
    write(IMAIN,*)
  endif

  end subroutine write_PNM_initialize


!=============================================================

  subroutine write_PNM_background()

  use image_PNM_par

  use specfem_par, only: myrank,NPROC

  implicit none
  ! local parameters
  real(kind=CUSTOM_REAL) :: bg_val
  integer :: i,j,k,iglob,ispec,iproc

  ! getting velocity for each local pixels
  image_color_background_display(:,:) = 0.d0

  do k = 1, nb_pixel_loc
    j = ceiling(real(num_pixel_loc(k)) / real(NX_IMAGE_color))
    i = num_pixel_loc(k) - (j-1)*NX_IMAGE_color

    iglob = iglob_image_color(i,j)
    ispec = ispec_image_color(i,j)

    ! select background value
    select case(BACKGROUND_TYPE)
    case (1)
      call get_iglob_vp(iglob,ispec,bg_val)
    case (2)
      call get_iglob_vs(iglob,ispec,bg_val)
    case (3)
      call get_iglob_rho(iglob,ispec,bg_val)
    end select

    data_pixel_send(k) = bg_val
    image_color_background_display(i,j) = bg_val
  enddo

  ! MPI assembling array image_color_background_display on process zero for color output
  if (NPROC > 1) then
    ! main collects
    if (myrank == 0) then
      do iproc = 1, NPROC-1
        if (nb_pixel_per_proc(iproc) > 0) then
          call recvv_cr(data_pixel_recv(1),nb_pixel_per_proc(iproc),iproc,43)
          ! fills vp display array
          do k = 1, nb_pixel_per_proc(iproc)
            j = ceiling(real(num_pixel_recv(k,iproc)) / real(NX_IMAGE_color))
            i = num_pixel_recv(k,iproc) - (j-1)*NX_IMAGE_color
            image_color_background_display(i,j) = data_pixel_recv(k)
          enddo
        endif
      enddo
    else
      if (nb_pixel_loc > 0) then
        ! secondary processes send
        call sendv_cr(data_pixel_send,nb_pixel_loc,0,43)
      endif
    endif
  endif

  end subroutine write_PNM_background


!================================================================

  subroutine write_PNM_create_image()

! creates color PNM image

!! DK DK Jan 2013: here for performance and to reduce the size of the files, one day
!! DK DK Jan 2013: we should switch to using the JPEG library directly, as already implemented in SPECFEM2D

  use image_PNM_par

  use constants, only: NDIM
  use specfem_par, only: NPROC,it,myrank

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL),dimension(NDIM) :: val_vector
  real(kind=CUSTOM_REAL):: temp
  integer :: i,j,k,iglob,ispec,iproc

  ! initializes color data
  image_color_data(:,:) = 0.d0

  ! reads/retrieves color data
  do k = 1, nb_pixel_loc
    j = ceiling(real(num_pixel_loc(k)) / real(NX_IMAGE_color))
    i = num_pixel_loc(k) - (j-1)*NX_IMAGE_color

    ! global point and element indices of GLL point in this pixel
    iglob = iglob_image_color(i,j)
    ispec = ispec_image_color(i,j)

    ! data type
    select case (PNM_IMAGE_TYPE)
    case (1)
      ! x-velocity component
      ! gets velocity for point iglob
      call get_iglob_veloc(iglob,ispec,val_vector)
      temp = val_vector(1)
    case (2)
      ! y-velocity component
      ! gets velocity for point iglob
      call get_iglob_veloc(iglob,ispec,val_vector)
      temp = val_vector(2)
    case (3)
      ! z-velocity component
      ! gets velocity for point iglob
      call get_iglob_veloc(iglob,ispec,val_vector)
      temp = val_vector(3)
    case (4)
      ! velocity norm
      ! gets velocity for point iglob
      call get_iglob_veloc(iglob,ispec,val_vector)
      temp = sqrt( val_vector(1)**2 + val_vector(2)**2 + val_vector(3)**2)
    case (5)
      ! model vp
      ! gets model value for point iglob
      call get_iglob_vp(iglob,ispec,temp)
    case (6)
      ! model vs
      ! gets model value for point iglob
      call get_iglob_vs(iglob,ispec,temp)
    case (7)
      ! model rho
      ! gets model value for point iglob
      call get_iglob_rho(iglob,ispec,temp)
    case default
      call exit_MPI(myrank,'Error invalid PNM_IMAGE_TYPE for image data selection')
    end select

    ! stores data
    image_color_data(i,j) = temp
    data_pixel_send(k) = temp
  enddo

  ! MPI assembling array image_color_data on process zero for color output
  if (NPROC > 1) then
    if (myrank == 0) then
      do iproc = 1, NPROC-1
        if (nb_pixel_per_proc(iproc) > 0) then
          call recvv_cr(data_pixel_recv(1),nb_pixel_per_proc(iproc),iproc,43)
          ! distributes on image pixels
          do k = 1, nb_pixel_per_proc(iproc)
            j = ceiling(real(num_pixel_recv(k,iproc)) / real(NX_IMAGE_color))
            i = num_pixel_recv(k,iproc) - (j-1)*NX_IMAGE_color
            image_color_data(i,j) = data_pixel_recv(k)
          enddo
        endif
      enddo
    else
      if (nb_pixel_loc > 0) then
        ! secondary processes send
        call sendv_cr(data_pixel_send(1),nb_pixel_loc,0,43)
      endif
    endif
  endif
  ! synchronizes processes
  call synchronize_all()

  ! main process writes out file
  if (myrank == 0) then
    ! writes output file
    call write_PNM_data(image_color_data,iglob_image_color, &
                        NX_IMAGE_color,NZ_IMAGE_color,it,image_cutsnaps,image_color_background_display)
  endif

  end subroutine write_PNM_create_image


!================================================================


  subroutine write_PNM_data(color_image_2D_data,iglob_image_color_2D, &
                            NX,NY,it,cutsnaps,image_color_background_display)

! display a given field as a red and blue color image
! to display the snapshots : display image*.gif
! when compiling with Intel ifort, use " -assume byterecl " option to create binary PNM images

  use constants, only: HUGEVAL,TINYVAL,VERYSMALLVAL,CUSTOM_REAL,OUTPUT_FILES,MAX_STRING_LEN,IOUT

  use image_PNM_par, only: BINARY_FILE,POWER_DISPLAY_COLOR,BACKGROUND_TYPE,PNM_IMAGE_TYPE

  implicit none

  integer,intent(in) :: NX,NY,it
  real(kind=CUSTOM_REAL) :: cutsnaps

  integer, dimension(NX,NY) :: iglob_image_color_2D

  real(kind=CUSTOM_REAL), dimension(NX,NY) :: color_image_2D_data
  real(kind=CUSTOM_REAL), dimension(NX,NY) :: image_color_background_display

  ! local parameter
  integer :: ix,iy,R,G,B,tenthousands,thousands,hundreds,tens,units,remainder,current_rec
  real(kind=CUSTOM_REAL) :: amplitude_max,amplitude_min,normalized_value
  real(kind=CUSTOM_REAL) :: bg_min,bg_max,x1
  character(len=MAX_STRING_LEN) :: file_name
  ! ASCII code of character '0' and of carriage return character
  integer, parameter :: ascii_code_of_zero = 48, ascii_code_of_carriage_return = 10

  ! open the image file
  select case (PNM_IMAGE_TYPE)
  case (1,2,3,4)
    ! wavefield images
    write(file_name,"(a,'/image',i7.7,'.pnm')") OUTPUT_FILES(1:len_trim(OUTPUT_FILES)),it
    ! first delete the file, just in case it was previously bigger
    open(unit=IOUT,file=trim(file_name),status='unknown')
    close(unit=IOUT,status='delete')
  case (5)
    ! vp
    write(file_name,"(a,'/image_model_vp',i5.5,'.pnm')") OUTPUT_FILES(1:len_trim(OUTPUT_FILES)),it
  case (6)
    ! vs
    write(file_name,"(a,'/image_model_vs',i5.5,'.pnm')") OUTPUT_FILES(1:len_trim(OUTPUT_FILES)),it
  case (7)
    ! rho
    write(file_name,"(a,'/image_model_rho',i5.5,'.pnm')") OUTPUT_FILES(1:len_trim(OUTPUT_FILES)),it
  end select

  if (BINARY_FILE) then
    open(unit=IOUT,file=trim(file_name),status='unknown',access='direct',recl=1)
    write(IOUT,rec=1) 'P'
    write(IOUT,rec=2) '6' ! write P6 = binary PNM image format
    write(IOUT,rec=3) char(ascii_code_of_carriage_return)

    ! compute and write horizontal size
    remainder = NX

    tenthousands = remainder / 10000
    remainder = remainder - 10000 * tenthousands

    thousands = remainder / 1000
    remainder = remainder - 1000 * thousands

    hundreds = remainder / 100
    remainder = remainder - 100 * hundreds

    tens = remainder / 10
    remainder = remainder - 10 * tens

    units = remainder

    write(IOUT,rec=4) char(tenthousands + ascii_code_of_zero)
    write(IOUT,rec=5) char(thousands + ascii_code_of_zero)
    write(IOUT,rec=6) char(hundreds + ascii_code_of_zero)
    write(IOUT,rec=7) char(tens + ascii_code_of_zero)
    write(IOUT,rec=8) char(units + ascii_code_of_zero)
    write(IOUT,rec=9) ' '

    ! compute and write vertical size
    remainder = NY

    tenthousands = remainder / 10000
    remainder = remainder - 10000 * tenthousands

    thousands = remainder / 1000
    remainder = remainder - 1000 * thousands

    hundreds = remainder / 100
    remainder = remainder - 100 * hundreds

    tens = remainder / 10
    remainder = remainder - 10 * tens

    units = remainder

    write(IOUT,rec=10) char(tenthousands + ascii_code_of_zero)
    write(IOUT,rec=11) char(thousands + ascii_code_of_zero)
    write(IOUT,rec=12) char(hundreds + ascii_code_of_zero)
    write(IOUT,rec=13) char(tens + ascii_code_of_zero)
    write(IOUT,rec=14) char(units + ascii_code_of_zero)
    write(IOUT,rec=15) char(ascii_code_of_carriage_return)

    ! number of shades
    write(IOUT,rec=16) '2'
    write(IOUT,rec=17) '5'
    write(IOUT,rec=18) '5'
    write(IOUT,rec=19) char(ascii_code_of_carriage_return)

    ! block of image data starts at sixteenth character
    current_rec = 20
  else
    open(unit=IOUT,file=trim(file_name),status='unknown')
    write(IOUT,"('P3')") ! write P3 = ASCII PNM image format
    write(IOUT,*) NX,NY  ! write image size
    write(IOUT,*) '255'  ! number of shades
  endif

  ! compute background maximum amplitude
  bg_min = HUGEVAL
  bg_max = TINYVAL
  do iy=1,NY
    do ix=1,NX
      if (iglob_image_color_2D(ix,iy) > -1) then
        bg_min = min(bg_min,image_color_background_display(ix,iy))
        bg_max = max(bg_max,image_color_background_display(ix,iy))
      endif
    enddo
  enddo

  ! color data maximum
  amplitude_min = HUGEVAL
  amplitude_max = -HUGEVAL
  do iy=1,NY
    do ix=1,NX
      if (iglob_image_color_2D(ix,iy) > -1) then
        amplitude_min = min(amplitude_min,color_image_2D_data(ix,iy))
        amplitude_max = max(amplitude_max,color_image_2D_data(ix,iy))
      endif
    enddo
  enddo
  ! wavefields
  if (PNM_IMAGE_TYPE >= 1 .and. PNM_IMAGE_TYPE <= 4) then
    ! sets absolute maximum, to have symmetric color limits
    amplitude_max = maxval(abs(color_image_2D_data))
    if (amplitude_max < VERYSMALLVAL) amplitude_max = HUGEVAL
  endif

  ! debug
  !print *,'PNM image: file ',trim(file_name),' min/max = ',minval(color_image_2D_data),'/',maxval(color_image_2D_data), &
  !         'amplitude min/max = ',amplitude_min,amplitude_max

  ! in the PNM format, the image starts in the upper-left corner
  do iy=NY,1,-1
    do ix=1,NX
      ! check if pixel is defined or not (can be above topography for instance)
      if (iglob_image_color_2D(ix,iy) == -1) then
        ! use white (/ black /light blue) to display undefined region above topography
        R = 255 ! 0 !204
        G = 255 ! 0 !255
        B = 255 ! 0 !255
      else
        ! determine pixel color
        select case (PNM_IMAGE_TYPE)
        case (1,2,3,4)
          ! wavefields
          ! suppress small amplitudes considered as noise
          if (abs(color_image_2D_data(ix,iy)) < amplitude_max * cutsnaps) then
            ! determines background color
            if (BACKGROUND_TYPE > 0) then
              ! use velocity model as background where amplitude is negligible
              if (abs(bg_min) > TINYVAL) then
                if ((bg_max-bg_min)/bg_min > 0.02d0) then
                  if (abs(bg_max - bg_max) > TINYVAL) then
                    x1 = (image_color_background_display(ix,iy)-bg_min)/(bg_max-bg_min)
                  else
                    x1 = 0.5_CUSTOM_REAL
                  endif
                else
                  x1 = 0.5_CUSTOM_REAL
                endif
              else
                x1 = 0.5_CUSTOM_REAL
              endif

              ! rescale to avoid very dark gray levels
              x1 = x1*0.7 + 0.2
              if (x1 > 1._CUSTOM_REAL) x1 = 1._CUSTOM_REAL

              ! invert scale: white = bg_min, dark gray = bg_max
              x1 = 1._CUSTOM_REAL - x1

              ! map to [0,255]
              x1 = x1 * 255._CUSTOM_REAL

              R = nint(x1)
              if (R < 0) R = 0
              if (R > 255) R = 255
              G = R
              B = R
            else
              ! white
              R = 255
              G = 255
              B = 255
            endif

          else
            ! define normalized image data in [-1:1] and convert to nearest integer
            ! keeping in mind that data values can be negative
            if (abs(amplitude_max) > TINYVAL) then
              normalized_value = color_image_2D_data(ix,iy) / amplitude_max
            else
              normalized_value = 0.d0
            endif

            ! suppress values outside of [-1:+1]
            if (normalized_value < -1.d0) normalized_value = -1.d0
            if (normalized_value > 1.d0) normalized_value = 1.d0

            ! sets color value range [0,1]
            x1 = 2._CUSTOM_REAL * normalized_value - 1._CUSTOM_REAL

            ! sets color
            call get_color_rgb(x1,R,G,B)
          endif

        case(5,6,7)
          ! model values
          if (abs(amplitude_max - amplitude_min) > TINYVAL) then
            x1 = (color_image_2D_data(ix,iy)-amplitude_min)/(amplitude_max-amplitude_min)
          else
            ! uniform model, centers value for range [0,1]
            x1 = 0.5_CUSTOM_REAL
          endif

          ! rescale to saturate colors
          x1 = 1.5 * (x1 - 0.5) + 0.5

          ! limits x range [0,1]
          if (x1 > 1._CUSTOM_REAL) x1 = 1._CUSTOM_REAL
          if (x1 < 0._CUSTOM_REAL) x1 = 0._CUSTOM_REAL

          ! invert scale: white = bg_min, dark gray = bg_max
          x1 = 1._CUSTOM_REAL - x1

          ! sets colors
          call get_color_rgb(x1,R,G,B)
        end select
      endif

      ! write color image
      if (BINARY_FILE) then
        ! first write red
        write(IOUT,rec=current_rec) char(R)
        current_rec = current_rec + 1
        ! then write green
        write(IOUT,rec=current_rec) char(G)
        current_rec = current_rec + 1
        ! then write blue
        write(IOUT,rec=current_rec) char(B)
        current_rec = current_rec + 1
      else
        write(IOUT,"(i3,' ',i3,' ',i3)") R,G,B
      endif
    enddo
  enddo

  ! close the file
  close(IOUT)

  end subroutine write_PNM_data

!=============================================================

  subroutine get_iglob_vp(iglob,ispec,vp)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,FOUR_THIRDS

  use specfem_par, only: mustore,kappastore,rhostore,ibool,myrank
  use shared_parameters, only: ACOUSTIC_SIMULATION,ELASTIC_SIMULATION
  use specfem_par_elastic, only: rho_vp

  implicit none

  integer,intent(in) :: iglob,ispec
  real(kind=CUSTOM_REAL),intent(out):: vp

  !local parameters
  integer :: i,j,k

  ! returns first vp encountered for iglob index
  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        if (ibool(i,j,k,ispec) == iglob) then
          ! calculates vp
          if (ELASTIC_SIMULATION) then
            vp =  (FOUR_THIRDS * mustore(i,j,k,ispec) + kappastore(i,j,k,ispec)) / rho_vp(i,j,k,ispec)
          else if (ACOUSTIC_SIMULATION) then
            vp = sqrt( kappastore(i,j,k,ispec) / rhostore(i,j,k,ispec) )
          else
            call exit_mpi(myrank,'error vp not implemented')
          endif
          return
        endif
      enddo
    enddo
  enddo

  end subroutine get_iglob_vp

!=============================================================

  subroutine get_iglob_vs(iglob,ispec,vs)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,FOUR_THIRDS

  use specfem_par, only: mustore,rhostore,ibool

  implicit none

  integer,intent(in) :: iglob,ispec
  real(kind=CUSTOM_REAL),intent(out):: vs

  !local parameters
  integer :: i,j,k

  ! returns first vs encountered for iglob index
  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        if (ibool(i,j,k,ispec) == iglob) then
          ! calculates vs
          vs = sqrt( mustore(i,j,k,ispec) / rhostore(i,j,k,ispec) )
          return
        endif
      enddo
    enddo
  enddo

  end subroutine get_iglob_vs

!=============================================================

  subroutine get_iglob_rho(iglob,ispec,rho)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,FOUR_THIRDS

  use specfem_par, only: rhostore,ibool

  implicit none

  integer,intent(in) :: iglob,ispec
  real(kind=CUSTOM_REAL),intent(out):: rho

  !local parameters
  integer :: i,j,k

  ! returns first vs encountered for iglob index
  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        if (ibool(i,j,k,ispec) == iglob) then
          ! returns rho
          rho = rhostore(i,j,k,ispec)
          return
        endif
      enddo
    enddo
  enddo

  end subroutine get_iglob_rho

!=============================================================

  subroutine get_iglob_veloc(iglob,ispec,val_vector)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM
  use shared_parameters, only: ACOUSTIC_SIMULATION,ELASTIC_SIMULATION
  use specfem_par_acoustic, only: potential_acoustic,potential_dot_acoustic, &
                                  ispec_is_acoustic,b_potential_acoustic,b_potential_dot_acoustic
  use specfem_par_elastic, only: displ,veloc,ispec_is_elastic
  use specfem_par, only: SIMULATION_TYPE,SAVE_DISPLACEMENT,ibool

  implicit none

  integer,intent(in) :: iglob,ispec
  real(kind=CUSTOM_REAL),dimension(NDIM),intent(out):: val_vector

  ! local parameters
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ):: val_element
  integer :: i,j,k

  ! returns first element encountered for iglob index
  if (ELASTIC_SIMULATION) then
    if (ispec_is_elastic(ispec)) then
      if (SAVE_DISPLACEMENT) then
        if (SIMULATION_TYPE == 3) then
          ! to display re-constructed wavefield
          !val_vector(:) = b_displ(:,iglob)
          ! to display adjoint wavefield
          val_vector(:) = displ(:,iglob)
        else
          val_vector(:) = displ(:,iglob)
        endif
      else
        if (SIMULATION_TYPE == 3) then
          ! to display re-constructed wavefield
          !val_vector(:) = b_veloc(:,iglob)
          ! to display adjoint wavefield
          val_vector(:) = veloc(:,iglob)
        else
          val_vector(:) = veloc(:,iglob)
        endif
      endif

      ! returns with this result
      return
    endif
  endif

  if (ACOUSTIC_SIMULATION) then
    if (ispec_is_acoustic(ispec)) then
      if (SAVE_DISPLACEMENT) then
        if (SIMULATION_TYPE == 3) then
          ! displacement vector from backward potential
          call compute_gradient_in_acoustic(ispec,b_potential_acoustic,val_element)
        else
          ! displacement vector
          call compute_gradient_in_acoustic(ispec,potential_acoustic,val_element)
        endif
      else
        if (SIMULATION_TYPE == 3) then
          ! velocity vector for backward/reconstructed wavefield
          call compute_gradient_in_acoustic(ispec,b_potential_dot_acoustic,val_element)
        else
          ! velocity vector
          call compute_gradient_in_acoustic(ispec,potential_dot_acoustic,val_element)
        endif
      endif

      ! returns corresponding iglob velocity entry
      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX
            if (ibool(i,j,k,ispec) == iglob) then
              val_vector(:) = val_element(:,i,j,k)
              return
            endif
          enddo
        enddo
      enddo

    endif
  endif

  ! should not reach this point
  call exit_mpi(0,'error image velocity not found')

  end subroutine get_iglob_veloc

!=============================================================

  subroutine get_color_rgb(val_in,R,G,B)

! determines RGB values for input in range [0,1]
  use constants, only: CUSTOM_REAL
  use image_PNM_par, only: PNM_COLOR_PALETTE,POWER_DISPLAY_COLOR

  implicit none
  real(kind=CUSTOM_REAL), intent(in) :: val_in
  integer, intent(inout) :: R,G,B

  ! local parameters
  real(kind=CUSTOM_REAL) :: val,val_tmp,fac,x1,x2
  real(kind=CUSTOM_REAL) :: h,s,l
  real(kind=CUSTOM_REAL) :: h0,h1,h2,h2a,h2b,h3,h4,h5,h6
  real(kind=CUSTOM_REAL) :: l0,l1,l2,l2a,l2b,l3,l4,l5,l6
  real(kind=CUSTOM_REAL) :: s0,s1,s2,s2a,s2b,s3,s4,s5,s6
  real(kind=CUSTOM_REAL) :: split,delta,thr1,thr2,v0
  real(kind=CUSTOM_REAL),dimension(3) :: array1_h,array1_l,array1_s
  real(kind=CUSTOM_REAL),dimension(5) :: array2_h,array2_l,array2_s
  real(kind=CUSTOM_REAL),dimension(5,3) :: array_hsl
  integer :: i,ii,jj,incr
  logical :: done_color
  ! tomo colormap
  real, dimension(4,9) :: tomo_keyrgb = reshape( (/ 0.0,0.1,0.0,0.0,     0.2,0.8,0.0,0.0,    0.3,1.0,0.7,0.0, &
                                                    0.48,0.92,0.92,0.92, 0.5,0.92,0.92,0.92, 0.52,0.92,0.92,0.92, &
                                                    0.7,0.0,0.6,0.7,     0.8,0.0,0.0,0.8,    1.0,0.0,0.0,0.1 /), &
                                                 (/4,9/) )

  ! initializes
  val = val_in
  R = 255
  G = 255
  B = 255

  select case(PNM_COLOR_PALETTE)
  case (1)
    ! blue-red
    ! use red if above average value, blue if below, no green
    ! scale between [-1,1]
    val = 2._CUSTOM_REAL * val - 1._CUSTOM_REAL
    if (val >= 0._CUSTOM_REAL) then
      R = nint(255.0 * (val)**POWER_DISPLAY_COLOR)
      G = 0
      B = 0
    else
      R = 0
      G = 0
      B = nint(255.0 * abs(val)**POWER_DISPLAY_COLOR)
    endif

  case (2)
    ! rainbow (blue-cyan-green-yellow-orange-red)
    ! scale between [-1,1]
    val = 2._CUSTOM_REAL * val - 1._CUSTOM_REAL

    val_tmp = 1.5 - abs(2.0 * val - 1.0)
    R = nint(255.0 * clamp(val_tmp))

    val_tmp = 1.5 - abs(2.0 * val)
    G = nint(255.0 * clamp(val_tmp))

    val_tmp = 1.5 - abs(2.0 * val + 1.0)
    B = nint(255.0 * clamp(val_tmp))

  case (3)
    ! tomo
    ! uses value range [0,1]
    do i = 1,8
      x1 = tomo_keyrgb(1,i)
      x2 = tomo_keyrgb(1,i+1)
      if (val >= x1 .and. val <= x2 ) then
        ! interpolation factor
        fac = (val - x1)/(x2 - x1)
        R = nint( 255.0 * ((1.0 - fac) * tomo_keyrgb(2,i) + fac * tomo_keyrgb(2,i+1)) )
        G = nint( 255.0 * ((1.0 - fac) * tomo_keyrgb(3,i) + fac * tomo_keyrgb(3,i+1)) )
        B = nint( 255.0 * ((1.0 - fac) * tomo_keyrgb(4,i) + fac * tomo_keyrgb(4,i+1)) )
        ! value found and set, all done
        return
      endif
    enddo
    ! no value set yet, uses end values
    R = nint(255.0 * tomo_keyrgb(2,9))
    G = nint(255.0 * tomo_keyrgb(3,9))
    B = nint(255.0 * tomo_keyrgb(4,9))

  case (4)
    ! terrain (blue-dark-blue-white-brown-red-yellow-black)
    ! color function: val in [0,1]
    ! HSL (hue-saturation-lightness)
    ! http://www.workwithcolor.com/hsl-color-picker-01.htm?cp=0000FF
    ! turquise
    h0 = 185./360.
    s0 = 1.0
    l0 = 0.61
    ! dark blue
    h1 = 225./360.
    s1 = 0.35
    l1 = 0.28
    ! white (from blue)
    h2a = 240./360.
    s2a = 1.0
    l2a = 1.0
    ! white
    h2b = 20./360.
    s2b = 0.2
    l2b = 1.0
    ! brown-green
    h3 = 45./360.
    s3 = 0.75
    l3 = 0.34
    ! red
    h4 = 0.0
    s4 = 1.0
    l4 = 0.5
    ! yellow
    h5 = 61./360.
    s5 = 1.0
    l5 = 0.5
    ! black
    h6 = 61./360.
    s6 = 0.0
    l6 = 0.1
    ! turquise-blue-white
    array1_h(1:3) = (/ h0,h1,h2a /)
    array1_s(1:3) = (/ s0,s1,s2a /)
    array1_l(1:3) = (/ l0,l1,l2a /)
    ! white-brown-red-yellow-black
    array2_h(1:5) = (/ h2b,h3,h4,h5,h6 /)
    array2_s(1:5) = (/ s2b,s3,s4,s5,s6 /)
    array2_l(1:5) = (/ l2b,l3,l4,l5,l6 /)
    ! sets the split value between the 2 color ranges (separated by white)
    split = 0.5

    array_hsl(:,:) = 0.0
    h = 0.0; s = 0.0; l = 0.0

    do ii = 1,2
      h = 0.0
      s = 1.0
      l = 0.0
      ! value scale and number of color increments
      if (ii == 1) then
        incr = 2
        delta = split/incr
        v0 = 0.0
        array_hsl(1:3,1) = array1_h   ! (array1_h,array1_s,array1_l)
        array_hsl(1:3,2) = array1_s
        array_hsl(1:3,3) = array1_l
      else
        incr = 4
        delta = (1.0 - split)/incr
        v0 = split
        array_hsl(1:5,1) = array2_h
        array_hsl(1:5,2) = array2_s
        array_hsl(1:5,3) = array2_l
      endif

      done_color = .false.
      do jj = 1,incr
        thr1 = v0 + (jj-1) * delta
        thr2 = v0 + jj * delta
        if (val >= thr1 .and. val <= thr2) then
          fac = (val - thr1)/delta
          h1 = array_hsl(jj  ,1)
          h2 = array_hsl(jj+1,1)
          s1 = array_hsl(jj  ,2)
          s2 = array_hsl(jj+1,2)
          l1 = array_hsl(jj  ,3)
          l2 = array_hsl(jj+1,3)
          h = h1 + (h2-h1) * fac
          s = s1 + (s2-s1) * fac
          l = l1 + (l2-l1) * fac
          !print *,"range ",jj,v,fac,h,s,l
          done_color = .true.
          exit
        endif
      enddo
      if (done_color) exit
    enddo
    if (val >= 1.0) then
      h = 0.0
      s = 0.0
      l = 0.01
    endif
    !print *,"terrain hsl ",i,v,h,s,l
    call hsl_2_rgb(h, s, l, R, G, B)

  case default
    ! grayscale
    ! map to [0,255]
    val = val * 255._CUSTOM_REAL

    R = nint(val)
    if (R < 0) R = 0
    if (R > 255) R = 255
    G = R
    B = R
  end select

contains

    ! helper functions

    real(kind=CUSTOM_REAL) function clamp(x)

    implicit none
    real(kind=CUSTOM_REAL), intent(in) :: x
    real(kind=CUSTOM_REAL), parameter :: xmin = 0._CUSTOM_REAL, xmax = 1._CUSTOM_REAL

    if (x < xmin) then
      clamp = xmin
    else if (x > xmax) then
      clamp = xmax
    else
      clamp = x
    endif

    end function clamp

    !------------------------------------------------------

    subroutine hsl_2_rgb(h_in, s_in, l_in, R, G, B)

    implicit none
    real(kind=CUSTOM_REAL), intent(in) :: h_in,s_in,l_in
    integer, intent(out) :: R,G,B
    ! local
    real(kind=CUSTOM_REAL) :: h,s,l,v,m,sv,fract,vsf,mid1,mid2
    integer :: sextant
    ! HSL in range 0-1 to rgb in range 0-1 (H-hue, S-saturation, L-lightness)
    ! see: http://geekymonkey.com/Programming/CSharp/RGB2HSL_HSL2RGB.htm

    ! default (RGB values between [0,255])
    R = 255
    G = 255
    B = 255

    h = h_in
    s = s_in
    l = l_in

    if (l <= 0.5) then
      v = l * (1.0 + s)
    else
      v = l + s - l * s
    endif

    if (v > 0) then
      m = l + l - v
      sv = (v - m ) / v
      h = h * 6.0

      sextant = int(h)
      fract = h - sextant
      vsf = v * sv * fract
      mid1 = m + vsf
      mid2 = v - vsf

      if (sextant == 0) then
        R = nint(255.0 * v)
        G = nint(255.0 * mid1)
        B = nint(255.0 * m)
      else if (sextant == 1) then
        R = nint(255.0 * mid2)
        G = nint(255.0 * v)
        B = nint(255.0 * m)
      else if (sextant == 2) then
        R = nint(255.0 * m)
        G = nint(255.0 * v)
        B = nint(255.0 * mid1)
      else if (sextant == 3) then
        R = nint(255.0 * m)
        G = nint(255.0 * mid2)
        B = nint(255.0 * v)
      else if (sextant == 4) then
        R = nint(255.0 * mid1)
        G = nint(255.0 * m)
        B = nint(255.0 * v)
      else if (sextant == 5) then
        R = nint(255.0 * v)
        G = nint(255.0 * m)
        B = nint(255.0 * mid2)
      endif
    endif
    end subroutine hsl_2_rgb

  end subroutine get_color_rgb

!
!-------------------------------------------------------------------------------------------------
!

  subroutine detect_surface_PNM_image(NPROC,nglob,nspec,ibool, &
                                      ispec_is_image_surface, &
                                      iglob_is_image_surface, &
                                      num_iglob_image_surface, &
                                      num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                      nibool_interfaces_ext_mesh,my_neighbors_ext_mesh, &
                                      ibool_interfaces_ext_mesh, &
                                      xstore,ystore,zstore,myrank)

! this returns points on a cross-section surface through model
!
! returns: ispec_is_image_surface, iglob_is_image_surface & num_iglob_image_surface

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

  use image_PNM_par, only: PNM_section_xorg,PNM_section_yorg,PNM_section_zorg, &
                           PNM_section_nx,PNM_section_ny,PNM_section_nz

  implicit none

  ! global indexing
  integer,intent(in) :: NPROC,nglob,nspec,myrank
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool

  ! surface
  logical, dimension(nspec),intent(inout) :: ispec_is_image_surface
  logical, dimension(nglob),intent(inout) :: iglob_is_image_surface
  integer,intent(inout) :: num_iglob_image_surface

  ! MPI partitions
  integer,intent(in) :: num_interfaces_ext_mesh
  integer,intent(in) :: max_nibool_interfaces_ext_mesh
  integer,dimension(num_interfaces_ext_mesh),intent(in) :: nibool_interfaces_ext_mesh
  integer,dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),intent(in) :: ibool_interfaces_ext_mesh
  integer,dimension(num_interfaces_ext_mesh),intent(in) :: my_neighbors_ext_mesh

  ! mesh global point coordinates
  real(kind=CUSTOM_REAL), dimension(nglob),intent(in) :: xstore,ystore,zstore

  ! local parameters
  real(kind=CUSTOM_REAL) :: min_dist,distance
  integer, dimension(:), allocatable :: valence
  integer :: ispec,i,j,k,iglob,ier,countval
  real(kind=CUSTOM_REAL),parameter :: TOLERANCE_DISTANCE = 0.9

  ! detecting surface points/elements (based on valence check on NGLL points) for external mesh
  allocate(valence(nglob),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1225')
  if (ier /= 0) stop 'error allocate valence array'

  ! initialize surface indices
  ispec_is_image_surface(:) = .false.
  iglob_is_image_surface(:) = .false.
  valence(:) = 0
  num_iglob_image_surface = 0

  ! an estimation of the minimum distance between global points
  min_dist = minval( (xstore(ibool(1,1,1,:)) - xstore(ibool(2,1,1,:)))**2 &
                   + (ystore(ibool(1,1,1,:)) - ystore(ibool(2,1,1,:)))**2 &
                   + (zstore(ibool(1,1,1,:)) - zstore(ibool(2,1,1,:)))**2)
  min_dist = sqrt(min_dist)
  distance = TOLERANCE_DISTANCE*min_dist

  ! sets valence value to one corresponding to process rank  for points on cross-sections
  do ispec = 1, nspec
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool(i,j,k,ispec)

          ! chooses points close to cross-section
          if (abs((xstore(iglob)-PNM_section_xorg)*PNM_section_nx &
                + (ystore(iglob)-PNM_section_yorg)*PNM_section_ny &
                + (zstore(iglob)-PNM_section_zorg)*PNM_section_nz ) < distance) then
            ! sets valence to 1 for points on cross-sections
            valence(iglob) = myrank+1
          endif
        enddo
      enddo
    enddo
  enddo

  ! adds contributions from different partitions to valence
  call assemble_MPI_scalar_i_blocking(NPROC,nglob,valence, &
                                      num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                      nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                      my_neighbors_ext_mesh)


  ! determines spectral elements containing points on surface
  countval = 0
  do ispec = 1, nspec
    ! loops over GLL points not on edges or corners, but inside faces
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool(i,j,k,ispec)
          ! considers only points in same process or, if point is shared between two processes,
          ! only with higher process ranks than itself
          if (valence(iglob) == myrank+1 .or. valence(iglob) > 2*(myrank+1)) then
            if (iglob_is_image_surface(iglob) .eqv. .false. ) countval = countval + 1
            iglob_is_image_surface(iglob) = .true.
            ispec_is_image_surface(ispec) = .true.
          endif
        enddo
      enddo
    enddo
  enddo ! nspec
  num_iglob_image_surface = countval

  end subroutine detect_surface_PNM_image
