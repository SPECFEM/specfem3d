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

  module image_PNM_par

  use constants,only: CUSTOM_REAL,IMAIN
  use specfem_par,only: myrank,NPROC,it

  ! ----------------------------------------------
  ! USER PARAMETER

  ! image data output:
  !   type = 1 : displ/velocity x-component
  !   type = 2 : displ/velocity y-component
  !   type = 3 : displ/velocity z-component
  !   type = 4 : displ/velocity norm
  integer,parameter:: IMAGE_TYPE = 3 ! 4

  ! cross-section surface
  ! cross-section origin point
  real(kind=CUSTOM_REAL),parameter:: section_xorg = 0.0 ! 67000.0
  real(kind=CUSTOM_REAL),parameter:: section_yorg = 0.0 ! 0.0
  real(kind=CUSTOM_REAL),parameter:: section_zorg = -100.0 ! 0.0

  ! cross-section surface normal
  real(kind=CUSTOM_REAL),parameter:: section_nx = 0.0 !1.0
  real(kind=CUSTOM_REAL),parameter:: section_ny = 0.0 !0.0
  real(kind=CUSTOM_REAL),parameter:: section_nz = 1.0 !0.0

  ! cross-section (in-plane) horizontal-direction
  real(kind=CUSTOM_REAL),parameter:: section_hdirx = 1.0 ! 0.0
  real(kind=CUSTOM_REAL),parameter:: section_hdiry = 0.0 !1.0
  real(kind=CUSTOM_REAL),parameter:: section_hdirz = 0.0 ! 0.0

  ! cross-section (in-plane) vertical-direction
  real(kind=CUSTOM_REAL),parameter:: section_vdirx = 0.0 ! 0.0
  real(kind=CUSTOM_REAL),parameter:: section_vdiry = 1.0 ! 0.0
  real(kind=CUSTOM_REAL),parameter:: section_vdirz = 0.0 ! 1.0

  ! non linear display to enhance small amplitudes in color images
  real(kind=CUSTOM_REAL), parameter :: POWER_DISPLAY_COLOR = 0.30_CUSTOM_REAL

  ! amplitude threshold
  real(kind=CUSTOM_REAL),parameter :: image_cutsnaps  = 1.e-2

  ! use vp as gray background
  logical, parameter :: VP_BACKGROUND = .false.

  ! create temporary image files in binary PNM P6 format (smaller)
  ! or ASCII PNM P3 format (easier to edit)
  logical, parameter :: BINARY_FILE = .true.

  ! ----------------------------------------------

  ! image data
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable:: image_color_vp_display
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
  use specfem_par,only: NGLOB_AB,NSPEC_AB,ibool,xstore,ystore,zstore,&
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,my_neighbours_ext_mesh, &
                        ibool_interfaces_ext_mesh !,prname
  use constants,only: HUGEVAL,NGLLX,NGLLY,NGLLZ
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
  integer:: npgeo,npgeo_glob
  integer:: i,j,k,iproc,iglob,ispec,ier
  ! data from mesh
  real(kind=CUSTOM_REAL),dimension(:),allocatable:: xcoord,zcoord
  integer,dimension(:),allocatable :: iglob_coord,ispec_coord
  logical,dimension(:),allocatable:: ispec_is_image_surface,iglob_is_image_surface
  integer :: num_iglob_image_surface
  integer :: count,loc(1),irank
  !character(len=256) :: vtkfilename
  integer :: zoom_factor = 4
  logical :: zoom
  integer, dimension(1) :: tmp_pixel_loc
  integer, dimension(1,0:NPROC-1) :: tmp_pixel_per_proc

  ! checks image type
  if(IMAGE_TYPE > 4 .or. IMAGE_TYPE < 1) then
    call exit_mpi(myrank,'That type is not implemented for PNM images yet')
  endif

  ! user output
  if( myrank == 0 ) then
    write(IMAIN,*)
    write(IMAIN,*) '********'
    !   type = 1 : velocity V_x component
    if( IMAGE_TYPE == 1 ) write(IMAIN,*) 'PNM image: velocity V_x component'
    !   type = 2 : velocity V_y component
    if( IMAGE_TYPE == 2 ) write(IMAIN,*) 'PNM image: velocity V_y component'
    !   type = 3 : velocity V_z component
    if( IMAGE_TYPE == 3 ) write(IMAIN,*) 'PNM image: velocity V_z component'
    !   type = 4 : velocity V norm
    if( IMAGE_TYPE == 4 ) write(IMAIN,*) 'PNM image: velocity norm'
  endif

  ! finds global points on image surface
  allocate(ispec_is_image_surface(NSPEC_AB),iglob_is_image_surface(NGLOB_AB),stat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error allocating image ispec and iglob')

  call detect_surface_PNM_image(NPROC,NGLOB_AB,NSPEC_AB,ibool,&
                            ispec_is_image_surface, &
                            iglob_is_image_surface, &
                            num_iglob_image_surface, &
                            num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                            nibool_interfaces_ext_mesh,my_neighbours_ext_mesh, &
                            ibool_interfaces_ext_mesh,&
                            section_xorg,section_yorg,section_zorg, &
                            section_nx,section_ny,section_nz, &
                            xstore,ystore,zstore,myrank)

  ! extracts points on surface
  allocate( xcoord(num_iglob_image_surface),&
           zcoord(num_iglob_image_surface),&
           iglob_coord(num_iglob_image_surface),&
           ispec_coord(num_iglob_image_surface),stat=ier )
  if( ier /= 0 ) call exit_mpi(myrank,'error allocating xyz image coordinates')

  count=0
  do ispec=1,NSPEC_AB
    if( ispec_is_image_surface(ispec) ) then
      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX
            iglob = ibool(i,j,k,ispec)
            if( iglob_is_image_surface(iglob) ) then
              count = count+1
              ! coordinates with respect to horizontal and vertical direction
              xcoord(count)= xstore(iglob)*section_hdirx &
                                + ystore(iglob)*section_hdiry &
                                + zstore(iglob)*section_hdirz
              zcoord(count)= xstore(iglob)*section_vdirx &
                                + ystore(iglob)*section_vdiry &
                                + zstore(iglob)*section_vdirz
              iglob_coord(count) = iglob
              ispec_coord(count) = ispec

              ! reset iglob flag
              iglob_is_image_surface(iglob) = .false.
            endif
          enddo
        enddo
      enddo
    endif
  enddo

  if( count /= num_iglob_image_surface) call exit_mpi(myrank,'error image point number')

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
  if( NX_IMAGE_color > 4096 .or. NZ_IMAGE_color > 4096 ) then
    ! half of it
    NX_IMAGE_color = NX_IMAGE_color / 2
    NZ_IMAGE_color = NZ_IMAGE_color / 2
    ! even numbers
    NX_IMAGE_color = 2 * (NX_IMAGE_color / 2)
    NZ_IMAGE_color = 2 * (NZ_IMAGE_color / 2)
  endif

  ! ...and not too small
  zoom = .false.
  if( NX_IMAGE_color < 200 .or. NZ_IMAGE_color < 200 ) then
    ! increases it
    NX_IMAGE_color = NX_IMAGE_color * zoom_factor
    NZ_IMAGE_color = NZ_IMAGE_color * zoom_factor
    zoom = .true.
  endif

  ! create all the pixels
  if( NX_IMAGE_color /= 0 ) then
    size_pixel_horizontal = (xmax_color_image - xmin_color_image) / dble(NX_IMAGE_color)
  else
    size_pixel_horizontal = 0.0
  endif

  if( NZ_IMAGE_color /= 0 ) then
    size_pixel_vertical = (zmax_color_image - zmin_color_image) / dble(NZ_IMAGE_color)
  else
    size_pixel_vertical = 0.0
  endif

  if (myrank == 0) then
    write(IMAIN,*) '  image points: ',npgeo_glob
    write(IMAIN,*) '  xmin/xmax: ',xmin_color_image,'/',xmax_color_image
    write(IMAIN,*) '  zmin/zmax: ',zmin_color_image,'/',zmax_color_image
    write(IMAIN,*) '  pixel numbers: ',NX_IMAGE_color,' x ',NZ_IMAGE_color
    write(IMAIN,*) '  pixel sizes  : ',size_pixel_horizontal,' x ',size_pixel_vertical
  endif

  ! allocate an array for the grid point that corresponds to a given image data point
  allocate(iglob_image_color(NX_IMAGE_color,NZ_IMAGE_color),&
          ispec_image_color(NX_IMAGE_color,NZ_IMAGE_color),stat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error allocating iglob_image_color')

  allocate(dist_pixel_image(NX_IMAGE_color,NZ_IMAGE_color),stat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error allocating dist pixel image')

  iglob_image_color(:,:) = -1
  ispec_image_color(:,:) = 0
  dist_pixel_image(:,:) = HUGEVAL

  if( zoom ) then
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
      if( pixel_midpoint_x < xmin_color_image_loc &
        .or. pixel_midpoint_x > xmax_color_image_loc ) cycle
      if( pixel_midpoint_z < zmin_color_image_loc &
        .or. pixel_midpoint_z > zmax_color_image_loc ) cycle

      ! looks for closest point to midpoint of pixel
      dist_min_pixel = HUGEVAL
      do iglob=1,num_iglob_image_surface
        ! point location with respect to image surface
        x_loc = xcoord(iglob)
        z_loc = zcoord(iglob)

        ! checks if inside pixel range for larger numbers of points, minimizing computation time
        if( x_loc < xtmp - distance_x1 .or. x_loc > xtmp + distance_x2 ) cycle
        if( z_loc < ztmp - distance_z1 .or. z_loc > ztmp + distance_z2 ) cycle

        ! stores closest iglob
        x_loc = pixel_midpoint_x - x_loc
        z_loc = pixel_midpoint_z - z_loc
        dist_pixel = x_loc*x_loc + z_loc*z_loc
        if( dist_pixel < dist_min_pixel) then
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
  if(ier /= 0 ) call exit_mpi(myrank,'error allocating dist pixel recv')
  dist_pixel_recv(:,:) = HUGEVAL
  nb_pixel_loc = 0
  do j=1,NZ_IMAGE_color
    ! compares with other processes
    call gather_all_all_cr(dist_pixel_image(:,j),dist_pixel_recv,NX_IMAGE_color,NPROC)

    ! selects entries
    do i=1,NX_IMAGE_color
      ! note: minimum location will be between 1 and NPROC
      loc = minloc(dist_pixel_recv(i,:))
      irank = loc(1) - 1
      ! store only own best points
      if( irank == myrank .and. dist_pixel_recv(i,irank) < HUGEVAL) then
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
  if( nb_pixel_loc > 0 ) then
    allocate(num_pixel_loc(nb_pixel_loc),stat=ier)
    if( ier /= 0 ) stop 'error allocating array num_pixel_loc'
  endif
  nb_pixel_loc = 0
  do i = 1, NX_IMAGE_color
    do j = 1, NZ_IMAGE_color
      if ( iglob_image_color(i,j) /= -1 ) then
        nb_pixel_loc = nb_pixel_loc + 1
        num_pixel_loc(nb_pixel_loc) = (j-1)*NX_IMAGE_color + i
      endif
    enddo
  enddo
  ! checks if array is allocated
  if( nb_pixel_loc > 0 ) then
    if( .not. allocated(num_pixel_loc) ) call exit_MPI(myrank,'error num_pixel_loc allocation')
  endif

  ! filling array iglob_image_color, containing info on which process owns which pixels.
  allocate(nb_pixel_per_proc(0:NPROC-1),stat=ier)
  if( ier /= 0 ) stop 'error allocating array nb_pixel_per_proc'

  tmp_pixel_loc(1) = nb_pixel_loc
  call gather_all_i(tmp_pixel_loc,1,tmp_pixel_per_proc,1,NPROC)
  nb_pixel_per_proc(:) = tmp_pixel_per_proc(1,:)

  ! allocates receiving array
  if ( myrank == 0 ) then
    allocate( num_pixel_recv(maxval(nb_pixel_per_proc(:)),0:NPROC-1),stat=ier)
    if( ier /= 0 ) stop 'error allocating array num_pixel_recv'
  endif
  ! fills iglob_image_color index array
  if( NPROC > 1 ) then
    if (myrank == 0) then
      do iproc = 1, NPROC-1
        if( nb_pixel_per_proc(iproc) > 0 ) then
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
      if( nb_pixel_loc > 0 ) then
        call send_i(num_pixel_loc(:),nb_pixel_loc,0,42)
      endif
    endif
  endif

  ! allocate an array for image data
  allocate(image_color_data(NX_IMAGE_color,NZ_IMAGE_color), &
          image_color_vp_display(NX_IMAGE_color,NZ_IMAGE_color),stat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error allocating image data')

  image_color_data(:,:) = 0._CUSTOM_REAL
  image_color_vp_display(:,:) = 0._CUSTOM_REAL

  if ( myrank == 0 ) then
    allocate( data_pixel_recv(maxval(nb_pixel_per_proc(:))),stat=ier)
    if( ier /= 0 ) stop 'error allocating array data_pixel_recv'
    data_pixel_recv(:) = 0._CUSTOM_REAL
  endif
  if( nb_pixel_loc > 0 ) then
    allocate(data_pixel_send(nb_pixel_loc),stat=ier)
    if(ier /= 0 ) call exit_mpi(myrank,'error allocating image send data')
    data_pixel_send(:) = 0._CUSTOM_REAL
  endif

  ! handles vp background data
  call write_PNM_vp_background()

  ! user output
  if( myrank == 0 ) then
    write(IMAIN,*) '******** '
    write(IMAIN,*)
  endif


  end subroutine write_PNM_initialize


!=============================================================


  subroutine write_PNM_vp_background

  use image_PNM_par
  use specfem_par,only:myrank
  implicit none
  ! local parameters
  real(kind=CUSTOM_REAL) :: vp
  integer :: i,j,k,iglob,ispec,iproc

  ! getting velocity for each local pixels
  image_color_vp_display(:,:) = 0.d0

  do k = 1, nb_pixel_loc
    j = ceiling(real(num_pixel_loc(k)) / real(NX_IMAGE_color))
    i = num_pixel_loc(k) - (j-1)*NX_IMAGE_color

    iglob = iglob_image_color(i,j)
    ispec = ispec_image_color(i,j)
    call get_iglob_vp(iglob,ispec,vp)

    data_pixel_send(k) = vp
    image_color_vp_display(i,j) = vp
  enddo

  ! MPI assembling array image_color_vp_display on process zero for color output
  if (NPROC > 1) then
    ! master collects
    if (myrank == 0) then
      do iproc = 1, NPROC-1
        if( nb_pixel_per_proc(iproc) > 0 ) then
          call recvv_cr(data_pixel_recv(1),nb_pixel_per_proc(iproc),iproc,43)
          ! fills vp display array
          do k = 1, nb_pixel_per_proc(iproc)
            j = ceiling(real(num_pixel_recv(k,iproc)) / real(NX_IMAGE_color))
            i = num_pixel_recv(k,iproc) - (j-1)*NX_IMAGE_color
            image_color_vp_display(i,j) = data_pixel_recv(k)
          enddo
        endif
      enddo
    else
      if( nb_pixel_loc > 0 ) then
        ! slave processes send
        call sendv_cr(data_pixel_send,nb_pixel_loc,0,43)
      endif
    endif
  endif

  end subroutine write_PNM_vp_background


!================================================================

  subroutine write_PNM_create_image

! creates color PNM image

!! DK DK Jan 2013: here for performance and to reduce the size of the files, one day
!! DK DK Jan 2013: we should switch to using the JPEG library directly, as already implemented in SPECFEM2D

  use image_PNM_par
  use constants,only: NDIM
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

    ! global point and element indices of gll point in this pixel
    iglob = iglob_image_color(i,j)
    ispec = ispec_image_color(i,j)

    ! gets velocity for point iglob
    call get_iglob_veloc(iglob,ispec,val_vector)

    ! data type
    if( IMAGE_TYPE == 4 ) then
      ! velocity norm
      temp = sqrt( val_vector(1)**2 + val_vector(2)**2 + val_vector(3)**2 )
    else
      ! velocity component
      temp = val_vector(IMAGE_TYPE)
    endif

    ! stores data
    image_color_data(i,j) = temp
    data_pixel_send(k) = temp
  enddo

  ! MPI assembling array image_color_data on process zero for color output
  if (NPROC > 1) then
    if (myrank == 0) then
      do iproc = 1, NPROC-1
        if( nb_pixel_per_proc(iproc) > 0 ) then
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
      if( nb_pixel_loc > 0 ) then
        ! slave processes send
        call sendv_cr(data_pixel_send(1),nb_pixel_loc,0,43)
      endif
    endif
  endif

  ! master process writes out file
  if (myrank == 0) then
    ! writes output file
    call write_PNM_data(image_color_data,iglob_image_color,&
                            NX_IMAGE_color,NZ_IMAGE_color,it,image_cutsnaps,image_color_vp_display)
  endif


  end subroutine write_PNM_create_image


!================================================================


  subroutine write_PNM_data(color_image_2D_data,iglob_image_color_2D,&
                                NX,NY,it,cutsnaps,image_color_vp_display)

! display a given field as a red and blue color image
! to display the snapshots : display image*.gif
! when compiling with Intel ifort, use " -assume byterecl " option to create binary PNM images
  use constants,only: HUGEVAL,TINYVAL,CUSTOM_REAL,OUTPUT_FILES_PATH
  use image_PNM_par,only: BINARY_FILE,VP_BACKGROUND,&
                        POWER_DISPLAY_COLOR
  implicit none

  integer :: NX,NY,it
  real(kind=CUSTOM_REAL) :: cutsnaps

  integer, dimension(NX,NY) :: iglob_image_color_2D

  real(kind=CUSTOM_REAL), dimension(NX,NY) :: color_image_2D_data
  real(kind=CUSTOM_REAL), dimension(NX,NY) :: image_color_vp_display

  ! local parameter
  integer :: ix,iy,R,G,B,tenthousands,thousands,hundreds,tens,units,remainder,current_rec
  real(kind=CUSTOM_REAL) :: amplitude_max,normalized_value,vpmin,vpmax,x1
  character(len=256) :: file_name
  ! ASCII code of character '0' and of carriage return character
  integer, parameter :: ascii_code_of_zero = 48, ascii_code_of_carriage_return = 10

  ! open the image file
  write(file_name,"(a,'/image',i7.7,'.pnm')") OUTPUT_FILES_PATH(1:len_trim(OUTPUT_FILES_PATH)),it

! first delete the file, just in case it was previously bigger
  open(unit=27,file=file_name,status='unknown')
  close(unit=27,status='delete')

  if(BINARY_FILE) then
    open(unit=27,file=file_name,status='unknown',access='direct',recl=1)
    write(27,rec=1) 'P'
    write(27,rec=2) '6' ! write P6 = binary PNM image format
    write(27,rec=3) char(ascii_code_of_carriage_return)

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

    write(27,rec=4) char(tenthousands + ascii_code_of_zero)
    write(27,rec=5) char(thousands + ascii_code_of_zero)
    write(27,rec=6) char(hundreds + ascii_code_of_zero)
    write(27,rec=7) char(tens + ascii_code_of_zero)
    write(27,rec=8) char(units + ascii_code_of_zero)
    write(27,rec=9) ' '

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

    write(27,rec=10) char(tenthousands + ascii_code_of_zero)
    write(27,rec=11) char(thousands + ascii_code_of_zero)
    write(27,rec=12) char(hundreds + ascii_code_of_zero)
    write(27,rec=13) char(tens + ascii_code_of_zero)
    write(27,rec=14) char(units + ascii_code_of_zero)
    write(27,rec=15) char(ascii_code_of_carriage_return)

    ! number of shades
    write(27,rec=16) '2'
    write(27,rec=17) '5'
    write(27,rec=18) '5'
    write(27,rec=19) char(ascii_code_of_carriage_return)

    ! block of image data starts at sixteenth character
    current_rec = 20
  else
    open(unit=27,file=file_name,status='unknown')
    write(27,"('P3')") ! write P3 = ASCII PNM image format
    write(27,*) NX,NY  ! write image size
    write(27,*) '255'  ! number of shades
  endif

  ! compute maximum amplitude
  vpmin = HUGEVAL
  vpmax = TINYVAL
  do iy=1,NY
    do ix=1,NX
      if ( iglob_image_color_2D(ix,iy) > -1 ) then
        vpmin = min(vpmin,image_color_vp_display(ix,iy))
        vpmax = max(vpmax,image_color_vp_display(ix,iy))
      endif
    enddo
  enddo
  amplitude_max = maxval(abs(color_image_2D_data))
  if( amplitude_max < TINYVAL ) amplitude_max = HUGEVAL

  ! in the PNM format, the image starts in the upper-left corner
  do iy=NY,1,-1
    do ix=1,NX
      ! check if pixel is defined or not (can be above topography for instance)
      if(iglob_image_color_2D(ix,iy) == -1) then
        ! use black (/light blue) to display undefined region above topography
        R = 0 !204
        G = 0 !255
        B = 0 !255

      ! suppress small amplitudes considered as noise
      else if (abs(color_image_2D_data(ix,iy)) < amplitude_max * cutsnaps) then

        if( VP_BACKGROUND ) then
          ! use P velocity model as background where amplitude is negligible
          if((vpmax-vpmin)/vpmin > 0.02d0) then
            x1 = (image_color_vp_display(ix,iy)-vpmin)/(vpmax-vpmin)
          else
            x1 = 0.5d0
          endif

          ! rescale to avoid very dark gray levels
          x1 = x1*0.7 + 0.2
          if(x1 > 1.d0) x1=1.d0

          ! invert scale: white = vpmin, dark gray = vpmax
          x1 = 1.d0 - x1

          ! map to [0,255]
          x1 = x1 * 255.d0

          R = nint(x1)
          if(R < 0) R = 0
          if(R > 255) R = 255
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
        normalized_value = color_image_2D_data(ix,iy) / amplitude_max

        ! suppress values outside of [-1:+1]
        if(normalized_value < -1.d0) normalized_value = -1.d0
        if(normalized_value > 1.d0) normalized_value = 1.d0

        ! use red if positive value, blue if negative, no green
        if(normalized_value >= 0.d0) then
          R = nint(255.d0*normalized_value**POWER_DISPLAY_COLOR)
          G = 0
          B = 0
        else
          R = 0
          G = 0
          B = nint(255.d0*abs(normalized_value)**POWER_DISPLAY_COLOR)
        endif
      endif

      ! write color image
      if(BINARY_FILE) then
        ! first write red
        write(27,rec=current_rec) char(R)
        current_rec = current_rec + 1
        ! then write green
        write(27,rec=current_rec) char(G)
        current_rec = current_rec + 1
        ! then write blue
        write(27,rec=current_rec) char(B)
        current_rec = current_rec + 1
      else
        write(27,"(i3,' ',i3,' ',i3)") R,G,B
      endif
    enddo
  enddo

  ! close the file
  close(27)

  end subroutine write_PNM_data

!=============================================================

  subroutine get_iglob_vp(iglob,ispec,vp)

  use constants,only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,FOUR_THIRDS
  use specfem_par,only: mustore,kappastore,rhostore,ibool,myrank
  use specfem_par_acoustic,only: ACOUSTIC_SIMULATION
  use specfem_par_elastic,only: ELASTIC_SIMULATION,rho_vp
  implicit none

  integer,intent(in) :: iglob,ispec
  real(kind=CUSTOM_REAL),intent(out):: vp

  !local parameters
  integer :: i,j,k

  ! returns first vp encountered for iglob index
  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        if( ibool(i,j,k,ispec) == iglob ) then
          ! calculates vp
          if( ELASTIC_SIMULATION ) then
            vp =  (FOUR_THIRDS * mustore(i,j,k,ispec) + kappastore(i,j,k,ispec)) / rho_vp(i,j,k,ispec)
          else if( ACOUSTIC_SIMULATION ) then
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

  subroutine get_iglob_veloc(iglob,ispec,val_vector)

  use constants,only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM
  use specfem_par_acoustic,only: ACOUSTIC_SIMULATION,potential_acoustic,potential_dot_acoustic, &
                                ispec_is_acoustic, &
                                b_potential_acoustic,b_potential_dot_acoustic
  use specfem_par_elastic,only: ELASTIC_SIMULATION,displ,veloc, &
                                ispec_is_elastic ! ,b_displ,b_veloc
  use specfem_par,only: NSPEC_AB,NGLOB_AB,hprime_xx,hprime_yy,hprime_zz, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        ibool,rhostore,SIMULATION_TYPE,GRAVITY
  use specfem_par_movie,only:SAVE_DISPLACEMENT
  implicit none

  integer,intent(in) :: iglob,ispec
  real(kind=CUSTOM_REAL),dimension(NDIM),intent(out):: val_vector

  ! local parameters
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ):: val_element
  integer :: i,j,k

  ! returns first element encountered for iglob index
  if( ELASTIC_SIMULATION ) then
    if( ispec_is_elastic(ispec) ) then
      if(SAVE_DISPLACEMENT) then
        if( SIMULATION_TYPE == 3 ) then
          ! to display re-constructed wavefield
          !val_vector(:) = b_displ(:,iglob)
          ! to display adjoint wavefield
          val_vector(:) = displ(:,iglob)
        else
          val_vector(:) = displ(:,iglob)
        endif
      else
        if( SIMULATION_TYPE == 3 ) then
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

  if( ACOUSTIC_SIMULATION ) then
    if( ispec_is_acoustic(ispec) ) then
      if(SAVE_DISPLACEMENT) then
        if( SIMULATION_TYPE == 3 ) then
          ! displacement vector from backward potential
          call compute_gradient(ispec,NSPEC_AB,NGLOB_AB, &
                          b_potential_acoustic, val_element,&
                          hprime_xx,hprime_yy,hprime_zz, &
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                          ibool,rhostore,GRAVITY)
        else
          ! displacement vector
          call compute_gradient(ispec,NSPEC_AB,NGLOB_AB, &
                          potential_acoustic, val_element,&
                          hprime_xx,hprime_yy,hprime_zz, &
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                          ibool,rhostore,GRAVITY)
        endif
      else
        if( SIMULATION_TYPE == 3 ) then
          ! velocity vector for backward/reconstructed wavefield
          call compute_gradient(ispec,NSPEC_AB,NGLOB_AB, &
                          b_potential_dot_acoustic, val_element,&
                          hprime_xx,hprime_yy,hprime_zz, &
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                          ibool,rhostore,GRAVITY)
        else
          ! velocity vector
          call compute_gradient(ispec,NSPEC_AB,NGLOB_AB, &
                          potential_dot_acoustic, val_element,&
                          hprime_xx,hprime_yy,hprime_zz, &
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                          ibool,rhostore,GRAVITY)
        endif
      endif

      ! returns corresponding iglob velocity entry
      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX
            if( ibool(i,j,k,ispec) == iglob ) then
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
