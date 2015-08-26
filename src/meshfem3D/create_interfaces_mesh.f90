!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
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

  subroutine create_interfaces_mesh()

  use meshfem3D_par, only: myrank,INTERFACES_FILE, &
    number_of_interfaces, &
    SUPPRESS_UTM_PROJECTION, &
    UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,UTM_PROJECTION_ZONE,Z_DEPTH_BLOCK, &
    max_npx_interface,max_npy_interface,npx,npy, &
    number_of_layers,ner_layer,iproc_xi_current,iproc_eta_current,NPROC_XI,NPROC_ETA, &
    xgrid,ygrid,zgrid

  use constants,only: IMAIN,IIN,MF_IN_DATA_FILES,DONT_IGNORE_JUNK,IUTM2LONGLAT,MAX_STRING_LEN,HUGEVAL

  implicit none

  ! local parameters
  integer :: ilayer
  integer :: ix,iy,ir
  integer :: ier,icount

  ! auxiliary variables to generate interface
  double precision :: xin,etan
  double precision :: x_current,y_current

  ! use integer array to store topography values
  integer :: icornerlat,icornerlong
  double precision :: lat,long
  double precision :: long_corner,lat_corner,ratio_xi,ratio_eta
  double precision, dimension(:,:),allocatable :: interface_bottom,interface_top
  double precision :: elevation
  double precision :: min_elevation,max_elevation

! to compute the coordinate transformation
  integer :: ioffset
  double precision :: gamma

  double precision :: z_interface_bottom,z_interface_top

  integer :: npx_interface_bottom,npy_interface_bottom
  integer :: npx_interface_top,npy_interface_top

  double precision :: orig_x_interface_bottom,orig_y_interface_bottom
  double precision :: orig_x_interface_top,orig_y_interface_top

  double precision :: spacing_x_interface_bottom,spacing_y_interface_bottom
  double precision :: spacing_x_interface_top,spacing_y_interface_top

  character(len=MAX_STRING_LEN) :: interface_top_file

  logical :: SUPPRESS_UTM_PROJECTION_BOTTOM,SUPPRESS_UTM_PROJECTION_TOP

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'Reading interface data from file ',trim(MF_IN_DATA_FILES)//trim(INTERFACES_FILE)
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  open(unit=IIN,file=trim(MF_IN_DATA_FILES)//trim(INTERFACES_FILE),status='old',iostat=ier)
  if (ier /= 0) stop 'Error opening interfaces file'

  ! allocates interface arrays
  allocate(interface_bottom(max_npx_interface,max_npy_interface),stat=ier)
  if (ier /= 0) stop 'Error allocating array interface_bottom'
  allocate(interface_top(max_npx_interface,max_npy_interface),stat=ier)
  if (ier /= 0) stop 'Error allocating array interface_top'

  ! read number of interfaces
  call read_value_integer_mesh(IIN,DONT_IGNORE_JUNK,number_of_interfaces,'NINTERFACES', ier)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'number of interfaces: ',number_of_interfaces
    call flush_IMAIN()
  endif

  ! initializes
  SUPPRESS_UTM_PROJECTION_BOTTOM = SUPPRESS_UTM_PROJECTION
  npx_interface_bottom = 2
  npy_interface_bottom = 2

  orig_x_interface_bottom = UTM_X_MIN
  orig_y_interface_bottom = UTM_Y_MIN

  spacing_x_interface_bottom = UTM_X_MAX - UTM_X_MIN
  spacing_y_interface_bottom = UTM_Y_MAX - UTM_Y_MIN

  interface_bottom(:,:) = - dabs(Z_DEPTH_BLOCK)

  ! loop on all the layers
  ! note: number of layers and number of interfaces are equal
  do ilayer = 1,number_of_layers

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  reading interface ',ilayer
      call flush_IMAIN()
    endif

    ! read top interface
    call read_interface_parameters(IIN,SUPPRESS_UTM_PROJECTION_TOP,interface_top_file,&
                                   npx_interface_top,npy_interface_top,&
                                   orig_x_interface_top,orig_y_interface_top,&
                                   spacing_x_interface_top,spacing_y_interface_top,ier)
    if (ier /= 0) stop 'Error reading interface parameters'

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  interface file   : ',trim(interface_top_file)
      write(IMAIN,*)
      write(IMAIN,*) '  number of points x/y = ',npx_interface_top,npy_interface_top
      write(IMAIN,*) '  origin x/y           = ',sngl(orig_x_interface_top),sngl(orig_y_interface_top)
      write(IMAIN,*) '  spacing x/y          = ',sngl(spacing_x_interface_top),sngl(spacing_y_interface_top)
      call flush_IMAIN()
    endif

    ! loop on all the points describing this interface
    open(unit=45,file=trim(MF_IN_DATA_FILES)//trim(interface_top_file),status='old',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(MF_IN_DATA_FILES)//trim(interface_top_file)
      stop 'Error opening interface file'
    endif

    ! counts number of lines, i.e. elevation point entries
    icount = 0
    do while(ier == 0)
      read(45,*,iostat=ier) elevation
      if (ier == 0) icount = icount + 1
    enddo
    rewind(45)

    ! checks number of points npoints_interface_top = npx_interface_top * npy_interface
    if (icount /= npx_interface_top * npy_interface_top) then
      print *,'Error number of interface points ',icount,' should be ',npx_interface_top * npy_interface_top
      print *,'Please check interface definition in file: ',trim(MF_IN_DATA_FILES)//trim(interface_top_file)
      stop 'Error invalid number of interface file points'
    endif

    ! statistics
    max_elevation = - HUGEVAL
    min_elevation = HUGEVAL

    ! reads in interface points
    do iy = 1,npy_interface_top
      do ix = 1,npx_interface_top
        call read_value_dble_precision_mesh(45,DONT_IGNORE_JUNK,elevation,'Z_INTERFACE_TOP',ier)
        if (ier /= 0) stop 'Error reading interface value'

        ! stores values for interpolation
        interface_top(ix,iy) = elevation

        ! statstics
        if (elevation < min_elevation) min_elevation = elevation
        if (elevation > max_elevation) max_elevation = elevation
      enddo
    enddo
    close(45)

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  elevation min/max    = ',sngl(min_elevation),sngl(max_elevation)
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! compute the offset of this layer in terms of number of spectral elements below along Z
    if (ilayer > 1) then
       ioffset = sum(ner_layer(1:ilayer-1))
    else
       ioffset = 0
    endif

    !--- definition of the mesh

    do iy = 0,npy
      do ix = 0,npx

        !   define the mesh points on the top and the bottom
        xin = dble(ix)/dble(npx)
        x_current = UTM_X_MIN + (dble(iproc_xi_current)+xin)*(UTM_X_MAX-UTM_X_MIN)/dble(NPROC_XI)

        etan = dble(iy)/dble(npy)
        y_current = UTM_Y_MIN + (dble(iproc_eta_current)+etan)*(UTM_Y_MAX-UTM_Y_MIN)/dble(NPROC_ETA)

        ! get bottom interface value
        ! project x and y in UTM back to long/lat since topo file is in long/lat
        call utm_geo(long,lat,x_current,y_current,UTM_PROJECTION_ZONE,IUTM2LONGLAT,SUPPRESS_UTM_PROJECTION_BOTTOM)

        ! get coordinate of corner in bathy/topo model
        icornerlong = int((long - orig_x_interface_bottom) / spacing_x_interface_bottom) + 1
        icornerlat = int((lat - orig_y_interface_bottom) / spacing_y_interface_bottom) + 1

        ! avoid edge effects and extend with identical point if outside model
        if (icornerlong < 1) icornerlong = 1
        if (icornerlong > npx_interface_bottom-1) icornerlong = npx_interface_bottom-1
        if (icornerlat < 1) icornerlat = 1
        if (icornerlat > npy_interface_bottom-1) icornerlat = npy_interface_bottom-1

        ! compute coordinates of corner
        long_corner = orig_x_interface_bottom + (icornerlong-1)*spacing_x_interface_bottom
        lat_corner = orig_y_interface_bottom + (icornerlat-1)*spacing_y_interface_bottom

        ! compute ratio for interpolation
        ratio_xi = (long - long_corner) / spacing_x_interface_bottom
        ratio_eta = (lat - lat_corner) / spacing_y_interface_bottom

        ! avoid edge effects
        if (ratio_xi < 0.) ratio_xi = 0.
        if (ratio_xi > 1.) ratio_xi = 1.
        if (ratio_eta < 0.) ratio_eta = 0.
        if (ratio_eta > 1.) ratio_eta = 1.

        ! interpolate elevation at current point
        z_interface_bottom = &
              interface_bottom(icornerlong,icornerlat)*(1.-ratio_xi)*(1.-ratio_eta) + &
              interface_bottom(icornerlong+1,icornerlat)*ratio_xi*(1.-ratio_eta) + &
              interface_bottom(icornerlong+1,icornerlat+1)*ratio_xi*ratio_eta + &
              interface_bottom(icornerlong,icornerlat+1)*(1.-ratio_xi)*ratio_eta

        ! get top interface value
        ! project x and y in UTM back to long/lat since topo file is in long/lat
        call utm_geo(long,lat,x_current,y_current,UTM_PROJECTION_ZONE,IUTM2LONGLAT,SUPPRESS_UTM_PROJECTION_TOP)

        ! get coordinate of corner in bathy/topo model
        icornerlong = int((long - orig_x_interface_top) / spacing_x_interface_top) + 1
        icornerlat = int((lat - orig_y_interface_top) / spacing_y_interface_top) + 1

        ! avoid edge effects and extend with identical point if outside model
        if (icornerlong < 1) icornerlong = 1
        if (icornerlong > npx_interface_top-1) icornerlong = npx_interface_top-1
        if (icornerlat < 1) icornerlat = 1
        if (icornerlat > npy_interface_top-1) icornerlat = npy_interface_top-1

        ! compute coordinates of corner
        long_corner = orig_x_interface_top + (icornerlong-1)*spacing_x_interface_top
        lat_corner = orig_y_interface_top + (icornerlat-1)*spacing_y_interface_top

        ! compute ratio for interpolation
        ratio_xi = (long - long_corner) / spacing_x_interface_top
        ratio_eta = (lat - lat_corner) / spacing_y_interface_top

        ! avoid edge effects
        if (ratio_xi < 0.) ratio_xi = 0.
        if (ratio_xi > 1.) ratio_xi = 1.
        if (ratio_eta < 0.) ratio_eta = 0.
        if (ratio_eta > 1.) ratio_eta = 1.

        ! interpolate elevation at current point
        z_interface_top = &
             interface_top(icornerlong,icornerlat)*(1.-ratio_xi)*(1.-ratio_eta) + &
             interface_top(icornerlong+1,icornerlat)*ratio_xi*(1.-ratio_eta) + &
             interface_top(icornerlong+1,icornerlat+1)*ratio_xi*ratio_eta + &
             interface_top(icornerlong,icornerlat+1)*(1.-ratio_xi)*ratio_eta

        do ir = 0,ner_layer(ilayer)
          ! linear interpolation between bottom and top
          gamma = dble(ir) / dble(ner_layer(ilayer))

          ! coordinates of the grid points
          xgrid(ir + ioffset,ix,iy) = x_current
          ygrid(ir + ioffset,ix,iy) = y_current
          zgrid(ir + ioffset,ix,iy) = gamma*z_interface_top + (1.d0 - gamma)*z_interface_bottom
        enddo

      enddo
    enddo

    ! the top interface becomes the bottom interface before switching to the next layer
    SUPPRESS_UTM_PROJECTION_BOTTOM = SUPPRESS_UTM_PROJECTION_TOP

    npx_interface_bottom = npx_interface_top
    npy_interface_bottom = npy_interface_top

    orig_x_interface_bottom = orig_x_interface_top
    orig_y_interface_bottom = orig_y_interface_top

    spacing_x_interface_bottom = spacing_x_interface_top
    spacing_y_interface_bottom = spacing_y_interface_top

    interface_bottom(:,:) = interface_top(:,:)
  enddo

  ! free memory
  deallocate(interface_top,interface_bottom)

  ! make sure everybody is synchronized
  call synchronize_all()

  end subroutine create_interfaces_mesh

!
!----------------------------------------------------------------------------------
!

  subroutine get_interfaces_mesh_count()

  use meshfem3D_par, only: myrank,INTERFACES_FILE, &
    number_of_interfaces,max_npx_interface,max_npy_interface

  use constants,only: IMAIN,IIN,MF_IN_DATA_FILES,DONT_IGNORE_JUNK,IUTM2LONGLAT,MAX_STRING_LEN

  implicit none

  ! local parameters
  integer :: npx_interface,npy_interface
  integer :: interface_current
  integer :: ier

  ! dummy values
  double precision :: orig_x_dummy,orig_y_dummy
  double precision :: spacing_x_dummy,spacing_y_dummy
  character(len=MAX_STRING_LEN) :: dummy_file
  logical :: SUPPRESS_UTM_PROJECTION_DUMMY

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Reading interface data from file ',trim(MF_IN_DATA_FILES)//trim(INTERFACES_FILE)
    call flush_IMAIN()
  endif

  ! opens interfaces file
  open(unit=IIN,file=trim(MF_IN_DATA_FILES)//trim(INTERFACES_FILE),status='old',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening interface file: ',trim(MF_IN_DATA_FILES)//trim(INTERFACES_FILE)
    stop 'Error opening interface file'
  endif

  ! initializes maximum count
  max_npx_interface  = -1
  max_npy_interface  = -1

  ! read number of interfaces
  call read_value_integer_mesh(IIN,DONT_IGNORE_JUNK,number_of_interfaces,'NINTERFACES', ier)
  if (ier /= 0) stop 'Error reading interface parameter for NINTERFACES'

  if (number_of_interfaces < 1) stop 'Error not enough interfaces (minimum is 1, for topography)'

  ! loop on all the interfaces
  do interface_current = 1,number_of_interfaces
    call read_interface_parameters(IIN,SUPPRESS_UTM_PROJECTION_DUMMY,dummy_file, &
                                   npx_interface,npy_interface,&
                                   orig_x_dummy,orig_y_dummy,&
                                   spacing_x_dummy,spacing_y_dummy,ier)
    if (ier /= 0) then
      print *,'Error reading interface parameters: interface ',interface_current
      stop 'Error reading interface parameters for interfaces'
    endif

    max_npx_interface = max(npx_interface,max_npx_interface)
    max_npy_interface = max(npy_interface,max_npy_interface)

    if ((max_npx_interface < 2) .or.(max_npy_interface < 2)) then
      print *,'Error interface ',interface_current,': has not enough interface points (minimum is 2x2)'
      stop 'Error not enough interface points (minimum is 2x2)'
    endif
  enddo

  if (myrank == 0) then
    write(IMAIN,*) 'maximum interface points x/y = ',max_npx_interface,max_npy_interface
    call flush_IMAIN()
  endif

  end subroutine get_interfaces_mesh_count
