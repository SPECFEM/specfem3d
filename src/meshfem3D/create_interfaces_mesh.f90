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

  subroutine create_interfaces_mesh()

  use meshfem_par, only: myrank,INTERFACES_FILE, &
    number_of_interfaces, &
    SUPPRESS_UTM_PROJECTION, &
    UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK, &
    max_npx_interface,max_npy_interface, &
    npx_element_steps,npy_element_steps, &
    number_of_layers,ner_layer,iproc_xi_current,iproc_eta_current,NPROC_XI,NPROC_ETA, &
    xgrid,ygrid,zgrid

  use constants, only: IMAIN,IIN,MF_IN_DATA_FILES,DONT_IGNORE_JUNK,IUTM2LONGLAT,MAX_STRING_LEN,HUGEVAL

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
  double precision :: min_long,max_long
  double precision :: min_lat,max_lat
  double precision :: min_x,max_x
  double precision :: min_y,max_y

  double precision :: min_elevation_all,max_elevation_all
  double precision :: min_long_all,max_long_all
  double precision :: min_lat_all,max_lat_all
  double precision :: min_x_all,max_x_all
  double precision :: min_y_all,max_y_all

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
  character(len=12) :: str_unit
  character(len=12),parameter :: str_unit_m = "(m)",str_unit_deg = "(deg)"

  logical :: SUPPRESS_UTM_PROJECTION_BOTTOM,SUPPRESS_UTM_PROJECTION_TOP,SUPPRESS_UTM_PROJECTION_COPY

  ! conversion factor:
  ! at equator earth circumference 40,075,161.2 divided by 360.0 degree
  double precision, parameter :: DEGREE_TO_METERS = 111319.8922222222d0

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
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1281')
  if (ier /= 0) stop 'Error allocating array interface_bottom'
  allocate(interface_top(max_npx_interface,max_npy_interface),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1282')
  if (ier /= 0) stop 'Error allocating array interface_top'

  ! read number of interfaces
  call read_value_integer_mesh(IIN,DONT_IGNORE_JUNK,number_of_interfaces,'NINTERFACES', ier)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'number of interfaces: ',number_of_interfaces
      write(IMAIN,*)
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

  ! collects mesh statistics
  call min_all_dp(UTM_X_MIN,min_x_all)
  call max_all_dp(UTM_X_MAX,max_x_all)
  call min_all_dp(UTM_Y_MIN,min_y_all)
  call max_all_dp(UTM_Y_MAX,max_y_all)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'mesh:'
    write(IMAIN,*) '  origin UTM minimum x/y        (m) = ',sngl(min_x_all),sngl(min_y_all)
    if (.not. SUPPRESS_UTM_PROJECTION) then
      ! project x and y in UTM back to long/lat since topo file is in long/lat
      call utm_geo(long,lat,min_x_all,min_y_all,IUTM2LONGLAT)
      write(IMAIN,*) '                     lat/lon  (deg) = ',sngl(lat),sngl(long)
    endif
    write(IMAIN,*) '  origin UTM maximum x/y        (m) = ',sngl(max_x_all),sngl(max_y_all)
    if (.not. SUPPRESS_UTM_PROJECTION) then
      call utm_geo(long,lat,max_x_all,max_y_all,IUTM2LONGLAT)
      write(IMAIN,*) '                     lat/lon  (deg) = ',sngl(lat),sngl(long)
    endif
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  SUPPRESS_UTM_PROJECTION_COPY = SUPPRESS_UTM_PROJECTION

  ! loop on all the layers
  ! note: number of layers and number of interfaces are equal
  do ilayer = 1,number_of_layers

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) 'reading interface ',ilayer
      call flush_IMAIN()
    endif

    ! read top interface
    call read_interface_parameters(IIN,SUPPRESS_UTM_PROJECTION_TOP,interface_top_file, &
                                   npx_interface_top,npy_interface_top, &
                                   orig_x_interface_top,orig_y_interface_top, &
                                   spacing_x_interface_top,spacing_y_interface_top,ier)
    if (ier /= 0) stop 'Error reading interface parameters'

    ! user output
    if (myrank == 0) then
      if (SUPPRESS_UTM_PROJECTION_TOP) then
        str_unit = str_unit_m
      else
        str_unit = str_unit_deg
      endif
      write(IMAIN,*) '  interface file   : ',trim(interface_top_file)
      write(IMAIN,*)
      write(IMAIN,*) '  number of points x/y = ',npx_interface_top,npy_interface_top
      write(IMAIN,*) '  origin x/y     ',trim(str_unit),' = ', sngl(orig_x_interface_top),sngl(orig_y_interface_top)
      write(IMAIN,*) '  spacing x/y    ',trim(str_unit),' = ',sngl(spacing_x_interface_top),sngl(spacing_y_interface_top)
      if (.not. SUPPRESS_UTM_PROJECTION_TOP) then
        write(IMAIN,*) '                   ',trim(str_unit_m),' = ', &
             sngl(spacing_x_interface_top*DEGREE_TO_METERS),sngl(spacing_y_interface_top*DEGREE_TO_METERS)
      endif
      write(IMAIN,*)
      write(IMAIN,*) '  dimension x-direction ',trim(str_unit),' = ',sngl(orig_x_interface_top),'/', &
                            sngl(orig_x_interface_top + (npx_interface_top-1)*spacing_x_interface_top)
      write(IMAIN,*) '  dimension y-direction ',trim(str_unit),' = ',sngl(orig_y_interface_top),'/', &
                            sngl(orig_y_interface_top + (npy_interface_top-1)*spacing_y_interface_top)
      write(IMAIN,*)
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

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  total number of file points = ',icount,' should be ',npx_interface_top * npy_interface_top
      if (icount == npx_interface_top * npy_interface_top) then
        write(IMAIN,*) '  this point total is okay'
      else
        write(IMAIN,*) '  this point total is incorrect'
      endif
      write(IMAIN,*)
      call flush_IMAIN()
    endif

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
        ! reading without trailing/triming...
        read(45,*,iostat=ier) elevation
        !in case data file would have comment lines...
        !call read_value_dble_precision_mesh(45,DONT_IGNORE_JUNK,elevation,'Z_INTERFACE_TOP',ier)
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
      write(IMAIN,*) '  original elevation min/max             = ',sngl(min_elevation),sngl(max_elevation)
      call flush_IMAIN()
    endif

    ! compute the offset of this layer in terms of number of spectral elements below along Z
    if (ilayer > 1) then
       ioffset = sum(ner_layer(1:ilayer-1))
    else
       ioffset = 0
    endif

    !--- definition of the mesh
    ! statistics
    max_elevation = - HUGEVAL
    min_elevation = HUGEVAL
    min_long = HUGEVAL
    max_long = - HUGEVAL
    min_lat = HUGEVAL
    max_lat = - HUGEVAL
    min_x = HUGEVAL
    max_x = - HUGEVAL
    min_y = HUGEVAL
    max_y = - HUGEVAL

    do iy = 0,npy_element_steps
      do ix = 0,npx_element_steps

        !   define the mesh points on the top and the bottom
        xin = dble(ix)/dble(npx_element_steps)
        x_current = UTM_X_MIN + (dble(iproc_xi_current)+xin)*(UTM_X_MAX-UTM_X_MIN)/dble(NPROC_XI)

        etan = dble(iy)/dble(npy_element_steps)
        y_current = UTM_Y_MIN + (dble(iproc_eta_current)+etan)*(UTM_Y_MAX-UTM_Y_MIN)/dble(NPROC_ETA)

        ! get bottom interface value
        ! project x and y in UTM back to long/lat since topo file is in long/lat
        SUPPRESS_UTM_PROJECTION = SUPPRESS_UTM_PROJECTION_BOTTOM
        call utm_geo(long,lat,x_current,y_current,IUTM2LONGLAT)

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
        if (ratio_xi < 0.d0) ratio_xi = 0.d0
        if (ratio_xi > 1.d0) ratio_xi = 1.d0
        if (ratio_eta < 0.d0) ratio_eta = 0.d0
        if (ratio_eta > 1.d0) ratio_eta = 1.d0

        ! interpolate elevation at current point
        z_interface_bottom = &
              interface_bottom(icornerlong,icornerlat)*(1.d0-ratio_xi)*(1.d0-ratio_eta) + &
              interface_bottom(icornerlong+1,icornerlat)*ratio_xi*(1.d0-ratio_eta) + &
              interface_bottom(icornerlong+1,icornerlat+1)*ratio_xi*ratio_eta + &
              interface_bottom(icornerlong,icornerlat+1)*(1.d0-ratio_xi)*ratio_eta

        ! get top interface value
        ! project x and y in UTM back to long/lat since topo file is in long/lat
        SUPPRESS_UTM_PROJECTION = SUPPRESS_UTM_PROJECTION_TOP
        call utm_geo(long,lat,x_current,y_current,IUTM2LONGLAT)

        ! debug
        !if (long < 138.d0) print *,'long:',long,lat,x_current,y_current,'rank',myrank,ix,iy,'xi/eta',xin,etan

        ! statstics
        if (long < min_long) min_long = long
        if (long > max_long) max_long = long
        if (lat < min_lat) min_lat = lat
        if (lat > max_lat) max_lat = lat
        if (x_current < min_x) min_x = x_current
        if (x_current > max_x) max_x = x_current
        if (y_current < min_y) min_y = y_current
        if (y_current > max_y) max_y = y_current

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
        if (ratio_xi < 0.d0) ratio_xi = 0.d0
        if (ratio_xi > 1.d0) ratio_xi = 1.d0
        if (ratio_eta < 0.d0) ratio_eta = 0.d0
        if (ratio_eta > 1.d0) ratio_eta = 1.d0

        ! interpolate elevation at current point
        z_interface_top = &
             interface_top(icornerlong,icornerlat)*(1.d0-ratio_xi)*(1.d0-ratio_eta) + &
             interface_top(icornerlong+1,icornerlat)*ratio_xi*(1.d0-ratio_eta) + &
             interface_top(icornerlong+1,icornerlat+1)*ratio_xi*ratio_eta + &
             interface_top(icornerlong,icornerlat+1)*(1.d0-ratio_xi)*ratio_eta

        ! statstics
        if (z_interface_top < min_elevation) min_elevation = z_interface_top
        if (z_interface_top > max_elevation) max_elevation = z_interface_top

        ! stores location of grid points
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

    ! collects statistics
    call min_all_dp(min_elevation,min_elevation_all)
    call max_all_dp(max_elevation,max_elevation_all)
    call min_all_dp(min_long,min_long_all)
    call max_all_dp(max_long,max_long_all)
    call min_all_dp(min_lat,min_lat_all)
    call max_all_dp(max_lat,max_lat_all)
    call min_all_dp(min_x,min_x_all)
    call max_all_dp(max_x,max_x_all)
    call min_all_dp(min_y,min_y_all)
    call max_all_dp(max_y,max_y_all)

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  interpolated mesh elevation min/max    = ',sngl(min_elevation_all),sngl(max_elevation_all)
      write(IMAIN,*)
      if (.not. SUPPRESS_UTM_PROJECTION_TOP) then
        write(IMAIN,*) '  interpolated mesh longitude min/max ',trim(str_unit_deg),' = ',sngl(min_long_all),'/',sngl(max_long_all)
        write(IMAIN,*) '  interpolated mesh latitude  min/max ',trim(str_unit_deg),' = ',sngl(min_lat_all),'/',sngl(max_lat_all)
        write(IMAIN,*)
      endif
      write(IMAIN,*) '  interpolated mesh UTM minimum x/y ',trim(str_unit_m),' = ',sngl(min_x_all),sngl(min_y_all)
      write(IMAIN,*) '  interpolated mesh UTM maximum x/y ',trim(str_unit_m),' = ',sngl(max_x_all),sngl(max_y_all)
      write(IMAIN,*)
      call flush_IMAIN()
    endif

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

  SUPPRESS_UTM_PROJECTION = SUPPRESS_UTM_PROJECTION_COPY

  ! free memory
  deallocate(interface_top,interface_bottom)

  ! make sure everybody is synchronized
  call synchronize_all()

  end subroutine create_interfaces_mesh

!
!----------------------------------------------------------------------------------
!

  subroutine get_interfaces_mesh_count()

  use meshfem_par, only: INTERFACES_FILE, &
    number_of_interfaces,max_npx_interface,max_npy_interface, &
    number_of_layers,ner_layer

  use constants, only: IMAIN,IIN,MF_IN_DATA_FILES,DONT_IGNORE_JUNK,MAX_STRING_LEN,myrank

  implicit none

  ! local parameters
  integer :: npx_interface,npy_interface
  integer :: interface_current,ilayer
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
                                   npx_interface,npy_interface, &
                                   orig_x_dummy,orig_y_dummy, &
                                   spacing_x_dummy,spacing_y_dummy,ier)
    if (ier /= 0) then
      print *,'Error reading interface parameters: interface ',interface_current
      stop 'Error reading interface parameters for interfaces'
    endif

    max_npx_interface = max(npx_interface,max_npx_interface)
    max_npy_interface = max(npy_interface,max_npy_interface)

    if ((max_npx_interface < 2) .or. (max_npy_interface < 2)) then
      print *,'Error interface ',interface_current,': has not enough interface points (minimum is 2x2)'
      stop 'Error not enough interface points (minimum is 2x2)'
    endif
  enddo

  if (myrank == 0) then
    write(IMAIN,*) '  maximum interface points x/y = ',max_npx_interface,max_npy_interface
    call flush_IMAIN()
  endif

  ! define number of layers
  number_of_layers = number_of_interfaces

  allocate(ner_layer(number_of_layers),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1283')
  if (ier /= 0) stop 'Error allocating array ner_layer'
  ner_layer(:) = 0

  ! loop on all the layers
  do ilayer = 1,number_of_layers

    ! read number of spectral elements in vertical direction in this layer
    call read_value_integer_mesh(IIN,DONT_IGNORE_JUNK,ner_layer(ilayer),'NER_LAYER', ier)
    if (ier /= 0) stop 'Error reading interface parameter for NER_LAYER'

    ! checks
    if (ner_layer(ilayer) < 1) then
      print *,'Error invalid layering number ',ner_layer(ilayer),' for layer ',ilayer
      print *,'Please use a minimum element layer number of 1 in interface file'
      stop 'Error not enough spectral elements along Z in layer (minimum is 1)'
    endif
  enddo

  close(IIN)

  if (myrank == 0) then
    write(IMAIN,*) '  interfaces done'
    call flush_IMAIN()
  endif

  end subroutine get_interfaces_mesh_count
