!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

!----
!---- locate_receivers finds the correct position of the receivers
!----
  subroutine locate_receivers(ibool,myrank,NSPEC_AB,NGLOB_AB, &
                 xstore,ystore,zstore,xigll,yigll,zigll,rec_filename, &
                 nrec,islice_selected_rec,ispec_selected_rec, &
                 xi_receiver,eta_receiver,gamma_receiver,station_name,network_name,nu, &
                 NPROC,utm_x_source,utm_y_source, &
                 TOPOGRAPHY,itopo_bathy_basin,UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION, &
                 NX_TOPO,NY_TOPO,ORIG_LAT_TOPO,ORIG_LONG_TOPO,DEGREES_PER_CELL_TOPO)

  implicit none

#ifdef USE_MPI
! standard include of the MPI library
  include 'mpif.h'
#endif

  include "constants.h"
#ifdef USE_MPI
  include "precision.h"
#endif

  integer :: NPROC,UTM_PROJECTION_ZONE,NX_TOPO,NY_TOPO

  logical :: TOPOGRAPHY,SUPPRESS_UTM_PROJECTION

  integer :: nrec,myrank

  integer :: NSPEC_AB,NGLOB_AB

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

  double precision :: ORIG_LAT_TOPO,ORIG_LONG_TOPO,DEGREES_PER_CELL_TOPO

! arrays containing coordinates of the points
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: xstore,ystore,zstore

! use integer array to store topography values
  integer :: itopo_bathy_basin(NX_TOPO,NY_TOPO)
  double precision :: long_corner,lat_corner,ratio_xi,ratio_eta

  integer, allocatable, dimension(:) :: ix_initial_guess,iy_initial_guess,iz_initial_guess

  integer :: iprocloop
  integer :: nrec_dummy
  integer :: ios

  double precision, allocatable, dimension(:) :: x_target,y_target,z_target
  double precision, allocatable, dimension(:) :: horiz_dist,elevation
  double precision, allocatable, dimension(:) :: x_found,y_found,z_found
  double precision, allocatable, dimension(:,:) :: x_found_all,y_found_all,z_found_all

  integer :: irec
  integer :: i,j,k,ispec,iglob
#ifdef USE_MPI
  integer :: ier
#endif

  integer :: icornerlong,icornerlat
  double precision :: utm_x_source,utm_y_source
  double precision :: dist
  double precision :: xi,eta,gamma,dx,dy,dz,dxi,deta,dgamma

! Gauss-Lobatto-Legendre points of integration
  double precision :: xigll(NGLLX)
  double precision :: yigll(NGLLY)
  double precision :: zigll(NGLLZ)

! input receiver file name
  character(len=*) :: rec_filename

! topology of the control points of the surface element
  integer :: iax,iay,iaz
  integer :: iaddx(NGNOD),iaddy(NGNOD),iaddz(NGNOD)

! coordinates of the control points of the surface element
  double precision :: xelm(NGNOD),yelm(NGNOD),zelm(NGNOD)

  integer :: iter_loop,ispec_iterate

  integer :: ia
  double precision :: x,y,z
  double precision :: xix,xiy,xiz
  double precision :: etax,etay,etaz
  double precision :: gammax,gammay,gammaz

! timer MPI
  double precision :: time_start,tCPU

! use dynamic allocation
  double precision, dimension(:), allocatable :: final_distance
  double precision, dimension(:,:), allocatable :: final_distance_all
  double precision distmin,final_distance_max

! receiver information
! timing information for the stations
! station information for writing the seismograms

  integer, dimension(nrec) :: islice_selected_rec,ispec_selected_rec
  double precision, dimension(nrec) :: xi_receiver,eta_receiver,gamma_receiver
  double precision, dimension(3,3,nrec) :: nu
  character(len=MAX_LENGTH_STATION_NAME), dimension(nrec) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(nrec) :: network_name

  integer, allocatable, dimension(:,:) :: ispec_selected_rec_all
  double precision, allocatable, dimension(:) :: stlat,stlon,stele,stbur,stutm_x,stutm_y
  double precision, allocatable, dimension(:,:) :: xi_receiver_all,eta_receiver_all,gamma_receiver_all

  character(len=150) :: OUTPUT_FILES

! **************

! get MPI starting time
#ifdef USE_MPI
  time_start = MPI_WTIME()
#else
  time_start = 0.d0
#endif

  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '********************'
    write(IMAIN,*) ' locating receivers'
    write(IMAIN,*) '********************'
    write(IMAIN,*)
  endif

! define topology of the control element
  call usual_hex_nodes(iaddx,iaddy,iaddz)

  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '*****************************************************************'
    write(IMAIN,'(1x,a,a,a)') 'reading receiver information from ', trim(rec_filename), ' file'
    write(IMAIN,*) '*****************************************************************'
  endif

! get number of stations from receiver file
  open(unit=1,file=trim(rec_filename),status='old',action='read',iostat=ios)
  if (ios /= 0) call exit_mpi(myrank,'error opening file '//trim(rec_filename))
  read(1,*) nrec_dummy

  if (nrec_dummy /= nrec) call exit_MPI(myrank,'problem with number of receivers')

! allocate memory for arrays using number of stations
  allocate(stlat(nrec))
  allocate(stlon(nrec))
  allocate(stele(nrec))
  allocate(stbur(nrec))
  allocate(stutm_x(nrec))
  allocate(stutm_y(nrec))
  allocate(horiz_dist(nrec))
  allocate(elevation(nrec))

  allocate(ix_initial_guess(nrec))
  allocate(iy_initial_guess(nrec))
  allocate(iz_initial_guess(nrec))
  allocate(x_target(nrec))
  allocate(y_target(nrec))
  allocate(z_target(nrec))
  allocate(x_found(nrec))
  allocate(y_found(nrec))
  allocate(z_found(nrec))
  allocate(final_distance(nrec))

  allocate(ispec_selected_rec_all(nrec,0:NPROC-1))
  allocate(xi_receiver_all(nrec,0:NPROC-1))
  allocate(eta_receiver_all(nrec,0:NPROC-1))
  allocate(gamma_receiver_all(nrec,0:NPROC-1))
  allocate(x_found_all(nrec,0:NPROC-1))
  allocate(y_found_all(nrec,0:NPROC-1))
  allocate(z_found_all(nrec,0:NPROC-1))
  allocate(final_distance_all(nrec,0:NPROC-1))

! loop on all the stations
  do irec=1,nrec

! set distance to huge initial value
  distmin=HUGEVAL

    read(1,*,iostat=ios) station_name(irec),network_name(irec),stlat(irec),stlon(irec),stele(irec),stbur(irec)
    if (ios /= 0) call exit_mpi(myrank, 'Error reading station file '//trim(rec_filename))

! convert station location to UTM
    call utm_geo(stlon(irec),stlat(irec),stutm_x(irec),stutm_y(irec),UTM_PROJECTION_ZONE,ILONGLAT2UTM,SUPPRESS_UTM_PROJECTION)

! compute horizontal distance between source and receiver in km
    horiz_dist(irec) = dsqrt((stutm_y(irec)-utm_y_source)**2 + (stutm_x(irec)-utm_x_source)**2) / 1000.

! print some information about stations
    if (myrank == 0) &
        write(IMAIN,*) 'Station #',irec,': ',station_name(irec)(1:len_trim(station_name(irec))), &
                       '.',network_name(irec)(1:len_trim(network_name(irec))), &
                       '    horizontal distance:  ',sngl(horiz_dist(irec)),' km'

!     get the Cartesian components of n in the model: nu

! orientation consistent with the UTM projection

!     East
      nu(1,1,irec) = 1.d0
      nu(1,2,irec) = 0.d0
      nu(1,3,irec) = 0.d0

!     North
      nu(2,1,irec) = 0.d0
      nu(2,2,irec) = 1.d0
      nu(2,3,irec) = 0.d0

!     Vertical
      nu(3,1,irec) = 0.d0
      nu(3,2,irec) = 0.d0
      nu(3,3,irec) = 1.d0

! compute elevation of topography at the receiver location
! we assume that receivers are always at the surface i.e. not buried
  if (TOPOGRAPHY) then

! get coordinate of corner in bathy/topo model
    icornerlong = int((stlon(irec) - ORIG_LONG_TOPO) / DEGREES_PER_CELL_TOPO) + 1
    icornerlat = int((stlat(irec) - ORIG_LAT_TOPO) / DEGREES_PER_CELL_TOPO) + 1

! avoid edge effects and extend with identical point if outside model
    if (icornerlong < 1) icornerlong = 1
    if (icornerlong > NX_TOPO-1) icornerlong = NX_TOPO-1
    if (icornerlat < 1) icornerlat = 1
    if (icornerlat > NY_TOPO-1) icornerlat = NY_TOPO-1

! compute coordinates of corner
    long_corner = ORIG_LONG_TOPO + (icornerlong-1)*DEGREES_PER_CELL_TOPO
    lat_corner = ORIG_LAT_TOPO + (icornerlat-1)*DEGREES_PER_CELL_TOPO

! compute ratio for interpolation
    ratio_xi = (stlon(irec) - long_corner) / DEGREES_PER_CELL_TOPO
    ratio_eta = (stlat(irec) - lat_corner) / DEGREES_PER_CELL_TOPO

! avoid edge effects
    if (ratio_xi < 0.) ratio_xi = 0.
    if (ratio_xi > 1.) ratio_xi = 1.
    if (ratio_eta < 0.) ratio_eta = 0.
    if (ratio_eta > 1.) ratio_eta = 1.

! interpolate elevation at current point
    elevation = &
      itopo_bathy_basin(icornerlong,icornerlat)*(1.-ratio_xi)*(1.-ratio_eta) + &
      itopo_bathy_basin(icornerlong+1,icornerlat)*ratio_xi*(1.-ratio_eta) + &
      itopo_bathy_basin(icornerlong+1,icornerlat+1)*ratio_xi*ratio_eta + &
      itopo_bathy_basin(icornerlong,icornerlat+1)*(1.-ratio_xi)*ratio_eta

  else
    elevation(irec) = 0.d0
  endif

! compute the Cartesian position of the receiver
      x_target(irec) = stutm_x(irec)
      y_target(irec) = stutm_y(irec)
      z_target(irec) = elevation(irec) - stbur(irec)
      if (myrank == 0) write(IOVTK,*) x_target(irec), y_target(irec), z_target(irec)

! examine top of the elements only (receivers always at the surface)
!      k = NGLLZ

      do ispec=1,NSPEC_AB

! modification by Qinya Liu: idoubling is no longer used because receivers can now be in depth
! (in previous versions, receivers were always assumed to be at the surface)

! loop only on points inside the element
! exclude edges to ensure this point is not shared with other elements
       do k=2,NGLLZ-1
        do j=2,NGLLY-1
          do i=2,NGLLX-1

            iglob = ibool(i,j,k,ispec)
            dist = dsqrt((x_target(irec)-dble(xstore(iglob)))**2 &
                        +(y_target(irec)-dble(ystore(iglob)))**2 &
                        +(z_target(irec)-dble(zstore(iglob)))**2)

!           keep this point if it is closer to the receiver
            if (dist < distmin) then
              distmin = dist
              ispec_selected_rec(irec) = ispec
              ix_initial_guess(irec) = i
              iy_initial_guess(irec) = j
              iz_initial_guess(irec) = k
            endif

          enddo
        enddo
       enddo
!      endif

! end of loop on all the spectral elements in current slice
      enddo

! end of loop on all the stations
  enddo

! close receiver file
  close(1)

! ****************************************
! find the best (xi,eta,gamma) for each receiver
! ****************************************

! loop on all the receivers to iterate in that slice
    do irec = 1,nrec

        ispec_iterate = ispec_selected_rec(irec)

! use initial guess in xi and eta
        xi = xigll(ix_initial_guess(irec))
        eta = yigll(iy_initial_guess(irec))
        gamma = zigll(iz_initial_guess(irec))

! define coordinates of the control points of the element

  do ia=1,NGNOD

    if (iaddx(ia) == 0) then
      iax = 1
    else if (iaddx(ia) == 1) then
      iax = (NGLLX+1)/2
    else if (iaddx(ia) == 2) then
      iax = NGLLX
    else
      call exit_MPI(myrank,'incorrect value of iaddx')
    endif

    if (iaddy(ia) == 0) then
      iay = 1
    else if (iaddy(ia) == 1) then
      iay = (NGLLY+1)/2
    else if (iaddy(ia) == 2) then
      iay = NGLLY
    else
      call exit_MPI(myrank,'incorrect value of iaddy')
    endif

    if (iaddz(ia) == 0) then
      iaz = 1
    else if (iaddz(ia) == 1) then
      iaz = (NGLLZ+1)/2
    else if (iaddz(ia) == 2) then
      iaz = NGLLZ
    else
      call exit_MPI(myrank,'incorrect value of iaddz')
    endif

    iglob = ibool(iax,iay,iaz,ispec_iterate)
    xelm(ia) = dble(xstore(iglob))
    yelm(ia) = dble(ystore(iglob))
    zelm(ia) = dble(zstore(iglob))

  enddo

! iterate to solve the non linear system
  do iter_loop = 1,NUM_ITER

! impose receiver exactly at the surface
!    gamma = 1.d0

! recompute jacobian for the new point
    call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
           xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

! compute distance to target location
  dx = - (x - x_target(irec))
  dy = - (y - y_target(irec))
  dz = - (z - z_target(irec))

! compute increments
! gamma does not change since we know the receiver is exactly on the surface
  dxi  = xix*dx + xiy*dy + xiz*dz
  deta = etax*dx + etay*dy + etaz*dz
  dgamma = gammax*dx + gammay*dy + gammaz*dz

! update values
  xi = xi + dxi
  eta = eta + deta
  gamma = gamma + dgamma

! impose that we stay in that element
! (useful if user gives a receiver outside the mesh for instance)
! we can go slightly outside the [1,1] segment since with finite elements
! the polynomial solution is defined everywhere
! this can be useful for convergence of itertive scheme with distorted elements
  if (xi > 1.10d0) xi = 1.10d0
  if (xi < -1.10d0) xi = -1.10d0
  if (eta > 1.10d0) eta = 1.10d0
  if (eta < -1.10d0) eta = -1.10d0
  if (gamma > 1.10d0) gamma = 1.10d0
  if (gamma < -1.10d0) gamma = -1.10d0

! end of non linear iterations
  enddo

! impose receiver exactly at the surface after final iteration
!  gamma = 1.d0

! compute final coordinates of point found
  call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
         xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

! store xi,eta and x,y,z of point found
  xi_receiver(irec) = xi
  eta_receiver(irec) = eta
  gamma_receiver(irec) = gamma
  x_found(irec) = x
  y_found(irec) = y
  z_found(irec) = z

! compute final distance between asked and found (converted to km)
  final_distance(irec) = dsqrt((x_target(irec)-x_found(irec))**2 + &
    (y_target(irec)-y_found(irec))**2 + (z_target(irec)-z_found(irec))**2)

    enddo

! synchronize all the processes to make sure all the estimates are available
#ifdef USE_MPI
  call MPI_BARRIER(MPI_COMM_WORLD,ier)
#endif

! for MPI version, gather information from all the nodes
  ispec_selected_rec_all(:,:) = -1
#ifdef USE_MPI
  call MPI_GATHER(ispec_selected_rec,nrec,MPI_INTEGER,ispec_selected_rec_all,nrec,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  call MPI_GATHER(xi_receiver,nrec,MPI_DOUBLE_PRECISION,xi_receiver_all,nrec,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(eta_receiver,nrec,MPI_DOUBLE_PRECISION,eta_receiver_all,nrec,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(gamma_receiver,nrec,MPI_DOUBLE_PRECISION,gamma_receiver_all,nrec,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(final_distance,nrec,MPI_DOUBLE_PRECISION,final_distance_all,nrec,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(x_found,nrec,MPI_DOUBLE_PRECISION,x_found_all,nrec,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(y_found,nrec,MPI_DOUBLE_PRECISION,y_found_all,nrec,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(z_found,nrec,MPI_DOUBLE_PRECISION,z_found_all,nrec,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
#else
  ispec_selected_rec_all(:,0) = ispec_selected_rec

  xi_receiver_all(:,0) = xi_receiver(:)
  eta_receiver_all(:,0) = eta_receiver(:)
  gamma_receiver_all(:,0) = gamma_receiver(:)
  final_distance_all(:,0) = final_distance(:)
  x_found_all(:,0) = x_found(:)
  y_found_all(:,0) = y_found(:)
  z_found_all(:,0) = z_found(:)
#endif

! this is executed by main process only
  if (myrank == 0) then

! check that the gather operation went well
  if (any(ispec_selected_rec_all(:,:) == -1)) call exit_MPI(myrank,'gather operation failed for receivers')

! MPI loop on all the results to determine the best slice
  islice_selected_rec(:) = -1
  do irec = 1,nrec
  distmin = HUGEVAL
  do iprocloop = 0,NPROC-1
    if (final_distance_all(irec,iprocloop) < distmin) then
      distmin = final_distance_all(irec,iprocloop)
      islice_selected_rec(irec) = iprocloop
      ispec_selected_rec(irec) = ispec_selected_rec_all(irec,iprocloop)
      xi_receiver(irec) = xi_receiver_all(irec,iprocloop)
      eta_receiver(irec) = eta_receiver_all(irec,iprocloop)
      gamma_receiver(irec) = gamma_receiver_all(irec,iprocloop)
      x_found(irec) = x_found_all(irec,iprocloop)
      y_found(irec) = y_found_all(irec,iprocloop)
      z_found(irec) = z_found_all(irec,iprocloop)
    endif
  enddo
  final_distance(irec) = distmin
  enddo

  do irec=1,nrec

    write(IMAIN,*)
    write(IMAIN,*) 'station # ',irec,'    ',station_name(irec),network_name(irec)

    if (final_distance(irec) == HUGEVAL) call exit_MPI(myrank,'error locating receiver')

    write(IMAIN,*) '     original latitude: ',sngl(stlat(irec))
    write(IMAIN,*) '    original longitude: ',sngl(stlon(irec))
    write(IMAIN,*) '        original UTM x: ',sngl(stutm_x(irec))
    write(IMAIN,*) '        original UTM y: ',sngl(stutm_y(irec))
    write(IMAIN,*) '   horizontal distance: ',sngl(horiz_dist(irec))
    if (TOPOGRAPHY) write(IMAIN,*) '  topography elevation: ',sngl(elevation(irec))
    write(IMAIN,*) '   target x, y, z: ',sngl(x_target(irec)),sngl(y_target(irec)),sngl(z_target(irec))

    write(IMAIN,*) 'closest estimate found: ',sngl(final_distance(irec)),' m away'
    write(IMAIN,*) ' in slice ',islice_selected_rec(irec),' in element ',ispec_selected_rec(irec)
    write(IMAIN,*) ' at xi,eta,gamma coordinates = ',xi_receiver(irec),eta_receiver(irec),gamma_receiver(irec)

! add warning if estimate is poor
! (usually means receiver outside the mesh given by the user)
    if (final_distance(irec) > 3000.d0) then
      write(IMAIN,*) '*******************************************************'
      write(IMAIN,*) '***** WARNING: receiver location estimate is poor *****'
      write(IMAIN,*) '*******************************************************'
    endif

    write(IMAIN,*)

  enddo

! compute maximal distance for all the receivers
    final_distance_max = maxval(final_distance(:))

! display maximum error for all the receivers
    write(IMAIN,*) 'maximum error in location of all the receivers: ',sngl(final_distance_max),' m'

! add warning if estimate is poor
! (usually means receiver outside the mesh given by the user)
    if (final_distance_max > 1000.d0) then
      write(IMAIN,*)
      write(IMAIN,*) '************************************************************'
      write(IMAIN,*) '************************************************************'
      write(IMAIN,*) '***** WARNING: at least one receiver is poorly located *****'
      write(IMAIN,*) '************************************************************'
      write(IMAIN,*) '************************************************************'
    endif

! get the base pathname for output files
  call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', 'OUTPUT_FILES')

! write the list of stations and associated epicentral distance
  open(unit=27,file=trim(OUTPUT_FILES)//'/output_list_stations.txt',status='unknown')
  do irec=1,nrec
    write(27,*) station_name(irec),'.',network_name(irec),' : ',horiz_dist(irec),' km horizontal distance'
  enddo
  close(27)

! elapsed time since beginning of mesh generation
#ifdef USE_MPI
  tCPU = MPI_WTIME() - time_start
#else
  tCPU = 0.d0
#endif
  write(IMAIN,*)
  write(IMAIN,*) 'Elapsed time for receiver detection in seconds = ',tCPU
  write(IMAIN,*)
  write(IMAIN,*) 'End of receiver detection - done'
  write(IMAIN,*)

  endif    ! end of section executed by main process only

#ifdef USE_MPI
! main process broadcasts the results to all the slices
  call MPI_BCAST(islice_selected_rec,nrec,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(ispec_selected_rec,nrec,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(xi_receiver,nrec,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(eta_receiver,nrec,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(gamma_receiver,nrec,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

! synchronize all the processes to make sure everybody has finished
  call MPI_BARRIER(MPI_COMM_WORLD,ier)
#endif

! deallocate arrays
  deallocate(stlat)
  deallocate(stlon)
  deallocate(stele)
  deallocate(stbur)
  deallocate(stutm_x)
  deallocate(stutm_y)
  deallocate(horiz_dist)
  deallocate(ix_initial_guess)
  deallocate(iy_initial_guess)
  deallocate(iz_initial_guess)
  deallocate(x_target)
  deallocate(y_target)
  deallocate(z_target)
  deallocate(x_found)
  deallocate(y_found)
  deallocate(z_found)
  deallocate(final_distance)
  deallocate(ispec_selected_rec_all)
  deallocate(xi_receiver_all)
  deallocate(eta_receiver_all)
  deallocate(gamma_receiver_all)
  deallocate(x_found_all)
  deallocate(y_found_all)
  deallocate(z_found_all)
  deallocate(final_distance_all)

  end subroutine locate_receivers

!===========================

  subroutine station_filter(myrank,filename,filtered_filename,nfilter, &
      LATITUDE_MIN, LATITUDE_MAX, LONGITUDE_MIN, LONGITUDE_MAX)

  implicit none

  include 'constants.h'

! input
  integer :: myrank
  character(len=*) filename,filtered_filename
  double precision LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX

! output
  integer nfilter

  integer irec, nrec, nrec_filtered, ios

  double precision stlat,stlon,stele,stbur
  character(len=MAX_LENGTH_STATION_NAME) station_name
  character(len=MAX_LENGTH_NETWORK_NAME) network_name

  nrec_filtered = 0

  open(unit=IIN, file=trim(filename), status = 'old', iostat = ios)
  if (ios /= 0) call exit_mpi(myrank, 'No file '//trim(filename)//', exit')
  read(IIN, *) nrec
  do irec = 1, nrec
    read(IIN, *) station_name, network_name, stlat, stlon, stele, stbur
    if (stlat >= LATITUDE_MIN .and. stlat <= LATITUDE_MAX .and. stlon >= LONGITUDE_MIN .and. stlon <= LONGITUDE_MAX) &
          nrec_filtered = nrec_filtered + 1
  enddo
  close(IIN)

  if (myrank == 0) then
    open(unit=IIN,file=trim(filename),status='old',action='read')
    open(unit=IOUT,file=trim(filtered_filename),status='unknown')
    read(IIN,*) nrec
    write(IOUT,*) nrec_filtered
    do irec = 1,nrec
      read(IIN,*) station_name,network_name,stlat,stlon,stele,stbur
      if (stlat >= LATITUDE_MIN .and. stlat <= LATITUDE_MAX .and. stlon >= LONGITUDE_MIN .and. stlon <= LONGITUDE_MAX) &
            write(IOUT,*) station_name,' ',network_name,' ',sngl(stlat),' ',sngl(stlon), ' ',sngl(stele), ' ',sngl(stbur)
    enddo
    close(IIN)
    close(IOUT)

    write(IMAIN,*)
    write(IMAIN,*) 'there are ',nrec,' stations in file ', trim(filename)
    write(IMAIN,*) 'saving ',nrec_filtered,' stations inside the model in file ', trim(filtered_filename)
    write(IMAIN,*) 'excluding ',nrec - nrec_filtered,' stations located outside the model'
    write(IMAIN,*)

  endif

  nfilter = nrec_filtered

  end subroutine station_filter

