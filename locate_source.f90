!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 1
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology October 2002
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

!----
!----  locate_source finds the correct position of the source
!----

  subroutine locate_source(ibool,isource,NSOURCES,myrank,NSPEC_AB,NGLOB_AB,xstore,ystore,zstore, &
                 xigll,yigll,zigll,NPROC, &
                 sec,t_cmt,yr,jda,ho,mi,utm_x_source,utm_y_source, &
                 NSTEP,DT,hdur,Mxx,Myy,Mzz,Mxy,Mxz,Myz, &
                 islice_selected_source,ispec_selected_source, &
                 xi_source,eta_source,gamma_source, &
                 LAT_MIN,LAT_MAX,LONG_MIN,LONG_MAX,Z_DEPTH_BLOCK, &
                 TOPOGRAPHY,itopo_bathy_basin,UTM_PROJECTION_ZONE)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

  integer NPROC,UTM_PROJECTION_ZONE
  integer NSTEP,NSPEC_AB,NGLOB_AB,NSOURCES

  logical TOPOGRAPHY

  double precision DT,LAT_MIN,LAT_MAX,LONG_MIN,LONG_MAX,Z_DEPTH_BLOCK

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! use integer array to store topography values
  integer itopo_bathy_basin(NX_TOPO,NY_TOPO)
  double precision long_corner,lat_corner,ratio_xi,ratio_eta

  integer myrank,isource

! arrays containing coordinates of the points
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: xstore,ystore,zstore

  integer yr,jda,ho,mi

  double precision sec,t_cmt

  integer iprocloop

  integer i,j,k,ispec,iglob
  integer ier

  double precision utm_x_source,utm_y_source
  double precision dist
  double precision xi,eta,gamma,dx,dy,dz,dxi,deta

! Gauss-Lobatto-Legendre points of integration
  double precision xigll(NGLLX)
  double precision yigll(NGLLY)
  double precision zigll(NGLLZ)

! topology of the control points of the surface element
  integer iax,iay,iaz
  integer iaddx(NGNOD),iaddy(NGNOD),iaddz(NGNOD)

! coordinates of the control points of the surface element
  double precision xelm(NGNOD),yelm(NGNOD),zelm(NGNOD)

  integer iter_loop

  integer ia
  double precision x,y,z
  double precision xix,xiy,xiz
  double precision etax,etay,etaz
  double precision gammax,gammay,gammaz
  double precision dgamma

  double precision final_distance_source

  double precision x_target_source,y_target_source,z_target_source

  integer islice_selected_source

! timer MPI
  double precision time_start,tCPU

  integer ispec_selected_source

  integer, dimension(0:NPROC-1) :: ispec_selected_source_all
  double precision, dimension(0:NPROC-1) :: xi_source_all,eta_source_all,gamma_source_all, &
     final_distance_source_all,x_found_source_all,y_found_source_all,z_found_source_all

  double precision hdur

  double precision Mxx,Myy,Mzz,Mxy,Mxz,Myz
  double precision xi_source,eta_source,gamma_source

  integer icornerlong,icornerlat
  double precision lat,long,depth,elevation
  double precision moment_tensor(6)

  character(len=150) cmt_file,plot_file

  double precision x_found_source,y_found_source,z_found_source
  double precision distmin

  integer ix_initial_guess_source,iy_initial_guess_source,iz_initial_guess_source

! for calculation of source time function
  integer it
  double precision time_source
  double precision, external :: comp_source_time_function

! number of points to plot the source time function
  integer, parameter :: NSAMP_PLOT_SOURCE = 1000

! receiver information
! timing information for the stations
! station information for writing the seismograms

! **************

! get MPI starting time
  time_start = MPI_WTIME()

  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '*****************'
    write(IMAIN,*) ' locating source ',isource
    write(IMAIN,*) '*****************'
    write(IMAIN,*)
  endif

! define topology of the control element
  call usual_hex_nodes(iaddx,iaddy,iaddz)

! get source information
! get source information
  if(NSOURCES == 1) then
    cmt_file='DATA/CMTSOLUTION'
  else
    if(isource < 10) then
      write(cmt_file,"('DATA/CMTSOLUTION',i1)") isource
    elseif(isource < 100) then
      write(cmt_file,"('DATA/CMTSOLUTION',i2)") isource
    elseif(isource < 1000) then
      write(cmt_file,"('DATA/CMTSOLUTION',i3)") isource
    else
      call exit_MPI(myrank,'too many sources')
    endif
  endif

  call get_cmt(cmt_file,yr,jda,ho,mi,sec,t_cmt,hdur,lat,long,depth,moment_tensor,DT)

  if(isource == 1) then
    if(t_cmt /= 0.) call exit_MPI(myrank,'t_cmt for the first source should be zero')
  else
    if(t_cmt < 0.) call exit_MPI(myrank,'t_cmt should not be less than zero')
  endif

! check that the source is inside the basin model
  if(lat <= LAT_MIN .or. lat >= LAT_MAX .or. long <= LONG_MIN .or. long >= LONG_MAX) &
    call exit_MPI(myrank,'the source is outside the model')

  if(depth >= dabs(Z_DEPTH_BLOCK/1000.d0)) &
    call exit_MPI(myrank,'the source is below the bottom of the model')

!
! r -> z, theta -> -y, phi -> x
!
!  Mrr =  Mzz
!  Mtt =  Myy
!  Mpp =  Mxx
!  Mrt = -Myz
!  Mrp =  Mxz
!  Mtp = -Mxy

! get the moment tensor
  Mzz = + moment_tensor(1)
  Mxx = + moment_tensor(3)
  Myy = + moment_tensor(2)
  Mxz = + moment_tensor(5)
  Myz = - moment_tensor(4)
  Mxy = - moment_tensor(6)

  call utm_geo(long,lat,utm_x_source,utm_y_source,UTM_PROJECTION_ZONE,ILONGLAT2UTM)

! compute elevation of topography at the epicenter
  if(TOPOGRAPHY) then

! get coordinate of corner in bathy/topo model
    icornerlong = int((long - ORIG_LONG_TOPO) / DEGREES_PER_CELL_TOPO) + 1
    icornerlat = int((lat - ORIG_LAT_TOPO) / DEGREES_PER_CELL_TOPO) + 1

! avoid edge effects and extend with identical point if outside model
    if(icornerlong < 1) icornerlong = 1
    if(icornerlong > NX_TOPO-1) icornerlong = NX_TOPO-1
    if(icornerlat < 1) icornerlat = 1
    if(icornerlat > NY_TOPO-1) icornerlat = NY_TOPO-1

! compute coordinates of corner
    long_corner = ORIG_LONG_TOPO + (icornerlong-1)*DEGREES_PER_CELL_TOPO
    lat_corner = ORIG_LAT_TOPO + (icornerlat-1)*DEGREES_PER_CELL_TOPO

! compute ratio for interpolation
    ratio_xi = (long - long_corner) / DEGREES_PER_CELL_TOPO
    ratio_eta = (lat - lat_corner) / DEGREES_PER_CELL_TOPO

! avoid edge effects
    if(ratio_xi < 0.) ratio_xi = 0.
    if(ratio_xi > 1.) ratio_xi = 1.
    if(ratio_eta < 0.) ratio_eta = 0.
    if(ratio_eta > 1.) ratio_eta = 1.

! interpolate elevation at current point
    elevation = &
      itopo_bathy_basin(icornerlong,icornerlat)*(1.-ratio_xi)*(1.-ratio_eta) + &
      itopo_bathy_basin(icornerlong+1,icornerlat)*ratio_xi*(1.-ratio_eta) + &
      itopo_bathy_basin(icornerlong+1,icornerlat+1)*ratio_xi*ratio_eta + &
      itopo_bathy_basin(icornerlong,icornerlat+1)*(1.-ratio_xi)*ratio_eta

  else
    elevation = 0.d0
  endif

! compute the Cartesian position of the source
! take elevation of the surface into account
  x_target_source = utm_x_source
  y_target_source = utm_y_source
  z_target_source = - depth*1000.0d0 + elevation

! set distance to huge initial value
  distmin = HUGEVAL

  do ispec=1,NSPEC_AB

! loop only on points inside the element
! exclude edges to ensure this point is not shared with other elements
  do k=2,NGLLZ-1
    do j=2,NGLLY-1
      do i=2,NGLLX-1

!       keep this point if it is closer to the receiver
        iglob = ibool(i,j,k,ispec)
        dist=dsqrt((x_target_source-dble(xstore(iglob)))**2 &
                  +(y_target_source-dble(ystore(iglob)))**2 &
                  +(z_target_source-dble(zstore(iglob)))**2)
        if(dist < distmin) then
          distmin=dist
          ispec_selected_source=ispec
          ix_initial_guess_source = i
          iy_initial_guess_source = j
          iz_initial_guess_source = k
        endif

      enddo
    enddo
  enddo

! end of loop on all the elements in current slice
  enddo

! *******************************************
! find the best (xi,eta,gamma) for the source
! *******************************************

! use initial guess in xi, eta and gamma
  xi = xigll(ix_initial_guess_source)
  eta = yigll(iy_initial_guess_source)
  gamma = zigll(iz_initial_guess_source)

! define coordinates of the control points of the element

  do ia=1,NGNOD

    if(iaddx(ia) == 0) then
      iax = 1
    else if(iaddx(ia) == 1) then
      iax = (NGLLX+1)/2
    else if(iaddx(ia) == 2) then
      iax = NGLLX
    else
      call exit_MPI(myrank,'incorrect value of iaddx')
    endif

    if(iaddy(ia) == 0) then
      iay = 1
    else if(iaddy(ia) == 1) then
      iay = (NGLLY+1)/2
    else if(iaddy(ia) == 2) then
      iay = NGLLY
    else
      call exit_MPI(myrank,'incorrect value of iaddy')
    endif

    if(iaddz(ia) == 0) then
      iaz = 1
    else if(iaddz(ia) == 1) then
      iaz = (NGLLZ+1)/2
    else if(iaddz(ia) == 2) then
      iaz = NGLLZ
    else
      call exit_MPI(myrank,'incorrect value of iaddz')
    endif

    iglob = ibool(iax,iay,iaz,ispec_selected_source)
    xelm(ia) = dble(xstore(iglob))
    yelm(ia) = dble(ystore(iglob))
    zelm(ia) = dble(zstore(iglob))

  enddo

! iterate to solve the non linear system
  do iter_loop = 1,NUM_ITER

! recompute jacobian for the new point
    call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
           xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

! compute distance to target location
  dx = - (x - x_target_source)
  dy = - (y - y_target_source)
  dz = - (z - z_target_source)

! compute increments
  dxi  = xix*dx + xiy*dy + xiz*dz
  deta = etax*dx + etay*dy + etaz*dz
  dgamma = gammax*dx + gammay*dy + gammaz*dz

! update values
  xi = xi + dxi
  eta = eta + deta
  gamma = gamma + dgamma

! impose that we stay in that element
! (useful if user gives a source outside the mesh for instance)
  if (xi > 1.d0) xi = 1.d0
  if (xi < -1.d0) xi = -1.d0
  if (eta > 1.d0) eta = 1.d0
  if (eta < -1.d0) eta = -1.d0
  if (gamma > 1.d0) gamma = 1.d0
  if (gamma < -1.d0) gamma = -1.d0

  enddo

! compute final coordinates of point found
  call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
         xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

! store xi,eta,gamma and x,y,z of point found
  xi_source = xi
  eta_source = eta
  gamma_source = gamma
  x_found_source = x
  y_found_source = y
  z_found_source = z

! compute final distance between asked and found (converted to km)
  final_distance_source = dsqrt((x_target_source-x_found_source)**2 + &
    (y_target_source-y_found_source)**2 + (z_target_source-z_found_source)**2)

! synchronize all the processes to make sure all the estimates are available
  call MPI_BARRIER(MPI_COMM_WORLD,ier)

! for MPI version, now gather information from all the nodes
  ispec_selected_source_all(:) = -1
  call MPI_GATHER(ispec_selected_source,1,MPI_INTEGER,ispec_selected_source_all,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  call MPI_GATHER(xi_source,1,MPI_DOUBLE_PRECISION,xi_source_all,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(eta_source,1,MPI_DOUBLE_PRECISION,eta_source_all,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(gamma_source,1,MPI_DOUBLE_PRECISION,gamma_source_all,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(final_distance_source,1,MPI_DOUBLE_PRECISION, &
    final_distance_source_all,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(x_found_source,1,MPI_DOUBLE_PRECISION, &
    x_found_source_all,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(y_found_source,1,MPI_DOUBLE_PRECISION, &
    y_found_source_all,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(z_found_source,1,MPI_DOUBLE_PRECISION, &
    z_found_source_all,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

! MPI this is executed by main process only
  if(myrank == 0) then

! check that the gather operation went well
  if(any(ispec_selected_source_all(:) == -1)) call exit_MPI(myrank,'gather operation failed for source')

! MPI loop on all the results to determine the best slice
  distmin = HUGEVAL
  do iprocloop = 0,NPROC-1
    if(final_distance_source_all(iprocloop) < distmin) then
      distmin = final_distance_source_all(iprocloop)
      islice_selected_source = iprocloop
      ispec_selected_source = ispec_selected_source_all(iprocloop)
      xi_source = xi_source_all(iprocloop)
      eta_source = eta_source_all(iprocloop)
      gamma_source = gamma_source_all(iprocloop)
      x_found_source = x_found_source_all(iprocloop)
      y_found_source = y_found_source_all(iprocloop)
      z_found_source = z_found_source_all(iprocloop)
    endif
  enddo
  final_distance_source = distmin

    write(IMAIN,*)
    write(IMAIN,*) 'source located in slice ',islice_selected_source
    write(IMAIN,*) '               in element ',ispec_selected_source
    write(IMAIN,*)
    write(IMAIN,*) '   xi coordinate of source in that element: ',xi_source
    write(IMAIN,*) '  eta coordinate of source in that element: ',eta_source
    write(IMAIN,*) 'gamma coordinate of source in that element: ',gamma_source

! add message if source is a Heaviside
    if(hdur < 6.*DT) then
      write(IMAIN,*)
      write(IMAIN,*) 'Source time function is a Heaviside, convolve later'
      write(IMAIN,*)
    endif

    write(IMAIN,*)
    write(IMAIN,*) 'original (requested) position of the source:'
    write(IMAIN,*)
    write(IMAIN,*) '      latitude: ',lat
    write(IMAIN,*) '     longitude: ',long
    write(IMAIN,*)
    write(IMAIN,*) '         UTM x: ',utm_x_source
    write(IMAIN,*) '         UTM y: ',utm_y_source
    write(IMAIN,*) '         depth: ',depth,' km'
    if(TOPOGRAPHY) write(IMAIN,*) 'topo elevation: ',elevation,' m'

    write(IMAIN,*)
    write(IMAIN,*) 'position of the source that will be used:'
    write(IMAIN,*)
    write(IMAIN,*) '         UTM x: ',x_found_source
    write(IMAIN,*) '         UTM y: ',y_found_source
    write(IMAIN,*) '         depth: ',dabs(z_found_source - elevation)/1000.,' km'
    write(IMAIN,*)

! display error in location estimate
    write(IMAIN,*) 'error in location of the source: ',sngl(final_distance_source),' m'

! add warning if estimate is poor
! (usually means source outside the mesh given by the user)
    if(final_distance_source > 3000.d0) then
      write(IMAIN,*)
      write(IMAIN,*) '*****************************************************'
      write(IMAIN,*) '*****************************************************'
      write(IMAIN,*) '***** WARNING: source location estimate is poor *****'
      write(IMAIN,*) '*****************************************************'
      write(IMAIN,*) '*****************************************************'
    endif

  write(IMAIN,*)
  write(IMAIN,*) 'half-duration of the source: ',hdur,' s'

  if(PRINT_SOURCE_TIME_FUNCTION) then

  write(IMAIN,*)
  write(IMAIN,*) 'printing the source-time function'

! print the source-time function
  if(NSOURCES == 1) then
    plot_file = 'OUTPUT_FILES/plot_source_time_function.txt'
  else
   if(isource < 10) then
      write(plot_file,"('OUTPUT_FILES/plot_source_time_function',i1,'.txt')") isource
    else
      write(plot_file,"('OUTPUT_FILES/plot_source_time_function',i2,'.txt')") isource
    endif
  endif
  open(unit=27,file=plot_file(1:len_trim(plot_file)),status='unknown')

  do it=1,NSTEP
    time_source = dble(it-1)*DT
    write(27,*) sngl(time_source),sngl(comp_source_time_function(time_source,hdur))
  enddo
  close(27)

  endif

! elapsed time since beginning of mesh generation
  tCPU = MPI_WTIME() - time_start
  write(IMAIN,*)
  write(IMAIN,*) 'Elapsed time for source detection in seconds = ',tCPU
  write(IMAIN,*)
  write(IMAIN,*) 'End of source detection - done'
  write(IMAIN,*)

  endif     ! end of section executed by main process only

! main process broadcasts the results to all the slices
  call MPI_BCAST(islice_selected_source,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(ispec_selected_source,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(xi_source,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(eta_source,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(gamma_source,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

! synchronize all the processes to make sure everybody has finished
  call MPI_BARRIER(MPI_COMM_WORLD,ier)

  end subroutine locate_source

