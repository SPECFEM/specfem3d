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

!----
!----  locate_source finds the correct position of the source
!----

  subroutine locate_source(ibool,NSOURCES,myrank,NSPEC_AB,NGLOB_AB,xstore,ystore,zstore, &
                 xigll,yigll,zigll,NPROC, &
                 sec,t_cmt,yr,jda,ho,mi,utm_x_source,utm_y_source, &
                 NSTEP,DT,hdur,Mxx,Myy,Mzz,Mxy,Mxz,Myz, &
                 islice_selected_source,ispec_selected_source, &
                 xi_source,eta_source,gamma_source, &
                 LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX,Z_DEPTH_BLOCK, &
                 TOPOGRAPHY,itopo_bathy,UTM_PROJECTION_ZONE, &
                 PRINT_SOURCE_TIME_FUNCTION,SUPPRESS_UTM_PROJECTION, &
                 NX_TOPO,NY_TOPO,ORIG_LAT_TOPO,ORIG_LONG_TOPO,DEGREES_PER_CELL_TOPO)

  implicit none

  include "constants.h"

  integer NPROC,UTM_PROJECTION_ZONE
  integer NSTEP,NSPEC_AB,NGLOB_AB,NSOURCES,NX_TOPO,NY_TOPO

  logical TOPOGRAPHY,PRINT_SOURCE_TIME_FUNCTION,SUPPRESS_UTM_PROJECTION

  double precision DT,LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX,Z_DEPTH_BLOCK
  double precision ORIG_LAT_TOPO,ORIG_LONG_TOPO,DEGREES_PER_CELL_TOPO

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! use integer array to store topography values
  integer itopo_bathy(NX_TOPO,NY_TOPO)
  double precision long_corner,lat_corner,ratio_xi,ratio_eta

  integer myrank

! arrays containing coordinates of the points
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: xstore,ystore,zstore

  integer yr,jda,ho,mi

  double precision t_cmt(NSOURCES)
  double precision sec

  integer iprocloop

  integer i,j,k,ispec,iglob,isource

  double precision, dimension(NSOURCES) :: utm_x_source,utm_y_source
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

  double precision final_distance_source(NSOURCES)

  double precision x_target_source,y_target_source,z_target_source

  integer islice_selected_source(NSOURCES)

! timer MPI
  double precision, external :: wtime
  double precision time_start,tCPU

  integer ispec_selected_source(NSOURCES)

  integer ngather, ns, ne, ig, is, ng

  integer, dimension(NGATHER_SOURCES,0:NPROC-1) :: ispec_selected_source_all
  double precision, dimension(NGATHER_SOURCES,0:NPROC-1) :: xi_source_all,eta_source_all,gamma_source_all, &
     final_distance_source_all,x_found_source_all,y_found_source_all,z_found_source_all

  double precision hdur(NSOURCES), hdur_gaussian(NSOURCES), t0

  double precision, dimension(NSOURCES) :: Mxx,Myy,Mzz,Mxy,Mxz,Myz
  double precision, dimension(NSOURCES) :: xi_source,eta_source,gamma_source

  integer icornerlong,icornerlat
  double precision, dimension(NSOURCES) :: lat,long,depth,elevation
  double precision moment_tensor(6,NSOURCES)

  character(len=150) OUTPUT_FILES,plot_file

  double precision, dimension(NSOURCES) :: x_found_source,y_found_source,z_found_source
  double precision distmin

  integer ix_initial_guess_source,iy_initial_guess_source,iz_initial_guess_source

! for calculation of source time function
  integer it
  double precision time_source
  double precision, external :: comp_source_time_function

! **************

! get the base pathname for output files
  call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', 'OUTPUT_FILES')

! read all the sources
  call get_cmt(yr,jda,ho,mi,sec,t_cmt,hdur,lat,long,depth,moment_tensor,DT,NSOURCES)

! convert the half duration for triangle STF to the one for gaussian STF
  hdur_gaussian = hdur/SOURCE_DECAY_MIMIC_TRIANGLE

! define t0 as the earliest start time
  t0 = - 1.5d0 * minval(t_cmt-hdur)

! define topology of the control element
  call usual_hex_nodes(iaddx,iaddy,iaddz)

! get MPI starting time
  time_start = wtime()

! loop on all the sources
  do isource = 1,NSOURCES

! check that the current source is inside the model
  if(lat(isource) < LATITUDE_MIN .or. lat(isource) > LATITUDE_MAX .or. long(isource) < LONGITUDE_MIN &
       .or. long(isource) > LONGITUDE_MAX)  call exit_MPI(myrank,'the current source is outside the model')

  if(depth(isource) >= dabs(Z_DEPTH_BLOCK/1000.d0)) &
    call exit_MPI(myrank,'the current source is below the bottom of the model')

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
  Mzz(isource) = + moment_tensor(1,isource)
  Mxx(isource) = + moment_tensor(3,isource)
  Myy(isource) = + moment_tensor(2,isource)
  Mxz(isource) = + moment_tensor(5,isource)
  Myz(isource) = - moment_tensor(4,isource)
  Mxy(isource) = - moment_tensor(6,isource)

  call utm_geo(long(isource),lat(isource),utm_x_source(isource),utm_y_source(isource), &
                   UTM_PROJECTION_ZONE,ILONGLAT2UTM,SUPPRESS_UTM_PROJECTION)

! compute elevation of topography at the epicenter
  if(TOPOGRAPHY) then

! get coordinate of corner in bathy/topo model
    icornerlong = int((long(isource) - ORIG_LONG_TOPO) / DEGREES_PER_CELL_TOPO) + 1
    icornerlat = int((lat(isource) - ORIG_LAT_TOPO) / DEGREES_PER_CELL_TOPO) + 1

! avoid edge effects and extend with identical point if outside model
    if(icornerlong < 1) icornerlong = 1
    if(icornerlong > NX_TOPO-1) icornerlong = NX_TOPO-1
    if(icornerlat < 1) icornerlat = 1
    if(icornerlat > NY_TOPO-1) icornerlat = NY_TOPO-1

! compute coordinates of corner
    long_corner = ORIG_LONG_TOPO + (icornerlong-1)*DEGREES_PER_CELL_TOPO
    lat_corner = ORIG_LAT_TOPO + (icornerlat-1)*DEGREES_PER_CELL_TOPO

! compute ratio for interpolation
    ratio_xi = (long(isource) - long_corner) / DEGREES_PER_CELL_TOPO
    ratio_eta = (lat(isource) - lat_corner) / DEGREES_PER_CELL_TOPO

! avoid edge effects
    if(ratio_xi < 0.) ratio_xi = 0.
    if(ratio_xi > 1.) ratio_xi = 1.
    if(ratio_eta < 0.) ratio_eta = 0.
    if(ratio_eta > 1.) ratio_eta = 1.

! interpolate elevation at current point
    elevation(isource) = &
      itopo_bathy(icornerlong,icornerlat)*(1.-ratio_xi)*(1.-ratio_eta) + &
      itopo_bathy(icornerlong+1,icornerlat)*ratio_xi*(1.-ratio_eta) + &
      itopo_bathy(icornerlong+1,icornerlat+1)*ratio_xi*ratio_eta + &
      itopo_bathy(icornerlong,icornerlat+1)*(1.-ratio_xi)*ratio_eta

  else
    elevation = 0.d0
  endif

! compute the Cartesian position of the source
! take elevation of the surface into account
  x_target_source = utm_x_source(isource)
  y_target_source = utm_y_source(isource)
  z_target_source = - depth(isource)*1000.0d0 + elevation(isource)
  if(myrank == 0) write(IOVTK,*) x_target_source,y_target_source,z_target_source

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
          ispec_selected_source(isource)=ispec
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

    iglob = ibool(iax,iay,iaz,ispec_selected_source(isource))
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
  xi_source(isource) = xi
  eta_source(isource) = eta
  gamma_source(isource) = gamma
  x_found_source(isource) = x
  y_found_source(isource) = y
  z_found_source(isource) = z

! compute final distance between asked and found (converted to km)
  final_distance_source(isource) = dsqrt((x_target_source-x_found_source(isource))**2 + &
    (y_target_source-y_found_source(isource))**2 + (z_target_source-z_found_source(isource))**2)

! end of loop on all the sources
  enddo

! now gather information from all the nodes
  ngather = NSOURCES/NGATHER_SOURCES
  if (mod(NSOURCES,NGATHER_SOURCES)/= 0) ngather = ngather+1
  do ig = 1, ngather
    ns = (ig-1) * NGATHER_SOURCES + 1
    ne = min(ig*NGATHER_SOURCES, NSOURCES)
    ng = ne - ns + 1

    ispec_selected_source_all(:,:) = -1

  call gather_all_i(ispec_selected_source(ns:ne),ng,ispec_selected_source_all(1:ng,:),ng)

  call gather_all_dp(xi_source(ns:ne),ng,xi_source_all(1:ng,:),ng,NPROC)
  call gather_all_dp(eta_source(ns:ne),ng,eta_source_all(1:ng,:),ng,NPROC)
  call gather_all_dp(gamma_source(ns:ne),ng,gamma_source_all(1:ng,:),ng,NPROC)
  call gather_all_dp(final_distance_source(ns:ne),ng,final_distance_source_all(1:ng,:),ng,NPROC)
  call gather_all_dp(x_found_source(ns:ne),ng,x_found_source_all(1:ng,:),ng,NPROC)
  call gather_all_dp(y_found_source(ns:ne),ng,y_found_source_all(1:ng,:),ng,NPROC)
  call gather_all_dp(z_found_source(ns:ne),ng,z_found_source_all(1:ng,:),ng,NPROC)

! this is executed by main process only
  if(myrank == 0) then

! check that the gather operation went well
  if(any(ispec_selected_source_all(1:ng,:) == -1)) call exit_MPI(myrank,'gather operation failed for source')

! loop on all the sources
  do is = 1,ng
    isource = ns + is - 1

! loop on all the results to determine the best slice
  distmin = HUGEVAL
  do iprocloop = 0,NPROC-1
    if(final_distance_source_all(is,iprocloop) < distmin) then
      distmin = final_distance_source_all(is,iprocloop)
      islice_selected_source(isource) = iprocloop
      ispec_selected_source(isource) = ispec_selected_source_all(is,iprocloop)
      xi_source(isource) = xi_source_all(is,iprocloop)
      eta_source(isource) = eta_source_all(is,iprocloop)
      gamma_source(isource) = gamma_source_all(is,iprocloop)
      x_found_source(isource) = x_found_source_all(is,iprocloop)
      y_found_source(isource) = y_found_source_all(is,iprocloop)
      z_found_source(isource) = z_found_source_all(is,iprocloop)
    endif
  enddo
  final_distance_source(isource) = distmin

  enddo
  endif
  enddo

  if (myrank == 0) then

  do isource = 1,NSOURCES

  if(SHOW_DETAILS_LOCATE_SOURCE .or. NSOURCES == 1) then

    write(IMAIN,*)
    write(IMAIN,*) '*************************************'
    write(IMAIN,*) ' locating source ',isource
    write(IMAIN,*) '*************************************'
    write(IMAIN,*)
    write(IMAIN,*) 'source located in slice ',islice_selected_source(isource)
    write(IMAIN,*) '               in element ',ispec_selected_source(isource)
    write(IMAIN,*)
    write(IMAIN,*) '   xi coordinate of source in that element: ',xi_source(isource)
    write(IMAIN,*) '  eta coordinate of source in that element: ',eta_source(isource)
    write(IMAIN,*) 'gamma coordinate of source in that element: ',gamma_source(isource)

! add message if source is a Heaviside
    if(hdur(isource) < 5.*DT) then
      write(IMAIN,*)
      write(IMAIN,*) 'Source time function is a Heaviside, convolve later'
      write(IMAIN,*)
    endif

    write(IMAIN,*)
    write(IMAIN,*) ' half duration: ',hdur(isource),' seconds'
    write(IMAIN,*) '    time shift: ',t_cmt(isource),' seconds'

    write(IMAIN,*)
    write(IMAIN,*) 'original (requested) position of the source:'
    write(IMAIN,*)
    write(IMAIN,*) '      latitude: ',lat(isource)
    write(IMAIN,*) '     longitude: ',long(isource)
    write(IMAIN,*)
    write(IMAIN,*) '         UTM x: ',utm_x_source(isource)
    write(IMAIN,*) '         UTM y: ',utm_y_source(isource)
    write(IMAIN,*) '         depth: ',depth(isource),' km'
    if(TOPOGRAPHY) write(IMAIN,*) 'topo elevation: ',elevation(isource),' m'

    write(IMAIN,*)
    write(IMAIN,*) 'position of the source that will be used:'
    write(IMAIN,*)
    write(IMAIN,*) '         UTM x: ',x_found_source(isource)
    write(IMAIN,*) '         UTM y: ',y_found_source(isource)
    write(IMAIN,*) '         depth: ',dabs(z_found_source(isource) - elevation(isource))/1000.,' km'
    write(IMAIN,*)

! display error in location estimate
    write(IMAIN,*) 'error in location of the source: ',sngl(final_distance_source(isource)),' m'

! add warning if estimate is poor
! (usually means source outside the mesh given by the user)
    if(final_distance_source(isource) > 3000.d0) then
      write(IMAIN,*)
      write(IMAIN,*) '*****************************************************'
      write(IMAIN,*) '*****************************************************'
      write(IMAIN,*) '***** WARNING: source location estimate is poor *****'
      write(IMAIN,*) '*****************************************************'
      write(IMAIN,*) '*****************************************************'
    endif

  endif  ! end of detailed output to locate source

  if(PRINT_SOURCE_TIME_FUNCTION) then

  write(IMAIN,*)
  write(IMAIN,*) 'printing the source-time function'

! print the source-time function
  if(NSOURCES == 1) then
    plot_file = '/plot_source_time_function.txt'
  else
   if(isource < 10) then
      write(plot_file,"('/plot_source_time_function',i1,'.txt')") isource
    else
      write(plot_file,"('/plot_source_time_function',i2,'.txt')") isource
    endif
  endif
  open(unit=27,file=trim(OUTPUT_FILES)//plot_file,status='unknown')

  do it=1,NSTEP
    time_source = dble(it-1)*DT
    write(27,*) sngl(time_source-t0),sngl(comp_source_time_function(time_source-t0-t_cmt(isource),hdur_gaussian(isource)))
  enddo
  close(27)

  endif

! end of loop on all the sources
  enddo

! display maximum error in location estimate
    write(IMAIN,*)
    write(IMAIN,*) 'maximum error in location of the sources: ',sngl(maxval(final_distance_source)),' m'
    write(IMAIN,*)

  endif     ! end of section executed by main process only

! main process broadcasts the results to all the slices
  call bcast_all_i(islice_selected_source,NSOURCES)
  call bcast_all_i(ispec_selected_source,NSOURCES)
  call bcast_all_dp(xi_source,NSOURCES)
  call bcast_all_dp(eta_source,NSOURCES)
  call bcast_all_dp(gamma_source,NSOURCES)

! elapsed time since beginning of source detection
  if(myrank == 0) then
    tCPU = wtime() - time_start
    write(IMAIN,*)
    write(IMAIN,*) 'Elapsed time for detection of sources in seconds = ',tCPU
    write(IMAIN,*)
    write(IMAIN,*) 'End of source detection - done'
    write(IMAIN,*)
  endif

  end subroutine locate_source

