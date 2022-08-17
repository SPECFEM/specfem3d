module wave2d_sub2

  use wave2d_constants
  use wave2d_variables

  implicit none

! this module contains subroutines pertaining to filtering gridpoints
! and converting amoung UTM, mesh, and index coordinates

contains

!=====================================================================

  subroutine locate_targets(nrec, x_target, z_target, &
       iglob_selected_rec, ispec_selected_rec, xi_receiver, gamma_receiver)

    ! adapted from locate_receivers.f90 (3D code)
    ! This does NOT remove redundant closest gridpoints, like in set_glob.f90

    ! input
    integer, intent(inout) :: nrec
    double precision, intent(inout) :: x_target(nrec), z_target(nrec)

    ! output -- note that (ispec, xi, gamma) are optional output
    integer,          dimension(nrec), intent(out)           :: iglob_selected_rec
    integer,          dimension(nrec), intent(out), optional :: ispec_selected_rec
    double precision, dimension(nrec), intent(out), optional :: xi_receiver, gamma_receiver

    ! local
    double precision, dimension(nrec) :: x_found, z_found, final_distance, dmin_selected_rec
    integer, dimension(nrec) :: ix_initial_guess, iz_initial_guess
    double precision :: dmin, dist
    double precision :: xtemp, ztemp, xi, gamma, dx, dz, dxi, dgamma
    double precision :: xix,xiz,gammax,gammaz
    integer :: irec, ispec, i, j, iglob, iter_loop, ispec_iterate
    logical :: i_xi_gamma

    ! control points of the surface element
    !integer iax, iaz, ia
    !integer iaddx(NGNOD2D), iaddz(NGNOD2D)
    !double precision xelm(NGNOD2D),zelm(NGNOD2D)

  !----------------------------------

  print *, nrec,' input target points into locate_targets.f90'
  print *, ' all target points should now lie within the mesh (station_filter.f90)'

  ! boolean: whether the (ispec, xi, gamma) are in the subroutine call
  i_xi_gamma = present(ispec_selected_rec) .and. present(xi_receiver) .and. present(gamma_receiver)

  ! loop over target points to get the (xi,gamma) for the closest gridpoint
  do irec=1,nrec

      dmin = sqrt(LENGTH**2+HEIGHT**2)  ! max possible distance

      do ispec=1,NSPEC

        ! loop only on points inside the element
        ! exclude edges to ensure this point is not shared with other elements
        do j = 2,NGLLZ-1
          do i = 2,NGLLX-1

            iglob = ibool(i,j,ispec)
            dist = sqrt((x_target(irec) - x(iglob))**2 + (z_target(irec) - z(iglob))**2)

            ! keep this point if it is closer to the receiver
            if (dist < dmin) then
              dmin = dist
              dmin_selected_rec(irec) = dmin
              iglob_selected_rec(irec) = iglob  ! closest gridpoint
              if (i_xi_gamma) then
                ispec_selected_rec(irec) = ispec
                ix_initial_guess(irec) = i
                iz_initial_guess(irec) = j
              endif
            endif

          enddo
        enddo
      enddo  ! ispec
   enddo  ! irec

   ! TESTING ONLY (01-Feb-2006)
   ! If the GLL points are your target points, then this should give the same
   ! "interpolated" seismograms as if you simply saved the wavefield at the GLL point.
   !x_target(:) = x(iglob_selected_rec(:))
   !z_target(:) = z(iglob_selected_rec(:))

  ! ****************************************
  ! find the best (xi,gamma) for each target point
  ! ****************************************

  ! if the (xi, gamma) are desired
  if (i_xi_gamma) then

     do irec = 1,nrec

        ! element that contains the target point
        ispec_iterate = ispec_selected_rec(irec)

        ! use initial guess in xi and gamma
        xi    = xigll(ix_initial_guess(irec))
        gamma = zigll(iz_initial_guess(irec))

!!$    ! define coordinates of the control points of the element
!!$    ! (ordering of these four points shouldn't matter)
!!$    ! We could also access the vectors x1,x2,z1,z2 instead.
!!$    do ia = 1,NGNOD2D    ! (NGNOD2D = 4)
!!$
!!$       if (ia==1) then
!!$          iax = 1
!!$          iaz = 1
!!$       else if (ia==2) then
!!$          iax = NGLLX
!!$          iaz = 1
!!$       else if (ia==3) then
!!$          iax = NGLLX
!!$          iaz = NGLLZ
!!$       else if (ia==4) then
!!$          iax = 1
!!$          iaz = NGLLZ
!!$       endif
!!$       iglob = ibool(iax,iaz,ispec_iterate)
!!$       xelm(ia) = x(iglob)
!!$       zelm(ia) = z(iglob)
!!$    enddo  ! ia

        ! iterate to solve the non linear system (it seems we only need one iteration for 2D)
        do iter_loop = 1,NUM_ITER

           ! recompute jacobian for the new point
           ! INPUT:  xelm,zelm,xi,gamma
           ! OUTPUT: xtemp,ztemp,xix,xiz,gammax,gammaz

           !call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,xtemp,ytemp,ztemp, &
           !   xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)
           call recompute_jacobian_2d(ispec_iterate,xi,gamma,xtemp,ztemp,xix,xiz,gammax,gammaz)

           ! compute distance to target location
           dx = -(xtemp - x_target(irec))
           dz = -(ztemp - z_target(irec))

           ! compute increments
           ! gamma does not change since we know the receiver is exactly on the surface
           dxi  = xix*dx + xiz*dz
           dgamma = gammax*dx + gammaz*dz

           !print *
           !print *, 'iteration   : ', iter_loop
           !print *, ' (dx, dz)   : ', dx,dz
           !print *, ' (xi,gamma) : ', xi,gamma
           !print *, xix,xiz,gammax,gammaz

           ! update values
           xi = xi + dxi
           gamma = gamma + dgamma

           !print *, xi,gamma,dxi,dgamma

           ! impose that we stay in that element
           ! (useful if user gives a receiver outside the mesh for instance)
           ! we can go slightly outside the [1,1] segment since with finite elements
           ! the polynomial solution is defined everywhere
           ! can be useful for convergence of iterative scheme with distorted elements
           if (xi > 1.10) xi = 1.10
           if (xi < -1.10) xi = -1.10
           if (gamma > 1.10) gamma = 1.10
           if (gamma < -1.10) gamma = -1.10

        enddo  ! iter_loop

        ! compute final coordinates of point found
        xtemp = 0. ; ztemp = 0.  ! re-initialize
        call recompute_jacobian_2d(ispec_iterate,xi,gamma,xtemp,ztemp,xix,xiz,gammax,gammaz)

        ! store xi,gamma and x,z of point found
        xi_receiver(irec) = xi
        gamma_receiver(irec) = gamma
        x_found(irec) = xtemp
        z_found(irec) = ztemp

        ! compute final distance between asked and found
        final_distance(irec) = sqrt((x_target(irec) - x_found(irec))**2 + &
             (z_target(irec) - z_found(irec))**2)
     enddo  ! irec

     ! display information
     do irec=1,nrec

        if (final_distance(irec) == HUGEVAL) stop 'error locating receiver'

        if (0 == 1) then
           print *
           print *,                'target point # ', irec
           write(*,'(a,1f18.8)')   '  target x (km)            : ', x_target(irec)/1000.
           write(*,'(a,1f18.8)')   '  target z (km)            : ', z_target(irec)/1000.
           write(*,'(a,1e18.8,a)') '  closest gridpoint found  : ', dmin_selected_rec(irec)/1000.,' km away'
           write(*,'(a,1f18.8)')   '  closest gridpoint x (km) : ', x(iglob_selected_rec(irec))/1000.
           write(*,'(a,1f18.8)')   '  closest gridpoint z (km) : ', z(iglob_selected_rec(irec))/1000.
           print *,                '     element = ', ispec_selected_rec(irec)
           print *,                '          xi = ', xigll(ix_initial_guess(irec))
           print *,                '       gamma = ', zigll(iz_initial_guess(irec))
           write(*,'(a,1f18.8)')   '  estimate x (km)            : ', x_found(irec)/1000.
           write(*,'(a,1f18.8)')   '  estimate z (km)            : ', z_found(irec)/1000.
           print *,                '  closest estimate found   : ', final_distance(irec)/1000.,' km away'
           print *,                '     element = ', ispec_selected_rec(irec)
           print *,                '          xi = ', xi_receiver(irec)
           print *,                '       gamma = ', gamma_receiver(irec)
        endif

        ! add warning if estimate is poor
        ! (usually means receiver outside the mesh given by the user)
        if (final_distance(irec) > 5.) then
           print *, 'station # ',irec
           print *, '*******************************************************'
           print *, '***** WARNING: receiver location estimate is poor *****'
           print *, '*******************************************************'
        endif

     enddo

     ! replace the input target points with the (x,z) corresponding to the (xi,gamma) values
     x_target(:) = x_found(:)
     z_target(:) = z_found(:)

  else
     ! replace the input target points with the (x,z) of the closest gridpoint
     x_target(:) = x(iglob_selected_rec(:))
     z_target(:) = z(iglob_selected_rec(:))

  endif

  !print *
  !print *, nrec,' output gridpoint indices from locate_targets.f90'

  end subroutine locate_targets

!=====================================================================

  subroutine recompute_jacobian_2d(ispec, xi, gamma, xtemp, ztemp, xix, xiz, gammax, gammaz)

! subroutine recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
!                  xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

  ! recompute jacobian for any (xi,eta,gamma) point, not necessarily a GLL point
  ! adapted from recompute_jacobian.f90 and mesher.f90 (above)
  !
  ! INPUT:  ispec,xi,gamma
  ! OUTPUT: xtemp,ztemp,xix,xiz,gammax,gammaz

  ! input
  integer ispec
  double precision :: xi,gamma
  !double precision xelm(NGNOD2D),zelm(NGNOD2D)  ! coordinates of the control points

  ! output
  double precision :: xtemp,ztemp,xix,xiz,gammax,gammaz

  double precision :: jacob

  !---------------------------------------------
  !
  ! jacobian = | dx/dxi dx/dgamma | = (z2-z1)*(x2-x1)/4  as dx/dgamma=dz/dxi = 0
  !            | dz/dxi dz/dgamma |

  ! jacobian, integration weight
  xix = 2. / (x2(ispec)-x1(ispec))
  xiz = 0.
  gammax = 0.
  gammaz = 2. / (z2(ispec)-z1(ispec))

  jacob = (z2(ispec)-z1(ispec))*(x2(ispec)-x1(ispec)) / 4.

  ! find the (x,z) corresponding to (xi,gamma,ispec)
  xtemp = 0.5*(1.-    xi)*x1(ispec) + 0.5*(1.+   xi)*x2(ispec)
  ztemp = 0.5*(1.- gamma)*z1(ispec) + 0.5*(1.+gamma)*z2(ispec)

  if (jacob <= 0.) stop '2D Jacobian undefined'

  end subroutine recompute_jacobian_2d

!=====================================================================

  subroutine lagrange_poly(xi,NGLL,xigll,h,hprime)

  ! subroutine to compute the Lagrange interpolants based upon the GLL points
  ! and their first derivatives at any point xi in [-1,1]
  ! copied from lagrange_poly.f90 in Qinya's original codes

    integer,intent(in) :: NGLL
    double precision,intent(in) :: xi, xigll(NGLL)
    double precision,intent(out) :: h(NGLL), hprime(NGLL)

    integer dgr,i,j
    double precision prod1,prod2

  !---------------------------------------------

  do dgr=1,NGLL

     prod1 = 1.0
     prod2 = 1.0
     do i=1,NGLL
        if (i /= dgr) then
           prod1 = prod1*(xi-xigll(i))
           prod2 = prod2*(xigll(dgr)-xigll(i))
        endif
     enddo
     h(dgr)=prod1/prod2

     hprime(dgr)=0.0
     do i=1,NGLL
        if (i /= dgr) then
           prod1=1.0
           do j=1,NGLL
              if (j /= dgr .and. j /= i) prod1 = prod1*(xi-xigll(j))
           enddo
           hprime(dgr) = hprime(dgr)+prod1
        endif
     enddo
     hprime(dgr) = hprime(dgr)/prod2

  enddo

  end subroutine lagrange_poly

!=====================================================================

  double precision function lagrange_deriv_GLL(I,j,ZGLL,NZ)

  ! subroutine to compute the derivative of the Lagrange interpolants
  ! at the GLL points at any given GLL point
  !
  !------------------------------------------------------------------------
  !  Compute the value of the derivative of the I-th
  !  Lagrange interpolant through the
  !  NZ Gauss-Lobatto Legendre points ZGLL at point ZGLL(j)
  !------------------------------------------------------------------------

    integer i,j,nz
    double precision zgll(0:nz-1)
    integer degpoly
    double precision, external :: pnleg,pndleg   ! external functions

  !---------------------------------------------

  degpoly = nz - 1
  if (i == 0 .and. j == 0) then
    lagrange_deriv_GLL = - dble(degpoly)*(dble(degpoly)+1.d0) / 4.d0

  else if (i == degpoly .and. j == degpoly) then
    lagrange_deriv_GLL = dble(degpoly)*(dble(degpoly)+1.d0) / 4.d0

  else if (i == j) then
    lagrange_deriv_GLL = 0.d0

  else
    lagrange_deriv_GLL = pnleg(zgll(j),degpoly) / &
      (pnleg(zgll(i),degpoly)*(zgll(j)-zgll(i))) &
      + (1.d0-zgll(j)*zgll(j))*pndleg(zgll(j),degpoly) / (dble(degpoly)* &
      (dble(degpoly)+1.d0)*pnleg(zgll(i),degpoly)*(zgll(j)-zgll(i))*(zgll(j)-zgll(i)))
  endif

  end function lagrange_deriv_GLL

!=====================================================================

  subroutine set_glob(nrec, x_rec, z_rec, rglob)

    integer, intent(inout) :: nrec
    double precision, intent(in) :: x_rec(nrec), z_rec(nrec)
    integer, intent(out) :: rglob(nrec)

    double precision :: dmin, d
    integer, dimension(nrec) :: rglobtemp
    integer :: irec, ispec, i, j, k, iglob, itemp, iflag

    !----------------------------------

    print *, nrec,' input target points into set_glob.f90'

    if (nrec /= 0) then

       ! find the closest gridpoint to the target point
       do irec = 1, nrec
          dmin = sqrt(LENGTH**2+HEIGHT**2)  ! max possible distance
          do iglob = 1,NGLOB
             d = sqrt((x_rec(irec)-x(iglob))**2+(z_rec(irec)-z(iglob))**2)
             if (d < dmin) then
                dmin  = d
                rglob(irec) = iglob
             endif
          enddo
       enddo

       ! remove redundant gridpoint indices
       ! (is there a more elegant way to do this?)
       rglobtemp(1) = rglob(1)
       k = 1
       do i = 2,nrec
          itemp = rglob(i)
          iflag = 0
          do j = 1,i
             if (rglobtemp(j) == itemp) iflag = 1
          enddo
          if (iflag == 0) then
             k = k+1
             rglobtemp(k) = itemp
          endif
       enddo

       ! update nrec and rglob
       nrec = k
       rglob(1:nrec) = rglobtemp(1:nrec)
       print *, nrec,' output gridpoint indices from set_glob.f90'

       ! NO NEED TO LOOP OVER i,j,k -- only iglob
!!$    do irec = 1, nrec
!!$       d_min_rec = sqrt(LENGTH**2+HEIGHT**2)  ! max possible distance
!!$       do ispec = 1,NSPEC
!!$          do j = 1,NGLLZ
!!$             do i = 1,NGLLX
!!$                iglob = ibool(i,j,ispec)
!!$                d = sqrt((x_rec(irec)-x(iglob))**2+(z_rec(irec)-z(iglob))**2)
!!$                if (d < d_min_rec) then
!!$                  d_min_rec  = d
!!$                  rglob(irec) = ibool(i,j,ispec)
!!$                endif
!!$             enddo
!!$          enddo
!!$       enddo
!!$    enddo

  endif

  end subroutine set_glob

!=====================================================================

!  subroutine station_filter(nrec, x_rec, z_rec, dmin_trsh, ncoast, coast_x, coast_z)
  subroutine station_filter(nrec, x_rec, z_rec, ifilter, dmin_trsh)

    integer, intent(inout) :: nrec
    double precision, intent(in) :: dmin_trsh, x_rec(nrec), z_rec(nrec)   ! NOTE INPUT ONLY
    integer, intent(out) :: ifilter(nrec)

    !integer, intent(in) :: ncoast
    !double precision, intent(in) :: dmin_trsh
    !double precision, intent(in) :: coast_x(ncoast), coast_z(ncoast)

    double precision :: dmin, d, xtar, ztar
    integer :: ispec, irec, i, j, iglob

    !---------------------------------------

    print *, nrec,' input target points into station_filter.f90 (some may be zeros)'
    print *, 'removing target points that are outside mesh or'
    print *, sngl(dmin_trsh/1000.),' km from grid boundary'
    !print *, 'removing target points that are ',sngl(dmin_trsh/1000.),' km from the SoCal coast'

    ! (1) exclude stations outside, or close to the boundary, of the mesh
    ! (2) exclude stations within a certain distance of the SoCal coast

    j = 0
    do irec = 1,nrec

       ! target point
       xtar = x_rec(irec)
       ztar = z_rec(irec)

       ! if target point is not near grid boundary
       if (xtar >= 0.+dmin_trsh      &
          .and. xtar <= LENGTH-dmin_trsh  &
          .and. ztar >= 0.+dmin_trsh      &
          .and. ztar <= HEIGHT-dmin_trsh) then

          j = j+1
          !x_rec(j) = xtar
          !z_rec(j) = ztar
          ifilter(j) = irec

          ! obtain the coastal point closest to the target point
          !dmin = sqrt(LENGTH**2+HEIGHT**2)
          !do i=1,ncoast
          !   ! distance from target point to coastal points
          !   d = sqrt((xtar-coast_x(i))**2+(ztar-coast_z(i))**2)
          !   if (d < dmin) then
          !      dmin  = d
          !   endif
          !enddo

          !! if target point is at least dmin_trsh from the coast, then keep it
          !if ( dmin >= dmin_trsh ) then
          !  j = j+1
          !  x_rec(j) = xtar
          !  z_rec(j) = ztar
          !endif

       endif
    enddo
    nrec = j
    print *, nrec,' output target points from station_filter.f90'

!!$    ! modified to exclude stations within a designated circle (CHT)
!!$    r    = 300.d+03
!!$    xcen = 50.d+03
!!$    zcen = 50.d+03
!!$    j = 0
!!$    do irec = 1, nrec
!!$       d = sqrt((x(rglob(irec)) - xcen)**2+(z(rglob(irec)) - zcen)**2)
!!$       if (d > r) then
!!$          j = j+1
!!$          rglob(j) = rglob(irec)
!!$       endif
!!$    enddo
!!$    nrec = j
!!$    print *, nrec,' output gridpoints'
!!$    print *

  end subroutine station_filter

!=====================================================================

  subroutine station_filter_2(nrec, x_rec, z_rec, ieast)

    integer, intent(inout) :: nrec
    double precision, intent(inout) :: x_rec(nrec), z_rec(nrec)
    integer, intent(in) :: ieast

    double precision :: xtar,ztar,m0,b0,mp,bp,x_div,z_div
    integer :: irec, j

    !---------------------------------------

    ! dividing line for stations
    m0 = HEIGHT/LENGTH    ! slope
    b0 = 0.d+03           ! y-intercept
    mp = -1./m0           ! perpendicular slope
    !ieast = -1             ! (1) for east points or (-1) for west points

    print *, nrec,' input target points into station_filter_2.f90'
    print *, 'removing a fraction of the target points'

    j = 0
    do irec = 1,nrec

       ! target point
       xtar = x_rec(irec)
       ztar = z_rec(irec)

       ! obtain the closest point on the dividing line to the target point
       bp       = ztar - mp*xtar
       x_div    = (bp - b0)/(m0 - mp)
       z_div    = mp*x_div + bp
       !write(*,'(4f16.6)') xtar/1000., ztar/1000., x_div/1000., z_div/1000

       ! if target point is at least dmin_trsh from the coast, then keep it
       if ( 0 < dble(ieast)*(xtar - x_div) ) then
         j = j+1
         x_rec(j) = xtar
         z_rec(j) = ztar
       endif
    enddo
    nrec = j
    print *, nrec,' output target points from station_filter_2.f90'

  end subroutine station_filter_2

!=====================================================================

  subroutine mesh_geo(npt,rlon,rlat,rx,rz,UTM_PROJECTION_ZONE,iway)

  integer, intent(in) :: npt, UTM_PROJECTION_ZONE, iway
  double precision, intent(inout) :: rx(npt),rz(npt),rlon(npt),rlat(npt)

  double precision :: xtemp,ztemp
  integer i

  if (iway == ILONLAT2MESH) then  ! lon-lat to mesh coordinates
     do i=1,npt
        call utm_geo(rlon(i),rlat(i),xtemp,ztemp,UTM_PROJECTION_ZONE,ILONGLAT2UTM)
        rx(i) = xtemp - utm_xmin
        rz(i) = ztemp - utm_zmin
     enddo

  else                           ! mesh coordinates to lon-lat
     do i=1,npt
        xtemp = rx(i) + utm_xmin
        ztemp = rz(i) + utm_zmin
        call utm_geo(rlon(i),rlat(i),xtemp,ztemp,UTM_PROJECTION_ZONE,IUTM2LONGLAT)
     enddo
  endif

  end subroutine mesh_geo

  !------------------------------------------

  subroutine utm_geo(rlon,rlat,rx,ry,UTM_PROJECTION_ZONE,iway)

!=====================================================================
!
!  UTM (Universal Transverse Mercator) projection from the USGS
!
!=====================================================================
!
! convert geodetic longitude and latitude to UTM, and back
! use iway = ILONGLAT2UTM for long/lat to UTM, IUTM2LONGLAT for UTM to lat/long
!
! NOTE: THIS APPEARS TO BE ONLY GOOD TO SINGLE PRECISION (25-JAN-2006)
!
!
!-----CAMx v2.03
!
!     UTM_GEO performs UTM to geodetic (long/lat) translation, and back.
!
!     This is a Fortran version of the BASIC program "Transverse Mercator
!     Conversion", Copyright 1986, Norman J. Berls (Stefan Musarra, 2/94)
!     Based on algorithm taken from "Map Projections Used by the USGS"
!     by John P. Snyder, Geological Survey Bulletin 1532, USDI.
!
!     Input/Output arguments:
!
!        rlon                  Longitude (deg, negative for West)
!        rlat                  Latitude (deg)
!        rx                    UTM easting (m)
!        ry                    UTM northing (m)
!        UTM_PROJECTION_ZONE  UTM zone
!        iway                  Conversion type
!                              ILONGLAT2UTM = geodetic to UTM
!                              IUTM2LONGLAT = UTM to geodetic
!

  integer UTM_PROJECTION_ZONE,iway
  double precision rx,ry,rlon,rlat

  double precision, parameter :: degrad=PI/180., raddeg=180./PI
  double precision, parameter :: semimaj=6378206.4d0, semimin=6356583.8d0 ! Clarke 1866 original uses this one
  !double precision, parameter :: semimaj=6378137.0d0, semimin=6356752.314245d0 ! WGS-84

  double precision, parameter :: scfa=.9996
  double precision, parameter :: north=0., east=500000.

  double precision e2,e4,e6,e8,ep2,xx,yy,dlat,dlon,zone,cm,cmr,delam
  double precision f1,f2,f3,f4,rm,rn,t,c,a,e1,u,rlat1,dlat1,c1,t1,rn1,r1,d
  double precision rx_save,ry_save,rlon_save,rlat_save

!---------------------------------------------------

  if (SUPPRESS_UTM_PROJECTION) then
    if (iway == ILONGLAT2UTM) then
      rx = rlon
      ry = rlat
    else
      rlon = rx
      rlat = ry
    endif
    return
  endif

  ! save original parameters
  rlon_save = rlon
  rlat_save = rlat
  rx_save = rx
  ry_save = ry

  ! define parameters of reference ellipsoid
  e2  = 1.0-(semimin/semimaj)**2.0
  e4  = e2*e2
  e6  = e2*e4
  e8  = e4*e4
  ep2 = e2/(1.-e2)

  if (iway == IUTM2LONGLAT) then
    xx = rx
    yy = ry
  else
    dlon = rlon
    dlat = rlat
  endif

  ! set zone parameters
  zone = dble(UTM_PROJECTION_ZONE)
  cm   = zone*6.0 - 183.0
  cmr  = cm*degrad

  if (iway == ILONGLAT2UTM) then  ! lat-lon to UTM

     rlon = degrad*dlon
     rlat = degrad*dlat

     delam = dlon - cm
     if (delam < -180.) delam = delam + 360.
     if (delam > 180.) delam = delam - 360.
     delam = delam*degrad

     f1 = (1. - e2/4. - 3.*e4/64. - 5.*e6/256)*rlat
     f2 = 3.*e2/8. + 3.*e4/32. + 45.*e6/1024.
     f2 = f2*sin(2.*rlat)
     f3 = 15.*e4/256. * 45.*e6/1024.  ! note: original used here contains an error, formula should be corrected to .. + 45.0*e6/1024.
     !f3 = 15.*e4/256. + 45.*e6/1024.  ! correct version
     f3 = f3*sin(4.*rlat)
     f4 = 35.*e6/3072.
     f4 = f4*sin(6.*rlat)
     rm  = semimaj*(f1 - f2 + f3 - f4)  ! mm in other code

     if (dlat == 90. .or. dlat == -90.) then
        xx = 0.
        yy = scfa*rm

     else
        rn = semimaj/sqrt(1. - e2*sin(rlat)**2)
        t = tan(rlat)**2
        c = ep2*cos(rlat)**2
        a = cos(rlat)*delam

        f1 = (1. - t + c)*a**3 / 6.
        f2 = 5. - 18.*t + t**2 + 72.*c - 58.*ep2
        f2 = f2*a**5 / 120.
        xx = scfa*rn*(a + f1 + f2)

        f1 = a**2 / 2.
        f2 = 5. - t + 9.*c + 4.*c**2
        f2 = f2*a**4 / 24.
        f3 = 61. - 58.*t + t**2 + 600.*c - 330.*ep2
        f3 = f3*a**6/720.
        yy = scfa*(rm + rn*tan(rlat)*(f1 + f2 + f3))

     endif

     ! apply false easting and northing
     xx = xx + east
     yy = yy + north

  else  ! UTM to lat-lon

     ! remove false easting and northing
     xx = xx - east
     yy = yy - north

     ! calculate the footpoint latitude
     e1 = sqrt(1. - e2)
     e1 = (1. - e1)/(1. + e1)
     rm = yy/scfa
     u = 1. - e2/4. - 3.*e4/64. - 5.*e6/256.
     u = rm / (semimaj*u)

     f1 = 3.*e1/2. - 27.*e1**3./32.
     f1 = f1*sin(2.*u)
     f2 = 21.*e1**2 / 16. - 55.*e1**4 / 32.
     f2 = f2*sin(4.*u)
     f3 = 151.*e1**3. / 96.
     f3 = f3*sin(6.*u)
     f4 = 1097.*e1**4. / 512.  ! added term (CHT)
     f4 = f4*sin(8.*u)         ! added term (CHT)

     rlat1 = u + f1 + f2 + f3 + f4
     dlat1 = rlat1*raddeg

     if (dlat1 >= 90. .or. dlat1 <= -90.) then
        dlat1 = dmin1(dlat1,dble(90.) )
        dlat1 = dmax1(dlat1,dble(-90.) )
        dlon = cm

     else
        c1 = ep2*cos(rlat1)**2.
        t1 = tan(rlat1)**2.
        f1 = 1. - e2*sin(rlat1)**2.
        rn1 = semimaj/sqrt(f1)
        r1 = semimaj*(1. - e2)/sqrt(f1**3)
        d = xx / (rn1*scfa)

        ! latitude
        f1 = rn1*tan(rlat1) / r1
        f2 = d**2 / 2.
        f3 = 5.*3.*t1 + 10.*c1 - 4.*c1**2 - 9.*ep2 ! note: original contains an error, should be 5. + 3.t1 .. instead of 5. * ..
        !f3 = 5. + 3.*t1 + 10.*c1 - 4.*c1**2 - 9.*ep2 ! correct version
        f3 = f3*d**2*d**2 / 24.
        f4 = 61. + 90.*t1 + 298.*c1 + 45.*t1**2. - 252.*ep2 - 3.*c1**2
        f4 = f4*(d**2)**3 / 720.
        rlat = rlat1 - f1*(f2 - f3 + f4)
        dlat = rlat*raddeg

        ! longitude
        f1 = 1. + 2.*t1 + c1
        f1 = f1*d**2*d/6.
        f2 = 5. - 2.*c1 + 28.*t1 - 3.*c1**2 + 8.*ep2 + 24.*t1**2
        f2 = f2*(d**2)**2*d/120.
        rlon = cmr + (d - f1 + f2)/cos(rlat1)
        dlon = rlon*raddeg
        if (dlon < -180.) dlon = dlon + 360.
        if (dlon > 180.) dlon = dlon - 360.
      endif

    endif

    if (iway == IUTM2LONGLAT) then
      rlon = dlon
      rlat = dlat
      rx = rx_save
      ry = ry_save
    else
      rx = xx
      ry = yy
      rlon = rlon_save
      rlat = rlat_save
    endif

  end subroutine utm_geo

!=====================================================================

end module wave2d_sub2
