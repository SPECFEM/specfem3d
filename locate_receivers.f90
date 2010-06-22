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
!---- locate_receivers finds the correct position of the receivers
!----
  subroutine locate_receivers(ibool,myrank,NSPEC_AB,NGLOB_AB, &
                 xstore,ystore,zstore,xigll,yigll,zigll,rec_filename, &
                 nrec,islice_selected_rec,ispec_selected_rec, &
                 xi_receiver,eta_receiver,gamma_receiver,station_name,network_name,nu, &
                 NPROC,utm_x_source,utm_y_source, &
                 UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION, &
                 iglob_is_surface_external_mesh,ispec_is_surface_external_mesh, &
                 num_free_surface_faces,free_surface_ispec,free_surface_ijk)

  implicit none

  include "constants.h"

  logical SUPPRESS_UTM_PROJECTION

  integer NPROC,UTM_PROJECTION_ZONE

  integer nrec,myrank

  integer NSPEC_AB,NGLOB_AB

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! arrays containing coordinates of the points
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: xstore,ystore,zstore

! for surface locating and normal computing with external mesh
  integer :: pt0_ix,pt0_iy,pt0_iz,pt1_ix,pt1_iy,pt1_iz,pt2_ix,pt2_iy,pt2_iz
  integer :: num_free_surface_faces
  real(kind=CUSTOM_REAL), dimension(3) :: u_vector,v_vector,w_vector
  logical, dimension(NGLOB_AB) :: iglob_is_surface_external_mesh
  logical, dimension(NSPEC_AB) :: ispec_is_surface_external_mesh
  integer, dimension(num_free_surface_faces) :: free_surface_ispec
  integer, dimension(3,NGLLSQUARE,num_free_surface_faces) :: free_surface_ijk

  integer, allocatable, dimension(:) :: ix_initial_guess,iy_initial_guess,iz_initial_guess

  integer iprocloop
  integer ios

  double precision,dimension(1) :: altitude_rec,distmin_ele
  double precision,dimension(4) :: elevation_node,dist_node
  double precision,dimension(NPROC) :: distmin_ele_all,elevation_all
  double precision, allocatable, dimension(:) :: x_target,y_target,z_target
  double precision, allocatable, dimension(:) :: horiz_dist
  double precision, allocatable, dimension(:) :: x_found,y_found,z_found
  double precision, allocatable, dimension(:,:) :: x_found_all,y_found_all,z_found_all

  integer irec
  integer i,j,k,ispec,iglob,iface,inode,imin,imax,jmin,jmax,kmin,kmax,igll,jgll,kgll
  integer iselected,jselected,iface_selected,iadjust,jadjust
  integer iproc(1)

  double precision utm_x_source,utm_y_source
  double precision dist
  double precision xi,eta,gamma,dx,dy,dz,dxi,deta,dgamma

! Gauss-Lobatto-Legendre points of integration
  double precision xigll(NGLLX)
  double precision yigll(NGLLY)
  double precision zigll(NGLLZ)

! input receiver file name
  character(len=*) rec_filename

! topology of the control points of the surface element
  integer iax,iay,iaz
  integer iaddx(NGNOD),iaddy(NGNOD),iaddz(NGNOD)

! coordinates of the control points of the surface element
  double precision xelm(NGNOD),yelm(NGNOD),zelm(NGNOD)

  integer iter_loop,ispec_iterate

  integer ia
  double precision x,y,z
  double precision xix,xiy,xiz
  double precision etax,etay,etaz
  double precision gammax,gammay,gammaz

! timer MPI
  double precision, external :: wtime
  double precision time_start,tCPU

! use dynamic allocation
  double precision, dimension(:), allocatable :: final_distance
  double precision, dimension(:,:), allocatable :: final_distance_all
  double precision distmin,final_distance_max

! receiver information
! timing information for the stations
! station information for writing the seismograms

  integer :: iglob_selected
  integer, dimension(nrec) :: islice_selected_rec,ispec_selected_rec
  double precision, dimension(nrec) :: xi_receiver,eta_receiver,gamma_receiver
  double precision, dimension(3,3,nrec) :: nu
  character(len=MAX_LENGTH_STATION_NAME), dimension(nrec) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(nrec) :: network_name

  integer, allocatable, dimension(:,:) :: ispec_selected_rec_all
  double precision, allocatable, dimension(:) :: stlat,stlon,stele,stbur,stutm_x,stutm_y,elevation
  double precision, allocatable, dimension(:,:) :: xi_receiver_all,eta_receiver_all,gamma_receiver_all
  double precision, allocatable, dimension(:,:,:,:) :: nu_all


  character(len=256) OUTPUT_FILES

! **************


! get MPI starting time
  time_start = wtime()

  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '********************'
    write(IMAIN,*) ' locating receivers'
    write(IMAIN,*) '********************'
    write(IMAIN,*)
  endif

! define topology of the control element
  call usual_hex_nodes(iaddx,iaddy,iaddz)

  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '*****************************************************************'
    write(IMAIN,'(1x,a,a,a)') 'reading receiver information from ', trim(rec_filename), ' file'
    write(IMAIN,*) '*****************************************************************'
  endif

! get number of stations from receiver file
  open(unit=1,file=trim(rec_filename),status='old',action='read',iostat=ios)
  if (ios /= 0) call exit_mpi(myrank,'error opening file '//trim(rec_filename))

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
  allocate(nu_all(3,3,nrec,0:NPROC-1))

! loop on all the stations
  do irec=1,nrec

    read(1,*,iostat=ios) station_name(irec),network_name(irec),stlat(irec),stlon(irec),stele(irec),stbur(irec)
    if (ios /= 0) call exit_mpi(myrank, 'Error reading station file '//trim(rec_filename))

! convert station location to UTM 
    call utm_geo(stlon(irec),stlat(irec),stutm_x(irec),stutm_y(irec),&
                UTM_PROJECTION_ZONE,ILONGLAT2UTM,SUPPRESS_UTM_PROJECTION)

! compute horizontal distance between source and receiver in km
    horiz_dist(irec) = dsqrt((stutm_y(irec)-utm_y_source)**2 + (stutm_x(irec)-utm_x_source)**2) / 1000.

! print some information about stations
    if(myrank == 0) &
        write(IMAIN,*) 'Station #',irec,': ',station_name(irec)(1:len_trim(station_name(irec))), &
                       '.',network_name(irec)(1:len_trim(network_name(irec))), &
                       '    horizontal distance:  ',sngl(horiz_dist(irec)),' km'

! get approximate topography elevation at source long/lat coordinates
!   set distance to huge initial value
    distmin = HUGEVAL
    if(num_free_surface_faces > 0) then
    iglob_selected = 1
! loop only on points inside the element
! exclude edges to ensure this point is not shared with other elements
        imin = 2
        imax = NGLLX - 1

        jmin = 2
        jmax = NGLLY - 1
    do iface=1,num_free_surface_faces
          do j=jmin,jmax
             do i=imin,imax

                ispec = free_surface_ispec(iface)
                igll = free_surface_ijk(1,(j-1)*NGLLY+i,iface)
                jgll = free_surface_ijk(2,(j-1)*NGLLY+i,iface)
                kgll = free_surface_ijk(3,(j-1)*NGLLY+i,iface)
                iglob = ibool(igll,jgll,kgll,ispec)

 !           keep this point if it is closer to the receiver
                dist = dsqrt((stutm_x(irec)-dble(xstore(iglob)))**2 + &
                     (stutm_y(irec)-dble(ystore(iglob)))**2)
                if(dist < distmin) then
                   distmin = dist
                   iglob_selected = iglob
                   iface_selected = iface
                   iselected = i
                   jselected = j
                   altitude_rec(1) = zstore(iglob_selected)
                endif
             enddo
          enddo
          ! end of loop on all the elements on the free surface
       end do
!  weighted mean at current point of topography elevation of the four closest nodes     
!  set distance to huge initial value
       distmin = HUGEVAL
       do j=jselected,jselected+1
          do i=iselected,iselected+1
             inode = 1
             do jadjust=0,1
                do iadjust= 0,1
                   ispec = free_surface_ispec(iface_selected)
                   igll = free_surface_ijk(1,(j-jadjust-1)*NGLLY+i-iadjust,iface_selected)
                   jgll = free_surface_ijk(2,(j-jadjust-1)*NGLLY+i-iadjust,iface_selected)
                   kgll = free_surface_ijk(3,(j-jadjust-1)*NGLLY+i-iadjust,iface_selected)
                   iglob = ibool(igll,jgll,kgll,ispec)

                   elevation_node(inode) = zstore(iglob)
                   dist_node(inode) = dsqrt((stutm_x(irec)-dble(xstore(iglob)))**2 + &
                        (stutm_y(irec)-dble(ystore(iglob)))**2)
                   inode = inode + 1
                end do
             end do
             dist = sum(dist_node)
             if(dist < distmin) then
                distmin = dist
                altitude_rec(1) = (dist_node(1)/dist)*elevation_node(1) + &
                     (dist_node(2)/dist)*elevation_node(2) + &
                     (dist_node(3)/dist)*elevation_node(3) + &
                     (dist_node(4)/dist)*elevation_node(4) 
             endif
          end do
       end do
    end if
!  MPI communications to determine the best slice
    distmin_ele(1)= distmin
    call gather_all_dp(distmin_ele,1,distmin_ele_all,1,NPROC)
    call gather_all_dp(altitude_rec,1,elevation_all,1,NPROC)
    if(myrank == 0) then
       iproc = minloc(distmin_ele_all)
       altitude_rec(1) = elevation_all(iproc(1))         
    end if
    call bcast_all_dp(altitude_rec,1)  
    elevation(irec) = altitude_rec(1)

! reset distance to huge initial value
  distmin=HUGEVAL

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


      x_target(irec) = stutm_x(irec)
      y_target(irec) = stutm_y(irec)
      z_target(irec) = elevation(irec) - stbur(irec)
      !z_target(irec) = stbur(irec)
      !if (myrank == 0) write(IOVTK,*) x_target(irec), y_target(irec), z_target(irec)

! examine top of the elements only (receivers always at the surface)
!      k = NGLLZ

      ispec_selected_rec(irec) = 0

      do ispec=1,NSPEC_AB

! define the interval in which we look for points
      if(FASTER_RECEIVERS_POINTS_ONLY) then
        imin = 1
        imax = NGLLX

        jmin = 1
        jmax = NGLLY

        kmin = 1
        kmax = NGLLZ

      else
! loop only on points inside the element
! exclude edges to ensure this point is not shared with other elements
        imin = 2
        imax = NGLLX - 1

        jmin = 2
        jmax = NGLLY - 1

        kmin = 2
        kmax = NGLLZ - 1
      endif

        do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

            iglob = ibool(i,j,k,ispec)

            if (.not. RECVS_CAN_BE_BURIED_EXT_MESH) then
              if ((.not. iglob_is_surface_external_mesh(iglob)) .or. (.not. ispec_is_surface_external_mesh(ispec))) then
                cycle
              endif
            endif

            dist = dsqrt((x_target(irec)-dble(xstore(iglob)))**2 &
                        +(y_target(irec)-dble(ystore(iglob)))**2 &
                        +(z_target(irec)-dble(zstore(iglob)))**2)

!           keep this point if it is closer to the receiver
            if(dist < distmin) then
              distmin = dist
              ispec_selected_rec(irec) = ispec
              ix_initial_guess(irec) = i
              iy_initial_guess(irec) = j
              iz_initial_guess(irec) = k

              xi_receiver(irec) = dble(ix_initial_guess(irec))
              eta_receiver(irec) = dble(iy_initial_guess(irec))
              gamma_receiver(irec) = dble(iz_initial_guess(irec))
              x_found(irec) = xstore(iglob)
              y_found(irec) = ystore(iglob)
              z_found(irec) = zstore(iglob)
            endif

          enddo
        enddo
       enddo

! compute final distance between asked and found (converted to km)
  final_distance(irec) = dsqrt((x_target(irec)-x_found(irec))**2 + &
    (y_target(irec)-y_found(irec))**2 + (z_target(irec)-z_found(irec))**2)
!      endif

! end of loop on all the spectral elements in current slice
      enddo

  if (ispec_selected_rec(irec) == 0) then
    final_distance(irec) = HUGEVAL
  endif

! get normal to the face of the hexaedra if receiver is on the surface
  if ((.not. RECVS_CAN_BE_BURIED_EXT_MESH) .and. &
       .not. (ispec_selected_rec(irec) == 0)) then
    pt0_ix = -1
    pt0_iy = -1
    pt0_iz = -1
    pt1_ix = -1
    pt1_iy = -1
    pt1_iz = -1
    pt2_ix = -1
    pt2_iy = -1
    pt2_iz = -1
! we get two vectors of the face (three points) to compute the normal
    if (ix_initial_guess(irec) == 1 .and. &
         iglob_is_surface_external_mesh(ibool(1,2,2,ispec_selected_rec(irec)))) then
      pt0_ix = 1
      pt0_iy = NGLLY
      pt0_iz = 1
      pt1_ix = 1
      pt1_iy = 1
      pt1_iz = 1
      pt2_ix = 1
      pt2_iy = NGLLY
      pt2_iz = NGLLZ
    endif
    if (ix_initial_guess(irec) == NGLLX .and. &
         iglob_is_surface_external_mesh(ibool(NGLLX,2,2,ispec_selected_rec(irec)))) then
      pt0_ix = NGLLX
      pt0_iy = 1
      pt0_iz = 1
      pt1_ix = NGLLX
      pt1_iy = NGLLY
      pt1_iz = 1
      pt2_ix = NGLLX
      pt2_iy = 1
      pt2_iz = NGLLZ
    endif
    if (iy_initial_guess(irec) == 1 .and. &
         iglob_is_surface_external_mesh(ibool(2,1,2,ispec_selected_rec(irec)))) then
      pt0_ix = 1
      pt0_iy = 1
      pt0_iz = 1
      pt1_ix = NGLLX
      pt1_iy = 1
      pt1_iz = 1
      pt2_ix = 1
      pt2_iy = 1
      pt2_iz = NGLLZ
    endif
    if (iy_initial_guess(irec) == NGLLY .and. &
         iglob_is_surface_external_mesh(ibool(2,NGLLY,2,ispec_selected_rec(irec)))) then
      pt0_ix = NGLLX
      pt0_iy = NGLLY
      pt0_iz = 1
      pt1_ix = 1
      pt1_iy = NGLLY
      pt1_iz = 1
      pt2_ix = NGLLX
      pt2_iy = NGLLY
      pt2_iz = NGLLZ
    endif
    if (iz_initial_guess(irec) == 1 .and. &
         iglob_is_surface_external_mesh(ibool(2,2,1,ispec_selected_rec(irec)))) then
      pt0_ix = NGLLX
      pt0_iy = 1
      pt0_iz = 1
      pt1_ix = 1
      pt1_iy = 1
      pt1_iz = 1
      pt2_ix = NGLLX
      pt2_iy = NGLLY
      pt2_iz = 1
    endif
    if (iz_initial_guess(irec) == NGLLZ .and. &
         iglob_is_surface_external_mesh(ibool(2,2,NGLLZ,ispec_selected_rec(irec)))) then
      pt0_ix = 1
      pt0_iy = 1
      pt0_iz = NGLLZ
      pt1_ix = NGLLX
      pt1_iy = 1
      pt1_iz = NGLLZ
      pt2_ix = 1
      pt2_iy = NGLLY
      pt2_iz = NGLLZ
    endif

    if (pt0_ix<0 .or.pt0_iy<0 .or. pt0_iz<0 .or. &
         pt1_ix<0 .or. pt1_iy<0 .or. pt1_iz<0 .or. &
         pt2_ix<0 .or. pt2_iy<0 .or. pt2_iz<0) then
       stop 'error in computing normal for receivers.'
    endif

    u_vector(1) = xstore(ibool(pt1_ix,pt1_iy,pt1_iz,ispec_selected_rec(irec))) &
         - xstore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected_rec(irec)))
    u_vector(2) = ystore(ibool(pt1_ix,pt1_iy,pt1_iz,ispec_selected_rec(irec))) &
         - ystore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected_rec(irec)))
    u_vector(3) = zstore(ibool(pt1_ix,pt1_iy,pt1_iz,ispec_selected_rec(irec))) &
         - zstore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected_rec(irec)))
    v_vector(1) = xstore(ibool(pt2_ix,pt2_iy,pt2_iz,ispec_selected_rec(irec))) &
         - xstore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected_rec(irec)))
    v_vector(2) = ystore(ibool(pt2_ix,pt2_iy,pt2_iz,ispec_selected_rec(irec))) &
         - ystore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected_rec(irec)))
    v_vector(3) = zstore(ibool(pt2_ix,pt2_iy,pt2_iz,ispec_selected_rec(irec))) &
         - zstore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected_rec(irec)))

! cross product
    w_vector(1) = u_vector(2)*v_vector(3) - u_vector(3)*v_vector(2)
    w_vector(2) = u_vector(3)*v_vector(1) - u_vector(1)*v_vector(3)
    w_vector(3) = u_vector(1)*v_vector(2) - u_vector(2)*v_vector(1)

! normalize vector w
    w_vector(:) = w_vector(:)/sqrt(w_vector(1)**2+w_vector(2)**2+w_vector(3)**2)

! build the two other vectors for a direct base: we normalize u, and v=w^u
    u_vector(:) = u_vector(:)/sqrt(u_vector(1)**2+u_vector(2)**2+u_vector(3)**2)
    v_vector(1) = w_vector(2)*u_vector(3) - w_vector(3)*u_vector(2)
    v_vector(2) = w_vector(3)*u_vector(1) - w_vector(1)*u_vector(3)
    v_vector(3) = w_vector(1)*u_vector(2) - w_vector(2)*u_vector(1)

! build rotation matrice nu for seismograms
    if (EXT_MESH_RECV_NORMAL) then
!     East (u)
      nu(1,1,irec) = u_vector(1)
      nu(1,2,irec) = v_vector(1)
      nu(1,3,irec) = w_vector(1)

!     North (v)
      nu(2,1,irec) = u_vector(2)
      nu(2,2,irec) = v_vector(2)
      nu(2,3,irec) = w_vector(2)

!     Vertical (w)
      nu(3,1,irec) = u_vector(3)
      nu(3,2,irec) = v_vector(3)
      nu(3,3,irec) = w_vector(3)
      else
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
      endif

  endif ! of if (.not. RECVS_CAN_BE_BURIED_EXT_MESH)

! end of loop on all the stations
  enddo

! close receiver file
  close(1)

! ****************************************
! find the best (xi,eta,gamma) for each receiver
! ****************************************

  if(.not. FASTER_RECEIVERS_POINTS_ONLY) then

! loop on all the receivers to iterate in that slice
    do irec = 1,nrec

        ispec_iterate = ispec_selected_rec(irec)

! use initial guess in xi and eta

        xi = xigll(ix_initial_guess(irec))
        eta = yigll(iy_initial_guess(irec))
        gamma = zigll(iz_initial_guess(irec))

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

  endif ! of if (.not. FASTER_RECEIVERS_POINTS_ONLY)

! synchronize all the processes to make sure all the estimates are available
  call sync_all()

! for MPI version, gather information from all the nodes
  ispec_selected_rec_all(:,:) = -1
  call gather_all_i(ispec_selected_rec,nrec,ispec_selected_rec_all,nrec,NPROC)
  call gather_all_dp(xi_receiver,nrec,xi_receiver_all,nrec,NPROC)
  call gather_all_dp(eta_receiver,nrec,eta_receiver_all,nrec,NPROC)
  call gather_all_dp(gamma_receiver,nrec,gamma_receiver_all,nrec,NPROC)
  call gather_all_dp(final_distance,nrec,final_distance_all,nrec,NPROC)
  call gather_all_dp(x_found,nrec,x_found_all,nrec,NPROC)
  call gather_all_dp(y_found,nrec,y_found_all,nrec,NPROC)
  call gather_all_dp(z_found,nrec,z_found_all,nrec,NPROC)
  call gather_all_dp(nu,3*3*nrec,nu_all,3*3*nrec,NPROC)

! this is executed by main process only
  if(myrank == 0) then

! check that the gather operation went well
    if(any(ispec_selected_rec_all(:,:) == -1)) call exit_MPI(myrank,'gather operation failed for receivers')

! MPI loop on all the results to determine the best slice
    islice_selected_rec(:) = -1
    do irec = 1,nrec
    distmin = HUGEVAL
    do iprocloop = 0,NPROC-1
      if(final_distance_all(irec,iprocloop) < distmin) then
        distmin = final_distance_all(irec,iprocloop)
        islice_selected_rec(irec) = iprocloop
        ispec_selected_rec(irec) = ispec_selected_rec_all(irec,iprocloop)
        xi_receiver(irec) = xi_receiver_all(irec,iprocloop)
        eta_receiver(irec) = eta_receiver_all(irec,iprocloop)
        gamma_receiver(irec) = gamma_receiver_all(irec,iprocloop)
        x_found(irec) = x_found_all(irec,iprocloop)
        y_found(irec) = y_found_all(irec,iprocloop)
        z_found(irec) = z_found_all(irec,iprocloop)
        nu(:,:,irec) = nu_all(:,:,irec,iprocloop)
      endif
    enddo
    final_distance(irec) = distmin
    enddo

    do irec=1,nrec

      write(IMAIN,*)
      write(IMAIN,*) 'station # ',irec,'    ',station_name(irec),network_name(irec)

      if(final_distance(irec) == HUGEVAL) call exit_MPI(myrank,'error locating receiver')

      write(IMAIN,*) '     original latitude: ',sngl(stlat(irec))
      write(IMAIN,*) '     original longitude: ',sngl(stlon(irec))
      if( SUPPRESS_UTM_PROJECTION ) then
        write(IMAIN,*) '     original x: ',sngl(stutm_x(irec))
        write(IMAIN,*) '     original y: ',sngl(stutm_y(irec))
      else
        write(IMAIN,*) '     original UTM x: ',sngl(stutm_x(irec))
        write(IMAIN,*) '     original UTM y: ',sngl(stutm_y(irec))      
      endif
        write(IMAIN,*) '     original depth: ',sngl(stbur(irec)),' m'  
      write(IMAIN,*) '     horizontal distance: ',sngl(horiz_dist(irec))
      write(IMAIN,*) '     target x, y, z: ',sngl(x_target(irec)),sngl(y_target(irec)),sngl(z_target(irec))

      write(IMAIN,*) '     closest estimate found: ',sngl(final_distance(irec)),' m away'
      write(IMAIN,*) '     in slice ',islice_selected_rec(irec),' in element ',ispec_selected_rec(irec)
      if(FASTER_RECEIVERS_POINTS_ONLY) then
        write(IMAIN,*) '     in point i,j,k = ',nint(xi_receiver(irec)),nint(eta_receiver(irec)),nint(gamma_receiver(irec))
        write(IMAIN,*) '     nu1 = ',nu(1,:,irec)
        write(IMAIN,*) '     nu2 = ',nu(2,:,irec)
        write(IMAIN,*) '     nu3 = ',nu(3,:,irec)
      else
        write(IMAIN,*) '     at coordinates: '
        write(IMAIN,*) '       xi    = ',xi_receiver(irec)
        write(IMAIN,*) '       eta   = ',eta_receiver(irec)
        write(IMAIN,*) '       gamma = ',gamma_receiver(irec)
      endif
      if( SUPPRESS_UTM_PROJECTION ) then
         write(IMAIN,*) '         x: ',x_found(irec)
         write(IMAIN,*) '         y: ',y_found(irec)
      else
         write(IMAIN,*) '     UTM x: ',x_found(irec)
         write(IMAIN,*) '     UTM y: ',y_found(irec)        
        endif
        write(IMAIN,*) '     depth: ',dabs(z_found(irec) - elevation(irec)),' m'
        write(IMAIN,*) '         z: ',z_found(irec)
        write(IMAIN,*)
      

! add warning if estimate is poor
! (usually means receiver outside the mesh given by the user)
      if(final_distance(irec) > 3000.d0) then
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
    if(final_distance_max > 1000.d0) then
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
    tCPU = wtime() - time_start
    write(IMAIN,*)
    write(IMAIN,*) 'Elapsed time for receiver detection in seconds = ',tCPU
    write(IMAIN,*)
    write(IMAIN,*) 'End of receiver detection - done'
    write(IMAIN,*)

  endif    ! end of section executed by main process only

! main process broadcasts the results to all the slices
  call bcast_all_i(islice_selected_rec,nrec)
  call bcast_all_i(ispec_selected_rec,nrec)
  call bcast_all_dp(xi_receiver,nrec)
  call bcast_all_dp(eta_receiver,nrec)
  call bcast_all_dp(gamma_receiver,nrec)
! synchronize all the processes to make sure everybody has finished
  call sync_all()

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

!=====================================================================


  subroutine station_filter(SUPPRESS_UTM_PROJECTION,UTM_PROJECTION_ZONE,myrank,filename,filtered_filename,nfilter, &
      LATITUDE_MIN, LATITUDE_MAX, LONGITUDE_MIN, LONGITUDE_MAX)

  implicit none

  include 'constants.h'

! input
  logical :: SUPPRESS_UTM_PROJECTION
  integer :: UTM_PROJECTION_ZONE
  integer :: myrank
  character(len=*) :: filename,filtered_filename
  double precision :: LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX

! output
  integer :: nfilter

  integer :: nrec, nrec_filtered, ios !, irec

  double precision :: stlat,stlon,stele,stbur,stutm_x,stutm_y
  character(len=MAX_LENGTH_STATION_NAME) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME) :: network_name
  character(len=256) :: dummystring

  nrec = 0
  nrec_filtered = 0

  ! counts number of lines in stations file
  open(unit=IIN, file=trim(filename), status = 'old', iostat = ios)
  if (ios /= 0) call exit_mpi(myrank, 'No file '//trim(filename)//', exit')
  do while(ios == 0)
    read(IIN,"(a256)",iostat = ios) dummystring
    if(ios /= 0) exit

    if( len_trim(dummystring) > 0 ) nrec = nrec + 1
  enddo
  close(IIN)

  ! reads in station locations
  open(unit=IIN, file=trim(filename), status = 'old', iostat = ios)
  !do irec = 1,nrec
  !    read(IIN,*) station_name,network_name,stlat,stlon,stele,stbur
  do while(ios == 0)
    read(IIN,"(a256)",iostat = ios) dummystring
    if( ios /= 0 ) exit

    ! counts number of stations in min/max region
    if( len_trim(dummystring) > 0 ) then
        dummystring = trim(dummystring)
        read(dummystring, *) station_name, network_name, stlat, stlon, stele, stbur
    
        ! convert station location to UTM 
        call utm_geo(stlon,stlat,stutm_x,stutm_y,&
             UTM_PROJECTION_ZONE,ILONGLAT2UTM,SUPPRESS_UTM_PROJECTION)

        ! counts stations within lon/lat region
        if( stutm_y >= LATITUDE_MIN .and. stutm_y <= LATITUDE_MAX .and. &
           stutm_x >= LONGITUDE_MIN .and. stutm_x <= LONGITUDE_MAX) &
          nrec_filtered = nrec_filtered + 1
     endif
  enddo
  close(IIN)

  ! writes out filtered stations file
  if (myrank == 0) then
    open(unit=IIN,file=trim(filename),status='old',action='read',iostat=ios)
    open(unit=IOUT,file=trim(filtered_filename),status='unknown')
    !write(IOUT,*) nrec_filtered
    !do irec = 1,nrec
    do while(ios == 0)
      read(IIN,"(a256)",iostat = ios) dummystring
      if( ios /= 0 ) exit

      !read(IIN,*) station_name,network_name,stlat,stlon,stele,stbur
      if( len_trim(dummystring) > 0 ) then
        dummystring = trim(dummystring)
        read(dummystring, *) station_name, network_name, stlat, stlon, stele, stbur
        
        ! convert station location to UTM 
        call utm_geo(stlon,stlat,stutm_x,stutm_y,&
             UTM_PROJECTION_ZONE,ILONGLAT2UTM,SUPPRESS_UTM_PROJECTION)

        if( stutm_y >= LATITUDE_MIN .and. stutm_y <= LATITUDE_MAX .and. &
           stutm_x >= LONGITUDE_MIN .and. stutm_x <= LONGITUDE_MAX) then
          write(IOUT,*) trim(station_name),' ',trim(network_name),' ',sngl(stlat), &
                       ' ',sngl(stlon), ' ',sngl(stele), ' ',sngl(stbur)
       endif
    end if
 enddo
    close(IIN)
    close(IOUT)

    write(IMAIN,*)
    write(IMAIN,*) 'there are ',nrec,' stations in file ', trim(filename)
    write(IMAIN,*) 'saving ',nrec_filtered,' stations inside the model in file ', trim(filtered_filename)
    write(IMAIN,*) 'excluding ',nrec - nrec_filtered,' stations located outside the model'
    write(IMAIN,*)

    if( nrec_filtered < 1 ) then    
      write(IMAIN,*) 'error filtered stations:'
      write(IMAIN,*) '  simulation needs at least 1 station but got ',nrec_filtered
      write(IMAIN,*) 
      write(IMAIN,*) '  check that stations in file '//trim(filename)//' are within'
      write(IMAIN,*) '    latitude min/max : ',LATITUDE_MIN,LATITUDE_MAX
      write(IMAIN,*) '    longitude min/max: ',LONGITUDE_MIN,LONGITUDE_MAX
      write(IMAIN,*) 
    endif

  endif

  nfilter = nrec_filtered
  
  end subroutine station_filter

