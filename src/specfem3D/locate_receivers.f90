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

!----
!---- locate_receivers finds the correct position of the receivers
!----
  subroutine locate_receivers(rec_filename,nrec,islice_selected_rec,ispec_selected_rec, &
                              xi_receiver,eta_receiver,gamma_receiver,station_name,network_name, &
                              utm_x_source,utm_y_source)

  use constants

  use specfem_par, only: USE_SOURCES_RECEIVERS_Z,ibool,myrank,NSPEC_AB,NGLOB_AB,NGNOD, &
                         xstore,ystore,zstore,xigll,yigll,zigll,NPROC,UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION, &
                         iglob_is_surface_external_mesh,ispec_is_surface_external_mesh, &
                         num_free_surface_faces,free_surface_ispec,free_surface_ijk,SU_FORMAT
  use specfem_par_acoustic, only: ispec_is_acoustic
  use specfem_par_elastic, only: ispec_is_elastic
  use specfem_par_poroelastic, only: ispec_is_poroelastic

  implicit none

  ! input receiver file name
  character(len=*) :: rec_filename

  ! receivers
  integer :: nrec
  integer, dimension(nrec),intent(out) :: islice_selected_rec,ispec_selected_rec
  double precision, dimension(nrec),intent(out) :: xi_receiver,eta_receiver,gamma_receiver
  character(len=MAX_LENGTH_STATION_NAME), dimension(nrec),intent(out) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(nrec),intent(out) :: network_name

  double precision :: utm_x_source,utm_y_source

  ! local parameters

  integer, allocatable, dimension(:) :: ix_initial_guess,iy_initial_guess,iz_initial_guess
  integer :: iprocloop

  double precision,dimension(1) :: altitude_rec,distmin_ele
  double precision,dimension(NPROC) :: distmin_ele_all,elevation_all

  real(kind=CUSTOM_REAL) :: xloc,yloc,loc_ele,loc_distmin

  double precision, allocatable, dimension(:) :: x_target,y_target,z_target

  double precision, allocatable, dimension(:) :: horiz_dist
  double precision, allocatable, dimension(:) :: x_found,y_found,z_found
  double precision :: dist_squared

  double precision :: xi,eta,gamma,dx,dy,dz,dxi,deta,dgamma
  double precision :: x,y,z
  double precision :: xix,xiy,xiz
  double precision :: etax,etay,etaz
  double precision :: gammax,gammay,gammaz

  ! coordinates of the control points of the surface element
  double precision :: xelm(NGNOD),yelm(NGNOD),zelm(NGNOD)

  integer :: irec
  integer :: i,j,k,ispec,iglob

  integer :: iproc(1)

  ! topology of the control points of the surface element
  integer :: iax,iay,iaz
  integer :: iaddx(NGNOD),iaddy(NGNOD),iaddz(NGNOD)

  integer :: iter_loop,ispec_iterate
  integer :: ia

  ! timer MPI
  double precision, external :: wtime
  double precision :: time_start,tCPU

  ! use dynamic allocation
  double precision, dimension(:), allocatable :: final_distance
  double precision :: distmin_squared,final_distance_max

  ! receiver information
  ! station information for writing the seismograms
  double precision, allocatable, dimension(:) :: stlat,stlon,stele,stbur,stutm_x,stutm_y,elevation
  double precision, allocatable, dimension(:) :: x_found_all,y_found_all,z_found_all
  double precision, dimension(:), allocatable :: final_distance_all
  double precision, allocatable, dimension(:) :: xi_receiver_all,eta_receiver_all,gamma_receiver_all
  integer, allocatable, dimension(:) :: ispec_selected_rec_all

  integer :: ier

  real(kind=CUSTOM_REAL) :: xmin,xmax,ymin,ymax,zmin,zmax

  ! SU_FORMAT parameters
  double precision :: llat,llon,lele,lbur
  logical :: SU_station_file_exists

  ! domains
  integer, dimension(:),allocatable :: idomain,idomain_all

  ! location search
  double precision :: typical_size_squared
  real(kind=CUSTOM_REAL) :: distance_min_glob,distance_max_glob
  real(kind=CUSTOM_REAL) :: elemsize_min_glob,elemsize_max_glob
  real(kind=CUSTOM_REAL) :: x_min_glob,x_max_glob
  real(kind=CUSTOM_REAL) :: y_min_glob,y_max_glob
  real(kind=CUSTOM_REAL) :: z_min_glob,z_max_glob

  ! get MPI starting time
  time_start = wtime()

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '********************'
    write(IMAIN,*) ' locating receivers'
    write(IMAIN,*) '********************'
    write(IMAIN,*)
    write(IMAIN,'(1x,a,a,a)') 'reading receiver information from ', trim(rec_filename), ' file'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! dimension of model in current proc
  xmin=minval(xstore(:));          xmax=maxval(xstore(:))
  ymin=minval(ystore(:));          ymax=maxval(ystore(:))
  zmin=minval(zstore(:));          zmax=maxval(zstore(:))

  ! define topology of the control element
  call usual_hex_nodes(NGNOD,iaddx,iaddy,iaddz)

  ! compute typical size of elements
  if (USE_DISTANCE_CRITERION_RECEIVERS) then
    ! gets mesh dimensions
    call check_mesh_distances(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                              x_min_glob,x_max_glob,y_min_glob,y_max_glob,z_min_glob,z_max_glob, &
                              elemsize_min_glob,elemsize_max_glob, &
                              distance_min_glob,distance_max_glob)

    ! sets typical element size for search
    typical_size_squared =  elemsize_max_glob

    ! use 10 times the distance as a criterion for source detection
    typical_size_squared = (10. * typical_size_squared)**2
  endif

  ! opens STATIONS or STATIONS_ADJOINT file
  open(unit=IIN,file=trim(rec_filename),status='old',action='read',iostat=ier)
  if (ier /= 0) call exit_mpi(myrank,'error opening file '//trim(rec_filename))

  ! checks if station locations already available
  if (SU_FORMAT) then
    ! checks if file with station infos located from previous run exists
    inquire(file=trim(OUTPUT_FILES)//'/SU_stations_info.bin',exist=SU_station_file_exists)
    if (SU_station_file_exists) then
      ! all processes read in stations names from STATIONS file
      do irec=1,nrec
        read(IIN,*,iostat=ier) station_name(irec),network_name(irec),llat,llon,lele,lbur
        if (ier /= 0) call exit_mpi(myrank, 'Error reading station file '//trim(rec_filename))
      enddo
      close(IIN)
      ! master reads in available station information
      if (myrank == 0) then
        open(unit=IOUT_SU,file=trim(OUTPUT_FILES)//'/SU_stations_info.bin', &
              status='old',action='read',form='unformatted',iostat=ier)
        if (ier /= 0) call exit_mpi(myrank,'error opening file '//trim(rec_filename))

        write(IMAIN,*) 'station details from SU_stations_info.bin'
        call flush_IMAIN()

        allocate(x_found(nrec),y_found(nrec),z_found(nrec))
        ! reads in station infos
        read(IOUT_SU) islice_selected_rec,ispec_selected_rec
        read(IOUT_SU) xi_receiver,eta_receiver,gamma_receiver
        read(IOUT_SU) x_found,y_found,z_found
        close(IOUT_SU)
        ! write the locations of stations, so that we can load them and write them to SU headers later
        open(unit=IOUT_SU,file=trim(OUTPUT_FILES)//'/output_list_stations.txt', &
              status='unknown',action='write',iostat=ier)
        if (ier /= 0) &
          call exit_mpi(myrank,'error opening file '//trim(OUTPUT_FILES)//'/output_list_stations.txt')

        do irec=1,nrec
          write(IOUT_SU,*) station_name(irec),network_name(irec),x_found(irec),y_found(irec),z_found(irec)
        enddo

        close(IOUT_SU)
        deallocate(x_found,y_found,z_found)
      endif
      ! main process broadcasts the results to all the slices
      call bcast_all_i(islice_selected_rec,nrec)
      call bcast_all_i(ispec_selected_rec,nrec)
      call bcast_all_dp(xi_receiver,nrec)
      call bcast_all_dp(eta_receiver,nrec)
      call bcast_all_dp(gamma_receiver,nrec)
      call synchronize_all()
      ! user output
      if (myrank == 0) then
        ! elapsed time since beginning of mesh generation
        tCPU = wtime() - time_start
        write(IMAIN,*)
        write(IMAIN,*) 'Elapsed time for receiver detection in seconds = ',tCPU
        write(IMAIN,*)
        write(IMAIN,*) 'End of receiver detection - done'
        write(IMAIN,*)
        call flush_IMAIN()
      endif
      ! everything done
      return
    endif
  endif

  ! allocate memory for arrays using number of stations
  allocate(stlat(nrec), &
           stlon(nrec), &
           stele(nrec), &
           stbur(nrec), &
           stutm_x(nrec), &
           stutm_y(nrec), &
           horiz_dist(nrec), &
           elevation(nrec), &
           ix_initial_guess(nrec), &
           iy_initial_guess(nrec), &
           iz_initial_guess(nrec), &
           x_target(nrec), &
           y_target(nrec), &
           z_target(nrec), &
           x_found(nrec), &
           y_found(nrec), &
           z_found(nrec), &
           final_distance(nrec),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays for locating receivers'

  ! for MPI collection
  allocate(ispec_selected_rec_all(nrec), &
           xi_receiver_all(nrec), &
           eta_receiver_all(nrec), &
           gamma_receiver_all(nrec), &
           x_found_all(nrec), &
           y_found_all(nrec), &
           z_found_all(nrec), &
           final_distance_all(nrec),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays for MPI locating receivers'

  allocate(idomain(nrec),idomain_all(nrec),stat=ier)
  if (ier /= 0) stop 'Error allocating idomain arrays'

  ! loop on all the stations
  do irec = 1,nrec

    read(IIN,*,iostat=ier) station_name(irec),network_name(irec), &
                           stlat(irec),stlon(irec),stele(irec),stbur(irec)

    if (ier /= 0) call exit_mpi(myrank, 'Error reading station file '//trim(rec_filename))

    ! convert station location to UTM
    call utm_geo(stlon(irec),stlat(irec),stutm_x(irec),stutm_y(irec), &
                UTM_PROJECTION_ZONE,ILONGLAT2UTM,SUPPRESS_UTM_PROJECTION)

    ! compute horizontal distance between source and receiver in km
    horiz_dist(irec) = dsqrt((stutm_y(irec)-utm_y_source)**2 &
                           + (stutm_x(irec)-utm_x_source)**2) / 1000.d0

    ! print some information about stations
    if (myrank == 0) then
      ! limits user output if too many receivers
      if (nrec < 1000 .and. (.not. SU_FORMAT)) then
        write(IMAIN,*) 'Station #',irec,': ', &
            network_name(irec)(1:len_trim(network_name(irec)))//'.'//station_name(irec)(1:len_trim(station_name(irec))), &
            '    horizontal distance:  ',sngl(horiz_dist(irec)),' km'
      endif
    endif

    ! get approximate topography elevation at source long/lat coordinates
    xloc = stutm_x(irec)
    yloc = stutm_y(irec)
    call get_topo_elevation_free(xloc,yloc,loc_ele,loc_distmin, &
                                 NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                                 num_free_surface_faces,free_surface_ispec,free_surface_ijk)
    altitude_rec(1) = loc_ele
    distmin_ele(1) = loc_distmin

    !  MPI communications to determine the best slice
    call gather_all_dp(distmin_ele,1,distmin_ele_all,1,NPROC)
    call gather_all_dp(altitude_rec,1,elevation_all,1,NPROC)

    if (myrank == 0) then
      iproc = minloc(distmin_ele_all)
      altitude_rec(1) = elevation_all(iproc(1))
    endif
    call bcast_all_dp(altitude_rec,1)
    elevation(irec) = altitude_rec(1)

    x_target(irec) = stutm_x(irec)
    y_target(irec) = stutm_y(irec)

    ! receiver's Z coordinate
    if (USE_SOURCES_RECEIVERS_Z) then
      ! alternative: burial depth is given as z value directly
      z_target(irec) = stbur(irec)
    else
      ! burial depth in STATIONS file given in m
      z_target(irec) = elevation(irec) - stbur(irec)
    endif

    ! reset distance to huge initial value
    distmin_squared = HUGEVAL

    ! determines closest GLL point
    ispec_selected_rec(irec) = 0

    do ispec = 1,NSPEC_AB

      ! exclude elements that are too far from target
      if (USE_DISTANCE_CRITERION_RECEIVERS) then
        iglob = ibool(MIDX,MIDY,MIDZ,ispec)
        dist_squared = (x_target(irec) - dble(xstore(iglob)))**2 &
                     + (y_target(irec) - dble(ystore(iglob)))**2 &
                     + (z_target(irec) - dble(zstore(iglob)))**2
        if (dist_squared > typical_size_squared) cycle
      endif

      do k = 2,NGLLZ - 1
        do j = 2,NGLLY - 1
          do i = 2,NGLLX - 1

            iglob = ibool(i,j,k,ispec)

            if (.not. RECEIVERS_CAN_BE_BURIED) then
              if ((.not. iglob_is_surface_external_mesh(iglob)) &
                .or. (.not. ispec_is_surface_external_mesh(ispec))) then
                cycle
              endif
            endif

            !  we compare squared distances instead of distances themselves to significantly speed up calculations
            dist_squared = (x_target(irec)-dble(xstore(iglob)))**2 &
                         + (y_target(irec)-dble(ystore(iglob)))**2 &
                         + (z_target(irec)-dble(zstore(iglob)))**2

            ! keep this point if it is closer to the receiver
            if (dist_squared < distmin_squared) then
              distmin_squared = dist_squared
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

    ! end of loop on all the spectral elements in current slice
    enddo

    if (ispec_selected_rec(irec) == 0) then
      ! receiver is NOT within this proc, assign trivial values
      ispec_selected_rec(irec) = 1
      ix_initial_guess(irec) = 1
      iy_initial_guess(irec) = 1
      iz_initial_guess(irec) = 1
    endif

    ! sets whether acoustic (1) or elastic (2)
    if (ispec_is_acoustic( ispec_selected_rec(irec) )) then
      idomain(irec) = IDOMAIN_ACOUSTIC
    else if (ispec_is_elastic( ispec_selected_rec(irec) )) then
      idomain(irec) = IDOMAIN_ELASTIC
    else if (ispec_is_poroelastic( ispec_selected_rec(irec) )) then
      idomain(irec) = IDOMAIN_POROELASTIC
    else
      idomain(irec) = 0
    endif

  ! *****************************
  ! find the best (xi,eta,gamma)
  ! *****************************

    ispec_iterate = ispec_selected_rec(irec)

    ! use initial guess in xi and eta
    xi = xigll(ix_initial_guess(irec))
    eta = yigll(iy_initial_guess(irec))
    gamma = zigll(iz_initial_guess(irec))

    ! define coordinates of the control points of the element
    do ia=1,NGNOD
      iax = 0
      iay = 0
      iaz = 0
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
              xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,NGNOD)

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
      ! this can be useful for convergence of itertive scheme with distorted elements.
      ! this is purposely set to 1.10, do *NOT* set it back to 1.00 because 1.10 gives more accurate locations
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
      xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,NGNOD)

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

  enddo ! loop over stations

  ! close receiver file
  close(IIN)

  ! synchronize all the processes to make sure all the estimates are available
  call synchronize_all()

  ! for MPI version, gather information from all the nodes
  if (myrank /= 0) then
    ! gather information from other processors (one at a time)
    call send_i(ispec_selected_rec, nrec,0,0)
    call send_dp(xi_receiver,       nrec,0,1)
    call send_dp(eta_receiver,      nrec,0,2)
    call send_dp(gamma_receiver,    nrec,0,3)
    call send_dp(final_distance,    nrec,0,4)
    call send_dp(x_found,           nrec,0,5)
    call send_dp(y_found,           nrec,0,6)
    call send_dp(z_found,           nrec,0,7)
    call send_i(idomain,            nrec,0,8)

  else
    ! master collects
    islice_selected_rec(:) = 0
    do iprocloop = 1,NPROC-1
      call recv_i(ispec_selected_rec_all, nrec,iprocloop,0)
      call recv_dp(xi_receiver_all,       nrec,iprocloop,1)
      call recv_dp(eta_receiver_all,      nrec,iprocloop,2)
      call recv_dp(gamma_receiver_all,    nrec,iprocloop,3)
      call recv_dp(final_distance_all,    nrec,iprocloop,4)
      call recv_dp(x_found_all,           nrec,iprocloop,5)
      call recv_dp(y_found_all,           nrec,iprocloop,6)
      call recv_dp(z_found_all,           nrec,iprocloop,7)
      call recv_i(idomain_all,            nrec,iprocloop,8)

      ! determines final locations
      do irec = 1,nrec
        if (final_distance_all(irec) < final_distance(irec)) then
          final_distance(irec) = final_distance_all(irec)
          islice_selected_rec(irec) = iprocloop
          ispec_selected_rec(irec) = ispec_selected_rec_all(irec)
          xi_receiver(irec) = xi_receiver_all(irec)
          eta_receiver(irec) = eta_receiver_all(irec)
          gamma_receiver(irec) = gamma_receiver_all(irec)
          x_found(irec) = x_found_all(irec)
          y_found(irec) = y_found_all(irec)
          z_found(irec) = z_found_all(irec)
          idomain(irec) = idomain_all(irec)
        endif
      enddo
    enddo
  endif

  ! this is executed by main process only
  if (myrank == 0) then

    do irec = 1,nrec

      ! checks stations location
      if (final_distance(irec) == HUGEVAL) then
        write(IMAIN,*) 'error locating station # ',irec,'    ',trim(network_name(irec)),'    ',trim(station_name(irec))
        call exit_MPI(myrank,'error locating receiver')
      endif

      ! limits user output if too many receivers
      if (nrec < 1000 .and. (.not. SU_FORMAT )) then

        ! receiver info
        write(IMAIN,*)
        write(IMAIN,*) 'station # ',irec,'    ',trim(network_name(irec)),'    ',trim(station_name(irec))

        ! location info
        write(IMAIN,*) '     original latitude: ',sngl(stlat(irec))
        write(IMAIN,*) '     original longitude: ',sngl(stlon(irec))
        if (SUPPRESS_UTM_PROJECTION) then
          write(IMAIN,*) '     original x: ',sngl(stutm_x(irec))
          write(IMAIN,*) '     original y: ',sngl(stutm_y(irec))
        else
          write(IMAIN,*) '     original UTM x: ',sngl(stutm_x(irec))
          write(IMAIN,*) '     original UTM y: ',sngl(stutm_y(irec))
        endif
        if (USE_SOURCES_RECEIVERS_Z) then
          write(IMAIN,*) '     original z: ',sngl(stbur(irec))
        else
          write(IMAIN,*) '     original depth: ',sngl(stbur(irec)),' m'
        endif
        write(IMAIN,*) '     horizontal distance: ',sngl(horiz_dist(irec))
        write(IMAIN,*) '     target x, y, z: ',sngl(x_target(irec)),sngl(y_target(irec)),sngl(z_target(irec))

        write(IMAIN,*) '     closest estimate found: ',sngl(final_distance(irec)),' m away'
        write(IMAIN,*)

        write(IMAIN,*) '     receiver located in slice ',islice_selected_rec(irec)
        write(IMAIN,*) '                      in element ',ispec_selected_rec(irec)
        if (idomain(irec) == IDOMAIN_ACOUSTIC) then
          write(IMAIN,*) '                      in acoustic domain'
        else if (idomain(irec) == IDOMAIN_ELASTIC) then
          write(IMAIN,*) '                      in elastic domain'
        else if (idomain(irec) == IDOMAIN_POROELASTIC) then
          write(IMAIN,*) '                      in poroelastic domain'
        else
          write(IMAIN,*) '                      in unknown domain'
        endif
        write(IMAIN,*) '     at coordinates: '
        write(IMAIN,*) '     xi    = ',xi_receiver(irec)
        write(IMAIN,*) '     eta   = ',eta_receiver(irec)
        write(IMAIN,*) '     gamma = ',gamma_receiver(irec)
        if (SUPPRESS_UTM_PROJECTION) then
          write(IMAIN,*) '     x: ',x_found(irec)
          write(IMAIN,*) '     y: ',y_found(irec)
        else
          write(IMAIN,*) '     UTM x: ',x_found(irec)
          write(IMAIN,*) '     UTM y: ',y_found(irec)
        endif
        if (USE_SOURCES_RECEIVERS_Z) then
          write(IMAIN,*) '     z: ',z_found(irec)
        else
          write(IMAIN,*) '     depth: ',dabs(z_found(irec) - elevation(irec)),' m'
          write(IMAIN,*) '     z: ',z_found(irec)
        endif
        write(IMAIN,*)

        ! add warning if estimate is poor
        ! (usually means receiver outside the mesh given by the user)
!! DK DK warning: this should be made a relative distance rather than absolute, now that we use the code at many different scales
!! DK DK warning: this should be made a relative distance rather than absolute, now that we use the code at many different scales
!! DK DK warning: this should be made a relative distance rather than absolute, now that we use the code at many different scales
        if (final_distance(irec) > 3000.d0) then
          write(IMAIN,*) '*******************************************************'
          write(IMAIN,*) '***** WARNING: receiver location estimate is poor *****'
          write(IMAIN,*) '*******************************************************'
        endif

        write(IMAIN,*)
      endif

    enddo

    ! compute maximal distance for all the receivers
    final_distance_max = maxval(final_distance(:))

    ! display maximum error for all the receivers
    write(IMAIN,*) 'maximum error in location of all the receivers: ',sngl(final_distance_max),' m'

    ! add warning if estimate is poor
    ! (usually means receiver outside the mesh given by the user)
!! DK DK warning: this should be made a relative distance rather than absolute, now that we use the code at many different scales
!! DK DK warning: this should be made a relative distance rather than absolute, now that we use the code at many different scales
!! DK DK warning: this should be made a relative distance rather than absolute, now that we use the code at many different scales
    if (final_distance_max > 1000.d0) then
      write(IMAIN,*)
      write(IMAIN,*) '************************************************************'
      write(IMAIN,*) '************************************************************'
      write(IMAIN,*) '***** WARNING: at least one receiver is poorly located *****'
      write(IMAIN,*) '************************************************************'
      write(IMAIN,*) '************************************************************'
    endif

    ! write the locations of stations, so that we can load them and write them to SU headers later
    open(unit=IOUT_SU,file=trim(OUTPUT_FILES)//'/output_list_stations.txt', &
         status='unknown',action='write',iostat=ier)
    if (ier /= 0) &
      call exit_mpi(myrank,'error opening file '//trim(OUTPUT_FILES)//'/output_list_stations.txt')

    do irec=1,nrec
      write(IOUT_SU,*) station_name(irec),network_name(irec),x_found(irec),y_found(irec),z_found(irec)
    enddo

    close(IOUT_SU)

    ! stores station infos for later runs
    if (SU_FORMAT) then
      open(unit=IOUT_SU,file=trim(OUTPUT_FILES)//'/SU_stations_info.bin', &
           status='unknown',action='write',form='unformatted',iostat=ier)
      if (ier == 0) then
        write(IOUT_SU) islice_selected_rec,ispec_selected_rec
        write(IOUT_SU) xi_receiver,eta_receiver,gamma_receiver
        write(IOUT_SU) x_found,y_found,z_found
        close(IOUT_SU)
      endif
    endif

    ! elapsed time since beginning of mesh generation
    tCPU = wtime() - time_start
    write(IMAIN,*)
    write(IMAIN,*) 'Elapsed time for receiver detection in seconds = ',tCPU
    write(IMAIN,*)
    write(IMAIN,*) 'End of receiver detection - done'
    write(IMAIN,*)
    call flush_IMAIN()
  endif    ! end of section executed by main process only

  ! main process broadcasts the results to all the slices
  call bcast_all_i(islice_selected_rec,nrec)
  call bcast_all_i(ispec_selected_rec,nrec)
  call bcast_all_dp(xi_receiver,nrec)
  call bcast_all_dp(eta_receiver,nrec)
  call bcast_all_dp(gamma_receiver,nrec)

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

  ! synchronize all the processes to make sure everybody has finished
  call synchronize_all()

  end subroutine locate_receivers

!
!-------------------------------------------------------------------------------------------------
!

  subroutine station_filter(SUPPRESS_UTM_PROJECTION,UTM_PROJECTION_ZONE, &
                            myrank,filename,filtered_filename,nfilter, &
                            LATITUDE_MIN, LATITUDE_MAX, LONGITUDE_MIN, LONGITUDE_MAX)

  use constants

  implicit none

! input
  logical :: SUPPRESS_UTM_PROJECTION
  integer :: UTM_PROJECTION_ZONE
  integer :: myrank
  character(len=*) :: filename,filtered_filename
  double precision :: LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX

! output
  integer :: nfilter

  integer :: nrec, nrec_filtered
  integer :: ier
  double precision :: stlat,stlon,stele,stbur,stutm_x,stutm_y
  double precision :: minlat,minlon,maxlat,maxlon
  character(len=MAX_LENGTH_STATION_NAME) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME) :: network_name
  character(len=MAX_STRING_LEN) :: dummystring

  nrec = 0
  nrec_filtered = 0

  ! counts number of lines in stations file
  open(unit=IIN, file=trim(filename), status = 'old', iostat = ier)
  if (ier /= 0) call exit_mpi(myrank, 'No file '//trim(filename)//', exit')
  do while (ier == 0)
    read(IIN,"(a)",iostat=ier) dummystring
    if (ier /= 0) exit

    if (len_trim(dummystring) > 0) nrec = nrec + 1
  enddo
  close(IIN)

  ! reads in station locations
  open(unit=IIN, file=trim(filename), status = 'old', iostat = ier)
  do while (ier == 0)
    read(IIN,"(a)",iostat=ier) dummystring
    if (ier /= 0) exit

    ! counts number of stations in min/max region
    if (len_trim(dummystring) > 0) then
      dummystring = trim(dummystring)
      read(dummystring, *) station_name, network_name, stlat, stlon, stele, stbur

      ! convert station location to UTM
      call utm_geo(stlon,stlat,stutm_x,stutm_y, &
           UTM_PROJECTION_ZONE,ILONGLAT2UTM,SUPPRESS_UTM_PROJECTION)

      ! counts stations within lon/lat region
      if (stutm_y >= LATITUDE_MIN .and. stutm_y <= LATITUDE_MAX .and. &
          stutm_x >= LONGITUDE_MIN .and. stutm_x <= LONGITUDE_MAX) &
        nrec_filtered = nrec_filtered + 1
    endif
  enddo
  close(IIN)

  ! writes out filtered stations file
  if (myrank == 0) then
    open(unit=IIN,file=trim(filename),status='old',action='read',iostat=ier)
    open(unit=IOUT,file=trim(filtered_filename),status='unknown')
    do while (ier == 0)
      read(IIN,"(a)",iostat=ier) dummystring
      if (ier /= 0) exit

      if (len_trim(dummystring) > 0) then
        dummystring = trim(dummystring)
        read(dummystring, *) station_name, network_name, stlat, stlon, stele, stbur

        ! convert station location to UTM
        call utm_geo(stlon,stlat,stutm_x,stutm_y, &
             UTM_PROJECTION_ZONE,ILONGLAT2UTM,SUPPRESS_UTM_PROJECTION)

        if (stutm_y >= LATITUDE_MIN .and. stutm_y <= LATITUDE_MAX .and. &
           stutm_x >= LONGITUDE_MIN .and. stutm_x <= LONGITUDE_MAX) then

          ! with specific format
          write(IOUT,'(a10,1x,a10,4e18.6)') &
                            trim(station_name),trim(network_name), &
                            sngl(stlat),sngl(stlon),sngl(stele),sngl(stbur)

        endif
      endif
    enddo
    close(IIN)
    close(IOUT)

    write(IMAIN,*)
    write(IMAIN,*) 'there are ',nrec,' stations in file ', trim(filename)
    write(IMAIN,*) 'saving ',nrec_filtered,' stations inside the model in file ', trim(filtered_filename)
    write(IMAIN,*) 'excluding ',nrec - nrec_filtered,' stations located outside the model'
    write(IMAIN,*)

    if (nrec_filtered < 1) then
      write(IMAIN,*) 'error filtered stations:'
      write(IMAIN,*) '  simulation needs at least 1 station but got ',nrec_filtered
      write(IMAIN,*)
      write(IMAIN,*) '  check that stations in file '//trim(filename)//' are within'

      if (SUPPRESS_UTM_PROJECTION) then
        write(IMAIN,*) '    latitude min/max : ',LATITUDE_MIN,LATITUDE_MAX
        write(IMAIN,*) '    longitude min/max: ',LONGITUDE_MIN,LONGITUDE_MAX
      else
        ! convert edge locations from UTM back to lat/lon
        call utm_geo(minlon,minlat,LONGITUDE_MIN,LATITUDE_MIN, &
             UTM_PROJECTION_ZONE,IUTM2LONGLAT,SUPPRESS_UTM_PROJECTION)
        call utm_geo(maxlon,maxlat,LONGITUDE_MAX,LATITUDE_MAX, &
             UTM_PROJECTION_ZONE,IUTM2LONGLAT,SUPPRESS_UTM_PROJECTION)
        write(IMAIN,*) '    longitude min/max: ',minlon,maxlon
        write(IMAIN,*) '    latitude min/max : ',minlat,maxlat
        write(IMAIN,*) '    UTM x min/max: ',LONGITUDE_MIN,LONGITUDE_MAX
        write(IMAIN,*) '    UTM y min/max : ',LATITUDE_MIN,LATITUDE_MAX
      endif

      write(IMAIN,*)
    endif

  endif

  nfilter = nrec_filtered

  end subroutine station_filter

