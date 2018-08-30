!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
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

!----
!---- locate_receivers finds the correct position of the receivers
!----
  subroutine locate_receivers(rec_filename,nrec,islice_selected_rec,ispec_selected_rec, &
                              xi_receiver,eta_receiver,gamma_receiver,station_name,network_name,nu, &
                              utm_x_source,utm_y_source)

  use constants

  use specfem_par, only: USE_SOURCES_RECEIVERS_Z,ibool,myrank,NSPEC_AB,NGLOB_AB, &
                         xstore,ystore,zstore,SUPPRESS_UTM_PROJECTION,INVERSE_FWI_FULL_PROBLEM, &
                         SU_FORMAT

  implicit none

  ! input receiver file name
  character(len=*) :: rec_filename

  ! receivers
  integer :: nrec
  integer, dimension(nrec),intent(out) :: islice_selected_rec,ispec_selected_rec
  double precision, dimension(nrec),intent(out) :: xi_receiver,eta_receiver,gamma_receiver
  character(len=MAX_LENGTH_STATION_NAME), dimension(nrec),intent(out) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(nrec),intent(out) :: network_name
  double precision, dimension(NDIM,NDIM,nrec),intent(out) :: nu
  double precision,intent(in) :: utm_x_source,utm_y_source

  ! local parameters
  double precision, allocatable, dimension(:) :: x_target,y_target,z_target
  double precision, allocatable, dimension(:) :: x_found,y_found,z_found

  integer :: irec,ier,i

  ! timer MPI
  double precision, external :: wtime
  double precision :: tstart,tCPU

  ! use dynamic allocation
  double precision, dimension(:), allocatable :: final_distance
  double precision :: final_distance_max

  ! receiver information
  ! station information for writing the seismograms
  double precision, allocatable, dimension(:) :: stlat,stlon,stele,stbur,stutm_x,stutm_y,elevation

  ! SU_FORMAT parameters
  double precision :: llat,llon,lele,lbur
  logical :: SU_station_file_exists

  ! domains
  integer, dimension(:),allocatable :: idomain

  ! location search
  real(kind=CUSTOM_REAL) :: distance_min_glob,distance_max_glob
  real(kind=CUSTOM_REAL) :: elemsize_min_glob,elemsize_max_glob
  real(kind=CUSTOM_REAL) :: x_min_glob,x_max_glob
  real(kind=CUSTOM_REAL) :: y_min_glob,y_max_glob
  real(kind=CUSTOM_REAL) :: z_min_glob,z_max_glob

  double precision :: x,y,z,x_new,y_new,z_new
  double precision :: xi,eta,gamma,final_distance_squared
  double precision, dimension(NDIM,NDIM) :: nu_found
  double precision, dimension(NDIM) :: nu_tmp ! to avoid intel compiler warning about temporary array in i/o routine

  integer :: ispec_found,idomain_found

  ! subset arrays
  double precision, dimension(NREC_SUBSET_MAX) :: xi_receiver_subset,eta_receiver_subset,gamma_receiver_subset
  double precision, dimension(NREC_SUBSET_MAX) :: x_found_subset,y_found_subset,z_found_subset
  double precision, dimension(NREC_SUBSET_MAX) :: final_distance_subset
  double precision, dimension(NDIM,NDIM,NREC_SUBSET_MAX) :: nu_subset
  integer, dimension(NREC_SUBSET_MAX) :: ispec_selected_rec_subset,idomain_subset
  integer :: nrec_subset_current_size,irec_in_this_subset,irec_already_done
  integer :: length_station_name,length_network_name
  integer, allocatable, dimension(:) :: station_duplet

  logical :: is_done_stations

  ! get MPI starting time
  tstart = wtime()

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '********************'
    write(IMAIN,*) ' locating receivers'
    write(IMAIN,*) '********************'
    write(IMAIN,*)
    write(IMAIN,'(1x,a,a,a)') 'reading receiver information from ', trim(rec_filename), ' file'
    write(IMAIN,*)
    if (USE_SOURCES_RECEIVERS_Z) then
      write(IMAIN,*) 'using sources/receivers Z:'
      write(IMAIN,*) '  (depth) becomes directly (z) coordinate'
    endif
    call flush_IMAIN()
  endif

  ! compute typical size of elements
  ! gets mesh dimensions
  call check_mesh_distances(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                            x_min_glob,x_max_glob,y_min_glob,y_max_glob,z_min_glob,z_max_glob, &
                            elemsize_min_glob,elemsize_max_glob, &
                            distance_min_glob,distance_max_glob)

  ! checks if station locations already available
  if (SU_FORMAT .and. (.not. INVERSE_FWI_FULL_PROBLEM) ) then
    call read_stations_from_previous_run(is_done_stations)
    ! check if done
    if (is_done_stations) return
  endif

  ! allocate memory for arrays using number of stations
  allocate(stlat(nrec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1955')
  allocate(stlon(nrec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1956')
  allocate(stele(nrec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1957')
  allocate(stbur(nrec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1958')
  allocate(stutm_x(nrec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1959')
  allocate(stutm_y(nrec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1960')
  allocate(elevation(nrec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1961')
  allocate(x_target(nrec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1962')
  allocate(y_target(nrec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1963')
  allocate(z_target(nrec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1964')
  allocate(x_found(nrec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1965')
  allocate(y_found(nrec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1966')
  allocate(z_found(nrec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1967')
  allocate(final_distance(nrec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1968')
  allocate(idomain(nrec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1969')
  if (ier /= 0) stop 'Error allocating arrays for locating receivers'

  ! loop on all the stations to read the file
  if (myrank == 0) then
    ! opens STATIONS or STATIONS_ADJOINT file
    open(unit=IIN,file=trim(rec_filename),status='old',action='read',iostat=ier)
    if (ier /= 0) call exit_mpi(myrank,'error opening file '//trim(rec_filename))
    ! reads all stations
    do irec = 1,nrec
      read(IIN,*,iostat=ier) station_name(irec),network_name(irec),stlat(irec),stlon(irec),stele(irec),stbur(irec)
      if (ier /= 0) call exit_mpi(myrank, 'Error reading station file '//trim(rec_filename))
    enddo
    ! close receiver file
    close(IIN)

    ! In case that the same station and network name appear twice (or more times) in the STATIONS
    ! file, problems occur, as two (or more) seismograms are written (with mode
    ! "append") to a file with same name. The philosophy here is to accept multiple
    ! appearances and to just add a count to the station name in this case.
    allocate(station_duplet(nrec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1970')
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating station_duplet array')
    station_duplet(:) = 0
    do irec = 1,nrec
      do i = 1,irec-1
        if ((station_name(irec) == station_name(i)) .and. (network_name(irec) == network_name(i))) then
            station_duplet(i) = station_duplet(i) + 1
            if (len_trim(station_name(irec)) <= MAX_LENGTH_STATION_NAME-3) then
              write(station_name(irec),"(a,'_',i2.2)") trim(station_name(irec)),station_duplet(i)+1
            else
              call exit_MPI(myrank,'Please increase MAX_LENGTH_STATION_NAME by at least 3 to name station duplets')
            endif
        endif
      enddo

      ! checks name lengths
      length_station_name = len_trim(station_name(irec))
      length_network_name = len_trim(network_name(irec))
      ! check that length conforms to standard
      if (length_station_name < 1 .or. length_station_name > MAX_LENGTH_STATION_NAME) then
        print *, 'Error: invalid station name ',trim(station_name(irec))
        call exit_MPI(myrank,'wrong length of station name')
      endif
      if (length_network_name < 1 .or. length_network_name > MAX_LENGTH_NETWORK_NAME) then
        print *, 'Error: invalid network name ',trim(network_name(irec))
        call exit_MPI(myrank,'wrong length of network name')
      endif

    enddo
    deallocate(station_duplet)

  endif

  ! broadcast values to other slices
  call bcast_all_ch_array(station_name,nrec,MAX_LENGTH_STATION_NAME)
  call bcast_all_ch_array(network_name,nrec,MAX_LENGTH_NETWORK_NAME)
  call bcast_all_dp(stlat,nrec)
  call bcast_all_dp(stlon,nrec)
  call bcast_all_dp(stele,nrec)
  call bcast_all_dp(stbur,nrec)

  ! determines target point locations (need to locate z coordinate of all receivers)
  ! note: we first locate all the target positions in the mesh to reduces the need of MPI communication
  call get_elevation_and_z_coordinate_all(nrec,stlon,stlat,stbur,stutm_x,stutm_y,elevation,x_target,y_target,z_target)

  ! note: we loop over subsets of receivers to fill first MPI buffers, thus reducing the MPI communication for each receiver
  !
  ! loop on all the stations to locate the stations
  do irec_already_done = 0, nrec, NREC_SUBSET_MAX

    ! the size of the subset can be the maximum size, or less (if we are in the last subset,
    ! or if there are fewer sources than the maximum size of a subset)
    nrec_subset_current_size = min(NREC_SUBSET_MAX, nrec - irec_already_done)

    ! initializes search results
    final_distance_subset(:) = HUGEVAL

    ! loop over the stations within this subset
    do irec_in_this_subset = 1,nrec_subset_current_size

      ! mapping from station number in current subset to real station number in all the subsets
      irec = irec_in_this_subset + irec_already_done

      x = x_target(irec)
      y = y_target(irec)
      z = z_target(irec)

      ! locates point in mesh
      call locate_point_in_mesh(x, y, z, &
                                RECEIVERS_CAN_BE_BURIED, elemsize_max_glob, &
                                ispec_found, xi, eta, gamma, &
                                x_new, y_new, z_new, &
                                idomain_found, nu_found, final_distance_squared)

      ispec_selected_rec_subset(irec_in_this_subset) = ispec_found

      x_found_subset(irec_in_this_subset) = x_new
      y_found_subset(irec_in_this_subset) = y_new
      z_found_subset(irec_in_this_subset) = z_new

      xi_receiver_subset(irec_in_this_subset) = xi
      eta_receiver_subset(irec_in_this_subset) = eta
      gamma_receiver_subset(irec_in_this_subset) = gamma

      idomain_subset(irec_in_this_subset) = idomain_found
      nu_subset(:,:,irec_in_this_subset) = nu_found(:,:)
      final_distance_subset(irec_in_this_subset) = sqrt(final_distance_squared)

      ! user output progress
      if (myrank == 0 .and. nrec > 1000) then
        if (mod(irec,500) == 0) then
          write(IMAIN,*) '  located receivers ',irec,'out of',nrec
          call flush_IMAIN()
        endif
      endif
    enddo ! loop over subset

    ! master process locates best location in all slices
    call locate_MPI_slice(nrec_subset_current_size,irec_already_done, &
                          ispec_selected_rec_subset, &
                          x_found_subset, y_found_subset, z_found_subset, &
                          xi_receiver_subset,eta_receiver_subset,gamma_receiver_subset, &
                          idomain_subset,nu_subset,final_distance_subset, &
                          nrec,ispec_selected_rec, islice_selected_rec, &
                          x_found,y_found,z_found, &
                          xi_receiver, eta_receiver, gamma_receiver, &
                          idomain,nu,final_distance)

  enddo ! loop over stations

  ! bcast from master process
  call bcast_all_i(islice_selected_rec,nrec)
  ! note: in principle, only islice must be updated on all slave processes, the ones containing the best location
  !       could have valid entries in all other arrays set before already.
  !       nevertheless, for convenience we broadcast all needed receiver arrays back to the slaves
  call bcast_all_i(idomain,nrec)
  call bcast_all_i(ispec_selected_rec,nrec)

  call bcast_all_dp(xi_receiver,nrec)
  call bcast_all_dp(eta_receiver,nrec)
  call bcast_all_dp(gamma_receiver,nrec)

  call bcast_all_dp(x_found,nrec)
  call bcast_all_dp(y_found,nrec)
  call bcast_all_dp(z_found,nrec)

  call bcast_all_dp(nu,NDIM*NDIM*nrec)
  call bcast_all_dp(final_distance,nrec)

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
        write(IMAIN,*) '     horizontal distance: ',sngl(dsqrt((stutm_y(irec)-utm_y_source)**2 &
                                                    + (stutm_x(irec)-utm_x_source)**2) / 1000.d0)
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

        write(IMAIN,*) '     rotation matrix: '
        nu_tmp(:) = nu(1,:,irec)
        write(IMAIN,*) '     nu1 = ',sngl(nu_tmp)
        nu_tmp(:) = nu(2,:,irec)
        write(IMAIN,*) '     nu2 = ',sngl(nu_tmp)
        nu_tmp(:) = nu(3,:,irec)
        write(IMAIN,*) '     nu3 = ',sngl(nu_tmp)

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
        if (final_distance(irec) > elemsize_max_glob) then
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
    if (final_distance_max > elemsize_max_glob) then
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
    ! writes station infos
    do irec=1,nrec
      write(IOUT_SU,'(a32,a8,3f24.12)') station_name(irec),network_name(irec),x_found(irec),y_found(irec),z_found(irec)
    enddo
    ! closes output file
    close(IOUT_SU)

    ! stores station infos for later runs
    if (SU_FORMAT) call write_stations_for_next_run()

    ! elapsed time since beginning of mesh generation
    tCPU = wtime() - tstart
    write(IMAIN,*)
    write(IMAIN,*) 'Elapsed time for receiver detection in seconds = ',tCPU
    write(IMAIN,*)
    write(IMAIN,*) 'End of receiver detection - done'
    write(IMAIN,*)
    call flush_IMAIN()
  endif    ! end of section executed by main process only

  ! deallocate arrays
  deallocate(stlat)
  deallocate(stlon)
  deallocate(stele)
  deallocate(stbur)
  deallocate(stutm_x)
  deallocate(stutm_y)
  deallocate(elevation)
  deallocate(x_target)
  deallocate(y_target)
  deallocate(z_target)
  deallocate(x_found)
  deallocate(y_found)
  deallocate(z_found)
  deallocate(final_distance)
  deallocate(idomain)

  ! synchronize all the processes to make sure everybody has finished
  call synchronize_all()

contains

!
!--------------------------------------------------------
!

    subroutine read_stations_from_previous_run(is_done_stations)

    implicit none
    logical, intent(out) :: is_done_stations

    ! initializes
    is_done_stations = .false.

    ! checks if file with station infos located from previous run exists
    inquire(file=trim(OUTPUT_FILES)//'/SU_stations_info.bin',exist=SU_station_file_exists)
    if (SU_station_file_exists) then
      if (myrank == 0) then
        ! opens STATIONS or STATIONS_ADJOINT file
        open(unit=IIN,file=trim(rec_filename),status='old',action='read',iostat=ier)
        if (ier /= 0) call exit_mpi(myrank,'error opening file '//trim(rec_filename))
        ! reads stations/network names
        do irec = 1,nrec
          read(IIN,*,iostat=ier) station_name(irec),network_name(irec),llat,llon,lele,lbur
          if (ier /= 0) call exit_mpi(myrank, 'Error reading station file '//trim(rec_filename))
        enddo
        close(IIN)
      endif
      call bcast_all_ch_array(station_name,nrec,MAX_LENGTH_STATION_NAME)
      call bcast_all_ch_array(network_name,nrec,MAX_LENGTH_NETWORK_NAME)

      ! master reads in available station information
      if (myrank == 0) then
        open(unit=IOUT_SU,file=trim(OUTPUT_FILES)//'/SU_stations_info.bin', &
              status='old',action='read',form='unformatted',iostat=ier)
        if (ier /= 0) call exit_mpi(myrank,'error opening file '//trim(rec_filename))

        write(IMAIN,*) 'station details from SU_stations_info.bin'
        call flush_IMAIN()

        allocate(x_found(nrec),y_found(nrec),z_found(nrec),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 1971')
        ! reads in station infos
        read(IOUT_SU) islice_selected_rec,ispec_selected_rec
        read(IOUT_SU) xi_receiver,eta_receiver,gamma_receiver
        read(IOUT_SU) x_found,y_found,z_found
        read(IOUT_SU) nu
        close(IOUT_SU)
        ! write the locations of stations, so that we can load them and write them to SU headers later
        open(unit=IOUT_SU,file=trim(OUTPUT_FILES)//'/output_list_stations.txt', &
              status='unknown',action='write',iostat=ier)
        if (ier /= 0) &
          call exit_mpi(myrank,'error opening file '//trim(OUTPUT_FILES)//'/output_list_stations.txt')

        do irec = 1,nrec
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
      call bcast_all_dp(nu,NDIM*NDIM*nrec)
      call synchronize_all()
      ! user output
      if (myrank == 0) then
        ! elapsed time since beginning of mesh generation
        tCPU = wtime() - tstart
        write(IMAIN,*)
        write(IMAIN,*) 'Elapsed time for receiver detection in seconds = ',tCPU
        write(IMAIN,*)
        write(IMAIN,*) 'End of receiver detection - done'
        write(IMAIN,*)
        call flush_IMAIN()
      endif

      ! everything done
      is_done_stations = .true.
    endif

    end subroutine read_stations_from_previous_run

!
!--------------------------------------------------------
!

    subroutine write_stations_for_next_run()

    implicit none

    open(unit=IOUT_SU,file=trim(OUTPUT_FILES)//'/SU_stations_info.bin', &
         status='unknown',action='write',form='unformatted',iostat=ier)
    if (ier == 0) then
      write(IOUT_SU) islice_selected_rec,ispec_selected_rec
      write(IOUT_SU) xi_receiver,eta_receiver,gamma_receiver
      write(IOUT_SU) x_found,y_found,z_found
      write(IOUT_SU) nu
      close(IOUT_SU)
    endif

    end subroutine write_stations_for_next_run


  end subroutine locate_receivers

