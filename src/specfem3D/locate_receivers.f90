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
  double precision :: utm_x_source,utm_y_source
  double precision, dimension(NDIM,NDIM,nrec),intent(out) :: nu

  ! local parameters
  double precision, allocatable, dimension(:) :: x_target,y_target,z_target
  double precision, allocatable, dimension(:) :: x_found,y_found,z_found
  integer :: irec

  ! timer MPI
  double precision, external :: wtime
  double precision :: time_start,tCPU

  ! use dynamic allocation
  double precision, dimension(:), allocatable :: final_distance
  double precision :: final_distance_max

  ! receiver information
  ! station information for writing the seismograms
  double precision, allocatable, dimension(:) :: stlat,stlon,stele,stbur,stutm_x,stutm_y,elevation

  integer :: ier

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

  ! compute typical size of elements
  ! gets mesh dimensions
  call check_mesh_distances(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                            x_min_glob,x_max_glob,y_min_glob,y_max_glob,z_min_glob,z_max_glob, &
                            elemsize_min_glob,elemsize_max_glob, &
                            distance_min_glob,distance_max_glob)


  ! opens STATIONS or STATIONS_ADJOINT file
  if (myrank == 0) then
    open(unit=IIN,file=trim(rec_filename),status='old',action='read',iostat=ier)
    if (ier /= 0) call exit_mpi(myrank,'error opening file '//trim(rec_filename))
  endif

  ! checks if station locations already available
  if (SU_FORMAT .and. (.not. INVERSE_FWI_FULL_PROBLEM) ) then
    ! checks if file with station infos located from previous run exists
    inquire(file=trim(OUTPUT_FILES)//'/SU_stations_info.bin',exist=SU_station_file_exists)
    if (SU_station_file_exists) then
      if (myrank == 0) then
        do irec=1,nrec
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

        allocate(x_found(nrec),y_found(nrec),z_found(nrec))
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
      call bcast_all_dp(nu,NDIM*NDIM*nrec)
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
           elevation(nrec), &
           x_target(nrec), &
           y_target(nrec), &
           z_target(nrec), &
           x_found(nrec), &
           y_found(nrec), &
           z_found(nrec), &
           final_distance(nrec), &
           idomain(nrec),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays for locating receivers'

  ! loop on all the stations to read the file
  if (myrank == 0) then
    do irec = 1,nrec
      read(IIN,*,iostat=ier) station_name(irec),network_name(irec),stlat(irec),stlon(irec),stele(irec),stbur(irec)
      if (ier /= 0) call exit_mpi(myrank, 'Error reading station file '//trim(rec_filename))
    enddo
  endif

  ! broadcast values to other slices
  call bcast_all_ch_array(station_name,nrec,MAX_LENGTH_STATION_NAME)
  call bcast_all_ch_array(network_name,nrec,MAX_LENGTH_NETWORK_NAME)
  call bcast_all_dp(stlat,nrec)
  call bcast_all_dp(stlon,nrec)
  call bcast_all_dp(stele,nrec)
  call bcast_all_dp(stbur,nrec)

  ! loop on all the stations to locate the stations
  do irec = 1,nrec

    ! get z target coordinate, depending on the topography
    call get_elevation_and_z_coordinate(stlon(irec),stlat(irec),stutm_x(irec),stutm_y(irec),z_target(irec), &
                                        elevation(irec),stbur(irec))
    x_target(irec) = stutm_x(irec)
    y_target(irec) = stutm_y(irec)

    call locate_point_in_mesh(x_target(irec), y_target(irec), z_target(irec), RECEIVERS_CAN_BE_BURIED, elemsize_max_glob, &
            ispec_selected_rec(irec), xi_receiver(irec), eta_receiver(irec), gamma_receiver(irec), &
            x_found(irec), y_found(irec), z_found(irec), idomain(irec),nu(:,:,irec))

    ! synchronize all the processes to make sure all the estimates are available
    call synchronize_all()

    call locate_MPI_slice_and_bcast_to_all(x_target(irec), y_target(irec), z_target(irec), &
                                           x_found(irec), y_found(irec), z_found(irec), &
                                           xi_receiver(irec), eta_receiver(irec), gamma_receiver(irec), &
                                           ispec_selected_rec(irec), islice_selected_rec(irec), &
                                           final_distance(irec), idomain(irec),nu(:,:,irec))

    ! user output progress
    if (myrank == 0 .and. nrec > 1000) then
      if (mod(irec,500) == 0) then
        write(IMAIN,*) '  located receivers ',irec,'out of',nrec
        call flush_IMAIN()
      endif
    endif

  enddo ! loop over stations

  ! close receiver file
  close(IIN)

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
        write(IMAIN,*) '     horizontal distance: ',dsqrt((stutm_y(irec)-utm_y_source)**2 &
                                                    + (stutm_x(irec)-utm_x_source)**2) / 1000.d0
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
        write(IMAIN,*) '     nu1 = ',nu(1,:,irec)
        write(IMAIN,*) '     nu2 = ',nu(2,:,irec)
        write(IMAIN,*) '     nu3 = ',nu(3,:,irec)

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

    do irec=1,nrec
      write(IOUT_SU,'(a32,a8,3f24.12)') station_name(irec),network_name(irec),x_found(irec),y_found(irec),z_found(irec)
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
        write(IOUT_SU) nu
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

  ! deallocate arrays
  deallocate(stlat)
  deallocate(stlon)
  deallocate(stele)
  deallocate(stbur)
  deallocate(stutm_x)
  deallocate(stutm_y)
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

  end subroutine locate_receivers

!-------------------------------------------------------------------------------------------------
! Remove stations located outside of the mesh
!-------------------------------------------------------------------------------------------------
  subroutine station_filter(filename,filtered_filename,nfilter)

  use constants
  use specfem_par, only: SUPPRESS_UTM_PROJECTION,myrank,xstore,ystore

  implicit none

  ! input
  character(len=*) :: filename,filtered_filename

  ! output
  integer :: nfilter

  ! local
  integer,dimension(1) :: nrec, nrec_filtered
  integer :: ier
  double precision :: stlat,stlon,stele,stbur,stutm_x,stutm_y
  double precision :: minlat,minlon,maxlat,maxlon
  character(len=MAX_LENGTH_STATION_NAME) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME) :: network_name
  character(len=MAX_STRING_LEN) :: dummystring
  real(kind=CUSTOM_REAL):: minl,maxl,min_all,max_all
  double precision :: LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX

  ! gets model dimensions
  minl = minval( xstore )
  maxl = maxval( xstore )
  call min_all_all_cr(minl,min_all)
  call max_all_all_cr(maxl,max_all)
  LONGITUDE_MIN = min_all
  LONGITUDE_MAX = max_all

  minl = minval( ystore )
  maxl = maxval( ystore )
  call min_all_all_cr(minl,min_all)
  call max_all_all_cr(maxl,max_all)
  LATITUDE_MIN = min_all
  LATITUDE_MAX = max_all

  ! initialization
  nrec = 0
  nrec_filtered = 0

  if (myrank == 0) then

    ! counts number of stations in stations file, filter them and output the list of active stations in STATIONS_FILTERED file
    open(unit=IIN, file=trim(filename), status = 'old', iostat = ier)
    if (ier /= 0) call exit_mpi(myrank, 'No file '//trim(filename)//', exit')
    open(unit=IOUT,file=trim(filtered_filename),status='unknown')
    do while (ier == 0)
      read(IIN,"(a)",iostat=ier) dummystring
      if (ier /= 0) exit

      if (len_trim(dummystring) > 0) then
        nrec(1) = nrec(1) + 1
        dummystring = trim(dummystring)
        read(dummystring, *) station_name, network_name, stlat, stlon, stele, stbur

        ! convert station location to UTM
        call utm_geo(stlon,stlat,stutm_x,stutm_y,ILONGLAT2UTM)

        ! counts stations within lon/lat region
        if (stutm_y >= LATITUDE_MIN .and. stutm_y <= LATITUDE_MAX .and. &
            stutm_x >= LONGITUDE_MIN .and. stutm_x <= LONGITUDE_MAX) then
          nrec_filtered(1) = nrec_filtered(1) + 1
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
    write(IMAIN,*) 'there are ',nrec(1),' stations in file ', trim(filename)
    write(IMAIN,*) 'saving ',nrec_filtered(1),' stations inside the model in file ', trim(filtered_filename)
    write(IMAIN,*) 'excluding ',nrec(1) - nrec_filtered(1),' stations located outside the model'
    write(IMAIN,*)

    if (nrec_filtered(1) < 1) then
      write(IMAIN,*) 'error filtered stations:'
      write(IMAIN,*) '  simulation needs at least 1 station but got ',nrec_filtered(1)
      write(IMAIN,*)
      write(IMAIN,*) '  check that stations in file '//trim(filename)//' are within'

      if (SUPPRESS_UTM_PROJECTION) then
        write(IMAIN,*) '    latitude min/max : ',LATITUDE_MIN,LATITUDE_MAX
        write(IMAIN,*) '    longitude min/max: ',LONGITUDE_MIN,LONGITUDE_MAX
      else
        ! convert edge locations from UTM back to lat/lon
        call utm_geo(minlon,minlat,LONGITUDE_MIN,LATITUDE_MIN,IUTM2LONGLAT)
        call utm_geo(maxlon,maxlat,LONGITUDE_MAX,LATITUDE_MAX,IUTM2LONGLAT)
        write(IMAIN,*) '    longitude min/max: ',minlon,maxlon
        write(IMAIN,*) '    latitude min/max : ',minlat,maxlat
        write(IMAIN,*) '    UTM x min/max: ',LONGITUDE_MIN,LONGITUDE_MAX
        write(IMAIN,*) '    UTM y min/max : ',LATITUDE_MIN,LATITUDE_MAX
      endif

      write(IMAIN,*)
    endif

  endif ! myrank == 0

  call bcast_all_i(nrec,1)
  call bcast_all_i(nrec_filtered,1)

  nfilter = nrec_filtered(1)

  end subroutine station_filter

!--------------------------------------------------------------------------------------------------------------------
! get z target coordinate, depending on the topography
!--------------------------------------------------------------------------------------------------------------------
  subroutine get_elevation_and_z_coordinate(lon,lat,utm_x,utm_y,z_target,elevation,bury)

  use constants
  use specfem_par, only: USE_SOURCES_RECEIVERS_Z,ibool,myrank,NSPEC_AB,NGLOB_AB, &
                         xstore,ystore,zstore,NPROC,num_free_surface_faces,free_surface_ispec,free_surface_ijk

  double precision,     intent(in)  :: lon,lat,utm_x,utm_y,bury
  double precision,     intent(out) :: z_target,elevation

  !local
  integer,dimension(1)              :: iproc
  double precision,dimension(1)     :: altitude_rec,distmin_ele
  double precision,dimension(NPROC) :: distmin_ele_all,elevation_all
  real(kind=CUSTOM_REAL)            :: xloc,yloc,loc_ele,loc_distmin

  ! convert station location to UTM
  call utm_geo(lon,lat,utm_x,utm_y,ILONGLAT2UTM)

  xloc = utm_x
  yloc = utm_y
  ! get approximate topography elevation at point long/lat coordinates
  call get_topo_elevation_free(xloc,yloc,loc_ele,loc_distmin, &
                               NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                               num_free_surface_faces,free_surface_ispec,free_surface_ijk)

  altitude_rec(1) = loc_ele
  distmin_ele(1)  = loc_distmin

  !  MPI communications to determine the best slice
  call gather_all_dp(distmin_ele,1,distmin_ele_all,1,NPROC)
  call gather_all_dp(altitude_rec,1,elevation_all,1,NPROC)

  if (myrank == 0) then
    iproc = minloc(distmin_ele_all)
    altitude_rec(1) = elevation_all(iproc(1))
  endif
  call bcast_all_dp(altitude_rec,1)
  elevation = altitude_rec(1)

  ! point's Z coordinate
  if (USE_SOURCES_RECEIVERS_Z) then
    ! alternative: burial depth is given as z value directly
    z_target = bury
  else
    ! burial depth read in file given in m
    z_target = elevation - bury
  endif

  end subroutine get_elevation_and_z_coordinate


