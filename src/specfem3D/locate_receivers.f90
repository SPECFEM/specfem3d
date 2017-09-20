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

  real(kind=CUSTOM_REAL) :: xmin,xmax,ymin,ymax,zmin,zmax

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

  ! dimension of model in current proc
  xmin=minval(xstore(:));          xmax=maxval(xstore(:))
  ymin=minval(ystore(:));          ymax=maxval(ystore(:))
  zmin=minval(zstore(:));          zmax=maxval(zstore(:))

  ! compute typical size of elements
  ! gets mesh dimensions
  call check_mesh_distances(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                            x_min_glob,x_max_glob,y_min_glob,y_max_glob,z_min_glob,z_max_glob, &
                            elemsize_min_glob,elemsize_max_glob, &
                            distance_min_glob,distance_max_glob)


  ! opens STATIONS or STATIONS_ADJOINT file
  if (myrank==0) then
    open(unit=IIN,file=trim(rec_filename),status='old',action='read',iostat=ier)
    if (ier /= 0) call exit_mpi(myrank,'error opening file '//trim(rec_filename))
  endif

  ! checks if station locations already available
  if (SU_FORMAT .and. (.not. INVERSE_FWI_FULL_PROBLEM) ) then
    ! checks if file with station infos located from previous run exists
    inquire(file=trim(OUTPUT_FILES)//'/SU_stations_info.bin',exist=SU_station_file_exists)
    if (SU_station_file_exists) then
      if (myrank==0) then
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
  if (myrank==0) then
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
    call get_elevation_and_z_coordinate(stlon(irec),stlat(irec),stutm_x(irec),stutm_y(irec),z_target(irec),&
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

!--------------------------------------------------------------------------------------------------------------------
!  locate MPI slice which contains the point and bcast to all
!--------------------------------------------------------------------------------------------------------------------
  subroutine locate_MPI_slice_and_bcast_to_all(x_to_locate, y_to_locate, z_to_locate, x_found, y_found, z_found, &
       xi, eta, gamma, ispec_selected, islice_selected, distance_from_target, domain, nu)

  use constants, only: HUGEVAL
  use specfem_par, only: NPROC,myrank

    integer,                                        intent(inout)  :: ispec_selected, islice_selected, domain
    double precision,                               intent(in)     :: x_to_locate, y_to_locate, z_to_locate
    double precision,                               intent(inout)  :: x_found,  y_found,  z_found
    double precision,                               intent(inout)  :: xi, eta, gamma, distance_from_target
    double precision, dimension(3,3),               intent(inout)  :: nu

    double precision,   dimension(:,:), allocatable                :: distance_from_target_all
    double precision,   dimension(:,:), allocatable                :: xi_all, eta_all, gamma_all
    double precision,   dimension(:,:,:), allocatable              :: nu_all
    double precision,   dimension(:,:), allocatable                :: x_found_all, y_found_all, z_found_all
    integer,            dimension(:,:), allocatable                :: ispec_selected_all, domain_all
    integer                                                        :: iproc

    !! to avoid compler error when calling gather_all*
    double precision,  dimension(1)                                :: distance_from_target_dummy
    double precision,  dimension(1)                                :: xi_dummy, eta_dummy, gamma_dummy
    double precision,  dimension(1)                                :: x_found_dummy, y_found_dummy, z_found_dummy
    integer,           dimension(1)                                :: ispec_selected_dummy, islice_selected_dummy, domain_dummy

    allocate(distance_from_target_all(1,0:NPROC-1), &
             xi_all(1,0:NPROC-1), &
             eta_all(1,0:NPROC-1), &
             gamma_all(1,0:NPROC-1), &
             x_found_all(1,0:NPROC-1), &
             y_found_all(1,0:NPROC-1), &
             z_found_all(1,0:NPROC-1), &
             nu_all(3,3,0:NPROC-1))

    allocate(ispec_selected_all(1,0:NPROC-1),domain_all(1,0:NPROC-1))

    distance_from_target = sqrt( (x_to_locate - x_found)**2&
                                +(y_to_locate - y_found)**2&
                                +(z_to_locate - z_found)**2)

    !! it's just to avoid compiler error
    distance_from_target_dummy(1)=distance_from_target
    xi_dummy(1)=xi
    eta_dummy(1)=eta
    gamma_dummy(1)=gamma
    ispec_selected_dummy(1)=ispec_selected
    x_found_dummy(1)=x_found
    y_found_dummy(1)=y_found
    z_found_dummy(1)=z_found
    domain_dummy(1)=domain

    ! gather all on myrank=0
    call gather_all_dp(distance_from_target_dummy, 1, distance_from_target_all, 1, NPROC)
    call gather_all_dp(xi_dummy,    1,  xi_all,    1,  NPROC)
    call gather_all_dp(eta_dummy,   1,  eta_all,   1,  NPROC)
    call gather_all_dp(gamma_dummy, 1,  gamma_all, 1,  NPROC)
    call gather_all_dp(x_found_dummy, 1,  x_found_all, 1,  NPROC)
    call gather_all_dp(y_found_dummy, 1,  y_found_all, 1,  NPROC)
    call gather_all_dp(z_found_dummy, 1,  z_found_all, 1,  NPROC)
    call gather_all_dp(nu, 3*3,  nu_all, 3*3,  NPROC)
    call gather_all_i(ispec_selected_dummy, 1, ispec_selected_all, 1, NPROC)
    call gather_all_i(domain_dummy, 1, domain_all, 1, NPROC)

    ! find the slice and element to put the source
    if (myrank == 0) then

       distance_from_target = HUGEVAL

       do iproc=0, NPROC-1
         if (distance_from_target >= distance_from_target_all(1,iproc)) then
           distance_from_target =  distance_from_target_all(1,iproc)
           islice_selected_dummy(1) = iproc
           ispec_selected_dummy(1) = ispec_selected_all(1,iproc)
           domain_dummy(1) = domain_all(1,iproc)
           xi_dummy(1)    = xi_all(1,iproc)
           eta_dummy(1)   = eta_all(1,iproc)
           gamma_dummy(1) = gamma_all(1,iproc)
           distance_from_target_dummy(1)=distance_from_target
           x_found_dummy(1)=x_found_all(1,iproc)
           y_found_dummy(1)=y_found_all(1,iproc)
           z_found_dummy(1)=z_found_all(1,iproc)
           nu(:,:)=nu_all(:,:,iproc)
         endif
       enddo

    endif

    ! bcast from myrank=0
    call bcast_all_i(islice_selected_dummy,1)
    call bcast_all_i(domain_dummy,1)
    call bcast_all_i(ispec_selected_dummy,1)
    call bcast_all_dp(xi_dummy,1)
    call bcast_all_dp(eta_dummy,1)
    call bcast_all_dp(gamma_dummy,1)
    call bcast_all_dp(distance_from_target_dummy,1)
    call bcast_all_dp(nu,3*3)
    call bcast_all_dp(x_found_dummy,1)
    call bcast_all_dp(y_found_dummy,1)
    call bcast_all_dp(z_found_dummy,1)

    !! it was just to avoid compler error
    islice_selected=islice_selected_dummy(1)
    domain=domain_dummy(1)
    ispec_selected=ispec_selected_dummy(1)
    xi=xi_dummy(1)
    eta=eta_dummy(1)
    gamma=gamma_dummy(1)
    x_found=x_found_dummy(1)
    y_found=y_found_dummy(1)
    z_found=z_found_dummy(1)
    distance_from_target=distance_from_target_dummy(1)

    deallocate(distance_from_target_all, xi_all, eta_all, gamma_all, x_found_all, y_found_all, z_found_all, nu_all)
    deallocate(ispec_selected_all,domain_all)

  end subroutine locate_MPI_slice_and_bcast_to_all

!--------------------------------------------------------------------------------------------------------------------
!  locate point in mesh.
!--------------------------------------------------------------------------------------------------------------------
  subroutine locate_point_in_mesh(x_target, y_target, z_target, POINT_CAN_BE_BURIED, elemsize_max_glob, &
                                  ispec_selected, xi_found, eta_found, gamma_found, x_found, y_found, z_found, domain, nu)

  use constants
  use specfem_par, only: ibool,myrank,NSPEC_AB,NGNOD, &
                         xstore,ystore,zstore,xigll,yigll,zigll,ispec_is_surface_external_mesh,iglob_is_surface_external_mesh
  use specfem_par_acoustic, only: ispec_is_acoustic
  use specfem_par_elastic, only: ispec_is_elastic
  use specfem_par_poroelastic, only: ispec_is_poroelastic


  double precision,                      intent(in)     :: x_target, y_target, z_target
  logical,                               intent(in)     :: POINT_CAN_BE_BURIED
  real(kind=CUSTOM_REAL),                intent(in)     :: elemsize_max_glob
  double precision,                      intent(out)    :: x_found,  y_found,  z_found
  double precision,                      intent(out)    :: xi_found, eta_found, gamma_found
  integer,                               intent(out)    :: ispec_selected, domain
  double precision, dimension(NDIM,NDIM),intent(out)    :: nu

  ! locals
  integer                                               :: iter_loop , ispec, iglob, i, j, k
  ! location search
  double precision                                      :: maximal_elem_size_squared, dist_squared
  double precision                                      :: distmin_squared
  double precision                                      :: x,y,z
  double precision                                      :: xi,eta,gamma,dx,dy,dz,dxi,deta
  double precision                                      :: xixs,xiys,xizs
  double precision                                      :: etaxs,etays,etazs
  double precision                                      :: gammaxs,gammays,gammazs, dgamma
  ! coordinates of the control points of the surface element
  double precision, dimension(NGNOD)                    :: xelm,yelm,zelm
  integer                                               :: ia,iax,iay,iaz
  integer                                               :: ix_initial_guess, iy_initial_guess, iz_initial_guess
  integer, dimension(NGNOD)                             :: iaddx,iaddy,iaddz
  integer                                               :: imin,imax,jmin,jmax,kmin,kmax

  ! sets maximal element size for search
  ! use 10 times the distance as a criterion for source detection
  maximal_elem_size_squared = (10. * elemsize_max_glob)**2

  ! INITIALIZE LOCATION --------
  ispec_selected   = 1   !! first element by default
  ix_initial_guess = 1
  iy_initial_guess = 1
  iz_initial_guess = 1
  ! set distance to huge initial value
  distmin_squared = HUGEVAL

  ! set which GLL points will be considered during research, depending on the possibility to bury the point or not
  if (.not. POINT_CAN_BE_BURIED) then
    imin  = 1
    imax  = NGLLX
    jmin  = 1
    jmax  = NGLLY
    kmin  = 1
    kmax  = NGLLZ
  else
    imin  = 2
    imax  = NGLLX-1
    jmin  = 2
    jmax  = NGLLY-1
    kmin  = 2
    kmax  = NGLLZ-1
  endif

  ! find the element candidate that may contain the target point
  do ispec = 1, NSPEC_AB

    if (.not. POINT_CAN_BE_BURIED .and. .not. ispec_is_surface_external_mesh(ispec)) cycle

    iglob = ibool(MIDX,MIDY,MIDZ,ispec)
    dist_squared = (x_target - dble(xstore(iglob)))**2 &
                 + (y_target - dble(ystore(iglob)))**2 &
                 + (z_target - dble(zstore(iglob)))**2
    if (dist_squared > maximal_elem_size_squared) cycle ! exclude elements that are too far from target

    ! find closest GLL point form target
    do k = kmin, kmax
      do j = jmin, jmax
        do i = imin, imax

          iglob=ibool(i,j,k,ispec)
          if (.not. POINT_CAN_BE_BURIED .and. .not. iglob_is_surface_external_mesh(iglob)) cycle
          dist_squared = (x_target - dble(xstore(iglob)))**2 &
                       + (y_target - dble(ystore(iglob)))**2 &
                       + (z_target - dble(zstore(iglob)))**2

          if (dist_squared < distmin_squared) then

            distmin_squared = dist_squared
            ispec_selected  = ispec
            ix_initial_guess = i
            iy_initial_guess = j
            iz_initial_guess = k

            x_found = xstore(iglob)
            y_found = ystore(iglob)
            z_found = zstore(iglob)

          endif

        enddo
      enddo
    enddo
  
  enddo

  ! get the rotation matrix that will be used to rotate-- source force vector/receiver seismogram --if the point is on the surface
  call define_rotation_matrix(POINT_CAN_BE_BURIED,ix_initial_guess,iy_initial_guess,iz_initial_guess,ispec_selected,nu)

  ! sets whether acoustic (1) or elastic (2)
  if (ispec_is_acoustic( ispec_selected )) then
    domain = IDOMAIN_ACOUSTIC
  else if (ispec_is_elastic( ispec_selected )) then
    domain = IDOMAIN_ELASTIC
  else if (ispec_is_poroelastic( ispec_selected )) then
    domain = IDOMAIN_POROELASTIC
  else
    domain = 0
  endif

  ! general coordinate of initial guess
  xi    = xigll(ix_initial_guess)
  eta   = yigll(iy_initial_guess)
  gamma = zigll(iz_initial_guess)

  ! define topology of the control element
  call usual_hex_nodes(NGNOD,iaddx,iaddy,iaddz)

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

     iglob = ibool(iax,iay,iaz,ispec_selected)
     xelm(ia) = dble(xstore(iglob))
     yelm(ia) = dble(ystore(iglob))
     zelm(ia) = dble(zstore(iglob))

  enddo

  ! iterate to solve the non linear system
  do iter_loop = 1, NUM_ITER

   ! recompute jacobian for the new point
     call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
          xixs,xiys,xizs,etaxs,etays,etazs,gammaxs,gammays,gammazs,NGNOD)

     ! compute distance to target location
     dx = - (x - x_target)
     dy = - (y - y_target)
     dz = - (z - z_target)

     ! compute increments
     dxi  = xixs*dx + xiys*dy + xizs*dz
     deta = etaxs*dx + etays*dy + etazs*dz
     dgamma = gammaxs*dx + gammays*dy + gammazs*dz

     ! update values
     xi = xi + dxi
     eta = eta + deta
     gamma = gamma + dgamma

     ! impose that we stay in that element
     ! (useful if user gives a point outside the mesh for instance)
     if (xi > 1.d0) xi     =  1.d0
     if (xi < -1.d0) xi     = -1.d0
     if (eta > 1.d0) eta    =  1.d0
     if (eta < -1.d0) eta    = -1.d0
     if (gamma > 1.d0) gamma  =  1.d0
     if (gamma < -1.d0) gamma  = -1.d0

  enddo

  ! compute final coordinates of point found
  call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
       xixs,xiys,xizs,etaxs,etays,etazs,gammaxs,gammays,gammazs,NGNOD)

  ! store xi,eta,gamma and x,y,z of point found
  ! note: xi/eta/gamma will be in range [-1,1]
  xi_found = xi
  eta_found = eta
  gamma_found = gamma

  x_found = x
  y_found = y
  z_found = z

  end subroutine locate_point_in_mesh

!--------------------------------------------------------------------------------------------------------------------
!  Define the rotation matrix in the selected point 
!--------------------------------------------------------------------------------------------------------------------
  subroutine define_rotation_matrix(POINT_CAN_BE_BURIED,ix_initial_guess,iy_initial_guess,iz_initial_guess,ispec_selected,nu)

  use constants
  use specfem_par, only: ibool,xstore,ystore,zstore,iglob_is_surface_external_mesh

  logical,                                intent(in)  :: POINT_CAN_BE_BURIED
  integer,                                intent(in)  :: ix_initial_guess,iy_initial_guess,iz_initial_guess,ispec_selected
  double precision, dimension(NDIM,NDIM), intent(out) :: nu

  ! local parameters
  ! for surface locating and normal computing with external mesh
  real(kind=CUSTOM_REAL), dimension(NDIM) :: u_vector,v_vector,w_vector
  integer                                 :: pt0_ix,pt0_iy,pt0_iz,pt1_ix,pt1_iy,pt1_iz,pt2_ix,pt2_iy,pt2_iz

  ! Rotation matrix is diagonal if the point is inside the mesh, or if decided by user
  if ( POINT_CAN_BE_BURIED .or. (.not. EXTERNAL_MESH_RECEIVERS_NORMAL) .or. (ispec_selected == 1) ) then

    !     East
    nu(1,1) = 1.d0
    nu(1,2) = 0.d0
    nu(1,3) = 0.d0
    !     North
    nu(2,1) = 0.d0
    nu(2,2) = 1.d0
    nu(2,3) = 0.d0
    !     Vertical
    nu(3,1) = 0.d0
    nu(3,2) = 0.d0
    nu(3,3) = 1.d0

  else

    ! get normal to the face of the hexaedra if receiver is on the surface
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
    if (ix_initial_guess == 1 .and. &
       iglob_is_surface_external_mesh(ibool(1,2,2,ispec_selected))) then
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
    if (ix_initial_guess == NGLLX .and. &
       iglob_is_surface_external_mesh(ibool(NGLLX,2,2,ispec_selected))) then
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
    if (iy_initial_guess == 1 .and. &
       iglob_is_surface_external_mesh(ibool(2,1,2,ispec_selected))) then
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
    if (iy_initial_guess == NGLLY .and. &
       iglob_is_surface_external_mesh(ibool(2,NGLLY,2,ispec_selected))) then
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
    if (iz_initial_guess == 1 .and. &
       iglob_is_surface_external_mesh(ibool(2,2,1,ispec_selected))) then
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
    if (iz_initial_guess == NGLLZ .and. &
       iglob_is_surface_external_mesh(ibool(2,2,NGLLZ,ispec_selected))) then
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

    if (pt0_ix < 0 .or. pt0_iy < 0 .or. pt0_iz < 0 .or. &
        pt1_ix < 0 .or. pt1_iy < 0 .or. pt1_iz < 0 .or. &
        pt2_ix < 0 .or. pt2_iy < 0 .or. pt2_iz < 0) then
      stop 'error in computing normal for receivers.'
    endif

    u_vector(1) = xstore(ibool(pt1_ix,pt1_iy,pt1_iz,ispec_selected)) &
                - xstore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected))
    u_vector(2) = ystore(ibool(pt1_ix,pt1_iy,pt1_iz,ispec_selected)) &
                - ystore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected))
    u_vector(3) = zstore(ibool(pt1_ix,pt1_iy,pt1_iz,ispec_selected)) &
                - zstore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected))
    v_vector(1) = xstore(ibool(pt2_ix,pt2_iy,pt2_iz,ispec_selected)) &
                - xstore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected))
    v_vector(2) = ystore(ibool(pt2_ix,pt2_iy,pt2_iz,ispec_selected)) &
                - ystore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected))
    v_vector(3) = zstore(ibool(pt2_ix,pt2_iy,pt2_iz,ispec_selected)) &
                - zstore(ibool(pt0_ix,pt0_iy,pt0_iz,ispec_selected))

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

    ! build rotation matrice nu 
    !     East (u)
    nu(1,1) = u_vector(1)
    nu(1,2) = v_vector(1)
    nu(1,3) = w_vector(1)
    !     North (v)
    nu(2,1) = u_vector(2)
    nu(2,2) = v_vector(2)
    nu(2,3) = w_vector(2)
    !     Vertical (w)
    nu(3,1) = u_vector(3)
    nu(3,2) = v_vector(3)
    nu(3,3) = w_vector(3)

  endif 

  end subroutine define_rotation_matrix
