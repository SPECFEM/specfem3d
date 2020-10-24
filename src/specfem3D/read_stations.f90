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


  subroutine read_stations(rec_filename,nrec,station_name,network_name,stlat,stlon,stele,stbur)

  use constants, only: MAX_LENGTH_STATION_NAME,MAX_LENGTH_NETWORK_NAME,IIN,MAX_STRING_LEN

  use specfem_par, only: myrank

  implicit none

  ! input receiver file name
  character(len=*),intent(in) :: rec_filename

  integer,intent(in) :: nrec
  character(len=MAX_LENGTH_STATION_NAME), dimension(nrec),intent(out) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(nrec),intent(out) :: network_name
  double precision, dimension(nrec),intent(out) :: stlat,stlon,stele,stbur

  ! local parameters
  integer :: i,irec,ier
  character(len=MAX_STRING_LEN) :: line

  integer, allocatable, dimension(:) :: station_duplet
  integer :: length_station_name,length_network_name

  ! loop on all the stations to read the file
  if (myrank == 0) then

    ! opens STATIONS or STATIONS_ADJOINT file
    open(unit=IIN,file=trim(rec_filename),status='old',action='read',iostat=ier)
    if (ier /= 0) call exit_mpi(myrank,'error opening file '//trim(rec_filename))

    ! reads all stations
    do irec = 1,nrec
      read(IIN,"(a)",iostat=ier) line
      if (ier /= 0) call exit_mpi(myrank, 'Error reading station file '//trim(rec_filename))

      ! skip comment lines
      if (len_trim(line) > 0 .and. line(1:1) /= '#') then
        line = trim(line)
        read(line,*,iostat=ier) station_name(irec),network_name(irec),stlat(irec),stlon(irec),stele(irec),stbur(irec)
        if (ier /= 0) call exit_mpi(myrank, 'Error reading station file line '//trim(line))
      endif
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
      ! checks if duplicate of another station read in so far
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

  end subroutine read_stations

!
!-------------------------------------------------------------------------------------------
!

  subroutine read_stations_SU_from_previous_run(nrec,station_name,network_name, &
                                                islice_selected_rec,ispec_selected_rec, &
                                                xi_receiver,eta_receiver,gamma_receiver, &
                                                nu_rec,is_done_stations)

  use constants, only: NDIM,IOUT_SU, &
    MAX_LENGTH_STATION_NAME,MAX_LENGTH_NETWORK_NAME,IMAIN,OUTPUT_FILES,MAX_STRING_LEN

  use specfem_par, only: myrank

  implicit none

  integer,intent(in) :: nrec
  character(len=MAX_LENGTH_STATION_NAME), dimension(nrec),intent(in) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(nrec),intent(in) :: network_name

  integer, dimension(nrec),intent(out) :: islice_selected_rec,ispec_selected_rec
  double precision, dimension(nrec),intent(out) :: xi_receiver,eta_receiver,gamma_receiver
  double precision, dimension(NDIM,NDIM,nrec),intent(out) :: nu_rec

  logical, intent(out) :: is_done_stations

  ! local parameters
  integer :: ier,irec

  ! SU_FORMAT parameters
  logical :: SU_station_file_exists

  double precision, allocatable, dimension(:) :: x_found,y_found,z_found
  character(len=MAX_STRING_LEN) :: filename

  ! initializes
  is_done_stations = .false.

  filename = trim(OUTPUT_FILES)//'/SU_stations_info.bin'

  ! checks if file with station infos located from previous run exists
  inquire(file=trim(filename),exist=SU_station_file_exists)

  if (SU_station_file_exists) then

    ! main reads in available station information
    if (myrank == 0) then
      open(unit=IOUT_SU,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
      if (ier /= 0) call exit_mpi(myrank,'error opening file '//trim(filename))

      write(IMAIN,*) 'station details from SU_stations_info.bin'
      call flush_IMAIN()

      allocate(x_found(nrec),y_found(nrec),z_found(nrec),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1971')

      ! reads in station infos
      read(IOUT_SU) islice_selected_rec,ispec_selected_rec
      read(IOUT_SU) xi_receiver,eta_receiver,gamma_receiver
      read(IOUT_SU) x_found,y_found,z_found
      read(IOUT_SU) nu_rec
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
    call bcast_all_dp(nu_rec,NDIM*NDIM*nrec)
    call synchronize_all()

    ! everything done
    is_done_stations = .true.
  endif

  end subroutine read_stations_SU_from_previous_run

!
!-------------------------------------------------------------------------------------------
!

  subroutine write_stations_SU_for_next_run(nrec,islice_selected_rec,ispec_selected_rec, &
                                            xi_receiver,eta_receiver,gamma_receiver, &
                                            x_found,y_found,z_found,nu_rec)


  use constants, only: NDIM,IOUT_SU,OUTPUT_FILES,MAX_STRING_LEN

  implicit none

  integer,intent(in) :: nrec
  integer, dimension(nrec),intent(in) :: islice_selected_rec,ispec_selected_rec
  double precision, dimension(nrec),intent(in) :: xi_receiver,eta_receiver,gamma_receiver
  double precision, dimension(nrec),intent(in) :: x_found,y_found,z_found
  double precision, dimension(NDIM,NDIM,nrec),intent(in) :: nu_rec

  ! local parameters
  integer :: ier
  character(len=MAX_STRING_LEN) :: filename

  filename = trim(OUTPUT_FILES)//'/SU_stations_info.bin'

  open(unit=IOUT_SU,file=trim(filename),status='unknown',action='write',form='unformatted',iostat=ier)
  if (ier == 0) then
    write(IOUT_SU) islice_selected_rec,ispec_selected_rec
    write(IOUT_SU) xi_receiver,eta_receiver,gamma_receiver
    write(IOUT_SU) x_found,y_found,z_found
    write(IOUT_SU) nu_rec
    close(IOUT_SU)
  endif

  end subroutine write_stations_SU_for_next_run
