
  program detect_duplicates_stations_file

! Author: Dimitri Komatitsch, University of Pau and INRIA, France, around 2009.

! This Fortran code also checks if two stations have the same latitude and longitude even if they do not have the same name
! (this can happen in STATIONS files: the same station appearing twice but under different names,
! in particular when the network name is different even if the station name itself is the same).
!
! Otherwise if the station appears twice but with the same name a simple Unix command such as:
!
!    sort STATIONS | uniq
!
! should remove multiples.

  use constants

  implicit none

! input station file to filter
  character(len=150), parameter :: STATIONS_FILE = '../../DATA/STATIONS'
! character(len=150), parameter :: STATIONS_FILE = 'STATIONS_all_20June2008'
! character(len=150), parameter :: STATIONS_FILE = 'STATIONS_SUBSET_35'
! character(len=150), parameter :: STATIONS_FILE = 'STATIONS_FULL_758'

  integer :: irec,irec2,nrec,ios

  character(len=MAX_LENGTH_STATION_NAME), dimension(:), allocatable :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(:), allocatable :: network_name

  logical, dimension(:), allocatable :: is_a_duplicate

  double precision, dimension(:), allocatable :: stlat,stlon,stele,stbur

  character(len=150) :: dummystring

! read the number of receivers
  open(unit=IIN,file=STATIONS_FILE,iostat=ios,status='old',action='read')
  nrec = 0
  do while(ios == 0)
    read(IIN,"(a)",iostat=ios) dummystring
    if(ios == 0) nrec = nrec + 1
  enddo
  close(IIN)

  print *
  print *,'the input station file contains ',nrec,' stations'
  print *

  allocate(stlat(nrec))
  allocate(stlon(nrec))
  allocate(stele(nrec))
  allocate(stbur(nrec))
  allocate(station_name(nrec))
  allocate(network_name(nrec))
  allocate(is_a_duplicate(nrec))

  is_a_duplicate(:) = .false.

  open(unit=IIN,file=STATIONS_FILE,status='old',action='read')
! loop on all the stations to read station information
  do irec = 1,nrec
    read(IIN,*) station_name(irec),network_name(irec),stlat(irec),stlon(irec),stele(irec),stbur(irec)
  enddo
! close receiver file
  close(IIN)

! look for duplicates in the station file in terms of station name
  do irec = 1,nrec-1
    if(is_a_duplicate(irec)) cycle
    do irec2 = irec+1,nrec
      if(is_a_duplicate(irec2)) cycle
      if(station_name(irec) == station_name(irec2) .and. network_name(irec) == network_name(irec2)) then
        print *, &
          network_name(irec2)(1:len_trim(network_name(irec2))),'.',station_name(irec2)(1:len_trim(station_name(irec2))), &
          ' is a duplicate of ',&
          network_name(irec)(1:len_trim(network_name(irec))),'.',station_name(irec)(1:len_trim(station_name(irec))), &
          ' (same name)'
        is_a_duplicate(irec2) = .true.
      endif
    enddo
  enddo

  print *

! look for duplicates in the station file in terms of position
  do irec = 1,nrec-1
    if(is_a_duplicate(irec)) cycle
    do irec2 = irec+1,nrec
      if(is_a_duplicate(irec2)) cycle
      if(stlat(irec) == stlat(irec2) .and. stlon(irec) == stlon(irec2)) then
        print *, &
          network_name(irec2)(1:len_trim(network_name(irec2))),'.',station_name(irec2)(1:len_trim(station_name(irec2))), &
          ' is a duplicate of ', &
          network_name(irec)(1:len_trim(network_name(irec))),'.',station_name(irec)(1:len_trim(station_name(irec))), &
          ' (same lat/long)'
        is_a_duplicate(irec2) = .true.
      endif
    enddo
  enddo

  print *
  print *,'found and removed a total of ',count(is_a_duplicate),' duplicates'
  print *

! write the new filtered station file
  open(unit=IOUT,file=STATIONS_FILE(1:len_trim(STATIONS_FILE))//'_cleaned',status='unknown',action='write')
! loop on all the stations to write station information
  do irec = 1,nrec
    if(.not. is_a_duplicate(irec)) write(IOUT,*) trim(station_name(irec)),' ', &
       trim(network_name(irec)),' ',sngl(stlat(irec)),' ',sngl(stlon(irec)),' ',sngl(stele(irec)),' ',sngl(stbur(irec))
  enddo
! close receiver file
  close(IOUT)

  end program detect_duplicates_stations_file

