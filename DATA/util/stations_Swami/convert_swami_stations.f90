
  program convert_swami_stations

  implicit none

  integer, parameter :: N_STATIONS = 392

  integer i,idummy

  double precision lat,long
  double precision dummy1,dummy2,dummy3,dummy4

  character(len=30) name

!! DK DK convert Swami's new stations to the right STATIONS format for SPECFEM3D, i.e.,
! ELKS  AZ   33.5813  -116.4496    0.0    0.0

! ignore first line, which is a text comment
  read(*,*)

! write total number of stations
  write(*,*) N_STATIONS

  do i=1,N_STATIONS
    read(*,*) idummy,name,dummy1,dummy2,lat,dummy3,dummy4,long
    write(*,*) name(1:len_trim(name)),' SWAMI ',sngl(lat),sngl(-long),' 0 0'
  enddo

  end program convert_swami_stations

