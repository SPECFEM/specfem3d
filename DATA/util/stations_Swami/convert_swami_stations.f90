
  program convert_swami_stations

  implicit none

  integer, parameter :: N_STATIONS = 25

  integer i,idummy
  double precision lat,long
  character(len=30) name

!! DK DK convert Swami's new stations to the right STATIONS format for SPECFEM3D, i.e.,
! ELKS  AZ   33.5813  -116.4496    0.0    0.0

  do i=1,N_STATIONS
    read(*,*) idummy,name,lat,long
    write(*,*) name(1:len_trim(name)),' SWAMI ',sngl(lat),sngl(-long),'    0.0    0.0'
  enddo

  end program convert_swami_stations

