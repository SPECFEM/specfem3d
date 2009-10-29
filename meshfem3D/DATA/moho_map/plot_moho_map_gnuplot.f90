
  program plot_moho_map

! plot Lupei Zhu's Moho map using gnuplot

  implicit none

  include "../../constants.h"

! use integer array to store Moho depth
  integer imoho_depth(NX_MOHO,NY_MOHO)

  integer ix,iy

  double precision long,lat
  double precision, dimension(NX_MOHO,NY_MOHO) :: utmx,utmy
  double precision depth_km

  imoho_depth(:,:) = 0

! read file
  open(unit=13,file='moho_lupei_zhu.dat',status='old')
! file starts from North-West corner
  do iy=NY_MOHO,1,-1
    do ix=1,NX_MOHO
      read(13,*) long,lat,depth_km
! convert depth to meters
      imoho_depth(ix,iy) = nint(depth_km * 1000.d0)
! convert to UTM
      call utm_geo(long,lat,utmx(ix,iy),utmy(ix,iy),11,ILONGLAT2UTM)
    enddo
  enddo
  close(13)

! create gnuplot file
  open(unit=3,file='moho_gnuplot.dat',status='unknown')
  do iy=1,NY_MOHO
    do ix=1,NX_MOHO-1
      write(3,*) sngl(utmx(ix,iy)),sngl(utmy(ix,iy)),- imoho_depth(ix,iy)
      write(3,*) sngl(utmx(ix+1,iy)),sngl(utmy(ix+1,iy)),- imoho_depth(ix+1,iy)
    enddo
    write(3,*)
  enddo
  close(3)

! create script for gnuplot
  open(unit=3,file='plot_all.gnu',status='unknown')
  write(3,*) 'set term x11'
  write(3,*) 'set view 40,60'
  write(3,*) 'splot "moho_gnuplot.dat" w linesp 1'
  write(3,*) 'pause -1 "hit key"'
  close(3)

  end program plot_moho_map

  include "../../utm_geo.f90"

