
  program plot_moho_map

! plot Lupei Zhu's Moho map using OpenDX

  implicit none

  include "../../constants.h"

! use integer array to store Moho depth
  integer imoho_depth(NX_MOHO,NY_MOHO)

  integer ix,iy,iglob1,iglob2,iglob3,iglob4

  double precision long,lat
  double precision, dimension(NX_MOHO,NY_MOHO) :: utmx,utmy
  double precision depth_km,sumval

  imoho_depth(:,:) = 0
  sumval = 0.d0

! read file
  open(unit=13,file='moho_lupei_zhu.dat',status='old')
! file starts from North-West corner
  do iy=NY_MOHO,1,-1
    do ix=1,NX_MOHO
      read(13,*) long,lat,depth_km
! convert depth to meters
      imoho_depth(ix,iy) = nint(depth_km * 1000.d0)
      sumval = sumval + dble(imoho_depth(ix,iy))
! convert to UTM
      call utm_geo(long,lat,utmx(ix,iy),utmy(ix,iy),11,ILONGLAT2UTM)
    enddo
  enddo
  close(13)

  print *,'Moho depth min max mean = ',minval(imoho_depth),maxval(imoho_depth),sumval/dble(NX_MOHO*NY_MOHO)

! create OpenDX file
  open(unit=3,file='moho_map.dx',status='unknown')
  write(3,*) 'object 1 class array type float rank 1 shape 3 items ',NX_MOHO*NY_MOHO,' data follows'

  do iy=1,NY_MOHO
    do ix=1,NX_MOHO
      write(3,*) sngl(utmx(ix,iy)),sngl(utmy(ix,iy)),- imoho_depth(ix,iy)
    enddo
  enddo

  write(3,*) 'object 2 class array type int rank 1 shape 4 items ',(NX_MOHO-1)*(NY_MOHO-1),' data follows'

  do iy=1,NY_MOHO-1
    do ix=1,NX_MOHO-1
      iglob1 = (iy-1)*NX_MOHO + ix
      iglob2 = (iy-1)*NX_MOHO + ix+1
      iglob3 = (iy+1-1)*NX_MOHO + ix+1
      iglob4 = (iy+1-1)*NX_MOHO + ix
      write(3,210) iglob1-1,iglob4-1,iglob2-1,iglob3-1
    enddo
  enddo

  write(3,*) 'attribute "element type" string "quads"'
  write(3,*) 'attribute "ref" string "positions"'
  write(3,*) 'object 3 class array type float rank 0 items ',NX_MOHO*NY_MOHO,' data follows'

  do iy=1,NY_MOHO
    do ix=1,NX_MOHO
      write(3,*) - imoho_depth(ix,iy)
    enddo
  enddo

  write(3,*) 'attribute "dep" string "positions"'
  write(3,*) 'object "irregular positions irregular connections" class field'
  write(3,*) 'component "positions" value 1'
  write(3,*) 'component "connections" value 2'
  write(3,*) 'component "data" value 3'
  write(3,*) 'end'

210 format(i6,1x,i6,1x,i6,1x,i6)

  close(3)

  end program plot_moho_map

  include "../../utm_geo.f90"

