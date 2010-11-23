
  program display_topo_bathy_opendx

!! DK DK display topo and bathy of L.A. region using OpenDX

  implicit none

  include "../../constants.h"

!! DK DK subsampled file for OpenDX display
  integer, parameter :: NX = 141,NY = 101

! amplify Z values if needed to see mountains
!!!  integer, parameter :: factor_Z = 30
  integer, parameter :: factor_Z = 1

  integer itopo(NX,NY)

  integer ix,iy

  double precision rlon,rlat,rx,ry

  open(unit=13,file='topo_bathy_subsampled_OpenDX.dat',status='old')
  do iy = 1,NY
    do ix = 1,NX
      read(13,*) itopo(ix,iy)
    enddo
  enddo
  close(13)

! write OpenDX file
  open(unit=13,file='topo_bathy_LA.dx',status='unknown')
  write(13,*) '# min max topo file = ',minval(itopo),maxval(itopo)
  write(13,*) '# number of points = ',NX*NY
  write(13,*) '# number of elements = ',(NX-1)*(NY-1)

! write list of points
  write(13,*) 'object 1 class array type float rank 1 shape 3 items ',NX*NY,' data follows '
  do iy = 1,NY
  do ix = 1,NX
    rlon = ORIG_LONG + (ix-1)*DEGREES_PER_CELL*SUBSAMP_FACTOR_OPENDX
    rlat = ORIG_LAT + (iy-1)*DEGREES_PER_CELL*SUBSAMP_FACTOR_OPENDX
    call utm_geo(rlon,rlat,rx,ry,IZONE_UTM_LA,ILONGLAT2UTM)
    write(13,*) sngl(rx),sngl(ry),itopo(ix,iy)*factor_Z
  enddo
  enddo

! write list of elements
  write(13,*) 'object 2 class array type int rank 1 shape 4 items ',(NX-1)*(NY-1),' data follows'
  do iy = 1,NY-1
  do ix = 1,NX-1
    write(13,*) (iy-1)*NX + ix - 1, (iy-1+1)*NX + ix - 1, (iy-1)*NX + ix+1 - 1, (iy-1+1)*NX + ix+1 - 1
  enddo
  enddo

! write data for points
  write(13,*) 'attribute "element type" string "quads"'
  write(13,*) 'attribute "ref" string "positions"'
  write(13,*) 'object 3 class array type float rank 0 items ',NX*NY,' data follows'
  do iy = 1,NY
  do ix = 1,NX
    write(13,*) itopo(ix,iy)
  enddo
  enddo

  write(13,*) 'attribute "dep" string "positions"'
  write(13,*) 'object "irregular connections  irregular positions" class field'
  write(13,*) 'component "positions" value 1'
  write(13,*) 'component "connections" value 2'
  write(13,*) 'component "data" value 3'
  write(13,*) 'end'
  close(13)

  end program display_topo_bathy_opendx

!! DK DK add UTM projection routine
  include "../../utm_geo.f90"

