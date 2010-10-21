
  program test_opendx_hexahedra

!-------------------------------------------------------------------------
! creates an OpenDX file showing the edges of the SEM grid
!-------------------------------------------------------------------------

  implicit none

  include "../../../constants.h"

  double precision, parameter :: LONGMIN = -120.5d0,LONGMAX = -114.5d0
  double precision, parameter :: LATMIN = 32.5d0,LATMAX = 36.5d0

  double precision, dimension(4) :: x,y

  write(*,*) 'creating OpenDX grid'

  call utm_geo(LONGMIN,LATMIN,x(1),y(1),IZONE_UTM_LA,ILONGLAT2UTM)
  call utm_geo(LONGMAX,LATMIN,x(2),y(2),IZONE_UTM_LA,ILONGLAT2UTM)
  call utm_geo(LONGMAX,LATMAX,x(3),y(3),IZONE_UTM_LA,ILONGLAT2UTM)
  call utm_geo(LONGMIN,LATMAX,x(4),y(4),IZONE_UTM_LA,ILONGLAT2UTM)

!---
!--- write OpenDX file
!---

  open(unit=10,file='basin_grid_edges.dx', status='unknown')

!--- write nodal coordinates
  write(10,*) 'object 1 class array type float rank 1 shape 3  items 4 data follows'
! use fictitious height of 1. to place on top of the rest on top view
  write(10,*) sngl(x(1)),sngl(y(1)),' 1'
  write(10,*) sngl(x(2)),sngl(y(2)),' 1'
  write(10,*) sngl(x(3)),sngl(y(3)),' 1'
  write(10,*) sngl(x(4)),sngl(y(4)),' 1'

!--- write connectivity pattern
!--- OpenDX node numbers start at zero
  write(10,*) 'object 2 class array type int rank 1 shape 4 items 1 data follows'
  write(10,*) '0 3 1 2'
  write(10,*) 'attribute "element type" string "quads"'
  write(10,*) 'attribute "ref" string "positions"'

!--- write element data
  write(10,*) 'object 3 class array type float rank 0  items 1 data follows'
  write(10,*) '100'
  write(10,*) 'attribute "dep" string "connections"'

!--- define field
  write(10,*) 'object "irregular connections  irregular positions" class field'
  write(10,*) 'component "positions" value 1'
  write(10,*) 'component "connections" value 2'
  write(10,*) 'component "data" value 3'
  write(10,*) 'end'

  close(10)

  end program test_opendx_hexahedra

!! DK DK add UTM projection routine
  include "../../../utm_geo.f90"

