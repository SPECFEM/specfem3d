
  program test_opendx_hexahedra

!-------------------------------------------------------------------------
! creates an OpenDX file showing the edges of the Gocad HR grid
!-------------------------------------------------------------------------

  implicit none

  include "../../../src/shared/constants.h"
  include "constants_gocad.h"

  double precision, dimension(4) :: x,y

  write(*,*) 'creating OpenDX grid'

  x(1) = ORIG_X_GOCAD_MR
  y(1) = ORIG_Y_GOCAD_MR

  x(2) = ORIG_X_GOCAD_MR + SPACING_X_GOCAD_MR * (NX_GOCAD_MR-1)
  y(2) = ORIG_Y_GOCAD_MR

  x(3) = ORIG_X_GOCAD_MR + SPACING_X_GOCAD_MR * (NX_GOCAD_MR-1)
  y(3) = ORIG_Y_GOCAD_MR + SPACING_Y_GOCAD_MR * (NY_GOCAD_MR-1)

  x(4) = ORIG_X_GOCAD_MR
  y(4) = ORIG_Y_GOCAD_MR + SPACING_Y_GOCAD_MR * (NY_GOCAD_MR-1)

!---
!--- write OpenDX file
!---

  open(unit=10,file='gocad_MR_grid_edges.dx', status='unknown')

!--- write nodal coordinates
  write(10,*) 'object 1 class array type float rank 1 shape 3  items 4 data follows'
! use fictitious height of 3000. to place on top of the rest on top view
  write(10,*) sngl(x(1)),sngl(y(1)),' 3000'
  write(10,*) sngl(x(2)),sngl(y(2)),' 3000'
  write(10,*) sngl(x(3)),sngl(y(3)),' 3000'
  write(10,*) sngl(x(4)),sngl(y(4)),' 3000'

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

