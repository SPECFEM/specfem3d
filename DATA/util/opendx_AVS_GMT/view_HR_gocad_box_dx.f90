
  program test_opendx_hexahedra

!-------------------------------------------------------------------------
! creates an OpenDX file showing the edges of the Gocad HR grid
!-------------------------------------------------------------------------

  implicit none

  double precision, dimension(4) :: x,y

  write(*,*) 'creating OpenDX grid'

  x(1) = 371052.25
  y(1) = 3725250.

  x(2) = 417052.250000
  y(2) = 3725250.

  x(3) = 417052.250000
  y(3) = 3774000.000

  x(4) = 371052.25
  y(4) = 3774000.000

!---
!--- write OpenDX file
!---

  open(unit=10,file='gocad_HR_grid_edges.dx', status='unknown')

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

