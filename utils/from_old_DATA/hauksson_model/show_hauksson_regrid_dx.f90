
  program jfdkfd

  implicit none

  include "../../constants.h"

  integer i,j,k

  double precision, dimension(NGRID_NEW_HAUKSSON,NGRID_NEW_HAUKSSON) :: utm_y,utm_x
  double precision, dimension(18,NGRID_NEW_HAUKSSON,NGRID_NEW_HAUKSSON) :: value

  do j=1,NGRID_NEW_HAUKSSON
  do i=1,NGRID_NEW_HAUKSSON
    read(*,*) (value(k,i,j),k=1,18)
    utm_x(i,j) = UTM_X_ORIG_HAUKSSON + dble(i-1) * SPACING_UTM_X_HAUKSSON
    utm_y(i,j) = UTM_Y_ORIG_HAUKSSON + dble(j-1) * SPACING_UTM_Y_HAUKSSON
  enddo
  enddo

print *,'object 1 class array type float rank 1 shape 3 items ',NGRID_NEW_HAUKSSON**2,' data follows'
  do j=1,NGRID_NEW_HAUKSSON
  do i=1,NGRID_NEW_HAUKSSON
    print *,utm_x(i,j),utm_y(i,j),' 0'
  enddo
  enddo

print *,'object 2 class array type int rank 1 shape 4 items ',(NGRID_NEW_HAUKSSON-1)**2,' data follows'
  do j=1,NGRID_NEW_HAUKSSON-1
  do i=1,NGRID_NEW_HAUKSSON-1
   print *,(j-1)*NGRID_NEW_HAUKSSON + i - 1, (j+1-1)*NGRID_NEW_HAUKSSON + i - 1, (j-1)*NGRID_NEW_HAUKSSON + i+1 - 1, (j+1-1)*NGRID_NEW_HAUKSSON + i+1 - 1
  enddo
  enddo

print *,'attribute "element type" string "quads"'
print *,'attribute "ref" string "positions"'
print *,'object 3 class array type float rank 0 items ',NGRID_NEW_HAUKSSON**2,' data follows'

  do j=1,NGRID_NEW_HAUKSSON
  do i=1,NGRID_NEW_HAUKSSON
   print *,value(1,i,j)
  enddo
  enddo

print *,'attribute "dep" string "positions"'
print *,'object "irregular connections  irregular positions" class field'
print *,'component "positions" value 1'
print *,'component "connections" value 2'
print *,'component "data" value 3'
print *,'end'

  end

