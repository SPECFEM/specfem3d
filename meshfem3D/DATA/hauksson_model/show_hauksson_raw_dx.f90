
  program jfdkfd

  implicit none

  include "../../constants.h"

  ! UTM projection zone
  integer, parameter :: UTM_PROJECTION_ZONE = 11
  logical, parameter :: SUPPRESS_UTM_PROJECTION = .false.

  integer i,j
  double precision utm_x,utm_y

  double precision, dimension(27,41) :: rlat1,rlat2,rlon1,rlon2,rlat,rlon,value

  do j=1,41
  do i=1,27
   read(*,*) rlat1(i,j)
  enddo
  enddo
  do j=1,41
  do i=1,27
   read(*,*) rlat2(i,j)
  enddo
  enddo

  do j=1,41
  do i=1,27
   read(*,*) rlon1(i,j)
  enddo
  enddo
  do j=1,41
  do i=1,27
   read(*,*) rlon2(i,j)
   rlat(i,j) = rlat1(i,j) + rlat2(i,j)/60.d0
   rlon(i,j) = rlon1(i,j) + rlon2(i,j)/60.d0
   rlon(i,j) = - rlon(i,j)
  enddo
  enddo
  do j=1,41
  do i=1,27
   read(*,*) value(i,j)
  enddo
  enddo

print *,'object 1 class array type float rank 1 shape 3 items 1107 data follows'
  do j=1,41
  do i=1,27
  call utm_geo(rlon(i,j),rlat(i,j),utm_x,utm_y,UTM_PROJECTION_ZONE,ILONGLAT2UTM,SUPPRESS_UTM_PROJECTION)
!!!  print *,rlon(i,j),rlat(i,j),' 0'
    print *,utm_x,utm_y,' 0'
  enddo
  enddo

print *,'object 2 class array type int rank 1 shape 4  items ',40*26,' data follows'
  do j=1,41-1
  do i=1,27-1
   print *,(j-1)*27 + i - 1, (j+1-1)*27 + i - 1, (j-1)*27 + i+1 - 1, (j+1-1)*27 + i+1 - 1
  enddo
  enddo

print *,'attribute "element type" string "quads"'
print *,'attribute "ref" string "positions"'
print *,'object 3 class array type float rank 0 items 1107 data follows'

  do j=1,41
  do i=1,27
   print *,value(i,j)
  enddo
  enddo

print *,'attribute "dep" string "positions"'
print *,'object "irregular connections  irregular positions" class field'
print *,'component "positions" value 1'
print *,'component "connections" value 2'
print *,'component "data" value 3'
print *,'end'

  end

!! DK DK add UTM projection routine
  include "../../utm_geo.f90"

