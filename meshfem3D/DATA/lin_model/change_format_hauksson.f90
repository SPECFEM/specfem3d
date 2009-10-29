
  program jfdkfd

  implicit none

  include "../../constants.h"

  integer i,j,ilayer

  double precision, dimension(27,41) :: rlat1,rlat2,rlon1,rlon2,rlat,rlon
  double precision, dimension(9,27,41) :: vp,vs

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

  do ilayer = 1,9
  do j=1,41
  do i=1,27
   read(*,*) vp(ilayer,i,j)
  enddo
  enddo
  enddo

  do ilayer = 1,9
  do j=1,41
  do i=1,27
   read(*,*) vs(ilayer,i,j)
  enddo
  enddo
  enddo

!! DK DK write in increasing order for lat and long
  do j=1,41
  do i=27,1,-1
!  call utm_geo(rlon(i,j),rlat(i,j),utm_x,utm_y,IZONE_UTM_LA,ILONGLAT2UTM)
  print *,rlon(i,j),rlat(i,j),vp(1,i,j),vp(2,i,j),vp(3,i,j), &
       vp(4,i,j),vp(5,i,j),vp(6,i,j),vp(7,i,j),vp(8,i,j),vp(9,i,j), &
       vs(1,i,j),vs(2,i,j),vs(3,i,j),vs(4,i,j),vs(5,i,j),vs(6,i,j),vs(7,i,j),vs(8,i,j),vs(9,i,j)
!    print *,utm_x,utm_y,' 0'
  enddo
  enddo

  end

!! DK DK add UTM projection routine
  include "../../utm_geo.f90"

