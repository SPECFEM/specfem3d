
 program rescalenonlinear

! rescale GMT data based upon a nonlinear function

! to use it:
! ./xscale < gmt_shaking_hollywood_accel.irreg > gmt_shaking_yorba_accel.irreg2

 implicit none

 double precision, parameter :: POWERVAL = 0.25d0

 integer, parameter :: NLINES = 623091

 integer i

 double precision, dimension(NLINES) :: a,b,c
 double precision cmax

 do i = 1,NLINES
  read(*,*) a(i),b(i),c(i)
 enddo

 cmax = maxval(c)

 do i = 1,NLINES
  write(*,*) a(i),b(i),cmax*(c(i)/cmax)**POWERVAL
 enddo

 end

