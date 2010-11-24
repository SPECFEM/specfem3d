
 program extract_curve_Sieh

 implicit none

 integer, parameter :: NX = 2638
 integer, parameter :: NY = 665

 integer image(NX,NY)

 integer ix,iy,imax,index_imax,numval

 double precision xval,yval,sum_values

 do iy=1,NY
 do ix=1,NX
! image is stored from top to bottom
   read(5,*) image(ix,NY-iy+1)
 enddo
 enddo

 sum_values = 0.
 numval = 0

! exclude first column, which has a problem
 do ix=2,NX
  imax = +1000000
  index_imax = -1
! detect minimum value for this column, corresponding to darkest pixel
 do iy=1,NY
   if(image(ix,iy) < imax) then
     imax = image(ix,iy)
     index_imax = iy
   endif
 enddo
! take into account the fact that we excluded the first column
 xval = dble(ix-2)*(327.95-(-34.45))/dble(NX-2) -34.45
 yval = (3.518-0.742)*dble(index_imax-5)/201.d0 + 0.742
 print *,xval,yval
 sum_values = sum_values + yval
 numval = numval + 1
 enddo

 print *,'mean value of slip (m) = ',sum_values/numval

 end program extract_curve_Sieh

! these values from DXF file provided by Andreas Plesch
!
! -34.45 3.518
! 327.95 0.742

