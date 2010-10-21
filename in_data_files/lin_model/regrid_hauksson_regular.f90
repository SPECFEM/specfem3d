
  program jfdkfd

!! DK DK regrid Hauksson on regular grid in lat-long for So-Cal

  implicit none

  include "../../constants.h"

  include "header.h"

  double precision, dimension(NLINES_HAUKSSON_DENSER) :: utm_y_ori,utm_x_ori
  double precision, dimension(16,NLINES_HAUKSSON_DENSER) :: value_ori
  double precision, dimension(NGRID_NEW_HAUKSSON,NGRID_NEW_HAUKSSON) :: utm_y_new,utm_x_new
  double precision, dimension(16,NGRID_NEW_HAUKSSON,NGRID_NEW_HAUKSSON) :: value_new

  integer i,ic,iold,j,k
  double precision dist,distmin

! read original titled non-evenly spaced Hauksson grid
  do i=1,NLINES_HAUKSSON_DENSER
   read(*,*) utm_x_ori(i),utm_y_ori(i),(value_ori(k,i),k=1,16)
  enddo

! generate coordinates of regular grid
  do j=1,NGRID_NEW_HAUKSSON
  do i=1,NGRID_NEW_HAUKSSON
    utm_x_new(i,j) = UTM_X_ORIG_HAUKSSON + dble(i-1) * SPACING_UTM_X_HAUKSSON
    utm_y_new(i,j) = UTM_Y_ORIG_HAUKSSON + dble(j-1) * SPACING_UTM_Y_HAUKSSON
  enddo
  enddo

! interpolate tilted grid onto regular grid
  do j=1,NGRID_NEW_HAUKSSON
  do i=1,NGRID_NEW_HAUKSSON

! look for closest point in denser grid
    distmin = +100000000000.d0
    do iold=1,NLINES_HAUKSSON_DENSER
      dist = dsqrt((utm_x_new(i,j) - utm_x_ori(iold))**2 + (utm_y_new(i,j) - utm_y_ori(iold))**2)
      if(dist < distmin) then
        ic = iold
        distmin = dist
      endif
    enddo

   value_new(:,i,j) = value_ori(:,ic)

  enddo
  enddo

! write new interpolated regular grid
  do j=1,NGRID_NEW_HAUKSSON
  do i=1,NGRID_NEW_HAUKSSON
    print *,sngl(value_new(:,i,j))
  enddo
  enddo

  end

