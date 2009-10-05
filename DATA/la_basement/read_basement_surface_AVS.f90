
  program read_basement_surf

  implicit none

  include "../../constants.h"
  include "../../constants_gocad.h"

! size of Gocad t-surf surface file for basement
 integer, parameter :: NX=144,NY=161

! threshold depth at which we truncate the surface to honor bottom only
 logical, parameter :: APPLY_THRESHOLD_BASEMENT = .false.

 double precision, dimension(NX,NY) :: utm_x,utm_y,z_value

 double precision zdummy1,zdummy2

 integer iline_x,iline_y,idummy,ielem

 do iline_y = 1,NY
  do iline_x = 1,NX
    read(*,*)  idummy,utm_x(iline_x,iline_y),utm_y(iline_x,iline_y),z_value(iline_x,iline_y),zdummy1,zdummy2

! apply threshold to get bottom of basin only
  if(APPLY_THRESHOLD_BASEMENT) then
    if(z_value(iline_x,iline_y) > Z_THRESHOLD_HONOR_BASEMENT) &
              z_value(iline_x,iline_y) = Z_THRESHOLD_HONOR_BASEMENT
  endif

  enddo
  enddo

! AVS header
 write(*,*) NX*NY,' ',(NX-1)*(NY-1),' 1 0 0'

! AVS points
 do iline_y = 1,NY
  do iline_x = 1,NX
    write(*,*)  (iline_y-1)*NX + iline_x,' ',sngl(utm_x(iline_x,iline_y)), &
       ' ',sngl(utm_y(iline_x,iline_y)),' ',sngl(z_value(iline_x,iline_y))
  enddo
  enddo

! AVS elements
 ielem = 0
 do iline_y = 1,NY-1
  do iline_x = 1,NX-1
    ielem = ielem + 1
    write(*,*)  ielem,' 1 quad ',(iline_y-1)*NX + iline_x,' ', &
      (iline_y-1)*NX + iline_x + 1,' ',(iline_y+1-1)*NX + iline_x + 1,' ',(iline_y+1-1)*NX + iline_x
  enddo
  enddo

! output AVS header for data
  write(*,*) '1 1'
  write(*,*) 'Zcoord, meters'

! AVS data for points
 do iline_y = 1,NY
  do iline_x = 1,NX
    write(*,*)  (iline_y-1)*NX + iline_x,' ',sngl(z_value(iline_x,iline_y))
  enddo
  enddo

 end program read_basement_surf

