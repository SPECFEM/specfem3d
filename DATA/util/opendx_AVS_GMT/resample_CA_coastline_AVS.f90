
  program dfdfd

! create new DX files for Science LA basin simulations
! remove elements outside of clipping box, and remove highways and county lines

  implicit none

  integer, parameter :: npoin = 63098
  integer, parameter :: nelem = 60412

! clipping box (UTM of region under study in Southern California)
  real, parameter :: xmin = 160000
  real, parameter :: xmax = 770000
  real, parameter :: ymin = 3560000
  real, parameter :: ymax = 4100000

! value for Z to be above topography in display
  integer, parameter :: IZ_VALUE = 2700

  real x(0:npoin-1),y(0:npoin-1),dataval(nelem)
  integer i1(nelem),i2(nelem)

  integer nelemnew,ipoin,ielem,ielemreal
  logical p1outside,p2outside
  real xval,yval

! read points
  read(5,*)

  do ipoin=0,npoin-1
! ignore Z, which we will set to high value to be above topography on display
    read(5,*) x(ipoin),y(ipoin)
    xval = x(ipoin)
    yval = y(ipoin)
    if(xval < xmin) xval = xmin
    if(xval > xmax) xval = xmax
    if(yval < ymin) yval = ymin
    if(yval > ymax) yval = ymax
    x(ipoin) = xval
    y(ipoin) = yval
  enddo

! read elements
  read(5,*)
  do ielem=1,nelem
    read(5,*) i1(ielem),i2(ielem)
  enddo

  read(5,*)
  read(5,*)
  read(5,*)

! read data
  do ielem=1,nelem
    read(5,*) dataval(ielem)
  enddo

! first count number of elements to keep
! exclude elements that are outside of clipping box, or wrong data value
 nelemnew = 0
 do ielem=1,nelem
   if(x(i1(ielem)) < xmin .or. x(i1(ielem)) > xmax .or. y(i1(ielem)) < ymin .or. y(i1(ielem)) > ymax) then
     p1outside = .true.
   else
     p1outside = .false.
   endif

   if(x(i2(ielem)) < xmin .or. x(i2(ielem)) > xmax .or. y(i2(ielem)) < ymin .or. y(i2(ielem)) > ymax) then
     p2outside = .true.
   else
     p2outside = .false.
   endif

!!! DK DK coastline only, no faults
  if(.not. p1outside .and. .not. p2outside .and. dataval(ielem) > 254.) nelemnew = nelemnew + 1
 enddo

! write points
  write(*,*) npoin,nelemnew,' 0 1 0'

  do ipoin=0,npoin-1
    write(*,*) ipoin+1,x(ipoin),y(ipoin),IZ_VALUE
  enddo

! then write elements kept
! exclude elements that are outside of clipping box, or wrong data value
 ielemreal = 0
 do ielem=1,nelem
   if(x(i1(ielem)) < xmin .or. x(i1(ielem)) > xmax .or. y(i1(ielem)) < ymin .or. y(i1(ielem)) > ymax) then
     p1outside = .true.
   else
     p1outside = .false.
   endif

   if(x(i2(ielem)) < xmin .or. x(i2(ielem)) > xmax .or. y(i2(ielem)) < ymin .or. y(i2(ielem)) > ymax) then
     p2outside = .true.
   else
     p2outside = .false.
   endif

!!! DK DK coastline only, no faults
  if(.not. p1outside .and. .not. p2outside .and. dataval(ielem) > 254.) then
    ielemreal = ielemreal + 1
    write(*,*) ielemreal,' 1 line ',i1(ielem)+1,i2(ielem)+1
  endif
 enddo


 write(*,*) '1 1'
 write(*,*) 'a, b'

! then write element data for elements kept
! exclude elements that are outside of clipping box, or wrong data value
 ielemreal = 0
 do ielem=1,nelem
   if(x(i1(ielem)) < xmin .or. x(i1(ielem)) > xmax .or. y(i1(ielem)) < ymin .or. y(i1(ielem)) > ymax) then
     p1outside = .true.
   else
     p1outside = .false.
   endif

   if(x(i2(ielem)) < xmin .or. x(i2(ielem)) > xmax .or. y(i2(ielem)) < ymin .or. y(i2(ielem)) > ymax) then
     p2outside = .true.
   else
     p2outside = .false.
   endif

!!! DK DK coastline only, no faults
  if(.not. p1outside .and. .not. p2outside .and. dataval(ielem) > 254.) then
    ielemreal = ielemreal + 1
    write(*,*) dataval(ielem)
  endif
 enddo

 end

