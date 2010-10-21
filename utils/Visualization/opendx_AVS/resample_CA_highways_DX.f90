
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
  integer, parameter :: izval = 2700

  real x(0:npoin-1),y(0:npoin-1),dataval(nelem)
  integer i1(nelem),i2(nelem)

  integer nelemnew,ipoin,ielem
  logical p1outside,p2outside
  real xval,yval

! read points
  read(5,*)
  write(*,*) 'object 1 class array type float rank 1 shape 3 items ',npoin,' data follows'

  do ipoin=0,npoin-1
! ignore Z, which we will set to high value to be above topography on display
    read(5,*) x(ipoin),y(ipoin)
    xval = x(ipoin)
    yval = y(ipoin)
    if(xval < xmin) xval = xmin
    if(xval > xmax) xval = xmax
    if(yval < ymin) yval = ymin
    if(yval > ymax) yval = ymax
    write(*,*) xval,yval,izval
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

  if(.not. p1outside .and. .not. p2outside .and. (dataval(ielem) < 0.1 .or. dataval(ielem) > 254.5)) nelemnew = nelemnew + 1
 enddo

 write(*,*) 'object 2 class array type int rank 1 shape 2 items ',nelemnew,' data follows'

! then write elements kept
! exclude elements that are outside of clipping box, or wrong data value
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

  if(.not. p1outside .and. .not. p2outside .and. (dataval(ielem) < 0.1 .or. dataval(ielem) > 254.5)) &
    write(*,*) i1(ielem),i2(ielem)
 enddo


 write(*,*) 'attribute "element type" string "lines"'
 write(*,*) 'attribute "ref" string "positions"'
 write(*,*) 'object 3 class array type float rank 0 items ',nelemnew,' data follows'

! then write element data for elements kept
! exclude elements that are outside of clipping box, or wrong data value
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

  if(.not. p1outside .and. .not. p2outside .and. (dataval(ielem) < 0.1 .or. dataval(ielem) > 254.5)) &
    write(*,*) dataval(ielem)
 enddo

write(*,*) 'attribute "dep" string "connections"'
write(*,*) 'object "irregular connections irregular positions" class field'
write(*,*) 'component "positions" value 1'
write(*,*) 'component "connections" value 2'
write(*,*) 'component "data" value 3'
write(*,*) 'end'

 end


