
  program convert_avs2dx_quads

! convert AVS files composed of quadrangles ("quads") to OpenDX file format

! Dimitri Komatitsch, University of Pau, January 2007

  implicit none

! maximum number of points in the AVS file
  integer, parameter :: maxpoints = 10000000

  real x(maxpoints)
  real y(maxpoints)
  real z(maxpoints)

  integer icorresp(maxpoints)

  integer :: npoin,nspec,ipoin,ispec,idum1,idum2,i1,i2,i3,i4

! filter "quad" from AVS file
  call system('sed -e ''1,$s/quad//g'' < cubed_sphere_surface.inp > ____tutu')

  open(unit=55,file='____tutu',status='old')

  read(55,*) npoin,nspec

  print *,'object 1 class array type float rank 1 shape 3 items ',npoin,' data follows'
  do ipoin=1,npoin
    read(55,*) i1,x(ipoin),y(ipoin),z(ipoin)
    icorresp(i1) = ipoin
    print *,x(ipoin),y(ipoin),z(ipoin)
  enddo

  print *,'object 2 class array type int rank 1 shape 4 items ',nspec,' data follows'
  do ispec=1,nspec
    read(55,*) idum1,idum2,i1,i2,i3,i4
! point order in OpenDX is 1,4,2,3 *not* 1,2,3,4 as in AVS
    print *,icorresp(i1)-1,icorresp(i4)-1,icorresp(i2)-1,icorresp(i3)-1
  enddo

! skip labels
  read(55,*)
  read(55,*)

  close(55)

  print *,'attribute "element type" string "quads"'
  print *,'attribute "ref" string "positions"'
  print *,'object 3 class array type float rank 0 items ',npoin,' data follows'
  do ipoin=1,npoin
    print *,ipoin
  enddo

  print *,'attribute "dep" string "positions"'
  print *,'object "irregular positions irregular connections" class field'
  print *,'component "positions" value 1'
  print *,'component "connections" value 2'
  print *,'component "data" value 3'
  print *,'end'

! remove temporary file
  call system('rm -f ____tutu')

  end program convert_avs2dx_quads

