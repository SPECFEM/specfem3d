
  program convert_DX_to_AVS_slices

  implicit none

  integer, parameter :: npoin = 8036, nelem = 7680

  integer :: i1,i2,i3,i4,datavalue,ipoin,ielem

  real x,y,z

! skip DX header
  read(*,*)

! write AVS header
  write(*,*) npoin,' ',nelem,' 1 0 0'

! read and write the points
  do ipoin = 1,npoin
    read(*,*) x,y,z
    write(*,*) ipoin,x,y,z
  enddo

! skip DX header
  read(*,*)

! read and write the elements
! in the case of OpenDX, node numbers start at zero
! in the case of AVS, node numbers start at one
! point order in OpenDX is 1,4,2,3 *not* 1,2,3,4 as in AVS
  do ielem = 1,nelem
    read(*,*) i1,i4,i2,i3
    write(*,*) ielem,' 1 quad ',i1+1,i2+1,i3+1,i4+1
  enddo

! skip DX header
  read(*,*)
  read(*,*)
  read(*,*)

! write AVS header
  write(*,*) '1 1'
  write(*,*) 'Zcoord, meters'

! read and write the point data
  do ipoin = 1,npoin
    read(*,*) datavalue
    write(*,*) ipoin,datavalue
  enddo

  end program convert_DX_to_AVS_slices

