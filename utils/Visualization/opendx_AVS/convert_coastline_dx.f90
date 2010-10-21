
  program jdjfdf

! convert AVS L.A. coastline to OpenDX

real, dimension(70000) :: x,y,z
integer, dimension(70000) :: i1,i2,idata,ibool

nelem = 7823
npoin = 8391

do ipoin=1,npoin
read(*,*) ibool(ipoin),x(ipoin),y(ipoin),z(ipoin)
enddo

do ielem=1,nelem
read(*,*) idummy1,idummy2,i1(ielem),i2(ielem)
enddo

do ielem=1,nelem
read(*,*) idummy1,idata(ielem)
enddo

print *,'aaaa'
do ipoin=1,npoin
!!!!!!!! write(*,*) x(ipoin),y(ipoin),z(ipoin)
!!! DK DK no z for flat surface
write(*,*) x(ipoin),y(ipoin),' 0'
enddo

print *,'aaaa'
do ielem=1,nelem
! locate point in list
do ipoin = 1,npoin
if(i1(ielem) == ibool(ipoin)) goto 700
enddo
700 i1val = ipoin
do ipoin = 1,npoin
if(i2(ielem) == ibool(ipoin)) goto 710
enddo
710 i2val = ipoin
! DK DK point number start at 0 in OpenDX
write(*,*) i1val-1,i2val-1
enddo

print *,'aaaa'
do ielem=1,nelem
write(*,*) idata(ielem)
enddo

end

