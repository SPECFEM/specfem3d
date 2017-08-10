program cree_frq
  implicit none
  integer i,myrank,imin,imax
  integer nfrqs,nproc,ifq
  write(*,*) 'nproc ?'
  read(*,*) nproc
  imin=0;imax=0
  open(20,file='frqs_out.txt')
  write(20,*) '2'
  do i=0,nproc-1
    write(20,*) i,i+imin,i+imax
  enddo
end program cree_frq
