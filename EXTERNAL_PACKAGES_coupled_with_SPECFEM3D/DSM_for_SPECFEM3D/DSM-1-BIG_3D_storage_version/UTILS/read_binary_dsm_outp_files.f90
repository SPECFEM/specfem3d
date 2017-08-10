program  read_binary_disp
  implicit none
  integer nsamples,nsta,irec,i,k
  double precision, allocatable :: X(:,:)
  character(len=120) fichier
  write(*,*) 'nsamples, nsta ?'
  read(*,*) nsamples, nsta
  write(*,*) 'irec ? '
  read(*,*) irec
  write(*,*) 'fichier ? '
  read(*,'(a)') fichier
  allocate(X(3,nsamples))
  open(10,file=trim(fichier),form='unformatted')
  open(20,file='sismo.txt')
  do i=1,nsta
    read(10) X(1,:)
    read(10) X(2,:)
    read(10) X(3,:)

    if (i == irec) then
      do k=1,nsamples
        write(20,*) X(1,k),X(2,k),X(3,k)
       enddo
    endif
  enddo
  close(10)
  close(20)
end program read_binary_disp
