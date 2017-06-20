program LectureSismoBin
  implicit none
  character*80 sismobinfile,sismoinfofile
  integer ibloc,imin,imax,irec,ns,komp,i,k,kk
  integer nbrc,Nbloc,Nbrec,nseis,ins,iirec,i1,i2,i3
  double precision, dimension(:,:,:,:),allocatable :: Sei
  double precision, dimension(:,:,:),allocatable :: Rec
  write(*,*) 'sismo file ?'
  read(*,'(a)') sismobinfile
  !write(*,*) 'info  file?'
  !read(*,'(a)') sismoinfofile
  open(20,file=sismobinfile,form='unformatted')
  !open(30,file=sismoinfofile)
  read(20) nbrc,Nbloc,Nbrec,nseis
  open(40,file='verifSismo')
  allocate(Sei(nbrc,3,nseis,Nbrec))
  allocate(Rec(3,nseis,Nbrec))
  write(*,*) '--'
  write(*,*) nbrc,Nbloc,Nbrec,nseis
  write(*,*)
  write(*,*) 'coor recepteur',nseis,Nbrec
  read(*,*) ins,iirec

  !read(30,*) i1,i2,i3
  !read(30,*) i1,i2,i3
!
  !do irec=1,Nbrec
  !   do ns=1, nseis
  !      read(30,*) Rec(1,ns,irec),Rec(2,ns,irec),Rec(3,ns,irec)
  !   enddo
  !enddo
!
  do ibloc=1,Nbloc
     IMIN = (ibloc - 1)*nbrc + 1
     IMAX = (ibloc - 1)*nbrc + nbrc
     write(*,*) ibloc, imin, imax
     do irec = 1,Nbrec
        do ns = 1, nseis
           do komp = 1, 3
              read(20) (Sei(i,komp,ns,irec),i=1,nbrc)
           enddo
           do komp =1, 3
           do i=1,nbrc
             if (isnan(Sei(i,komp,ns,irec)) ) then
                 write(*,*) i,komp,ns,irec,'NaN'
                 stop
             endif
           enddo
           enddo
           kk = imin
           if (ns == ins .and. irec == iirec) then
              do k= 1,nbrc
                 kk=kk+1
                 !write(*,*) ibloc,imin,kk,Sei(k,1,ns,irec)
                 write(40,'(2i10,3x,3f30.10)') ibloc,kk,1.d10*Sei(k,1,ns,irec),1.d10*Sei(k,2,ns,irec),1.d10*Sei(k,3,ns,irec)
              enddo
           endif
        enddo
     enddo

  enddo
end program LectureSismoBin


