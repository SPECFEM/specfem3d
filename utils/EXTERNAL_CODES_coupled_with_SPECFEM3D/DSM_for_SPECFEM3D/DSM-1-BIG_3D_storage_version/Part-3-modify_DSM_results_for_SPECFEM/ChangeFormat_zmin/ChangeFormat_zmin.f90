program ChangeFormat
   implicit none
   include 'mpif.h'
   INTEGER myrank,nbproc,ierr,iproc
   character(120) :: coutfile,outputDir,psvmodel,modelname,stationsinf,parentDir
   integer imin_global,imax_global,iprofmin,iprofmax,itranslat,r_n,r_n_global,i
   real(kind(0d0)) r0min,r0max,r0delta,r0lat,r0lon,tlen,f0,f1,dt
   real(kind(0d0)) :: mt(6)
   integer ir,nt,nbrc,ista,nsta,nsta_global,ng
   real(kind(0d0)), allocatable :: SeiLoc(:,:,:,:),SeiGlob(:,:,:,:),SeiLoc1(:,:,:,:)
   real(kind(0d0)), allocatable :: sig(:),sigf(:),gauss(:)
   character(len=11) debug
   integer   , allocatable :: Ind(:)
   integer   , allocatable :: NbInd(:)
   integer Coef,c8,d8,kk
   integer nbloc,i1,i2,i3,i4
   integer norder,irek
   integer ibloc,ns,irec,imin,imax,istamin,istamax
   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,nbproc,ierr)

   nbrc = 100

   !write(debug,'(a6,i5.5)') 'dg_Cgh',myrank
   !open(100,file=debug)

   if (myrank==0) then
   call pinputTra(outputDir,psvmodel,modelname,stationsinf,tlen,imin_global,imax_global,r0min,r0max,&
        r0delta,r0lat,r0lon,itranslat,mt,dt,f0,f1,myrank)

   open(2,file='recdepth')
   read(2,*) r_n_global
   close(2)

   open (1,file=stationsinf,status='old',action='read',position='rewind')
   read(1,*) nsta_global
   close(1)

   open(1,file='NbSmples.txt')
   read(1,*) nt
   close(1)
   endif

   call MPI_Bcast(outputDir,120,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(psvmodel,120,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(modelname,120,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)


   call MPI_Bcast(r_n_global,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(nsta_global,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(nt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(dt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(f0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(f1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
 !!
   norder = 4
   irek = 1
   write(*,*) myrank,nsta_global

   !call DistribDepth(r_n,r_n_global,myrank,nbproc,iprofmin,iprofmax)
   call DistribDepth(nsta,nsta_global,myrank,nbproc,istamin,istamax)
   call MPI_Barrier(MPI_COMM_WORLD,ierr)
   r_n = r_n_global
   r_n = 1
   r_n_global = 1
   iprofmin=1
   iprofmax=1

   if (myrank==0) then
    write(*,*) 'Nb sta', nsta, '/', nsta_global
    write(*,*) 'Istamin, istamax :',istamin,istamax
    write(*,*) 'Nb Prof', r_n , '/ ', r_n_global
    write(*,*) 'Iprof Min, Iprof Max :',iprofmin,iprofmax
   endif


   parentdir = trim(outputDir)//"/out/"

   coef  =  3 * r_n * nbrc
   ng = int(10/f0)
   allocate(SeiLoc(nt,3,nsta,r_n))
   allocate(SeiGlob(nbrc,3,nsta_global,r_n_global))
   allocate(SeiLoc1(nbrc,3,nsta,r_n))
   allocate(Ind(nbproc))
   allocate(NbInd(nbproc))
   allocate(sig(nt),sigf(nt),gauss(ng))
   call mpi_allgather(nsta,1,MPI_INTEGER,NbInd,1,MPI_INTEGER,MPI_COMM_WORLD, ierr )
   call define_gaussian(gauss,f0,dt,ng)

   Ind(1) = 0
   do ir = 2, nbproc
     c8 = NbInd(ir-1)
     Ind(ir) = Coef*c8 +   Ind(ir-1) !NbInd(ir-1)
   enddo

   do ir=1,nbproc
     NbInd(ir) = NbInd(ir) * Coef
   enddo


   !write(100,*) Ind
   !write(100,*)
   !write(100,*) NbInd
   ir = 0


   do i=iprofmin, iprofmax

    ir = ir + 1
    write(coutfile, '("Trac_",i8.8,"_",i8.8)') istamin,istamax
    coutfile = trim(parentDir)//"/"//trim(coutfile)
    !write(*,*) trim(coutfile)
    open(10,file=coutfile,form='unformatted')
    !write(100,'(a)') coutfile
    !read(10) i1,i2,i3,i4
    do ista = 1, nsta

      read(10) sig
      !call convolve_src_function(dt,sig,sigf,gauss,nt,ng)
      call  bwfilt (sig, sigf, dt, nt, irek, norder, f0, f1)
      SeiLoc(1:nt ,1,ista,ir)=sigf
!      SeiLoc(1:nt ,1,ista,ir)=sig

      read(10) sig
      !call convolve_src_function(dt,sig,sigf,gauss,nt,ng)
      call  bwfilt (sig, sigf, dt, nt, irek, norder, f0, f1)
      SeiLoc(1:nt ,2,ista,ir)=sigf
!      SeiLoc(1:nt ,2,ista,ir)=sig

      read(10) sig
      !call convolve_src_function(dt,sig,sigf,gauss,nt,ng)
      call  bwfilt (sig, sigf, dt, nt, irek, norder, f0, f1)
      SeiLoc(1:nt ,3,ista,ir)=sigf
!      SeiLoc(1:nt ,3,ista,ir)=sig

    enddo

    close(10)

   enddo

   Coef = 3*nbrc*r_n*nsta
   nbloc = int(nt/nbrc)
   if (myrank==0) then
     open(20,file=trim(outputdir)//'/Trac.bin',form='unformatted')
     write(20) nbrc,nbloc,r_n_global,nsta_global
   endif
   do ibloc = 1 , nbloc
!

      imin = (ibloc - 1)*nbrc  + 1
      imax = (ibloc - 1)*nbrc  + nbrc
      !write(100,*) 'Bloc ', ibloc,' : ', imin,imax
!
      !write(*,*) myrank,imin,imax,ibloc,nbrc,nt

      !write(*,*) SeiLoc(imin:imax,:,:,:)
      ! faire une convolution par une gaussienne
      ! au lieu de faire cette simple copie
      SeiLoc1(:,:,:,:) = SeiLoc(imin:imax,:,:,:)

      !write(*,*) 'avant gaterv', myrank,coef,Ind
      call mpi_gatherv( SeiLoc1,Coef ,MPI_DOUBLE_PRECISION,SeiGlob,NbInd,Ind,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr )
      !write(*,*) 'apres gatherv ', myrank
!
      do irec = 1, r_n_global
!
         if (myrank==0) write(*,*) irec,'/',r_n_global
         do ns = 1, nsta_global
!
           if (myrank==0) then

             write(20) SeiGlob(1:nbrc,1,ns,irec)
             write(20) SeiGlob(1:nbrc,2,ns,irec)
             write(20) SeiGlob(1:nbrc,3,ns,irec)

             !if (ns==1.and.irec==r_n_global) then
             !do kk=1,nbrc
               !write(100,*) kk, SeiGlob(kk,1,ns,irec)
             !enddo
            !endif


           endif
         enddo
!
      enddo
!    !
!
   enddo



  if (myrank==0) close(20)

   call MPI_FINALIZE(ierr)

end program ChangeFormat




 subroutine DistribDepth(n,n_global,myrank,nbproc,irmin,irmax)
   implicit none
   integer coef,n,n_global,myrank,nbproc,irmin,irmax




   n = int(n_global/nbproc)

   coef = mod(n_global,nbproc)


   if (myrank<mod(n_global,nbproc).and.mod(n_global,nbproc)/=0) then
      n = n + 1
      coef = 0
   endif


   irmin = (myrank ) * n + 1 + coef
   irmax = (myrank + 1 ) * n + coef



   !write(100,*) myrank, irmin, irmax
   !if (myrank<=mod(n_global,nbproc).and.mod(n_global,nbproc)/=0) then
   !   n = n + 1
   !   irmin = irmin + myrank
   !   irmax = irmax + myrank
   !endif



 end subroutine DistribDepth





