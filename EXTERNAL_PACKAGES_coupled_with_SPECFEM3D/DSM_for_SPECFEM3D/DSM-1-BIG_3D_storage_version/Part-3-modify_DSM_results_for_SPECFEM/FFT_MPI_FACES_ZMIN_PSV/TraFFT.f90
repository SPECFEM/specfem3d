program TraFFT


!-------------------------------------------------------------
!
!   TraFFT
!
!
!
!
!
!
!                                    2010.08. FUJI Nobuaki
!                                    2010.10. FUJI Nobuaki
!                                    2011.10  MONTEILLER Vadim
!
!-------------------------------------------------------------

  implicit none
  ! MPI  ---- VM VM
   include 'mpif.h'
   INTEGER myrank,nbproc,ierr,iproc! MPI  ---- VM VM


  integer :: iPS = 2 !3 ! 1:SH, 2:PSV, 3:ALL
  character(120) :: outputDir,psvmodel,modelname,stationsinf,parentDir
  integer :: itranslat
  character(120) :: coutfile
  real(kind(0d0)) :: tlen,r0min,r0max,r0delta,r0lat,r0lon,r_dum

  real(kind(0d0)), allocatable :: stla(:),stlo(:),stla_g(:),stlo_g(:)
  real(kind(0d0)), allocatable :: r_(:),r0(:),theta(:),phi(:),theta_g(:),phi_g(:)
  integer, allocatable :: updown(:)
  integer :: r_n,r0_n,ir0,imt,theta_n,nsta,r_n_global
  integer :: k,i,j,jj,imin,imax,np,updown_dum

  !for FFT, imin should be 0
  real(kind(0d0)) :: omegai,dt
  real(kind(0d0)), parameter :: pi=3.1415926535897932d0
  integer,parameter :: maxradsamples = 5
  ! FOR FOURIER TRANSFORM & BUTTERWORTH FILTERING


  real(kind(0d0)) :: mt(6)

  integer :: iwindowStart, iwindowEnd
  integer, parameter:: ibwfilt = 0 ! if =0, we don't perform filtering
  real(kind(0d0)) :: samplingHz,start,ending
  integer :: lsmooth,np0,np1
  !real(kind(0d0)),parameter :: start =  350.d0 ! Time window
  !real(kind(0d0)),parameter :: ending = 700.d0  ! in seconds

  ! RSGT & TSGT

  ! single !!!!!!
  complex(kind(0e0)),allocatable :: stresssngl(:,:,:,:), displacementsngl(:,:,:,:),tmpsngl(:,:,:),tmpsngl1(:,:,:)
  real(kind(0d0)), allocatable :: ygt(:,:)
  !!!!!

  complex(kind(0d0)), allocatable :: gt(:,:)!,normal(:,:),normal_g(:,:)
  real(kind(0d0)), allocatable :: tmpygt(:)

!----
  integer max_rank,irec,irank,iprof,r_n_loc,imin_global,imax_global,iprofmin,iprofmax,idum&
  ,imax_glob,imin_glob
  integer, allocatable :: Ifrq(:),Ifrq2(:,:)
! ---- VM VM
  real(kind(0d0)) SourceLat,SourceLong,epidis,azimuth,bazi,stla_dum,stlo_dum
  real(kind(0d0)) Orign_Lon_Chunk,Origin_Lat_Chunk,Lat_current,Lon_current,azi_chunk,delta_lat,delta_lon
  character(len=5) face
  !complex(kind(0d0)) normal(3)
  character(len=11) debug
  integer norder,icomp,iidum,iidum1,iidum2,para
  real(kind(0d0)) fgauss,nrx,nry,nrz
  integer ibloc,nbrc,Nbloc,i0,i1,nsamp_out
  real :: tstart, tfinish
  integer istamin, istamax, nsta_global
!-----------------------------------------------------------------------
 call MPI_INIT(ierr)
 call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
 call MPI_COMM_SIZE(MPI_COMM_WORLD,nbproc,ierr)


!! ---
   Nbrc = 100
   norder = 6
   !fl  =  0.004d0
   !fh = 0.04d0
!!----
  write(debug,'(a6,i5.5)') 'debug_',myrank
  open(100,file=debug)
! lecture des frequences
if (myrank == 0) then
    open(25,file='FrqsMpi.txt')
    read(25,*) para
    if (para == 1) then
      max_rank = 0
      do
        read(25,*,end=100) iidum
        max_rank = max_rank + 1
      enddo
100   close(25)
    else
      max_rank = 0
      do
        read(25,*,end=101) iidum,iidum1,iidum2
        max_rank = max_rank + 1
      enddo
101 close(25)
    endif
  endif
  call MPI_Bcast(max_rank,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  allocate(Ifrq(max_rank))
  allocate(Ifrq2(0:max_rank-1,2))

  if (myrank == 0) then
    open(25,file='FrqsMpi.txt')
    read(25,*) para
    if (para == 1) then
      do irank =1, max_rank
        read(25,*) Ifrq(irank)
      enddo
      Ifrq(1) = -1
      close(25)
     else
      do iidum=0,max_rank-1
         read(25,*)  irank,Ifrq2(irank,1),Ifrq2(irank,2)
      enddo
     endif
  endif
  call MPI_Bcast(Ifrq,max_rank,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(Ifrq2,2*max_rank,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(para,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)


! lecture des inputs
  face=stationsinf(3:6)
  if (myrank == 0) then
    call InputForTraction(Orign_Lon_Chunk,Origin_Lat_Chunk,azi_chunk,delta_lat,delta_lon)
    call pinputTra(outputDir,psvmodel,modelname,stationsinf,tlen,imin_glob,imax_glob,r0min,r0max, &
    r0delta,r0lat,r0lon,itranslat,mt,myrank,samplingHz,start,ending,fgauss)
  endif
  call MPI_Bcast(Orign_Lon_Chunk,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(Origin_Lat_Chunk,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(azi_chunk,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(delta_lat,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(delta_lon,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(outputDir,120,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(psvmodel,120,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(modelname,120,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(stationsinf,120,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(tlen,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(r0min,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(r0max,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(r0delta,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(r0lat,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(r0lon,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(imin_glob,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(imax_glob,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(mt,6,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(samplingHz,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(start,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ending,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(fgauss,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  !call MPI_Bcast(fh,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  itranslat=0

  ! bornes frequences
  if (para == 1) then
    imin_global = Ifrq(1) + 1
    imax_global =Ifrq(max_rank)
  else
     imin_global=minval(Ifrq2)
     imax_global=maxval(Ifrq2)
  endif

  np = imax_global

  parentdir = trim(outputDir)//"/out/"
  omegai = - dlog(1.d-2) / tlen



  if (myrank == 0) then
    open (1,file=stationsinf,status='old',action='read',position='rewind')
    !open(2,file='recdepth')
    !read(2,*) r_n_global
    read(1,*)nsta_global
    !r_n_global = 1
  endif
  !call MPI_Bcast(r_n_global,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(nsta_global,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  r_n_global = 1 ! on est sur zmin : une seule profondeur
  iprofmin=1
  iprofmax=r_n_global
  r_n=r_n_global
  ! parallelisation selon les distances epicentrales
  call DistribDepth(nsta,nsta_global,myrank,nbproc,istamin,istamax)
  theta_n = nsta
  allocate(r_(1:r_n))
  allocate(theta(1:theta_n))
  allocate(phi(1:theta_n))
  allocate(stla(1:nsta))
  allocate(stlo(1:nsta))
  allocate(updown(1:nsta))
  !allocate(normal(3,nsta))

  allocate(theta_g(1:nsta_global))
  allocate(phi_g(1:nsta_global))
  allocate(stla_g(1:nsta_global))
  allocate(stlo_g(1:nsta_global))
  !allocate(normal_g(3,nsta_global))


  if (myrank == 0) then
    !open(3,file='normal_stzmin')
    k = 0
    ! normale dans le repere local
    nrx=0.d0
    nry=0.d0
    nrz=1.d0
    do i = 1,nsta_global
       read(1,*) stlo_dum,stla_dum
      !read(3,*) nrx,nry,nrz
       !if (i >= istamin .and. i <= istamax) then
       !    k = k + 1
         stla_g(i) = stla_dum
         stlo_g(i) = stlo_dum
         !normal_g(2,i)=dcmplx(nrx)
         !normal_g(3,i)=dcmplx(nry)
         !normal_g(1,i)=dcmplx(nrz)
         !r_(i) = 6371.d0 -r_(i)
         !if (itranslat==1) call translat(stla(k),stla(k))
         !call calthetaphi(r0lat,r0lon,stla(k),stlo(k),theta(k),phi(k))
         call epitra1(r0lat,r0lon,stla_g(i),stlo_g(i),theta_g(i),phi_g(i),bazi)
      !enddo
    enddo

    k  =  0
    !do i=1,r_n_global
      !read(2,*) r_dum,idum,updown_dum
      !if (i >= iprofmin .and. i <= iprofmax) then
      !   k  = k + 1
       !r_(i) = 6371.d0 -r_dum
       !updown(i) = updown_dum
       !write(100,*) k,r_(k)
    !endif
    !enddo
    close(1)
    !close(2)
    close(3)
  endif

  !call MPI_Bcast(normal_g,3*nsta_global,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(theta_g,nsta_global,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(phi_g,nsta_global,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(stla_g,nsta_global,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(stlo_g,nsta_global,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(r_,r_n_global,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  !call MPI_Bcast(updown,r_n_global,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  !normal(:,:)=normal_g(:,istamin:istamax)
  theta(:)=theta_g(istamin:istamax)
  phi(:)=phi_g(istamin:istamax)
  stla(:)=stla_g(istamin:istamax)
  stlo(:)=stlo_g(istamin:istamax)

  r0_n =  int((r0max-r0min)/r0delta)+1
  allocate(r0(1:r0_n))
  do i = 1, r0_n
     r0(i) = r0min + dble(i-1)*r0delta
  enddo



  ! lsmoothfinder for FFT

  np0 = np
  call lsmoothfinder(tlen, np0, samplingHz, lsmooth)

  i=1
  do while (i < lsmooth)
     i = i*2
  enddo
  lsmooth = i
  i = 0

  np1 = 1
  do while (np1 < np0)
     np1 = np1*2
  enddo

  np1 = np1*lsmooth

  samplingHz = dble(2*np1)/tlen
  iWindowStart = int(start*samplingHz)
  iWindowEnd   = int(ending*samplingHz)
  dt = 1.d0 / samplingHz
  if (myrank == 0) then
  write(*,*) 'np1 :', np1
  write(*,*) 'sampling,dt : ', samplingHz,dt
  write(*,*)
  write(*,*) 'Window :', iWindowStart, iWindowEnd
  write(*,*) 'time w:', (iWindowStart-1)*dt,(iWindowEnd-1)*dt
  !if (myrank==0) then
     open(11,file='NbSmples.txt')
     write(11,*) iWindowEnd -iWindowStart + 1
     close(11)
  endif

  allocate(gt(1:9,0:2*np1-1))
  allocate(ygt(1:9,0:2*np1-1))
  allocate(tmpygt(0:2*np1-1))
  allocate(stresssngl(1:6,1:6,1:nsta,imin_global:imax_global))
  allocate(displacementsngl(1:3,1:6,1:nsta,imin_global:imax_global))
  allocate(tmpsngl(1:6,1:6,1:nsta_global))
  allocate(tmpsngl1(1:3,1:6,1:nsta_global))

  stresssngl=cmplx(0.e0)
  displacementsngl=cmplx(0.e0)
  ir0 = r0_n



  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  write(100,*) 'lecture des spectres '
  write(100,*)
  k=1
  do irank=1,imax_global + 1 !---------------------------------- FRQS
   !k  = 0
   write(100,*)
   write(100,*) '-----'
   write(100,*) ' lecture spectre calcule par myrank=',irank
   write(100,*) 'nsta global ',nsta_global
   write(100,*) 'r_n_global ', r_n_global
   i = irank - 1

     if (iPS /= 1) then ! PSV calculation
        write(coutfile, '(I5.5,".Stress_PSV_",i5.5)') int(r0(ir0)*10.d0),irank-1

        coutfile = trim(modelname)//"."//coutfile
        coutfile = trim(outputDir)//"/Stress/"//coutfile
        !! VM VM chagne the open statement to be consistent with TraPSV_MPI_read_Zmin
        !open(1,file=coutfile,status='unknown',form='unformatted', &
        !     access = 'direct', recl=2*6*6*kind(0e0)*nsta_global*r_n_global)
        !read(1,rec=k) tmpsngl(1:6,1:6,1:nsta_global)
        if (myrank == 0) then
        open(1,file=coutfile,form='unformatted',action='read')
        read(1) tmpsngl
        endif
        call mpi_bcast(tmpsngl,36*nsta_global,MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)
        !write(100,*) k,istamin,tmpsngl(1:6,1:6,istamin)
        stresssngl(1:6,1:6,1:nsta,i)=stresssngl(1:6,1:6,1:nsta,i)+tmpsngl(1:6,1:6,istamin:istamax)
        !write(100,*) '                suite>'
        !write(100,*) i,stresssngl(1:6,1:6,1,i)
        !close(1)
!! VM VM only one file
        !write(coutfile, '(I5.5,".Displ_PSV_",i5.5)') int(r0(ir0)*10.d0),irank-1

        !coutfile = trim(modelname)//"."//coutfile
        !coutfile = trim(outputDir)//"/Displacement/"//coutfile
        !write(100,*) trim(coutfile)
        !! VM VM chagne the open statement to be consistent with TraPSV_MPI_read_Zmin
        !open(1,file=coutfile,status='unknown',form='unformatted', &
        !     access = 'direct', recl=2*3*6*kind(0e0)*nsta_global*r_n_global)
        !read (1,rec=k) tmpsngl1(1:3,1:6,1:nsta_global)
        !open(1,file=coutfile,form='unformatted',action='read')
        if (myrank == 0) read(1) tmpsngl1
        call mpi_bcast(tmpsngl1,18*nsta_global,MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)
        displacementsngl(1:3,1:6,1:nsta,i)=displacementsngl(1:3,1:6,1:nsta,i)+tmpsngl1(1:3,1:6,istamin:istamax)

        !write(100,*) i,displacementsngl(1:3,1:6,1:1,i)
        close(1)
     endif
     if (iPS /= 2) then ! SH calculation   !! TO DO
        write(coutfile, '(I5.5,".Stress_SH_",i5.5)') int(r0(ir0)*10.d0),irank-1
        !do j = 1,7
        !   if (coutfile(j:j)==' ')coutfile(j:j) = '0'
        !enddo

        coutfile = trim(modelname)//"."//coutfile
        coutfile = trim(outputDir)//"/Stress/"//coutfile

        open(1,file=coutfile,status='unknown',form='unformatted', &
             access = 'direct', recl=2*6*6*kind(0e0)*nsta_global*r_n_global)
        read(1,rec=k) tmpsngl(1:6,1:6,1:nsta_global)
        stresssngl(1:6,1:6,1:nsta,i)=stresssngl(1:6,1:6,1:nsta,i)+tmpsngl(1:6,1:6,istamin:istamax)
        !print *, i,stresssngl(1:6,1:1,i)
        close(1)

        write(coutfile, '(I5.5,".Displ_SH_",i5.5)') int(r0(ir0)*10.d0),irank-1

        coutfile = trim(modelname)//"."//coutfile
        coutfile = trim(outputDir)//"/Displacement/"//coutfile

        open(1,file=coutfile,status='unknown',form='unformatted', &
             access = 'direct', recl=2*3*6*kind(0e0)*nsta_global*r_n_global)
        read (1,rec=k) tmpsngl1(1:3,1:6,1:nsta_global)

        ! il ne faut pas tout stocker de 1:r_n, on doit paralleliser cette boule !!!!!!!!!!!
        displacementsngl(1:3,1:6,1:nsta,i)=displacementsngl(1:3,1:6,1:nsta,i)+tmpsngl1(1:3,1:6,istamin:istamax)

        !write(100,*) i,displacementsngl(1:3,1:6,1:1,i)
        close(1)
     endif
  enddo

  !write(100,*) 'Fin lecture'

  imin = imin_global
  imax = imax_global

  irec = 0


  ! ouverture des fichiers

   nsamp_out = iWindowEnd - iWindowStart + 1
   write(coutfile, '("Trac_",i8.8,"_",i8.8)') istamin,istamax
   coutfile = trim(parentDir)//"/"//trim(coutfile)
   open(10,file=coutfile,form='unformatted')

   write(coutfile, '("Disp_",i8.8,"_",i8.8)') istamin,istamax
   coutfile = trim(parentDir)//"/"//trim(coutfile)
   open(20,file=coutfile,form='unformatted')

   do i = 1 , nsta
     gt = cmplx(0.d0)
     !call cpu_time(tstart)
     do imt=1,6
         !write(100,*) stresssngl(1:6,imt,i,iprof,imin:imax)
        gt(1:6,imin:imax)=gt(1:6,imin:imax)+dcmplx(stresssngl(1:6,imt,i,imin:imax))*dcmplx(mt(imt))
        gt(7:9,imin:imax)=gt(7:9,imin:imax)+dcmplx(displacementsngl(1:3,imt,i,imin:imax))*dcmplx(mt(imt))
     enddo
     call cpu_time(tfinish)
        !write(100,*)
        !write(100,*) ' tenseur 1:',tfinish - tstart
       ! rotation de gt de imin a imax
     Lon_current=stlo(i)
     Lat_current=stla(i)
     SourceLat=r0lat
     SourceLong=r0lon
     call epitra(SourceLat,SourceLong,Lat_current,Lon_current,epidis,azimuth,bazi)
     !write(100,*) '---ROTATION', i, Lat_current,Lon_current
     !call cpu_time(tstart)
     call Local2Chunk(gt,9,np1,Orign_Lon_Chunk,Origin_Lat_Chunk,azi_chunk,Lat_current,Lon_current,bazi,imin,imax)
     !call cpu_time(tfinish)

     ygt(:,:) = 0.d0
     !call cpu_time(tstart)
     call tensorFFT_real(9,imin,np1,gt,ygt,omegai,tlen)
     !call cpu_time(tfinish)

     ! FILTRE
     !call cpu_time(tstart)
     !do icomp=1,9
     !  call filtbutter(ygt(icomp,:),2*np1,dt,fl,fh,norder)
     !enddo
     !call cpu_time(tfinish)
     i0 = iWindowStart
     i1 = iWindowEnd!(ibloc - 1)*nbrc + nbrc
     !call cpu_time(tstart)
     write(10) ygt(2,i0:i1)  ! debut iWindowStart)
     write(10) ygt(3,i0:i1)
     write(10) ygt(1,i0:i1)
     write(20) ygt(8,i0:i1)
     write(20) ygt(9,i0:i1)
     write(20) ygt(7,i0:i1)
     !call cpu_time(tfinish)
     !call cpu_time(tfinish)
   !enddo
   !close(10)
   !close(20)
   enddo
   close(10)
   close(20)
   close(100)
   call MPI_FINALIZE(ierr)
end program TraFFT
