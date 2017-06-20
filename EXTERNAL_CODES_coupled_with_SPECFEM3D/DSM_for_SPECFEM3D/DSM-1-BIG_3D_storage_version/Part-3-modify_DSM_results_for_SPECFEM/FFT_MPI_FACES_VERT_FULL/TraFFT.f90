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
!                                    2013.03  MONTEILLER Vadim
!-----------------------------------------------------------------------

  implicit none
  ! MPI  ---- VM VM
   include 'mpif.h'
   INTEGER myrank,nbproc,ierr,iproc! MPI  ---- VM VM


  integer :: iPS = 3 !3 ! 1:SH, 2:PSV, 3:ALL
  character(250) :: outputDir,psvmodel,modelname,stationsinf,parentDir
  integer :: itranslat
  character(250) :: coutfile
  real(kind(0d0)) :: tlen,r0min,r0max,r0delta,delta,r0lat,r0lon,r_dum

  real(kind(0d0)), allocatable :: stla(:),stlo(:)
  real(kind(0d0)), allocatable :: r_(:),r0(:),theta(:),phi(:)
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
!  real(kind(0d0)) :: samplingHz = 10.d0
  integer :: lsmooth,np0,np1

   real(kind(0d0)) samplingHz,start,ending
!  real(kind(0d0)) :: samplingHz = 10.d0
!  real(kind(0d0)),parameter :: start =  0.d0 ! Time window
!  real(kind(0d0)),parameter :: ending = 1600.d0  ! in seconds

  ! RSGT & TSGT

  ! single !!!!!!
  complex(kind(0e0)),allocatable :: stresssngl(:,:,:,:,:), displacementsngl(:,:,:,:,:),tmpsngl(:,:,:,:),tmpsngl1(:,:,:,:)
  real(kind(0d0)), allocatable :: ygt(:,:)
  !!!!!

  complex(kind(0d0)), allocatable :: gt(:,:)
  real(kind(0d0)), allocatable :: tmpygt(:)

!----
  integer max_rank,irec,irank,iprof,r_n_loc,imin_global,imax_global,iprofmin,iprofmax,idum,iidum,iidum1,iidum2,para
  integer, allocatable :: Ifrq(:),Ifrq2(:,:)
! ---- VM VM
  real(kind(0d0)) SourceLat,SourceLong,epidis,azimuth,azi,bazi,bazii,axe_rotation(3)
  real(kind(0d0)) azi_chunk,Orign_Lon_Chunk,Origin_Lat_Chunk,Lat_current,Lon_current,delta_lat,delta_lon
  complex(kind(0d0)) normal(3),normal_ch1(3),Rot_azi_chunk(3,3),Rot_Lon_chunk(3,3),Rot_Lat_chunk(3,3)
  complex(kind(0.d0)) Rot_ch0_ch1(3,3),Rot_inv_azi_chunk(3,3),RL2L(3,3),RL2Li(3,3)
  character(len=11) debug
  integer norder,icomp
  real(kind(0d0)) fgauss
  integer ibloc,nbrc,Nbloc,i0,i1,nsamp_out
  real :: tstart, tfinish
  integer ista,iproff,kkk,size_of_array
  character(len=4) face
!-----------------------------------------------------------------------
 call MPI_INIT(ierr)
 call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
 call MPI_COMM_SIZE(MPI_COMM_WORLD,nbproc,ierr)


!! ---parametres
   Nbrc = 100
   norder = 6
!   fl  =  0.05d0
!   fh = 0.9d0
!!----to do ici
  write(debug,'(a6,i5.5)') 'debug_',myrank
  if (myrank == 0) open(100,file=debug)

 ! lecture du fichier split en frequences
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
  if (myrank == 0) then
    call pinputTra(outputDir,psvmodel,modelname,stationsinf,tlen,imin_global,imax_global,r0min,r0max, &
    r0delta,r0lat,r0lon,itranslat,mt,myrank,samplingHz,start,ending,fgauss)
  endif
  call MPI_Bcast(outputDir,250,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(psvmodel,250,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(modelname,250,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(stationsinf,250,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(tlen,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(r0min,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(r0max,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(r0delta,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(r0lat,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(r0lon,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(imin_global,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(imax_global,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(mt,6,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(samplingHz,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(start,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ending,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(fgauss,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!  call MPI_Bcast(fh,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  itranslat=0
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


  face=stationsinf(3:6)   ! a verifier
  call InputForTraction(Orign_Lon_Chunk,Origin_Lat_Chunk,azi_chunk,normal,delta_lat,delta_lon,face)

  open(2,file='recdepth')
  read(2,*) r_n_global
  ! // isation des profondeurs
  call DistribDepth(r_n,r_n_global,myrank,nbproc,iprofmin,iprofmax)
  if (myrank == 0) then
   write(100,*) 'r_n ', r_n
   write(100,*) 'irpofmin, iprofmax ',iprofmin,iprofmax
   write(100,*) 'imin_global, imax_global ',imin_global, imax_global
  endif
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  ! on lit les stations
  open (1,file=stationsinf,status='old',action='read',position='rewind')
  read(1,*)nsta
  theta_n = nsta
  allocate(r_(1:r_n))
  allocate(theta(1:theta_n))
  allocate(phi(1:theta_n))
  allocate(stla(1:nsta))
  allocate(stlo(1:nsta))
  allocate(updown(1:r_n))

  do i = 1,nsta
     read(1,*) stlo(i),stla(i)
     !if (itranslat==1) call translat(stla(i),stla(i))
     call epitra1(r0lat,r0lon,stla(i),stlo(i),theta(i),phi(i),bazi)
  enddo
  write(*,*) 'PROFONDEURS TRAITEES PAR MYRANK'
  k  =  0
  do i=1,r_n_global
    read(2,*) r_dum,idum,updown_dum
    if (i >= iprofmin .and. i <= iprofmax) then
       k  = k + 1
       r_(k) = 6371.d0 -r_dum
       updown(k) = updown_dum
       write(*,*) k,r_(k)
    endif
  enddo
  close(1)
  close(2)

  ! je ne sais pas ce que c'est VM VM
  r0_n =  int((r0max-r0min)/r0delta)+1
  allocate(r0(1:r0_n))
  do i = 1, r0_n
     r0(i) = r0min + dble(i-1)*r0delta
  enddo

  ! lsmoothfinder for FFT  !! ca non plus VM VM

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
    write(100,*) 'np1 :', np1
    write(100,*) 'sampling,dt : ', samplingHz,dt
    write(100,*)
    write(100,*) 'Window :', iWindowStart, iWindowEnd
    write(100,*) 'time w:', (iWindowStart-1)*dt,(iWindowEnd-1)*dt
    write(100,*)
  endif
  if (myrank == 0) then
    open(11,file='NbSmples.txt')
    write(11,*) iWindowEnd - iWindowStart + 1
    close(11)
  endif

  allocate(gt(1:9,0:2*np1-1))
  allocate(ygt(1:9,0:2*np1-1))
  allocate(tmpygt(0:2*np1-1))
  if (myrank == 0) write(100,*) 'ALLOCATION GROS TABLEAUX'
  allocate(stresssngl(1:6,1:6,1:nsta,1:r_n,imin_global:imax_global))
  allocate(displacementsngl(1:3,1:6,1:nsta,1:r_n,imin_global:imax_global))

  allocate(tmpsngl(1:6,1:6,1:r_n_global,1:nsta))
  allocate(tmpsngl1(1:3,1:6,1:r_n_global,1:nsta))
  if (myrank == 0) then
    write(100,*) 'nsta =',nsta
    write(100,*) 'r_n =',r_n
    write(100,*) 'imin_global, imax_global ',imin_global,imax_global
    write(100,*) 'r_n_global',r_n_global
    write(100,*)
  endif
  stresssngl=cmplx(0.e0)
  displacementsngl=cmplx(0.e0)
  ir0 = r0_n

  !calcul de la matrice de rotation pour la normale, cette matrice passe du chunk de reference
  ! au chunk positionne sur la zone d'etude
  ! rotation azimuth (axe Z)
  axe_rotation(1)=0.d0; axe_rotation(2)=0.d0; axe_rotation(3)=1.d0
  call rotation_matrix(Rot_azi_chunk,axe_rotation,Azi_Chunk)
  ! rotation inverse azimuth (axe Z)
  axe_rotation(1)=0.d0; axe_rotation(2)=0.d0; axe_rotation(3)=-1.d0
  call rotation_matrix(Rot_inv_azi_chunk,axe_rotation,Azi_Chunk)

  !! je crois que ce qui suit ne me servira pas
  ! rotation latitude (axe Y) (!!! il y un moins)
  !axe_rotation(1)=0.d0; axe_rotation(2)=-1.d0; axe_rotation(3)=0.d0
  !call rotation_matrix(Rot_lat_chunk,axe_rotation,Origin_Lat_Chunk)
  ! rotation longitude (axe Z)
  !axe_rotation(1)=0.d0; axe_rotation(2)=0.d0; axe_rotation(3)=1.d0
  !call rotation_matrix(Rot_lon_chunk,axe_rotation,Orign_Lon_Chunk)
  ! la matrice cherchee est la composition des 3 matrices precedentes
  !call compose3matrix(Rot_ch0_ch1,Rot_azi_chunk,Rot_lat_chunk,Rot_lon_chunk)
  ! normale exprimee dans le repere du chunk
  !!!!call matmulvect(normal_ch1,Rot_ch0_ch1,normal) ! ce truc est faux
  ! j'ai directement la normale dans le repre du chunk
  !!!!!!!!!!!!!!!!!!! FIN des choses inutiles

  ! initialisation
  irank = -1
  ! on synchronise les procs
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if (myrank == 0) then
   write(100,*) 'lecture des spectres '
   write(100,*)
  endif
  k=1
  do irank=1,imax_global + 1 ! max_rank-1 !------------------------------------------------------------------------------ FRQS : IRANK

   if (myrank == 0) then
     write(100,*)
     write(100,*) '-----'
     write(100,*) ' lecture spectre calcule par myrank=',irank
   endif
   i = irank - 1
     if (iPS /= 1) then ! PSV calculation
        write(coutfile, '(I5.5,".Stress_PSV_",i5.5)') int(r0(ir0)*10.d0),irank-1
!! VM VM supress one file and change read statement
        coutfile = trim(modelname)//"."//coutfile
        coutfile = trim(outputDir)//"/Stress/"//coutfile
        if (myrank == 0) then
          !open(1,file=coutfile,status='unknown',form='unformatted', &
          !    access = 'direct', recl=2*6*6*kind(0e0)*nsta*r_n_global)
          !read(1,rec=k) tmpsngl(1:6,1:6,1:r_n_global,1:nsta)
          !close(1)
           open(1,file=coutfile,form='unformatted',action='read')
           read(1)  tmpsngl
           write(100,*)
           write(100,*)  'rec ,frq, ',k,i
           write(100,*)  'nb rec,profmin profmax :', r_n,iprofmin,iprofmax
           write(100,'(a)') trim( coutfile)
        endif
        size_of_array=6*6*r_n_global*nsta  !  il faudrait faire un mpi_scatterv
        !call MPI_Bcast()
        call MPI_Bcast(tmpsngl,size_of_array,MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)
        do ista=1,nsta
           kkk=0
           do iproff=iprofmin,iprofmax
              kkk=kkk+1
              stresssngl(1:6,1:6,ista,kkk,i)=stresssngl(1:6,1:6,ista,kkk,i)+tmpsngl(1:6,1:6,iproff,ista)
           enddo
        enddo


        !write(coutfile, '(I5.5,".Displ_PSV_",i5.5)') int(r0(ir0)*10.d0),irank-1

        !coutfile = trim(modelname)//"."//coutfile
        !coutfile = trim(outputDir)//"/Displacement/"//coutfile

        if (myrank == 0) then
          !open(1,file=coutfile,status='unknown',form='unformatted', &
          !     access = 'direct', recl=2*3*6*kind(0e0)*nsta*r_n_global)
          !read (1,rec=k) tmpsngl1(1:3,1:6,1:r_n_global,1:nsta)
           read(1)  tmpsngl1
           close(1)
        endif
        size_of_array=3*6*r_n_global*nsta
        call MPI_Bcast(tmpsngl1,size_of_array,MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)
        kkk=0
        do iproff=iprofmin,iprofmax !! VM VM : It's seems to be useless to perfrom sumation because of initilisation of displacementsngl??
           kkk=kkk+1
           do ista=1,nsta
              displacementsngl(1:3,1:6,ista,kkk,i)=displacementsngl(1:3,1:6,ista,kkk,i)+tmpsngl1(1:3,1:6,iproff,ista)
           enddo
        enddo
        if (myrank == 0) write(100,'(a)') trim(coutfile)

     endif
     if (iPS /= 2) then ! SH calculation  ! TO DO changer l'ordre ista,irpof dans tmpsngl(1)
        write(coutfile, '(I5.5,".Stress_SH_",i5.5)') int(r0(ir0)*10.d0),irank-1
  coutfile = trim(modelname)//"."//coutfile
  coutfile = trim(outputDir)//"/Stress/"//coutfile
  if (myrank == 0) then

    open(1,file=coutfile,form='unformatted',action='read')
    read(1)  tmpsngl
!   write(100,*)
!   write(100,*)  'rec ,frq, ',k,i
!   write(100,*)  'nb rec,profmin profmax :', r_n,iprofmin,iprofmax
!   write(100,'(a)') trim( coutfile)
  endif
  size_of_array=6*6*r_n_global*nsta  !  il faudrait faire un mpi_scatterv
  !call MPI_Bcast()
  call MPI_Bcast(tmpsngl,size_of_array,MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)
  do ista=1,nsta
     kkk=0
     do iproff=iprofmin,iprofmax
        kkk=kkk+1
        stresssngl(1:6,1:6,ista,kkk,i)=stresssngl(1:6,1:6,ista,kkk,i)+tmpsngl(1:6,1:6,iproff,ista)
     enddo
  enddo

  if (myrank == 0) then
     read(1)  tmpsngl1
     close(1)
  endif
  size_of_array=3*6*r_n_global*nsta
  call MPI_Bcast(tmpsngl1,size_of_array,MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)
  kkk=0
  do iproff=iprofmin,iprofmax !! VM VM : It's seems to be useless to perfrom sumation because of initilisation of displacementsngl??
     kkk=kkk+1
     do ista=1,nsta
        displacementsngl(1:3,1:6,ista,kkk,i)=displacementsngl(1:3,1:6,ista,kkk,i)+tmpsngl1(1:3,1:6,iproff,ista)
           enddo
        enddo
  if (myrank == 0) write(100,'(a)') trim(coutfile)

     endif
     call MPI_Barrier(MPI_COMM_WORLD,ierr) !! test
  enddo

   if (myrank == 0) write(100,*) 'Fin lecture'

  imin = imin_global
  imax = imax_global

  irec = 0

  !
  do iprof = 1, r_n ! valeur locale attribuee a myrank

     ! 1 fichier par profondeur et toutes les distances epicentrales
     write(coutfile, '("Trac_",i5.5)') iprofmin+iprof-1
     coutfile = trim(parentDir)//"/"//trim(coutfile)
     open(10,file=coutfile,form='unformatted')
     write(coutfile, '("Disp_",i5.5)') iprofmin+iprof-1
     coutfile = trim(parentDir)//"/"//trim(coutfile)
     open(20,file=coutfile,form='unformatted')
     write(*,*) 'myrank :',myrank,trim(coutfile)
     do i = 1,nsta ! toutes distances epicentrales

        gt = cmplx(0.d0)
        call cpu_time(tstart)
        do imt=1,6
           gt(1:6,imin:imax)=gt(1:6,imin:imax)+stresssngl(1:6,imt,i,iprof,imin:imax)*mt(imt)
           gt(7:9,imin:imax)=gt(7:9,imin:imax)+displacementsngl(1:3,imt,i,iprof,imin:imax)*mt(imt)
        enddo
        call cpu_time(tfinish)

        ! rotation de gt de imin a imax
        Lat_current=stla(i)
        Lon_current=stlo(i)
        SourceLat=r0lat
        SourceLong=r0lon


        ! calcul de la distance epicentral source-station
        call epitra(SourceLat,SourceLong,Lat_current,Lon_current,delta,azi,bazi)
        if (myrank == 0) write(100,*) '---ROTATION ', i, Lon_current,Lat_current,delta

        ! matrice de rotation point courrant -> repere local du chunk
        !call local2localMatrix(Lat_current,Lon_current,Origin_Lat_Chunk,Orign_Lon_Chunk,RL2L)
        ! matrice inverse
        !call local2localMatrix(Origin_Lat_Chunk,Orign_Lon_Chunk,Lat_current,Lon_current,RL2Li)
        !call cpu_time(tstart)
        !! on doit faire la rotation d'azimuth
        ! on calcule la traction et la vitesse dans le repere du chunk
        call Local2Chunk(gt,9,np1,Origin_Lat_Chunk,Orign_Lon_Chunk,Lat_current,lon_current,normal,imin,imax,bazi,azi_chunk)
        !call cpu_time(tfinish)
        ! tranformee de Fourier

        ygt(:,:) = 0.d0
        call cpu_time(tstart)
        call tensorFFT_real(9,imin,np1,gt,ygt,omegai,tlen)
        !call cpu_time(tfinish)

        ! FILTRE
        !call cpu_time(tstart)
        !do icomp=1,9
        !  call filtbutter(ygt(icomp,:),2*np1,dt,fl,fh,norder)
        !enddo
        !call cpu_time(tfinish)

        nsamp_out = iWindowEnd - iWindowStart + 1
        Nbloc = int(nsamp_out / nbrc)
        !call cpu_time(tfinish)

        ! ecriture sur le disque
        i0 = iWindowStart
        i1 = iWindowEnd
        write(10) ygt(2,i0:i1)
        write(10) ygt(3,i0:i1)
        write(10) ygt(1,i0:i1)
        write(20) ygt(8,i0:i1)
        write(20) ygt(9,i0:i1)
        write(20) ygt(7,i0:i1)

   enddo
   close(10)
   close(20)
  enddo
  if (myrank == 0) close(100)
  call MPI_FINALIZE(ierr)
end program TraFFT
