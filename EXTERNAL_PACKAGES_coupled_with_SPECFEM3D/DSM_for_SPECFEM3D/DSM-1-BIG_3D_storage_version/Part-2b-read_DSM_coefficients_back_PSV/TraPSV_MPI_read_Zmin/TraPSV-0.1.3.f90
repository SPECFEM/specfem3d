
  program TraPSV

  !-----------------------------------------------------------------------
  !      TraPSV
  !
  !    calculation of the Green function for P wave inversion
  !
  !          2004.10.KAWAI Kenji
  !          2009.6. FUJI Nobuaki
  !          2011.9. FUJI Nobuaki
  !          2011.10 MONTEILLER Vadim (MPI)
  !          2013.03 MONTEILLER Vadim (MPI_SPLIT_COMM)
  !          2013.03 KOMATITSCH Dimitri (full loop vectorization, and suppressed useless memory copies in subroutine calls)
  !
  !-----------------------------------------------------------------------

  implicit none

   ! MPI  ---- VM VM
   include 'mpif.h'
   include '../../shared/constants.h'
   INTEGER myrank,nbproc,ierr
   integer SubCommunicators,mybigrank,nbbigproc
   integer, allocatable :: key(:),color(:)

  !------------------------- <  < input matrix >>----------------------------------

  character(120) :: inputdir,outputDir,psvmodel,modelname,stationsinf
  character(120) :: list,list1
  character(40) :: datex,timex
  real(kind(0d0)), parameter :: pi=3.1415926535897932d0
  real(kind(0d0)), parameter :: re=1.d-2, ratc=1.d-10, ratl=1.d-4
  !integer, parameter :: maxlmax = 35000, maxlmax_g=1000 !! VM VM I moved this in constants.h

  real(kind(0d0)) :: tlen
  real(kind(0d0)) :: r0min, r0max, r0delta  !!! JUST FOR ONE DEPTH FOR THE MOMENT !!
  real(kind(0d0)) :: r0lat, r0lon, stla_curr,stlo_curr
  real(kind(0d0)), allocatable :: stla(:),stlo(:),stla_g(:),stlo_g(:)
  real(kind(0d0)), allocatable :: r_(:),r0(:),theta(:),phi(:),phi_g(:),theta_g(:)
  real(kind(0d0)), allocatable :: A0sta(:),C0sta(:),F0sta(:),L0sta(:),N0sta(:)
  integer, allocatable :: updown(:),idum(:)
  real(kind(0d0)), allocatable :: rrsta(:,:)
  integer, allocatable :: iista(:,:)
  integer :: r_n,r0_n,ir_,ir0,itheta,theta_n,nsta
  integer :: imt,ier

  character(120) :: cinfile,coutfile
  integer :: itranslat


  ! --------------------------- <  < variables >>---------------------------
  ! variable for the trial function
  integer:: nlay
  integer, allocatable :: nlayer(:)
  integer:: l,m
  ! variable for the structure
  integer:: nzone,isl,ill,nsl,nll
  integer,allocatable:: iphase(:)
  integer::ndc
  real(kind(0d0)):: rmin,rmax
  real(kind(0d0)),allocatable:: vrmin(:),vrmax(:),rrho(:,:),vpv(:,:),vph(:,:),vsv(:,:), &
  vsh(:,:),eta(:,:),qmu(:),qkappa(:)
  complex(kind(0d0)),allocatable:: coef1(:),coef2(:),coef(:)
  real(kind(0d0)):: omega,omegai

  ! variable for the matrix elements
  complex(kind(0d0)),allocatable:: tabg0(:,:,:,:),tabg0der(:,:,:,:),tmp_tabg0(:,:,:,:),tmp_tabg0der(:,:,:,:)
  complex(kind(0d0)):: g0tmp(2,6),g0dertmp(2,6) ! forward
  complex(kind(0d0)):: g0tmp1,g0dertmp1 ! forward
  ! variable for the stack point
  integer,allocatable:: isp(:),issp(:),ilsp(:),jssp(:),jsp(:), ksp(:),lsp(:),lmax_r(:)
  ! variables for the output stack point
  integer,allocatable:: istazone(:)
  integer,allocatable:: ksta(:)   ! output stack point for g
  integer,allocatable:: jsta(:)   ! output stack point for d
  ! variables for the gridding
  integer,allocatable:: jjdr(:),kkdr(:)
  real(kind(0d0)),allocatable:: vmin(:),gridpar(:),dzpar(:)
  ! variables for l cut off
  integer:: llog,index_l,lref,lmax_lu

  ! other variables
  integer:: i,ig2,ista
  real(kind(0d0)):: l2,lsq,bazi

  !-----------------------------------------------------------------------
  complex(kind(0d0)), allocatable :: dvec(:,:,:),dvecdt(:,:,:),dvecdp(:,:,:)
  complex(kind(0d0)), allocatable :: stress(:,:,:), displacement(:,:,:)
  complex(kind(0e0)), allocatable :: stresssngl(:,:,:), displacementsngl(:,:,:)
  complex(kind(0e0)), allocatable :: stresssngl_global(:,:,:), displacementsngl_global(:,:,:)
  real(kind(0d0)), allocatable :: plm(:,:,:)
  integer r_n_dum

!! DK DK added this to precompute several things
  complex(kind(0d0)) :: inv_lsq
  real(kind(0d0)) :: Aval,Cval,Fval,Lval,Nval
  complex(kind(0d0)) :: u1,u2,u3,udr1,udr2,udr3,udt1,udt2,udt3,udp1,udp2,udp3
  complex(kind(0d0)) :: uder11 , uder12 , uder13 , uder21 , uder22 , uder23 , uder31 , uder32 , uder33
  real(kind(0d0)) :: rval,thetaval
  real(kind(0d0)) :: thetasin
  complex(kind(0d0)) :: thetacot_complex, inv_rval, inv_thetasin, inverse_of_omega_complex

!! DK DK added this conversion factor
  double precision, parameter :: convert_to_radians = pi/180.d0

  ! ajout de variables pour MPI --- VM VM
  integer i_color,nb_colors,current_color,imin_current,imax_current,nbproc_current
  integer iread,nb_tot_proc
  integer, allocatable :: Ifrq2(:,:),Ind(:),NbInd(:),IndSt(:),NbIndSt(:)
  integer istamin,istamax,k,ifq,c8,coeff
  integer imin_glob, imax_glob, irank,nsta_global
  integer ifrequ_min, ifrequ_max
  real time_total,time_to_read,time_to_project,time_to_plm
  real start_time_total,start_time_to_read,start_time_to_project,start_time_to_plm
  real finish_time_total,finish_time_to_read,finish_time_to_project,finish_time_to_plm
  real total_global_time
  character(len=250) debug_file
 !--------------------------------------------------------------------------

  ! MPI --- VM VM
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mybigrank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nbbigproc,ierr)

  if (USE_TIMER) then
     call cpu_time(start_time_total)
     time_to_read=0.
     time_to_project=0.
     time_to_plm=0.
     time_total=0.
  endif

  allocate(key(0:nbbigproc-1))
  allocate(color(0:nbbigproc-1))
  allocate(Ifrq2(2,0:nbbigproc-1))
  if (mybigrank == 0) then
    open(25,file='Double_para.txt')
    read(25,*) nb_colors,nb_tot_proc
    if (nbbigproc /= nb_tot_proc) then
       write(*,*) 'bad number MPI proc'
       stop
     endif
    do i_color=1,nb_colors
       read(25,*) current_color,imin_current,imax_current,nbproc_current
       do iread=1,nbproc_current
          read(25,*) irank
          color(irank)=current_color
          Ifrq2(1,irank)=imin_current
          Ifrq2(2,irank)=imax_current
          if (iread == 1) then
             key(irank)=-1
          else
             key(irank)=irank
          endif
       enddo
    enddo
  endif
  call MPI_Bcast(key,nbbigproc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(color,nbbigproc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(Ifrq2,2*nbbigproc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  if (mybigrank == 0) then
   call pinputTra_part2(ifrequ_min,ifrequ_max,inputdir,outputDir,psvmodel,modelname,stationsinf,tlen, &
               imin_glob,imax_glob,r0min,r0max,r0delta,r0lat,r0lon,itranslat,myrank)
  endif
  call MPI_Bcast(inputdir,120,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
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
  call MPI_Bcast(ifrequ_min,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ifrequ_max,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  itranslat=0
  ifrequ_min = Ifrq2(1,mybigrank)
  ifrequ_max = Ifrq2(2,mybigrank)

  if (mybigrank == 0) then
  open(20, file = psvmodel, status = 'old', action='read', position='rewind')
  read(20,*) nzone
  endif

  call MPI_Bcast(nzone,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  allocate(nlayer(1:nzone))
  allocate(iphase(1:nzone))
  allocate(vrmin(1:nzone))
  allocate(vrmax(1:nzone))
  allocate(rrho(1:4,1:nzone))
  allocate(vpv(1:4,1:nzone))
  allocate(vph(1:4,1:nzone))
  allocate(vsv(1:4,1:nzone))
  allocate(vsh(1:4,1:nzone))
  allocate(eta(1:4,1:nzone))
  allocate(qmu(1:nzone))
  allocate(qkappa(1:nzone))
  allocate(coef1(1:nzone))
  allocate(coef2(1:nzone))
  allocate(coef(1:nzone))
  allocate(jjdr(1:nzone))
  allocate(kkdr(1:nzone))
  allocate(vmin(1:nzone))
  allocate(gridpar(1:nzone))
  allocate(dzpar(1:nzone))
  allocate(isp(1:nzone))
  allocate(issp(1:nzone))
  allocate(ilsp(1:nzone))
  allocate(jssp(1:nzone))
  allocate(jsp(1:nzone))
  allocate(ksp(1:nzone))
  allocate(lsp(1:nzone))
  if (mybigrank == 0) then
  do i = 1, nzone
     read (20, *) vrmin(i), vrmax(i), rrho(1,i), rrho(2,i), rrho(3,i), rrho(4,i), vpv(1,i), vpv(2,i), vpv(3,i), &
     vpv(4,i), vph(1,i), vph(2,i), vph(3,i), vph(4,i), vsv(1,i), vsv(2,i), vsv(3,i), vsv(4,i), vsh(1,i), &
     vsh(2,i), vsh(3,i), vsh(4,i), eta(1,i), eta(2,i), eta(3,i), eta(4,i), qmu(i), qkappa(i)
  enddo
  close(20)
  endif

  call MPI_Bcast(vrmin,nzone,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(vrmax,nzone,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(rrho,4*nzone,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(vpv,4*nzone,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(vph,4*nzone,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(vsv,4*nzone,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(vsh,4*nzone,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(eta,4*nzone,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(qmu,nzone,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(qkappa,nzone,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  rmin = vrmin(1)
  rmax = vrmax(nzone)
  omegai = - dlog(1.d-2) / tlen

  ! --- computing the required parameters ---
  ! counting of the nsl and nll
  call calnl( nzone,vsv,iphase,nsl,nll )
  ndc = nzone - 1

  if (mybigrank == 0) then
   open (1,file=stationsinf,status='old',action='read',position='rewind')
   open(2,file='recdepth')
   read(2,*) r_n
   read(1,*)nsta_global
  endif
  call MPI_Bcast(r_n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(nsta_global,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  allocate(r_(1:r_n))
  allocate(lmax_r(r_n))
  allocate(A0sta(1:r_n))
  allocate(C0sta(1:r_n))
  allocate(F0sta(1:r_n))
  allocate(L0sta(1:r_n))
  allocate(N0sta(1:r_n))

  allocate(theta_g(1:nsta_global))
  allocate(phi_g(1:nsta_global))
  allocate(stla_g(1:nsta_global))
  allocate(stlo_g(1:nsta_global))

  allocate(updown(1:r_n))
  allocate(idum(1:r_n))

! it is sufficient for us to save the final result to disk in single precision
  allocate(stresssngl_global(1:6,1:6,1:nsta_global))
  allocate(displacementsngl_global(1:3,1:6,1:nsta_global))

  allocate(tabg0(1:maxlmax_g,2,1:6,-2:2))
  allocate(tabg0der(1:maxlmax_g,2,6,-2:2))
  allocate(tmp_tabg0(1:maxlmax_g,2,1:6,-2:2))
  allocate(tmp_tabg0der(1:maxlmax_g,2,6,-2:2))

  if (mybigrank == 0) then
     write(*,*)
     write(*,*) '   PROCESS ZMIN FACE :',r_n
     write(*,*)
     write(*,*)
     write(*,*) '**********************************'
     write(*,*)
     write(*,*) 'SOURCE :',r0lon,r0lat
     write(*,*)
     write(*,*) '**********************************'
     do i = 1,nsta_global
        read(1,*) stlo_curr,stla_curr
        stla_g(i)=stla_curr
        stlo_g(i)=stlo_curr
        call epitra(r0lat,r0lon,stla_g(i),stlo_g(i),theta_g(i),phi_g(i),bazi)
     enddo
  endif

  call MPI_Bcast(theta_g,nsta_global,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(phi_g,nsta_global,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(stla_g,nsta_global,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(stlo_g,nsta_global,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  if (mybigrank == 0) then
  do i=1,r_n
     read(2,*) r_(i),idum(i),updown(i)
     r_(i) = 6371.d0 -r_(i)
  enddo
  close(1)
  close(2)
  endif

  call MPI_Bcast(r_,r_n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(updown,r_n,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(idum,r_n,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  do i=1,r_n
     call calstg4onedepth(nlay,nzone,nzone,vrmin,vrmax,iphase,rrho,vpv,vph,vsv,vsh,eta,rmax,r_(i), &
     updown(i),A0sta(i),C0sta(i),F0sta(i),L0sta(i),N0sta(i),idum(i))
  enddo

! depths for stocking the Green function
  allocate(rrsta(1:3,1:r_n))
  allocate(iista(1:3,1:r_n))
  allocate(istazone(1:r_n))
  allocate(jsta(1:r_n))
  allocate(ksta(1:r_n))

  !
  ! ######################### sous communicateurs #####################"
  ! defition de l'ensemble des sous communicateurs
  call MPI_COMM_SPLIT(MPI_COMM_WORLD,color(mybigrank),mybigrank,subcommunicators,ierr)
  call MPI_COMM_RANK(SubCommunicators,myrank,ierr)
  call MPI_COMM_SIZE(SubCommunicators,nbproc,ierr)

! distribuer les stations sur les procs
  call DistribSta(nsta,nsta_global,myrank,nbproc,istamin,istamax)
  write(*,*) 'myrank  ista :', myrank, istamin,istamax

  call MPI_Barrier(SubCommunicators,ierr)

  allocate(Ind(nbproc))
  allocate(NbInd(nbproc))
  allocate(IndSt(nbproc))
  allocate(NbIndSt(nbproc))
  call mpi_allgather(nsta,1,MPI_INTEGER,NbInd,1,MPI_INTEGER,SubCommunicators, ierr )

  NbIndSt(:) = NbInd(:)

  ! indices et offsets pour mpi_gatherv
  r_n_dum=1  ! ici valeur de 1 car on ne stocke que la valeur du fond de la boite
  ! displ
  coeff  = 3*6*r_n_dum
  Ind(1) = 0
  do k = 2, nbproc
      c8 = NbInd(k-1)
      Ind(k) = Coeff*c8 +   Ind(k-1)
  enddo

  do k=1,nbproc
     NbInd(k) = NbInd(k) * Coeff
  enddo

  ! stress
  coeff  = 6*6*r_n_dum
  IndSt(1) = 0
  do k = 2, nbproc
     c8 = NbIndSt(k-1)
     IndSt(k) = Coeff*c8 +   IndSt(k-1)
  enddo

  do k=1,nbproc
     NbIndSt(k) = NbIndSt(k) * Coeff
  enddo


  theta_n = nsta
  allocate(theta(1:theta_n))
  allocate(phi(1:theta_n))
  allocate(stlo(theta_n))
  allocate(stla(theta_n))

  allocate(dvec(1:theta_n,1:3,-2:2))
  allocate(dvecdt(1:theta_n,1:3,-2:2))
  allocate(dvecdp(1:theta_n,1:3,-2:2))
  allocate(plm(1:theta_n,1:3,0:3),stat=ier)
  if (ier /= 0) stop 'error: not enough memory to allocate array'

  phi(:)=phi_g(istamin:istamax)
  theta(:)=theta_g(istamin:istamax)
  stla(:)=stla_g(istamin:istamax)
  stlo(:)=stlo_g(istamin:istamax)

  allocate(stress(1:nsta,1:6,1:6))
  allocate(displacement(1:nsta,1:3,1:6))

! it is sufficient for us to save the final result to disk in single precision
  allocate(stresssngl(1:6,1:6,1:nsta))
  allocate(displacementsngl(1:3,1:6,1:nsta))

  ! source depths
  ! currently we do not put any necessary allocation
  r0_n =  int((r0max-r0min)/r0delta)+1
  allocate(r0(1:r0_n))
  do i = 1, r0_n
     r0(i) = r0min + dble(i-1)*r0delta
  enddo
  ir0 = r0_n

  isl = 0
  ill = 0

  if (myrank == 0 .and. SLOW_DEBUG_MODE) then
    ! Record the date and time at the beginning of the job
    write(list, '(I6.6,".",I6.6)') ifrequ_min, ifrequ_max
    list = trim(outputDir)//"/log/calLog"//"."//trim(modelname)//"."//trim(list)
    open(25,file =list, status = 'unknown', form = 'formatted')
    call date_and_time(datex,timex)
    write(25,'(/a,a4,a,a2,a,a2,a,a2,a,a2,a,a4)') &
         '    Starting date and time:                     ', &
         datex(1:4),'-',datex(5:6),'-',datex(7:8),'.  ', &
         timex(1:2),':',timex(3:4),':',timex(5:8)


    write(list1, '(I6.6,".",I6.6)') ifrequ_min, ifrequ_max
    list1 = trim(outputDir)//"/log/list"//"."//trim(modelname)//"."//trim(list1)
    open(24, file = list1, status = 'unknown', form = 'formatted')
    write(24,*)
    write(24,*) 'nsta ', nsta
    write(24,*) 'phi  ',   phi(:)
    write(24,*) 'theta ', theta(:)
    write(24,*) 'stla ',   stla(:)
    write(24,*) 'stlo ',   stlo(:)
  endif
  if (SLOW_DEBUG_MODE) then
     write(debug_file,'("deb_",i5.5)') mybigrank
     open(36,file=trim(debug_file))
     !write(36,*) 'myrank  ista :', myrank, istamin,istamax
    ! write(36,*) 'nsta ', nsta
    do i=1,nsta
       write(36,'(i7,4f20.6)') i,phi(i),theta(i),stla(i),stlo(i)
    ! write(36,*) 'phi  ',   phi(:)
     ! write(36,*) 'theta ', theta(:)
     !write(36,*) 'stla ',   stla(:)
     !write(36,*) 'stlo ',   stlo(:)
    enddo
   close(36)
  endif
  llog = 0
  ! boucle sur les freqs
  do ifq=ifrequ_min+1, ifrequ_max+1
     i=ifq-1
     !imin=Ifrq(ifq)+1
     !imax=ifrq(ifq+1)
     time_to_read=0.
     time_to_project=0.
     time_to_plm=0.

!! DK DK now read displacement and stress from the same file to avoid opening too many files on the file system, which is usually slow
     write(cinfile, '(I5.5,".G0_PSV_",I5.5)') int(r0(ir0)*10.d0),ifq-1
     cinfile = trim(modelname)//"."//cinfile
     cinfile = trim(inputdir)//"/Stress/"//cinfile  ! TO DO : il faudrait mettre un inputdir
     if (myrank == 0) then
        if (SLOW_DEBUG_MODE) write(24,*) myrank,ifq,trim(cinfile)
        open(34,file=cinfile, form='unformatted',action='read')
     endif

     ir0 = r0_n

     ! initialisation
     stress = (0.d0,0.d0)
     displacement = (0.d0,0.d0)
     stresssngl = (0.e0,0.e0)
     displacementsngl = (0.e0,0.e0)
     stresssngl_global = (0.e0,0.e0)
     displacementsngl_global = (0.e0,0.e0)

     plm = 0.d0
     dvec = (0.d0,0.d0)
     dvecdp = (0.d0,0.d0)
     dvecdt = (0.d0,0.d0)

     lmax_r(:)=0
     omega = 2.d0 * pi * dble(i) / tlen

     if ( i /= 0 ) then
        lref=0
        index_l=0
        lmax_lu=0

! 33333333333333333333333333333333333333333333333333333333333333333333333333333333333333
        do l=0,maxlmax  ! on parcourt tous les l a priori
! 33333333333333333333333333333333333333333333333333333333333333333333333333333333333333

           call MPI_Barrier(SubCommunicators,ierr)  ! je veux que tous les procs aillent a la meme vitesse

           ! lecture des tableaux
           if (mod(l,maxlmax_g) == 0) then
              lref=l
              tabg0=(0.d0,0.d0)
              tabg0der=(0.d0,0.d0)
              if (myrank == 0) then
                if (SLOW_DEBUG_MODE) then
                  write(24,*)
                  write(24,*) 'reading l = ',l
                endif
                if (USE_TIMER) call cpu_time(start_time_to_read)

                do

                   read(34) ir_,llog

                   if (ir_ == -1) then ! this flag indicates the end of all the data to read
                      exit
                   else if (ir_ == r_n) then
! in the case of the Zmin face, we read but ignore all the data from the input file except the last one;
! in the other code (for vertical faces) we use all the data read, but in this code we use the last one only from the same file
                      read(34) tabg0(1:llog,:,:,:)
                      read(34) tabg0der(1:llog,:,:,:)
                      lmax_r(ir_) = llog + lref ! lmax courant pour ir_
                      lmax_lu=max(lmax_lu,lmax_r(r_n))  ! on ne prend le l que le zmin
                   else !!! VM VM : we need ir_=r_n which is not nessesary the last one thus I test it
                      read(34) tmp_tabg0(1:llog,:,:,:)
                      read(34) tmp_tabg0der(1:llog,:,:,:)
                      lmax_r(ir_) = llog + lref ! lmax courant pour ir_
                      lmax_lu=max(lmax_lu,lmax_r(r_n))
                   endif
                enddo

                if (USE_TIMER) then
                  call cpu_time(finish_time_to_read)
                  time_to_read = time_to_read + finish_time_to_read - start_time_to_read
                 endif
              endif

              index_l=0
              call MPI_Bcast(lmax_r,r_n,MPI_INTEGER,0,SubCommunicators,ierr)
              call MPI_Bcast(tabg0,60*maxlmax_g,MPI_DOUBLE_COMPLEX,0,SubCommunicators,ierr)
              call MPI_Bcast(tabg0der,60*maxlmax_g,MPI_DOUBLE_COMPLEX,0,SubCommunicators,ierr)
              call MPI_Bcast(ir_,1,MPI_INTEGER,0,SubCommunicators,ierr)
              call MPI_Bcast(llog,1,MPI_INTEGER,0,SubCommunicators,ierr)
              call MPI_Bcast(time_to_read,1,MPI_REAL,0,SubCommunicators,ierr)
           endif
           if (l >= lmax_r(r_n)) exit ! si jamais on atteind le lmax du dernier ir_
           if (ir_ == -1 .and. llog == -1) exit  !! VM VM: I add this test because of bugs
                                               !! VM VM: when mod(maxval(lmax_r),maxlmax_g)==0

           index_l=index_l+1
           l2 = dble(l)*dble(l+1)
           lsq = dsqrt( l2 )
           if (USE_TIMER) call cpu_time(start_time_to_plm)

! 33333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333
!! DK DK when l is greater than 5 we can get rid of all the "if" tests in this subroutine to make it much faster
           if (l <= 5) then
             call caldvec_for_l_less_than_5_no_store_l(l,theta,phi,plm,dvec,dvecdt,dvecdp,theta_n)
           else
             call caldvec_for_l_more_than_5_no_store_l(l,theta,phi,plm,dvec,dvecdt,dvecdp,theta_n)
           endif

! 33333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333

           if (USE_TIMER) then
              call cpu_time(finish_time_to_plm)
              time_to_plm =  time_to_plm + finish_time_to_plm - start_time_to_plm
              call cpu_time(start_time_to_project)
           endif

           do m=-2,2        ! m-loop start

              g0tmp(:,:) = tabg0(index_l,:,:,m)
              g0dertmp(:,:) = tabg0der(index_l,:,:,m)

              if ( iabs(m) <= iabs(l) ) then
                 ig2 = 0

                 if ( l == 0 ) then ! l-branch for calu (l=0)

! 33333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333
                    !  rearranging the matrix elements
                       do ir_ = r_n,r_n

                          if (l <= lmax_r(ir_)) then ! on continue si on a pas atteint lmax

!! VM VM I copied and pasted the case l>0 and then made some modifications because the
!! VM VM case (l=0) is simpler.
                             Aval = A0sta(ir_)
                             Cval = C0sta(ir_)
                             Fval = F0sta(ir_)
                             Lval = L0sta(ir_)
                             Nval = N0sta(ir_)

                             rval = r_(ir_)
                             inv_rval = 1.d0 / cmplx(rval)

                    do imt = 1,6 !! DK DK moved this loop here to try manual unrolling, but finally useless for this loop

                      g0tmp1 = tabg0(index_l,1,imt,m)
                      g0dertmp1 = tabg0der(index_l,1,imt,m)

!! DK DK added a compiler directive for the xlf compiler on IBM Blue Gene
!! from http://pic.dhe.ibm.com/infocenter/compbg/v121v141/index.jsp?topic=%2Fcom.ibm.xlf141.bg.doc%2Flanguage_ref%2Fassert.html
!IBM* ASSERT (MINITERCNT(1000))
!DIR$ loop count min(1000)
                             do itheta = 1, theta_n !! VM VM for l=0 the equations are simpler, only the first component is non zero

                                inv_thetasin = 1.d0 / cmplx(sin(theta(itheta) * convert_to_radians))

                                ! call calup0
                                u1 = g0tmp1*dvec(itheta,1,m)
                                udr1 = g0dertmp1*dvec(itheta,1,m)
                                udt1 = g0tmp1*dvecdt(itheta,1,m)
                                udp1 = g0tmp1*dvecdp(itheta,1,m)

                                uder11 = udr1
                                uder12 = udt1 * inv_rval
                                uder13 = udp1*inv_thetasin * inv_rval
                                uder22 = u1 * inv_rval
                                uder33 = u1 * inv_rval

                               ! call locallyCartesianDerivatives
                                stress(itheta,1,imt) = stress(itheta,1,imt)+Cval*uder11+Fval*uder22+Fval*uder33
                                stress(itheta,2,imt) = stress(itheta,2,imt)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder33
                                stress(itheta,3,imt) = stress(itheta,3,imt)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder22

                                stress(itheta,4,imt) = stress(itheta,4,imt)+Lval*uder12
                                stress(itheta,5,imt) = stress(itheta,5,imt)+Lval*uder13

                                ! call udertoStress
                                displacement(itheta,1,imt) = u1 + displacement(itheta,1,imt)

                             enddo ! of loop on itheta
                            enddo ! of loop on imt (mt-loop)
                          endif ! of if (l <= lmax_r(ir_)) then ... i.e. on the stack point
                       enddo ! of loop on ir_
! 33333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333

                 else ! for l /= 0

! 33333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333

!! DK DK added this to precompute the inverse
                    inv_lsq = 1.d0 / dcmplx(lsq)

                       do ir_ = r_n,r_n ! stack point

                          if (l <= lmax_r(ir_)) then ! on continue si on n'a pas atteint lmax

!! DK DK added this
                             Aval = A0sta(ir_)
                             Cval = C0sta(ir_)
                             Fval = F0sta(ir_)
                             Lval = L0sta(ir_)
                             Nval = N0sta(ir_)

!! DK DK added this
                             rval = r_(ir_)
!! DK DK it is always better to precompute an inverse and store it if we use it several times in the loop,
!! DK DK because on processors divisions are much more expensive than multiplications and thus it is much better
!! DK DK to compute the inverse once and for all and then later multiply by the inverse that has been stored
                             inv_rval = 1.d0 / cmplx(rval)

!! DK DK I put itheta as FIRST index in all arrays in order to be able to fully vectorize the loops
!! DK DK and optimize access to the processor cache; thus Vadim please make these permutations everywhere in the code
!! DK DK (for dvec, dvecdt, dvecdp, stress, and displacement)

!! DK DK to Vadim: finalement pas besoin de remonter cette boucle au dessus des deux autres car ca compliquerait le code
!! DK DK to Vadim: (surtout quand ir_ variera dans le cas des faces verticales) et ca n'apporterait quasiment rien
!! DK DK to Vadim: car la boucle actuelle ci-dessous est maintenant totalement vectorisee (suite a mes modifications)

!                   do imt = 1,6 !! DK DK moved this loop here to prepare for manual unrolling

!! DK DK added a compiler directive for the xlf compiler on IBM Blue Gene
!! from http://pic.dhe.ibm.com/infocenter/compbg/v121v141/index.jsp?topic=%2Fcom.ibm.xlf141.bg.doc%2Flanguage_ref%2Fassert.html
!IBM* ASSERT (MINITERCNT(1000))
!DIR$ loop count min(1000)
                             do itheta = 1, theta_n

                               thetaval = theta(itheta) * convert_to_radians
                               thetasin = sin(thetaval)
                               thetacot_complex = cmplx(cos(thetaval)/thetasin)
                               inv_thetasin = 1.d0 / cmplx(thetasin)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! case imt = 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!                               call calup
                                u1 = g0tmp(1,1)*dvec(itheta,1,m)
                                u2 = g0tmp(2,1)*dvec(itheta,2,m)*inv_lsq
                                u3 = g0tmp(2,1)*dvec(itheta,3,m)*inv_lsq

                                udr1 = g0dertmp(1,1)*dvec(itheta,1,m)
                                udr2 = g0dertmp(2,1)*dvec(itheta,2,m)*inv_lsq
                                udr3 = g0dertmp(2,1)*dvec(itheta,3,m)*inv_lsq

                                udt1 = g0tmp(1,1)*dvecdt(itheta,1,m)
                                udt2 = g0tmp(2,1)*dvecdt(itheta,2,m)*inv_lsq
                                udt3 = g0tmp(2,1)*dvecdt(itheta,3,m)*inv_lsq

                                udp1 = g0tmp(1,1)*dvecdp(itheta,1,m)
                                udp2 = g0tmp(2,1)*dvecdp(itheta,2,m)*inv_lsq
                                udp3 = g0tmp(2,1)*dvecdp(itheta,3,m)*inv_lsq

! call locallyCartesianDerivatives
! 1,2,3: r,theta,phi; , denotes the partial derivatives
  uder11 = udr1
  uder12 = (udt1-u2) * inv_rval
  uder13 = (udp1*inv_thetasin-u3) * inv_rval

  uder21 = udr2
  uder22 = (udt2+u1) * inv_rval
  uder23 = (udp2*inv_thetasin-u3*thetacot_complex) * inv_rval

  uder31 = udr3
  uder32 = udt3 * inv_rval
  uder33 = (udp3*inv_thetasin+u1+u2*thetacot_complex) * inv_rval

! call udertoStress
  stress(itheta,1,1) = stress(itheta,1,1)+Cval*uder11+Fval*uder22+Fval*uder33
  stress(itheta,2,1) = stress(itheta,2,1)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder33
  stress(itheta,3,1) = stress(itheta,3,1)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder22

  stress(itheta,4,1) = stress(itheta,4,1)+Lval*(uder12+uder21)
  stress(itheta,5,1) = stress(itheta,5,1)+Lval*(uder13+uder31)
  stress(itheta,6,1) = stress(itheta,6,1)+Nval*(uder23+uder32)

                                displacement(itheta,1,1) = u1 + displacement(itheta,1,1)
                                displacement(itheta,2,1) = u2 + displacement(itheta,2,1)
                                displacement(itheta,3,1) = u3 + displacement(itheta,3,1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! case imt = 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!                               call calup
                                u1 = g0tmp(1,2)*dvec(itheta,1,m)
                                u2 = g0tmp(2,2)*dvec(itheta,2,m)*inv_lsq
                                u3 = g0tmp(2,2)*dvec(itheta,3,m)*inv_lsq

                                udr1 = g0dertmp(1,2)*dvec(itheta,1,m)
                                udr2 = g0dertmp(2,2)*dvec(itheta,2,m)*inv_lsq
                                udr3 = g0dertmp(2,2)*dvec(itheta,3,m)*inv_lsq

                                udt1 = g0tmp(1,2)*dvecdt(itheta,1,m)
                                udt2 = g0tmp(2,2)*dvecdt(itheta,2,m)*inv_lsq
                                udt3 = g0tmp(2,2)*dvecdt(itheta,3,m)*inv_lsq

                                udp1 = g0tmp(1,2)*dvecdp(itheta,1,m)
                                udp2 = g0tmp(2,2)*dvecdp(itheta,2,m)*inv_lsq
                                udp3 = g0tmp(2,2)*dvecdp(itheta,3,m)*inv_lsq

! call locallyCartesianDerivatives
! 1,2,3: r,theta,phi; , denotes the partial derivatives
  uder11 = udr1
  uder12 = (udt1-u2) * inv_rval
  uder13 = (udp1*inv_thetasin-u3) * inv_rval

  uder21 = udr2
  uder22 = (udt2+u1) * inv_rval
  uder23 = (udp2*inv_thetasin-u3*thetacot_complex) * inv_rval

  uder31 = udr3
  uder32 = udt3 * inv_rval
  uder33 = (udp3*inv_thetasin+u1+u2*thetacot_complex) * inv_rval

! call udertoStress
  stress(itheta,1,2) = stress(itheta,1,2)+Cval*uder11+Fval*uder22+Fval*uder33
  stress(itheta,2,2) = stress(itheta,2,2)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder33
  stress(itheta,3,2) = stress(itheta,3,2)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder22

  stress(itheta,4,2) = stress(itheta,4,2)+Lval*(uder12+uder21)
  stress(itheta,5,2) = stress(itheta,5,2)+Lval*(uder13+uder31)
  stress(itheta,6,2) = stress(itheta,6,2)+Nval*(uder23+uder32)

                                displacement(itheta,1,2) = u1 + displacement(itheta,1,2)
                                displacement(itheta,2,2) = u2 + displacement(itheta,2,2)
                                displacement(itheta,3,2) = u3 + displacement(itheta,3,2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! case imt = 3 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!                               call calup
                                u1 = g0tmp(1,3)*dvec(itheta,1,m)
                                u2 = g0tmp(2,3)*dvec(itheta,2,m)*inv_lsq
                                u3 = g0tmp(2,3)*dvec(itheta,3,m)*inv_lsq

                                udr1 = g0dertmp(1,3)*dvec(itheta,1,m)
                                udr2 = g0dertmp(2,3)*dvec(itheta,2,m)*inv_lsq
                                udr3 = g0dertmp(2,3)*dvec(itheta,3,m)*inv_lsq

                                udt1 = g0tmp(1,3)*dvecdt(itheta,1,m)
                                udt2 = g0tmp(2,3)*dvecdt(itheta,2,m)*inv_lsq
                                udt3 = g0tmp(2,3)*dvecdt(itheta,3,m)*inv_lsq

                                udp1 = g0tmp(1,3)*dvecdp(itheta,1,m)
                                udp2 = g0tmp(2,3)*dvecdp(itheta,2,m)*inv_lsq
                                udp3 = g0tmp(2,3)*dvecdp(itheta,3,m)*inv_lsq

! call locallyCartesianDerivatives
! 1,2,3: r,theta,phi; , denotes the partial derivatives
  uder11 = udr1
  uder12 = (udt1-u2) * inv_rval
  uder13 = (udp1*inv_thetasin-u3) * inv_rval

  uder21 = udr2
  uder22 = (udt2+u1) * inv_rval
  uder23 = (udp2*inv_thetasin-u3*thetacot_complex) * inv_rval

  uder31 = udr3
  uder32 = udt3 * inv_rval
  uder33 = (udp3*inv_thetasin+u1+u2*thetacot_complex) * inv_rval

! call udertoStress
  stress(itheta,1,3) = stress(itheta,1,3)+Cval*uder11+Fval*uder22+Fval*uder33
  stress(itheta,2,3) = stress(itheta,2,3)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder33
  stress(itheta,3,3) = stress(itheta,3,3)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder22

  stress(itheta,4,3) = stress(itheta,4,3)+Lval*(uder12+uder21)
  stress(itheta,5,3) = stress(itheta,5,3)+Lval*(uder13+uder31)
  stress(itheta,6,3) = stress(itheta,6,3)+Nval*(uder23+uder32)

                                displacement(itheta,1,3) = u1 + displacement(itheta,1,3)
                                displacement(itheta,2,3) = u2 + displacement(itheta,2,3)
                                displacement(itheta,3,3) = u3 + displacement(itheta,3,3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! case imt = 4 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!                               call calup
                                u1 = g0tmp(1,4)*dvec(itheta,1,m)
                                u2 = g0tmp(2,4)*dvec(itheta,2,m)*inv_lsq
                                u3 = g0tmp(2,4)*dvec(itheta,3,m)*inv_lsq

                                udr1 = g0dertmp(1,4)*dvec(itheta,1,m)
                                udr2 = g0dertmp(2,4)*dvec(itheta,2,m)*inv_lsq
                                udr3 = g0dertmp(2,4)*dvec(itheta,3,m)*inv_lsq

                                udt1 = g0tmp(1,4)*dvecdt(itheta,1,m)
                                udt2 = g0tmp(2,4)*dvecdt(itheta,2,m)*inv_lsq
                                udt3 = g0tmp(2,4)*dvecdt(itheta,3,m)*inv_lsq

                                udp1 = g0tmp(1,4)*dvecdp(itheta,1,m)
                                udp2 = g0tmp(2,4)*dvecdp(itheta,2,m)*inv_lsq
                                udp3 = g0tmp(2,4)*dvecdp(itheta,3,m)*inv_lsq

! call locallyCartesianDerivatives
! 1,2,3: r,theta,phi; , denotes the partial derivatives
  uder11 = udr1
  uder12 = (udt1-u2) * inv_rval
  uder13 = (udp1*inv_thetasin-u3) * inv_rval

  uder21 = udr2
  uder22 = (udt2+u1) * inv_rval
  uder23 = (udp2*inv_thetasin-u3*thetacot_complex) * inv_rval

  uder31 = udr3
  uder32 = udt3 * inv_rval
  uder33 = (udp3*inv_thetasin+u1+u2*thetacot_complex) * inv_rval

! call udertoStress
  stress(itheta,1,4) = stress(itheta,1,4)+Cval*uder11+Fval*uder22+Fval*uder33
  stress(itheta,2,4) = stress(itheta,2,4)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder33
  stress(itheta,3,4) = stress(itheta,3,4)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder22

  stress(itheta,4,4) = stress(itheta,4,4)+Lval*(uder12+uder21)
  stress(itheta,5,4) = stress(itheta,5,4)+Lval*(uder13+uder31)
  stress(itheta,6,4) = stress(itheta,6,4)+Nval*(uder23+uder32)

                                displacement(itheta,1,4) = u1 + displacement(itheta,1,4)
                                displacement(itheta,2,4) = u2 + displacement(itheta,2,4)
                                displacement(itheta,3,4) = u3 + displacement(itheta,3,4)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! case imt = 5 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!                               call calup
                                u1 = g0tmp(1,5)*dvec(itheta,1,m)
                                u2 = g0tmp(2,5)*dvec(itheta,2,m)*inv_lsq
                                u3 = g0tmp(2,5)*dvec(itheta,3,m)*inv_lsq

                                udr1 = g0dertmp(1,5)*dvec(itheta,1,m)
                                udr2 = g0dertmp(2,5)*dvec(itheta,2,m)*inv_lsq
                                udr3 = g0dertmp(2,5)*dvec(itheta,3,m)*inv_lsq

                                udt1 = g0tmp(1,5)*dvecdt(itheta,1,m)
                                udt2 = g0tmp(2,5)*dvecdt(itheta,2,m)*inv_lsq
                                udt3 = g0tmp(2,5)*dvecdt(itheta,3,m)*inv_lsq

                                udp1 = g0tmp(1,5)*dvecdp(itheta,1,m)
                                udp2 = g0tmp(2,5)*dvecdp(itheta,2,m)*inv_lsq
                                udp3 = g0tmp(2,5)*dvecdp(itheta,3,m)*inv_lsq

! call locallyCartesianDerivatives
! 1,2,3: r,theta,phi; , denotes the partial derivatives
  uder11 = udr1
  uder12 = (udt1-u2) * inv_rval
  uder13 = (udp1*inv_thetasin-u3) * inv_rval

  uder21 = udr2
  uder22 = (udt2+u1) * inv_rval
  uder23 = (udp2*inv_thetasin-u3*thetacot_complex) * inv_rval

  uder31 = udr3
  uder32 = udt3 * inv_rval
  uder33 = (udp3*inv_thetasin+u1+u2*thetacot_complex) * inv_rval

! call udertoStress
  stress(itheta,1,5) = stress(itheta,1,5)+Cval*uder11+Fval*uder22+Fval*uder33
  stress(itheta,2,5) = stress(itheta,2,5)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder33
  stress(itheta,3,5) = stress(itheta,3,5)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder22

  stress(itheta,4,5) = stress(itheta,4,5)+Lval*(uder12+uder21)
  stress(itheta,5,5) = stress(itheta,5,5)+Lval*(uder13+uder31)
  stress(itheta,6,5) = stress(itheta,6,5)+Nval*(uder23+uder32)

                                displacement(itheta,1,5) = u1 + displacement(itheta,1,5)
                                displacement(itheta,2,5) = u2 + displacement(itheta,2,5)
                                displacement(itheta,3,5) = u3 + displacement(itheta,3,5)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! case imt = 6 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!                               call calup
                                u1 = g0tmp(1,6)*dvec(itheta,1,m)
                                u2 = g0tmp(2,6)*dvec(itheta,2,m)*inv_lsq
                                u3 = g0tmp(2,6)*dvec(itheta,3,m)*inv_lsq

                                udr1 = g0dertmp(1,6)*dvec(itheta,1,m)
                                udr2 = g0dertmp(2,6)*dvec(itheta,2,m)*inv_lsq
                                udr3 = g0dertmp(2,6)*dvec(itheta,3,m)*inv_lsq

                                udt1 = g0tmp(1,6)*dvecdt(itheta,1,m)
                                udt2 = g0tmp(2,6)*dvecdt(itheta,2,m)*inv_lsq
                                udt3 = g0tmp(2,6)*dvecdt(itheta,3,m)*inv_lsq

                                udp1 = g0tmp(1,6)*dvecdp(itheta,1,m)
                                udp2 = g0tmp(2,6)*dvecdp(itheta,2,m)*inv_lsq
                                udp3 = g0tmp(2,6)*dvecdp(itheta,3,m)*inv_lsq

! call locallyCartesianDerivatives
! 1,2,3: r,theta,phi; , denotes the partial derivatives
  uder11 = udr1
  uder12 = (udt1-u2) * inv_rval
  uder13 = (udp1*inv_thetasin-u3) * inv_rval

  uder21 = udr2
  uder22 = (udt2+u1) * inv_rval
  uder23 = (udp2*inv_thetasin-u3*thetacot_complex) * inv_rval

  uder31 = udr3
  uder32 = udt3 * inv_rval
  uder33 = (udp3*inv_thetasin+u1+u2*thetacot_complex) * inv_rval

! call udertoStress
  stress(itheta,1,6) = stress(itheta,1,6)+Cval*uder11+Fval*uder22+Fval*uder33
  stress(itheta,2,6) = stress(itheta,2,6)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder33
  stress(itheta,3,6) = stress(itheta,3,6)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder22

  stress(itheta,4,6) = stress(itheta,4,6)+Lval*(uder12+uder21)
  stress(itheta,5,6) = stress(itheta,5,6)+Lval*(uder13+uder31)
  stress(itheta,6,6) = stress(itheta,6,6)+Nval*(uder23+uder32)

                                displacement(itheta,1,6) = u1 + displacement(itheta,1,6)
                                displacement(itheta,2,6) = u2 + displacement(itheta,2,6)
                                displacement(itheta,3,6) = u3 + displacement(itheta,3,6)

                             enddo ! of loop on itheta
!                           enddo ! of loop on imt (mt-loop)
                          endif ! of if (l <= lmax_r(ir_)) then ... i.e. on the stack point
                       enddo ! of loop on ir_
! 33333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333

                 endif   ! l-branch for calu
              endif

             if (USE_TIMER) then
                call cpu_time(finish_time_to_project)
                time_to_project =  time_to_project + finish_time_to_project - start_time_to_project
             endif

           enddo            ! of loop on m
        enddo               ! of loop on l

        !********************
        if (myrank == 0 .and. SLOW_DEBUG_MODE) write(24,*) i, dble(i)/tlen, lmax_lu
        if (myrank == 0 .and. SLOW_DEBUG_MODE) write(24,*) 'TIME :',mybigrank,time_to_read, time_to_plm,time_to_project

! Sur le cas des data de Pyrope (sur CURIE), j'ai un nsta = ~4000 pour chaque
! proc. J'ai un nsta_global qui vaut ~200000 que je repartis en 48. Je fais
! par ailleurs aussi une parallelisation sur les frequences. Je me retrouve
! avec 10 sous communicateurs de 48 procs (480 procs au total sur CURIE).
! Chaque sous communicateur calcule sur une plage de frequence specifique et
! a l'interieur de chaque sous communicateur les procs se partagent le
! nsta_global. J'ai fait ce type de decoupage car c'etait le seul moyen d'arriver
! a faire passer le code avec le code non optimise.

! Maintenant ca doit valoir le coup de reflechir a optimiser les boucles
! de copies car on ne va pas toujours tourner sur autant de procs, je dirais
! que nsta va aller de 1000 (pour ~2000 procs) a 20000 (~100 procs). Si dans
! le futur quelqu'un fait tourner le code sur un centre regional, il va
! plutot utiliser 100 procs. Dans l'immediat pour nos deux prochains papiers, on
! va faire tous les calculs sur les centres nationaux, donc ce n'est pas crucial
! de faire l'optimization maintenant, mais il faudra y penser a terme.

!! DK DK not a good idea to put indices when working on the whole array
!! DK DK because some compilers will then perform loops instead of memory copies
! it is sufficient for us to save the final result to disk in single precision.
! Also precompute the inverse to speed the code up
  inverse_of_omega_complex = 1.d0 / dcmplx(0,omega)

!! DK DK here the order of the indices is purposely not the same between stress and stresssngl,
!! DK DK nor between displacement and displacementsngl: index 1..nsta is first in one case and last in the other.
!! DK DK we need it to remain first for stress() and displacement() in order to be able to efficiently
!! DK DK vectorize the big loops on that index above, and we need it to remain last in stresssngl()
!! DK DK and displacementsngl() for the mpi_gatherv() below to work fine.
!! DK DK Thus the loop below will generate cache misses and will not use memory prefetching in an optimal fashion,
!! DK DK but we have to do it.

!IBM* ASSERT (MINITERCNT(1000))
!DIR$ SIMD
!DIR$ loop count min(1000)
  do ista = 1,nsta

!   stresssngl(1:6,1:6,ista) = stress(ista,1:6,1:6) * inverse_of_omega_complex
!   displacementsngl(1:3,1:6,ista) = displacement(ista,1:3,1:6)

    stresssngl(1,1,ista) = stress(ista,1,1) * inverse_of_omega_complex
    stresssngl(2,1,ista) = stress(ista,2,1) * inverse_of_omega_complex
    stresssngl(3,1,ista) = stress(ista,3,1) * inverse_of_omega_complex
    stresssngl(4,1,ista) = stress(ista,4,1) * inverse_of_omega_complex
    stresssngl(5,1,ista) = stress(ista,5,1) * inverse_of_omega_complex
    stresssngl(6,1,ista) = stress(ista,6,1) * inverse_of_omega_complex
    stresssngl(1,2,ista) = stress(ista,1,2) * inverse_of_omega_complex
    stresssngl(2,2,ista) = stress(ista,2,2) * inverse_of_omega_complex
    stresssngl(3,2,ista) = stress(ista,3,2) * inverse_of_omega_complex
    stresssngl(4,2,ista) = stress(ista,4,2) * inverse_of_omega_complex
    stresssngl(5,2,ista) = stress(ista,5,2) * inverse_of_omega_complex
    stresssngl(6,2,ista) = stress(ista,6,2) * inverse_of_omega_complex
    stresssngl(1,3,ista) = stress(ista,1,3) * inverse_of_omega_complex
    stresssngl(2,3,ista) = stress(ista,2,3) * inverse_of_omega_complex
    stresssngl(3,3,ista) = stress(ista,3,3) * inverse_of_omega_complex
    stresssngl(4,3,ista) = stress(ista,4,3) * inverse_of_omega_complex
    stresssngl(5,3,ista) = stress(ista,5,3) * inverse_of_omega_complex
    stresssngl(6,3,ista) = stress(ista,6,3) * inverse_of_omega_complex
    stresssngl(1,4,ista) = stress(ista,1,4) * inverse_of_omega_complex
    stresssngl(2,4,ista) = stress(ista,2,4) * inverse_of_omega_complex
    stresssngl(3,4,ista) = stress(ista,3,4) * inverse_of_omega_complex
    stresssngl(4,4,ista) = stress(ista,4,4) * inverse_of_omega_complex
    stresssngl(5,4,ista) = stress(ista,5,4) * inverse_of_omega_complex
    stresssngl(6,4,ista) = stress(ista,6,4) * inverse_of_omega_complex
    stresssngl(1,5,ista) = stress(ista,1,5) * inverse_of_omega_complex
    stresssngl(2,5,ista) = stress(ista,2,5) * inverse_of_omega_complex
    stresssngl(3,5,ista) = stress(ista,3,5) * inverse_of_omega_complex
    stresssngl(4,5,ista) = stress(ista,4,5) * inverse_of_omega_complex
    stresssngl(5,5,ista) = stress(ista,5,5) * inverse_of_omega_complex
    stresssngl(6,5,ista) = stress(ista,6,5) * inverse_of_omega_complex
    stresssngl(1,6,ista) = stress(ista,1,6) * inverse_of_omega_complex
    stresssngl(2,6,ista) = stress(ista,2,6) * inverse_of_omega_complex
    stresssngl(3,6,ista) = stress(ista,3,6) * inverse_of_omega_complex
    stresssngl(4,6,ista) = stress(ista,4,6) * inverse_of_omega_complex
    stresssngl(5,6,ista) = stress(ista,5,6) * inverse_of_omega_complex
    stresssngl(6,6,ista) = stress(ista,6,6) * inverse_of_omega_complex
  enddo

!IBM* ASSERT (MINITERCNT(1000))
!DIR$ SIMD
!DIR$ loop count min(1000)
  do ista = 1,nsta
    displacementsngl(1,1,ista) = displacement(ista,1,1)
    displacementsngl(2,1,ista) = displacement(ista,2,1)
    displacementsngl(3,1,ista) = displacement(ista,3,1)
    displacementsngl(1,2,ista) = displacement(ista,1,2)
    displacementsngl(2,2,ista) = displacement(ista,2,2)
    displacementsngl(3,2,ista) = displacement(ista,3,2)
    displacementsngl(1,3,ista) = displacement(ista,1,3)
    displacementsngl(2,3,ista) = displacement(ista,2,3)
    displacementsngl(3,3,ista) = displacement(ista,3,3)
    displacementsngl(1,4,ista) = displacement(ista,1,4)
    displacementsngl(2,4,ista) = displacement(ista,2,4)
    displacementsngl(3,4,ista) = displacement(ista,3,4)
    displacementsngl(1,5,ista) = displacement(ista,1,5)
    displacementsngl(2,5,ista) = displacement(ista,2,5)
    displacementsngl(3,5,ista) = displacement(ista,3,5)
    displacementsngl(1,6,ista) = displacement(ista,1,6)
    displacementsngl(2,6,ista) = displacement(ista,2,6)
    displacementsngl(3,6,ista) = displacement(ista,3,6)
  enddo

!       il faut faire un allgatherv pour recuperer stresssngl et displacementsngl
        Coeff = 3*6*nsta
        call mpi_gatherv(displacementsngl,Coeff,MPI_COMPLEX,displacementsngl_global,NbInd,Ind,MPI_COMPLEX,0,SubCommunicators,ierr)
        Coeff = 6*6*nsta
        call mpi_gatherv(stresssngl,Coeff,MPI_COMPLEX,stresssngl_global,NbIndSt,IndSt,MPI_COMPLEX,0,SubCommunicators,ierr)

     endif
     if (SLOW_DEBUG_MODE) then
        !write(36,*) 'i ', i
       ! write(36,*)  '1 :',displacement(1,:,:)
        !write(36,*)  'nsta :',displacement(nsta,:,:)
     endif
     !*!******** ECRITURE SUR LE DISQUE
     if (myrank == 0) then

     write(coutfile, '(I5.5,".Stress_PSV_",I5.5)') int(r0(ir0)*10.d0),ifq-1
     coutfile = trim(modelname)//"."//coutfile
     coutfile = trim(outputDir)//"/Stress/"//coutfile

!! DK DK not a good idea to give indices when writing the whole array: some compilers could create copies or use loops in such a case
!! DK DK unformatted sequential I/Os (because sequential I/O is usually faster than direct-access I/O)
!! DK DK and writing a single big file is also much better
        open(44,file=coutfile,form='unformatted',action='write')
        write(44) stresssngl_global
        write(44) displacementsngl_global
        close(44)
     endif

     if (myrank == 0 .and. SLOW_DEBUG_MODE) then
        call date_and_time(datex,timex)
        write(25,'(/a,i5,a,a4,a,a2,a,a2,a,a2,a,a2,a,a4)') '    Frequency-index ', i, ' :', &
             datex(1:4),'-',datex(5:6),'-',datex(7:8),'.  ',timex(1:2),':',timex(3:4),':',timex(5:8)
     endif
    if (myrank == 0) close(34)

  enddo ! of frequency loop

  if (myrank == 0 .and. SLOW_DEBUG_MODE) then
   close(24)
   call date_and_time(datex,timex)
   write(25,'(/a,a4,a,a2,a,a2,a,a2,a,a2,a,a4)') '    Finishing date and time:                     ', &
       datex(1:4),'-',datex(5:6),'-',datex(7:8),'.  ',timex(1:2),':',timex(3:4),':',timex(5:8)
   close(25)
  endif
  if (USE_TIMER) then
    call cpu_time(finish_time_total)
     time_total = finish_time_total - start_time_total
     ! convert to hours
     time_total = time_total / 3600.
     time_to_read = time_to_read / 3600.
     time_to_plm = time_to_plm / 3600.
     time_to_project = time_to_project / 3600.
     call mpi_allreduce( time_total, total_global_time,1, MPI_REAL ,MPI_SUM, MPI_COMM_WORLD, ierr)
     time_total= total_global_time
     call mpi_allreduce( time_to_read, total_global_time,1, MPI_REAL ,MPI_SUM, MPI_COMM_WORLD, ierr)
     time_to_read= total_global_time
     call mpi_allreduce( time_to_plm, total_global_time,1, MPI_REAL ,MPI_SUM, MPI_COMM_WORLD, ierr)
     time_to_plm = total_global_time
     call mpi_allreduce( time_to_project, total_global_time,1, MPI_REAL ,MPI_SUM, MPI_COMM_WORLD, ierr)
     time_to_project = total_global_time
      if (myrank == 0) then
        open(25,file='timer_part3_zmin.txt')
        write(25,*) 'Total :',time_total
        write(25,*) 'plm calculation :',time_to_plm,100*time_to_plm/time_total,' %'
        write(25,*) 'reading :', time_to_read,100* time_to_read/time_total
        write(25,*) 'projection  :', time_to_project,100*time_to_project/time_total
        close(25)
     endif

  endif
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call MPI_COMM_FREE(SubCommunicators,ierr)
  call MPI_FINALIZE(ierr)

end program TraPSV

