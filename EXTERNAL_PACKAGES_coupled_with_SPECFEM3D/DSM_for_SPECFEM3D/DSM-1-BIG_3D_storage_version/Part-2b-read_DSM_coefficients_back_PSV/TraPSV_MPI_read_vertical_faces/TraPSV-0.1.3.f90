
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
  !integer, parameter :: maxlmax = 35000, maxlmax_g=1000  !! VM VM I moved this in constants.h

  real(kind(0d0)) :: tlen
  real(kind(0d0)) :: r0min, r0max, r0delta  !!! JUST FOR ONE DEPTH FOR THIS MOMENT !!
  real(kind(0d0)) :: r0lat, r0lon, stla_curr,stlo_curr
  real(kind(0d0)), allocatable :: stla(:),stlo(:),stla_g(:),stlo_g(:)
  real(kind(0d0)), allocatable :: r_(:),r0(:),theta(:),phi(:),phi_g(:),theta_g(:)
  real(kind(0d0)), allocatable :: A0sta(:),C0sta(:),F0sta(:),L0sta(:),N0sta(:)
  integer, allocatable :: updown(:),idum(:)
  real(kind(0d0)), allocatable :: rrsta(:,:)
  integer, allocatable :: iista(:,:)
  integer :: r_n,r0_n,ir_,ir0,imt,itheta,theta_n,nsta,ier

  character(120) :: coutfile
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
  complex(kind(0d0)),allocatable:: tabg0(:,:,:,:,:),tabg0der(:,:,:,:,:)
  complex(kind(0d0)), allocatable :: g0tmp(:,:,:),g0dertmp(:,:,:) ! forward
  complex(kind(0d0)) :: g0tmp1, g0dertmp1
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
  integer:: i,ig2
  real(kind(0d0)):: l2,lsq,bazi

  !-----------------------------------------------------------------------
  complex(kind(0d0)), allocatable :: dvec(:,:,:,:),dvecdt(:,:,:,:),dvecdp(:,:,:,:)
  complex(kind(0d0)), allocatable :: stress(:,:,:,:), displacement(:,:,:,:)
  complex(kind(0e0)), allocatable :: stresssngl(:,:,:,:), displacementsngl(:,:,:,:)
  complex(kind(0e0)), allocatable :: stresssngl_global(:,:,:,:), displacementsngl_global(:,:,:,:)
  real(kind(0d0)), allocatable :: plm(:,:,:)
  integer ifrequ_min, ifrequ_max

!! DK DK added this to precompute several things
  complex(kind(0d0)) :: inv_lsq
  real(kind(0d0)) :: Aval,Cval,Fval,Lval,Nval
  complex(kind(0d0)) :: u1,u2,u3,udr1,udr2,udr3,udt1,udt2,udt3,udp1,udp2,udp3
  complex(kind(0d0)) :: uder11 , uder12 , uder13 , uder21 , uder22 , uder23 , uder31 , uder32 , uder33
  real(kind(0d0)) :: rval,thetaval
  real(kind(0d0)) :: thetasin
  complex(kind(0d0)) :: thetacot_complex, inv_rval, inv_thetasin, inverse_of_omega_complex

  ! ajout de variables pour MPI --- VM VM
  integer i_color,nb_colors,current_color,imin_current,imax_current,nbproc_current
  integer iread,nb_tot_proc
  integer, allocatable :: Ifrq2(:,:),Ind(:),NbInd(:),IndSt(:),NbIndSt(:)
  integer istamin,istamax,k,ifq,c8,coeff
! integer irec
  integer imin_glob, imax_glob, irank,nsta_global,ista

!! DK DK added this conversion factor
  double precision, parameter :: convert_to_radians = pi/180.d0

  !--------------------------------------------------------------------------

  ! MPI --- VM VM
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mybigrank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nbbigproc,ierr)

  ! pour definir les sous communicateurs
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
  write(*,*) mybigrank,'1'
  call MPI_Bcast(key,nbbigproc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(color,nbbigproc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(Ifrq2,2*nbbigproc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  if (mybigrank == 0) then
    call pinputTra_part2(ifrequ_min,ifrequ_max,inputdir,outputDir,psvmodel,modelname,stationsinf,tlen, &
                    imin_glob,imax_glob,r0min,r0max,r0delta,r0lat,r0lon,itranslat,mybigrank)
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
   read(1,*) nsta_global
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

  allocate(stresssngl_global(1:6,1:6,1:r_n,1:nsta_global))
  allocate(displacementsngl_global(1:3,1:6,1:r_n,1:nsta_global))
  allocate(tabg0(1:maxlmax_g,2,1:6,r_n,-2:2))
  allocate(tabg0der(1:maxlmax_g,2,6,r_n,-2:2))
  allocate(g0tmp(2,6,r_n))
  allocate(g0dertmp(2,6,r_n))

  if (mybigrank == 0) then
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
  ! definition de l'ensemble des sous communicateurs
  call MPI_COMM_SPLIT(MPI_COMM_WORLD,color(mybigrank),key(mybigrank),subcommunicators,ierr)
  call MPI_COMM_RANK(SubCommunicators,myrank,ierr)
  call MPI_COMM_SIZE(SubCommunicators,nbproc,ierr)

! distribuer les stations sur les procs
  call DistribSta(nsta,nsta_global,myrank,nbproc,istamin,istamax)
  write(*,*) 'myrank  ista :', myrank, istamin,istamax

  call MPI_Barrier(SubCommunicators,ierr)

 ! parallelisation selon les distances epicentrales
  allocate(Ind(nbproc))
  allocate(NbInd(nbproc))
  allocate(IndSt(nbproc))
  allocate(NbIndSt(nbproc))
  call mpi_allgather(nsta,1,MPI_INTEGER,NbInd,1,MPI_INTEGER,SubCommunicators, ierr )

  NbIndSt(:) = NbInd(:)

  ! displ
  coeff  = 3*6*r_n
  Ind(1) = 0
  do k = 2, nbproc
      c8 = NbInd(k-1)
      Ind(k) = Coeff*c8 +   Ind(k-1)
  enddo

  do k=1,nbproc
     NbInd(k) = NbInd(k) * Coeff
  enddo

  ! stress
  coeff  = 6*6*r_n
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
  allocate(dvec(1:theta_n,1:3,-2:2,0:maxlmax))
  allocate(dvecdt(1:theta_n,1:3,-2:2,0:maxlmax))
  allocate(dvecdp(1:theta_n,1:3,-2:2,0:maxlmax),stat=ier)
  if (ier /= 0) stop 'error: not enough memory to allocate array'
!! DK DK moved index 1:theta_n to the first position in order to be able to vectorize some loops
!! DK DK and also got rid of the "l" index because we do not store them for all l any more
! allocate(plm(1:3,0:3,1:theta_n,0:maxlmax),stat=ier)
  allocate(plm(1:theta_n,1:3,0:3),stat=ier)
  if (ier /= 0) stop 'error: not enough memory to allocate array'

  phi(:)=phi_g(istamin:istamax)
  theta(:)=theta_g(istamin:istamax)
  stla(:)=stla_g(istamin:istamax)
  stlo(:)=stlo_g(istamin:istamax)

  ! source depths
  ! for this moment we don't put any necessary allocation
  r0_n =  int((r0max-r0min)/r0delta)+1
  allocate(r0(1:r0_n))
  do i = 1, r0_n
     r0(i) = r0min + dble(i-1)*r0delta
  enddo
  ir0 = r0_n

  isl = 0
  ill = 0

! From Vadim Monteiller, April 2013:
! TO DO (DK DK: NOW DONE): on pourrait gagner ici sans stocker les harmoniques spheriques.
! Ce commentaire date du moment ou j'essayais de gagner de la memoire. J'avais vu la possibilite de reduire la taille
! du tableau plm() sans le pre-calculer mais en calculant les harmoniques spheriques a l'interieur de la boucle sur l
! (comme on le fait pour la face zmin). Mais comme jusqu'a present pour les faces verticales j'arrive a stocker ce tableau,
! je l'ai laisse tel quel.
!! DK DK April 2013: I have now removed it and I recompute the plm() array inside, thus this optimization is now implemented.
! call clPLM(plm,maxlmax,theta,theta_n,myrank)
!! VM VM initialisation
  plm=0.d0
  dvec=(0.d0,0.d0)
  dvecdt=(0.d0,0.d0)
  dvecdp=(0.d0,0.d0)
  do l = 0,maxlmax
!    do itheta = 1,theta_n
!       call caldvec_already_plm(l,(theta(itheta)/180.d0*pi), (phi(itheta)/180.d0*pi),plm(1:3,0:3,itheta,l), &
!               dvec(itheta,1:3,-2:2,l),dvecdt(itheta,1:3,-2:2,l),dvecdp(itheta,1:3,-2:2,l))
!    enddo
!! DK DK created new, vectorized versions of these routines
!! DK DK when l is greater than 5 we can get rid of all the "if" tests in this subroutine to make it much faster
    if (l <= 5) then
      call caldvec_for_l_less_than_5_store_l(l,theta,phi,plm,dvec,dvecdt,dvecdp,theta_n,maxlmax)
    else
      call caldvec_for_l_more_than_5_store_l(l,theta,phi,plm,dvec,dvecdt,dvecdp,theta_n,maxlmax)
    endif
  enddo
  deallocate(plm)

  allocate(stress(1:nsta,1:r_n,1:6,1:6))
  allocate(displacement(1:nsta,1:r_n,1:3,1:6))
  allocate(stresssngl(1:6,1:6,1:r_n,1:nsta))
  allocate(displacementsngl(1:3,1:6,1:r_n,1:nsta))

  !******************** plm reading                 *********************

  if (myrank == 0 .and. SLOW_DEBUG_MODE) then
    ! record the date and time at the beginning of the job
    write(list, '(I6.6,".",I6.6)') ifrequ_min, ifrequ_max
     list = trim(outputDir)//"/log/calLog"//"."//trim(modelname)//"."//trim(list)

    open(25,file =list, status = 'unknown', form = 'formatted')
    call date_and_time(datex,timex)
    write(25,'(/a,a4,a,a2,a,a2,a,a2,a,a2,a,a4)') '    Starting date and time:                     ', &
        datex(1:4),'-',datex(5:6),'-',datex(7:8),'.  ',timex(1:2),':',timex(3:4),':',timex(5:8)
  endif

  !******************** Computing the displacement *********************

  llog = 0
  ir_=0

  if (myrank == 0) then
   write(list1, '(I6.6,".",I6.6)') ifrequ_min, ifrequ_max
   list1 = trim(outputDir)//"/log/list"//"."//trim(modelname)//"."//trim(list1)

   if (SLOW_DEBUG_MODE) then
     open(24, file = list1, status = 'unknown', form = 'formatted')
     write(24,*)


     call date_and_time(datex,timex)
     write(25,'(/a,a4,a,a2,a,a2,a,a2,a,a2,a,a4)') '    PLM calculation done:                     ', &
         datex(1:4),'-',datex(5:6),'-',datex(7:8),'.  ',timex(1:2),':',timex(3:4),':',timex(5:8)
   endif

  endif

!***************************************************************************************************
  do ifq= ifrequ_min+1, ifrequ_max+1   ! loop on frequencies, we start at 0
!***************************************************************************************************

     i=ifq-1

     write(coutfile, '(I5.5,".G0_PSV_",I5.5)') int(r0(ir0)*10.d0),ifq-1
     coutfile = trim(modelname)//"."//coutfile
     coutfile = trim(inputdir)//"/Stress/"//coutfile
     if (myrank == 0) then
        write(*,*) myrank,ifq,trim(coutfile)
        open(34,file=coutfile, form='unformatted',action='read')
     endif
!! VM VM the stress and displacement are now in the same file
     !write(coutfile, '(I5.5,".G0der_PSV_",I5.5)') int(r0(ir0)*10.d0),ifq-1
     !coutfile = trim(modelname)//"."//coutfile
     !coutfile = trim(inputdir)//"/Displacement/"//coutfile
     !if (myrank==0) then
     !   write(*,*) myrank,ifq,trim(coutfile)
     !   open(35,file=coutfile, form='unformatted')
     !endif

     write(coutfile, '(I5.5,".Stress_PSV_",I5.5)') int(r0(ir0)*10.d0),ifq-1
     coutfile = trim(modelname)//"."//coutfile
     coutfile = trim(outputDir)//"/Stress/"//coutfile
!    if (myrank==0) open(44,file=coutfile, access = 'direct',recl=2*6*6*kind(0e0)*nsta_global*r_n,action='write')
!! DK DK changed from direct access to sequential unformatted
     if (myrank == 0) open(44,file=coutfile,form='unformatted',action='write')

!! DK DK suppressed that file in order to create a single output file for both stress and displacement
!    write(coutfile, '(I5.5,".Displ_PSV_",I5.5)') int(r0(ir0)*10.d0),ifq-1
!    coutfile = trim(modelname)//"."//coutfile
!    coutfile = trim(outputDir)//"/Displacement/"//coutfile
!    if (myrank==0) open(45,file=coutfile, access = 'direct',recl=2*3*6*kind(0e0)*nsta_global*r_n,action='write')

!    irec=0

     ir0 = r0_n

     ! initialisation
     stress = dcmplx(0.d0)
     displacement = dcmplx(0.d0)
     stresssngl=cmplx(0.e0)
     displacementsngl=cmplx(0.e0)
     stresssngl_global=cmplx(0.e0)
     displacementsngl_global=cmplx(0.e0)
     lmax_r(:)=maxlmax !! VM VM je dois initialiser au max a priori
     omega = 2.d0 * pi * dble(i) / tlen

     if ( i /= 0 ) then

        lref=0
        index_l=0
        lmax_lu=0

!***************************************************************************************************
        do l = 0,maxlmax  ! on parcourt tous les l a priori
!***************************************************************************************************
           if (myrank == 0 .and. SLOW_DEBUG_MODE) write(*,*) 'l= ',l
           call MPI_Barrier(SubCommunicators,ierr)  ! je veux que tous les procs aillent a la meme vitesse

           ! lecture des tableaux
           if (mod(l,maxlmax_g) == 0) then
              lref=l
              tabg0=(0.d0,0.d0)
              tabg0der=(0.d0,0.d0)
              if (myrank == 0) then
                if (SLOW_DEBUG_MODE) then
                  write(24,*)
                  write(24,*) myrank,' reading l = ',l
                endif

                do
                   if (SLOW_DEBUG_MODE) then
                      write(24,*) 'I want to read for l=',l
                      write(24,*) 'max(lmax_r)=',maxval(lmax_r)
                   endif
                   read(34) ir_,llog
                   if (SLOW_DEBUG_MODE) then
                      write(24,*) 'I already read for l=',l,' and ir_= ',ir_,llog
                   endif


                   if (SLOW_DEBUG_MODE) then
                      write(24,*)
                      write(24,*) mybigrank,'is reading for ir,llog',ir_,llog
                   endif

                   if (ir_ == -1) then ! this flag indicates the end of all the data to read for the current block
                      if (SLOW_DEBUG_MODE) then
                        write(24,*) mybigrank,ir_,llog, ' exit '
                      endif

                      exit
                   else
                      read(34) tabg0(1:llog,:,:,ir_,:)
                      read(34) tabg0der(1:llog,:,:,ir_,:)
                      lmax_r(ir_) = llog + lref ! lmax courant pour ir_
                      lmax_lu=max(lmax_lu,lmax_r(ir_))
                      if (SLOW_DEBUG_MODE) then
                        write(24,*) mybigrank,'allready read TABS',ir_,llog,lmax_r(ir_),lmax_lu
                      endif

                   endif

                enddo
                if (SLOW_DEBUG_MODE) write(24,*) 'end reading block'
              endif
              index_l=0
              call MPI_Bcast(lmax_r,r_n,MPI_INTEGER,0,SubCommunicators,ierr)
              call MPI_Bcast(tabg0,60*r_n*maxlmax_g,MPI_DOUBLE_COMPLEX,0,SubCommunicators,ierr)
              call MPI_Bcast(tabg0der,60*r_n*maxlmax_g,MPI_DOUBLE_COMPLEX,0,SubCommunicators,ierr)
              call MPI_Bcast(ir_,1,MPI_INTEGER,0,SubCommunicators,ierr)
              call MPI_Bcast(llog,1,MPI_INTEGER,0,SubCommunicators,ierr)
           endif
           !! VM VM : bug ici quand mod(maxval(lmax_r),maxlmax_g)==0 : on sort a maxval(lmax_r)+1
           !! mais on passe dans la boucle de lecture juste au dessus, ca plante car il y a
           !! plus rien a lire.
           if (l > maxval(lmax_r)) then !!exit          !! VM VM : this mean we reach the big l that wee need
              if (SLOW_DEBUG_MODE) write(24,*) 'I exit because l too big'
              exit
           endif
           if (ir_ == -1 .and. llog == -1) then !!exit  !! VM VM I add this test because of bugs
               if (SLOW_DEBUG_MODE) write(24,*) 'I exit because end of file'
               exit
           endif                                  !! VM VM when mod(maxval(lmax_r),maxlmax_g) == 0
           index_l = index_l+1
           l2 = dble(l)*dble(l+1)
           lsq = dsqrt( l2 )

           do m=-2,2        ! m-loop start

              g0tmp(:,:,:) = tabg0(index_l,:,:,:,m)
              g0dertmp(:,:,:) = tabg0der(index_l,:,:,:,m)

              if ( iabs(m) <= iabs(l) ) then

                 ig2 = 0

!***************************************************************************************************
                 if ( l == 0 ) then ! l-branch for calu (l=0)
!***************************************************************************************************

                    !  rearranging the matrix elements
                    do ir_=1,r_n

                       if (l <= lmax_r(ir_)) then ! on continue si on a pas atteint lmax

                            Aval = A0sta(ir_)
                            Cval = C0sta(ir_)
                            Fval = F0sta(ir_)
                            Lval = L0sta(ir_)
                            Nval = N0sta(ir_)

                            rval = r_(ir_)
                            inv_rval = 1.d0 / cmplx(rval)

                            do imt = 1,6

                               g0tmp1 = tabg0(index_l,1,imt,ir_,m)
                               g0dertmp1 = tabg0der(index_l,1,imt,ir_,m)

!IBM* ASSERT (MINITERCNT(1000))
!DIR$ loop count min(1000)
                               do itheta = 1, theta_n

                                  inv_thetasin = 1.d0 / cmplx(sin(theta(itheta) * convert_to_radians))

                                  ! call calup0
                                  u1 = g0tmp1*dvec(itheta,1,m,l)
                                  udr1 = g0dertmp1*dvec(itheta,1,m,l)
                                  udt1 = g0tmp1*dvecdt(itheta,1,m,l)
                                  udp1 = g0tmp1*dvecdp(itheta,1,m,l)

                                  uder11 = udr1
                                  uder12 = udt1 * inv_rval
                                  uder13 = udp1*inv_thetasin * inv_rval
                                  uder22 = u1 * inv_rval
                                  uder33 = u1 * inv_rval

                  ! call locallyCartesianDerivatives
                  stress(itheta,ir_,1,imt) = stress(itheta,ir_,1,imt)+Cval*uder11+Fval*uder22+Fval*uder33
                  stress(itheta,ir_,2,imt) = stress(itheta,ir_,2,imt)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder33
                  stress(itheta,ir_,3,imt) = stress(itheta,ir_,3,imt)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder22

                  stress(itheta,ir_,4,imt) = stress(itheta,ir_,4,imt)+Lval*uder12
                  stress(itheta,ir_,5,imt) = stress(itheta,ir_,5,imt)+Lval*uder13

                  ! call udertoStress
                  displacement(itheta,ir_,1,imt) = u1 + displacement(itheta,ir_,1,imt)

                               enddo ! of loop on itheta
                            enddo ! of loop on imt (mt-loop)
                         endif ! of if (l <= lmax_r(ir_)) then ... i.e. on the stack point
                      enddo ! of loop on ir_

!***************************************************************************************************
                   else ! for l /= 0
!***************************************************************************************************

                      inv_lsq = 1.d0 / dcmplx(lsq)

                       do ir_=1,r_n ! stack point

                          if (l <= lmax_r(ir_)) then ! on continue si on a pas atteint lmax

                             Aval = A0sta(ir_)
                             Cval = C0sta(ir_)
                             Fval = F0sta(ir_)
                             Lval = L0sta(ir_)
                             Nval = N0sta(ir_)

                             rval = r_(ir_)

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
                                u1 = g0tmp(1,1,ir_)*dvec(itheta,1,m,l)
                                u2 = g0tmp(2,1,ir_)*dvec(itheta,2,m,l)*inv_lsq
                                u3 = g0tmp(2,1,ir_)*dvec(itheta,3,m,l)*inv_lsq

                                udr1 = g0dertmp(1,1,ir_)*dvec(itheta,1,m,l)
                                udr2 = g0dertmp(2,1,ir_)*dvec(itheta,2,m,l)*inv_lsq
                                udr3 = g0dertmp(2,1,ir_)*dvec(itheta,3,m,l)*inv_lsq

                                udt1 = g0tmp(1,1,ir_)*dvecdt(itheta,1,m,l)
                                udt2 = g0tmp(2,1,ir_)*dvecdt(itheta,2,m,l)*inv_lsq
                                udt3 = g0tmp(2,1,ir_)*dvecdt(itheta,3,m,l)*inv_lsq

                                udp1 = g0tmp(1,1,ir_)*dvecdp(itheta,1,m,l)
                                udp2 = g0tmp(2,1,ir_)*dvecdp(itheta,2,m,l)*inv_lsq
                                udp3 = g0tmp(2,1,ir_)*dvecdp(itheta,3,m,l)*inv_lsq

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
                                stress(itheta,ir_,1,1) = stress(itheta,ir_,1,1)+Cval*uder11+Fval*uder22+Fval*uder33
                                stress(itheta,ir_,2,1) = stress(itheta,ir_,2,1)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder33
                                stress(itheta,ir_,3,1) = stress(itheta,ir_,3,1)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder22

                                stress(itheta,ir_,4,1) = stress(itheta,ir_,4,1)+Lval*(uder12+uder21)
                                stress(itheta,ir_,5,1) = stress(itheta,ir_,5,1)+Lval*(uder13+uder31)
                                stress(itheta,ir_,6,1) = stress(itheta,ir_,6,1)+Nval*(uder23+uder32)

                                displacement(itheta,ir_,1,1) = u1 + displacement(itheta,ir_,1,1)
                                displacement(itheta,ir_,2,1) = u2 + displacement(itheta,ir_,2,1)
                                displacement(itheta,ir_,3,1) = u3 + displacement(itheta,ir_,3,1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! case imt = 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!                               call calup
                                u1 = g0tmp(1,2,ir_)*dvec(itheta,1,m,l)
                                u2 = g0tmp(2,2,ir_)*dvec(itheta,2,m,l)*inv_lsq
                                u3 = g0tmp(2,2,ir_)*dvec(itheta,3,m,l)*inv_lsq

                                udr1 = g0dertmp(1,2,ir_)*dvec(itheta,1,m,l)
                                udr2 = g0dertmp(2,2,ir_)*dvec(itheta,2,m,l)*inv_lsq
                                udr3 = g0dertmp(2,2,ir_)*dvec(itheta,3,m,l)*inv_lsq

                                udt1 = g0tmp(1,2,ir_)*dvecdt(itheta,1,m,l)
                                udt2 = g0tmp(2,2,ir_)*dvecdt(itheta,2,m,l)*inv_lsq
                                udt3 = g0tmp(2,2,ir_)*dvecdt(itheta,3,m,l)*inv_lsq

                                udp1 = g0tmp(1,2,ir_)*dvecdp(itheta,1,m,l)
                                udp2 = g0tmp(2,2,ir_)*dvecdp(itheta,2,m,l)*inv_lsq
                                udp3 = g0tmp(2,2,ir_)*dvecdp(itheta,3,m,l)*inv_lsq

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
                                stress(itheta,ir_,1,2) = stress(itheta,ir_,1,2)+Cval*uder11+Fval*uder22+Fval*uder33
                                stress(itheta,ir_,2,2) = stress(itheta,ir_,2,2)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder33
                                stress(itheta,ir_,3,2) = stress(itheta,ir_,3,2)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder22

                                stress(itheta,ir_,4,2) = stress(itheta,ir_,4,2)+Lval*(uder12+uder21)
                                stress(itheta,ir_,5,2) = stress(itheta,ir_,5,2)+Lval*(uder13+uder31)
                                stress(itheta,ir_,6,2) = stress(itheta,ir_,6,2)+Nval*(uder23+uder32)

                                displacement(itheta,ir_,1,2) = u1 + displacement(itheta,ir_,1,2)
                                displacement(itheta,ir_,2,2) = u2 + displacement(itheta,ir_,2,2)
                                displacement(itheta,ir_,3,2) = u3 + displacement(itheta,ir_,3,2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! case imt = 3 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!                               call calup
                                u1 = g0tmp(1,3,ir_)*dvec(itheta,1,m,l)
                                u2 = g0tmp(2,3,ir_)*dvec(itheta,2,m,l)*inv_lsq
                                u3 = g0tmp(2,3,ir_)*dvec(itheta,3,m,l)*inv_lsq

                                udr1 = g0dertmp(1,3,ir_)*dvec(itheta,1,m,l)
                                udr2 = g0dertmp(2,3,ir_)*dvec(itheta,2,m,l)*inv_lsq
                                udr3 = g0dertmp(2,3,ir_)*dvec(itheta,3,m,l)*inv_lsq

                                udt1 = g0tmp(1,3,ir_)*dvecdt(itheta,1,m,l)
                                udt2 = g0tmp(2,3,ir_)*dvecdt(itheta,2,m,l)*inv_lsq
                                udt3 = g0tmp(2,3,ir_)*dvecdt(itheta,3,m,l)*inv_lsq

                                udp1 = g0tmp(1,3,ir_)*dvecdp(itheta,1,m,l)
                                udp2 = g0tmp(2,3,ir_)*dvecdp(itheta,2,m,l)*inv_lsq
                                udp3 = g0tmp(2,3,ir_)*dvecdp(itheta,3,m,l)*inv_lsq

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
                                stress(itheta,ir_,1,3) = stress(itheta,ir_,1,3)+Cval*uder11+Fval*uder22+Fval*uder33
                                stress(itheta,ir_,2,3) = stress(itheta,ir_,2,3)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder33
                                stress(itheta,ir_,3,3) = stress(itheta,ir_,3,3)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder22

                                stress(itheta,ir_,4,3) = stress(itheta,ir_,4,3)+Lval*(uder12+uder21)
                                stress(itheta,ir_,5,3) = stress(itheta,ir_,5,3)+Lval*(uder13+uder31)
                                stress(itheta,ir_,6,3) = stress(itheta,ir_,6,3)+Nval*(uder23+uder32)

                                displacement(itheta,ir_,1,3) = u1 + displacement(itheta,ir_,1,3)
                                displacement(itheta,ir_,2,3) = u2 + displacement(itheta,ir_,2,3)
                                displacement(itheta,ir_,3,3) = u3 + displacement(itheta,ir_,3,3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! case imt = 4 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!                               call calup
                                u1 = g0tmp(1,4,ir_)*dvec(itheta,1,m,l)
                                u2 = g0tmp(2,4,ir_)*dvec(itheta,2,m,l)*inv_lsq
                                u3 = g0tmp(2,4,ir_)*dvec(itheta,3,m,l)*inv_lsq

                                udr1 = g0dertmp(1,4,ir_)*dvec(itheta,1,m,l)
                                udr2 = g0dertmp(2,4,ir_)*dvec(itheta,2,m,l)*inv_lsq
                                udr3 = g0dertmp(2,4,ir_)*dvec(itheta,3,m,l)*inv_lsq

                                udt1 = g0tmp(1,4,ir_)*dvecdt(itheta,1,m,l)
                                udt2 = g0tmp(2,4,ir_)*dvecdt(itheta,2,m,l)*inv_lsq
                                udt3 = g0tmp(2,4,ir_)*dvecdt(itheta,3,m,l)*inv_lsq

                                udp1 = g0tmp(1,4,ir_)*dvecdp(itheta,1,m,l)
                                udp2 = g0tmp(2,4,ir_)*dvecdp(itheta,2,m,l)*inv_lsq
                                udp3 = g0tmp(2,4,ir_)*dvecdp(itheta,3,m,l)*inv_lsq

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
                                stress(itheta,ir_,1,4) = stress(itheta,ir_,1,4)+Cval*uder11+Fval*uder22+Fval*uder33
                                stress(itheta,ir_,2,4) = stress(itheta,ir_,2,4)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder33
                                stress(itheta,ir_,3,4) = stress(itheta,ir_,3,4)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder22

                                stress(itheta,ir_,4,4) = stress(itheta,ir_,4,4)+Lval*(uder12+uder21)
                                stress(itheta,ir_,5,4) = stress(itheta,ir_,5,4)+Lval*(uder13+uder31)
                                stress(itheta,ir_,6,4) = stress(itheta,ir_,6,4)+Nval*(uder23+uder32)

                                displacement(itheta,ir_,1,4) = u1 + displacement(itheta,ir_,1,4)
                                displacement(itheta,ir_,2,4) = u2 + displacement(itheta,ir_,2,4)
                                displacement(itheta,ir_,3,4) = u3 + displacement(itheta,ir_,3,4)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! case imt = 5 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!                               call calup
                                u1 = g0tmp(1,5,ir_)*dvec(itheta,1,m,l)
                                u2 = g0tmp(2,5,ir_)*dvec(itheta,2,m,l)*inv_lsq
                                u3 = g0tmp(2,5,ir_)*dvec(itheta,3,m,l)*inv_lsq

                                udr1 = g0dertmp(1,5,ir_)*dvec(itheta,1,m,l)
                                udr2 = g0dertmp(2,5,ir_)*dvec(itheta,2,m,l)*inv_lsq
                                udr3 = g0dertmp(2,5,ir_)*dvec(itheta,3,m,l)*inv_lsq

                                udt1 = g0tmp(1,5,ir_)*dvecdt(itheta,1,m,l)
                                udt2 = g0tmp(2,5,ir_)*dvecdt(itheta,2,m,l)*inv_lsq
                                udt3 = g0tmp(2,5,ir_)*dvecdt(itheta,3,m,l)*inv_lsq

                                udp1 = g0tmp(1,5,ir_)*dvecdp(itheta,1,m,l)
                                udp2 = g0tmp(2,5,ir_)*dvecdp(itheta,2,m,l)*inv_lsq
                                udp3 = g0tmp(2,5,ir_)*dvecdp(itheta,3,m,l)*inv_lsq

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
                                stress(itheta,ir_,1,5) = stress(itheta,ir_,1,5)+Cval*uder11+Fval*uder22+Fval*uder33
                                stress(itheta,ir_,2,5) = stress(itheta,ir_,2,5)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder33
                                stress(itheta,ir_,3,5) = stress(itheta,ir_,3,5)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder22

                                stress(itheta,ir_,4,5) = stress(itheta,ir_,4,5)+Lval*(uder12+uder21)
                                stress(itheta,ir_,5,5) = stress(itheta,ir_,5,5)+Lval*(uder13+uder31)
                                stress(itheta,ir_,6,5) = stress(itheta,ir_,6,5)+Nval*(uder23+uder32)

                                displacement(itheta,ir_,1,5) = u1 + displacement(itheta,ir_,1,5)
                                displacement(itheta,ir_,2,5) = u2 + displacement(itheta,ir_,2,5)
                                displacement(itheta,ir_,3,5) = u3 + displacement(itheta,ir_,3,5)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! case imt = 6 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!                               call calup
                                u1 = g0tmp(1,6,ir_)*dvec(itheta,1,m,l)
                                u2 = g0tmp(2,6,ir_)*dvec(itheta,2,m,l)*inv_lsq
                                u3 = g0tmp(2,6,ir_)*dvec(itheta,3,m,l)*inv_lsq

                                udr1 = g0dertmp(1,6,ir_)*dvec(itheta,1,m,l)
                                udr2 = g0dertmp(2,6,ir_)*dvec(itheta,2,m,l)*inv_lsq
                                udr3 = g0dertmp(2,6,ir_)*dvec(itheta,3,m,l)*inv_lsq

                                udt1 = g0tmp(1,6,ir_)*dvecdt(itheta,1,m,l)
                                udt2 = g0tmp(2,6,ir_)*dvecdt(itheta,2,m,l)*inv_lsq
                                udt3 = g0tmp(2,6,ir_)*dvecdt(itheta,3,m,l)*inv_lsq

                                udp1 = g0tmp(1,6,ir_)*dvecdp(itheta,1,m,l)
                                udp2 = g0tmp(2,6,ir_)*dvecdp(itheta,2,m,l)*inv_lsq
                                udp3 = g0tmp(2,6,ir_)*dvecdp(itheta,3,m,l)*inv_lsq

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
                                stress(itheta,ir_,1,6) = stress(itheta,ir_,1,6)+Cval*uder11+Fval*uder22+Fval*uder33
                                stress(itheta,ir_,2,6) = stress(itheta,ir_,2,6)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder33
                                stress(itheta,ir_,3,6) = stress(itheta,ir_,3,6)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder22

                                stress(itheta,ir_,4,6) = stress(itheta,ir_,4,6)+Lval*(uder12+uder21)
                                stress(itheta,ir_,5,6) = stress(itheta,ir_,5,6)+Lval*(uder13+uder31)
                                stress(itheta,ir_,6,6) = stress(itheta,ir_,6,6)+Nval*(uder23+uder32)

                                displacement(itheta,ir_,1,6) = u1 + displacement(itheta,ir_,1,6)
                                displacement(itheta,ir_,2,6) = u2 + displacement(itheta,ir_,2,6)
                                displacement(itheta,ir_,3,6) = u3 + displacement(itheta,ir_,3,6)

                             enddo! of loop on itheta
                          endif  ! of if (l <= lmax_r(ir_)) then ... i.e. on the stack point

                       enddo   !  of loop on ir_ ie. stack point
                    !enddo ! mt-loop

                    endif   ! l-branch for calu
                 endif
              enddo            ! m-loop end
           enddo               ! l-loop end

        !********************

        if (myrank == 0 .and. SLOW_DEBUG_MODE) write(24,*) i, dble(i)/tlen, llog-1

!       il faut faire un allgatherv pour recuperer stresssgl et displacementsg
        ! stress(1:6,1:6,1:r_n,1:nsta) = stress(1:6,1:6,1:r_n,1:nsta)/dcmplx(0,omega)
        ! stresssngl(1:6,1:6,1:r_n,1:nsta) = stress(1:6,1:6,1:r_n,1:nsta)
        ! displacementsngl(1:3,1:6,1:r_n,1:nsta) = displacement(1:3,1:6,1:r_n,1:nsta)

        inverse_of_omega_complex = 1.d0 / dcmplx(0,omega)

        do ir_=1,r_n
!IBM* ASSERT (MINITERCNT(1000))
!DIR$ SIMD
!DIR$ loop count min(1000)
          do ista=1,nsta
              stresssngl(1,1,ir_,ista) = stress(ista,ir_,1,1) * inverse_of_omega_complex
              stresssngl(2,1,ir_,ista) = stress(ista,ir_,2,1) * inverse_of_omega_complex
              stresssngl(3,1,ir_,ista) = stress(ista,ir_,3,1) * inverse_of_omega_complex
              stresssngl(4,1,ir_,ista) = stress(ista,ir_,4,1) * inverse_of_omega_complex
              stresssngl(5,1,ir_,ista) = stress(ista,ir_,5,1) * inverse_of_omega_complex
              stresssngl(6,1,ir_,ista) = stress(ista,ir_,6,1) * inverse_of_omega_complex
              stresssngl(1,2,ir_,ista) = stress(ista,ir_,1,2) * inverse_of_omega_complex
              stresssngl(2,2,ir_,ista) = stress(ista,ir_,2,2) * inverse_of_omega_complex
              stresssngl(3,2,ir_,ista) = stress(ista,ir_,3,2) * inverse_of_omega_complex
              stresssngl(4,2,ir_,ista) = stress(ista,ir_,4,2) * inverse_of_omega_complex
              stresssngl(5,2,ir_,ista) = stress(ista,ir_,5,2) * inverse_of_omega_complex
              stresssngl(6,2,ir_,ista) = stress(ista,ir_,6,2) * inverse_of_omega_complex
              stresssngl(1,3,ir_,ista) = stress(ista,ir_,1,3) * inverse_of_omega_complex
              stresssngl(2,3,ir_,ista) = stress(ista,ir_,2,3) * inverse_of_omega_complex
              stresssngl(3,3,ir_,ista) = stress(ista,ir_,3,3) * inverse_of_omega_complex
              stresssngl(4,3,ir_,ista) = stress(ista,ir_,4,3) * inverse_of_omega_complex
              stresssngl(5,3,ir_,ista) = stress(ista,ir_,5,3) * inverse_of_omega_complex
              stresssngl(6,3,ir_,ista) = stress(ista,ir_,6,3) * inverse_of_omega_complex
              stresssngl(1,4,ir_,ista) = stress(ista,ir_,1,4) * inverse_of_omega_complex
              stresssngl(2,4,ir_,ista) = stress(ista,ir_,2,4) * inverse_of_omega_complex
              stresssngl(3,4,ir_,ista) = stress(ista,ir_,3,4) * inverse_of_omega_complex
              stresssngl(4,4,ir_,ista) = stress(ista,ir_,4,4) * inverse_of_omega_complex
              stresssngl(5,4,ir_,ista) = stress(ista,ir_,5,4) * inverse_of_omega_complex
              stresssngl(6,4,ir_,ista) = stress(ista,ir_,6,4) * inverse_of_omega_complex
              stresssngl(1,5,ir_,ista) = stress(ista,ir_,1,5) * inverse_of_omega_complex
              stresssngl(2,5,ir_,ista) = stress(ista,ir_,2,5) * inverse_of_omega_complex
              stresssngl(3,5,ir_,ista) = stress(ista,ir_,3,5) * inverse_of_omega_complex
              stresssngl(4,5,ir_,ista) = stress(ista,ir_,4,5) * inverse_of_omega_complex
              stresssngl(5,5,ir_,ista) = stress(ista,ir_,5,5) * inverse_of_omega_complex
              stresssngl(6,5,ir_,ista) = stress(ista,ir_,6,5) * inverse_of_omega_complex
              stresssngl(1,6,ir_,ista) = stress(ista,ir_,1,6) * inverse_of_omega_complex
              stresssngl(2,6,ir_,ista) = stress(ista,ir_,2,6) * inverse_of_omega_complex
              stresssngl(3,6,ir_,ista) = stress(ista,ir_,3,6) * inverse_of_omega_complex
              stresssngl(4,6,ir_,ista) = stress(ista,ir_,4,6) * inverse_of_omega_complex
              stresssngl(5,6,ir_,ista) = stress(ista,ir_,5,6) * inverse_of_omega_complex
              stresssngl(6,6,ir_,ista) = stress(ista,ir_,6,6) * inverse_of_omega_complex
           enddo
        enddo

        do ir_=1,r_n
!IBM* ASSERT (MINITERCNT(1000))
!DIR$ SIMD
!DIR$ loop count min(1000)
          do ista = 1,nsta
              displacementsngl(1,1,ir_,ista) = displacement(ista,ir_,1,1)
              displacementsngl(2,1,ir_,ista) = displacement(ista,ir_,2,1)
              displacementsngl(3,1,ir_,ista) = displacement(ista,ir_,3,1)
              displacementsngl(1,2,ir_,ista) = displacement(ista,ir_,1,2)
              displacementsngl(2,2,ir_,ista) = displacement(ista,ir_,2,2)
              displacementsngl(3,2,ir_,ista) = displacement(ista,ir_,3,2)
              displacementsngl(1,3,ir_,ista) = displacement(ista,ir_,1,3)
              displacementsngl(2,3,ir_,ista) = displacement(ista,ir_,2,3)
              displacementsngl(3,3,ir_,ista) = displacement(ista,ir_,3,3)
              displacementsngl(1,4,ir_,ista) = displacement(ista,ir_,1,4)
              displacementsngl(2,4,ir_,ista) = displacement(ista,ir_,2,4)
              displacementsngl(3,4,ir_,ista) = displacement(ista,ir_,3,4)
              displacementsngl(1,5,ir_,ista) = displacement(ista,ir_,1,5)
              displacementsngl(2,5,ir_,ista) = displacement(ista,ir_,2,5)
              displacementsngl(3,5,ir_,ista) = displacement(ista,ir_,3,5)
              displacementsngl(1,6,ir_,ista) = displacement(ista,ir_,1,6)
              displacementsngl(2,6,ir_,ista) = displacement(ista,ir_,2,6)
              displacementsngl(3,6,ir_,ista) = displacement(ista,ir_,3,6)
           enddo
        enddo

        Coeff = 3*6*r_n*nsta
        call mpi_gatherv( displacementsngl,Coeff ,MPI_COMPLEX,displacementsngl_global,NbInd,Ind,MPI_COMPLEX, &
                                         0,SubCommunicators, ierr )
        Coeff = 6*6*r_n*nsta
        call mpi_gatherv( stresssngl,Coeff ,MPI_COMPLEX,stresssngl_global,NbIndSt,IndSt,MPI_COMPLEX, &
                                         0,SubCommunicators, ierr )
     endif

!    irec = irec + 1

     !******** ECRITURE SUR LE DISQUE

!! DK DK not a good idea to give indices when writing the whole array: some compilers could create copies or use loops in such a case
!! DK DK unformatted sequential I/Os (because sequential I/O is usually faster than direct-access I/O)
!! DK DK and writing a single big file is also much better
     if (myrank == 0) then
        write(44) stresssngl_global
        write(44) displacementsngl_global
        close(44)
        close(34)
     endif

     if (myrank == 0 .and. SLOW_DEBUG_MODE) then
        call date_and_time(datex,timex)
        write(25,'(/a,i5,a,a4,a,a2,a,a2,a,a2,a,a2,a,a4)') '    Frequency-index ', i, ' :', &
             datex(1:4),'-',datex(5:6),'-',datex(7:8),'.  ',timex(1:2),':',timex(3:4),':',timex(5:8)
     endif

!**************************************************************************************
  enddo     ! of frequency loop
!**************************************************************************************

  if (myrank == 0) then
     !close(34)
!    close(35)
     !close(44)
!    close(45)
     if (SLOW_DEBUG_MODE) then
       close(24)
       call date_and_time(datex,timex)
       write(25,'(/a,a4,a,a2,a,a2,a,a2,a,a2,a,a4)') '    Finishing date and time:                     ', &
         datex(1:4),'-',datex(5:6),'-',datex(7:8),'.  ',timex(1:2),':',timex(3:4),':',timex(5:8)
       write(25,*) 'ir_ lmax'
       do ir_=1,r_n
          write(25,*) ir_,lmax_r(ir_)
       enddo
       close(25)
     endif
  endif

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call MPI_COMM_FREE(SubCommunicators,ierr)
  call MPI_FINALIZE(ierr)

end program TraPSV

