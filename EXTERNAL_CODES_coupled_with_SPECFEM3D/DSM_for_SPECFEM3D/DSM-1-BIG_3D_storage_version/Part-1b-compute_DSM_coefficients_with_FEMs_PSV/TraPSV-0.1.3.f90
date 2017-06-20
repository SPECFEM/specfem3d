
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
   INTEGER myrank,nbproc,ierr

   include "../shared/constants.h"

!! DK DK made these variables become parameters to allow for compiler optimization
  integer, parameter :: nsta = 1 ! on ne lit que la premiere station dans cette partie 1 du code
  integer, parameter :: theta_n = nsta
  integer :: nsta_ignored_in_this_code ! read from the input file, but then not used

  !------------------------- <  < input matrix >>----------------------------------

  character(120) :: outputDir,psvmodel,modelname,stationsinf
  character(120) :: list,list1
  character(40) :: datex,timex
  real(kind(0d0)), parameter :: pi=3.1415926535897932d0
  real(kind(0d0)), parameter :: re=1.d-2, ratc=1.d-10, ratl=1.d-4
  !integer, parameter :: maxlmax = 35000,maxlmax_g=1000 !! VM VM I moved this in constants.h

  real(kind(0d0)) :: tlen
  real(kind(0d0)) :: r0min, r0max, r0delta  !!! JUST FOR ONE DEPTH FOR THIS MOMENT !!
  real(kind(0d0)) :: r0lat, r0lon
  real(kind(0d0)), allocatable :: stla(:),stlo(:)
  real(kind(0d0)), allocatable :: r_(:),r0(:),theta(:),phi(:)
  real(kind(0d0)), allocatable :: A0sta(:),C0sta(:),F0sta(:),L0sta(:),N0sta(:)
  integer, allocatable :: updown(:),idum(:)
  real(kind(0d0)), allocatable :: rrsta(:,:)
  integer, allocatable :: iista(:,:)
  integer :: r_n,r0_n,ciista, ir_,ir0,imt

  character(120) :: coutfile
  integer :: imin, imax, imax_global
  integer :: itranslat


  ! --------------------------- <  < variables >>---------------------------
  ! variable for the trial function
  integer:: nnlayer,nlay
  integer, allocatable :: nlayer(:)
  integer:: nslay,nllay
  integer:: inlayer,jnlayer,jnslay,jnllay
  integer:: l,m
  real(kind(0d0)),allocatable:: ra(:)
  ! variable for the structure
  integer:: nzone,isl,ill,nsl,nll
  integer,allocatable:: iphase(:)
  integer::ndc,vnp
  real(kind(0d0)):: rmin,rmax
  real(kind(0d0)),allocatable:: vrmin(:),vrmax(:),rrho(:,:),vpv(:,:),vph(:,:),vsv(:,:), &
  vsh(:,:),eta(:,:),qmu(:),qkappa(:)
  real(kind(0d0)),allocatable::vra(:),rho(:),kappa(:)
  real(kind(0d0)),allocatable::ecKx(:) !3*Kx=3A-4N
  real(kind(0d0)),allocatable::ecKy(:) !3*Ky=3F+2N
  real(kind(0d0)),allocatable::ecKz(:) !3*Kz=2F+C
  real(kind(0d0)),allocatable::mu(:),ecL(:),ecN(:),rhoinv(:),kappainv(:)
  complex(kind(0d0)),allocatable:: coef1(:),coef2(:),coef(:)
  real(kind(0d0)):: omega,omegai
  ! variable for the source
  integer:: spn,ns
  real(kind(0d0)):: mt(3,3),spo
  real(kind(0d0)):: ecC0,ecF0,ecL0
  complex(kind(0d0)):: ya(4),yb(4),yc(4),yd(4)

  ! variable for the matrix elements
  complex(kind(0d0)),allocatable:: a0(:,:),a1(:,:),a2(:,:), a(:,:), c(:,:) !, ctmp(:,:)
  real(kind(0d0)), allocatable :: t(:)
  real(kind(0d0)), allocatable :: h1x(:), h1y(:), h1z(:), h2L(:), h2N(:), h3ax(:), h3ay(:), h3az(:), h4aL(:), &
  h4aN(:), h5ax(:), h5ay(:), h5az(:), h6aL(:), h6aN(:), h3x(:), h3y(:), h3z(:), h4L(:), h4N(:), &
  h5x(:), h5y(:), h5z(:), h6L(:), h6N(:), h7x(:), h7y(:), h7z(:), h8L(:), h8N(:), h3mx(:,:), &
  h3my(:,:), h3mz(:,:), h5mx(:,:), h5my(:,:), h5mz(:,:), h4m1L(:,:), h4m1N(:,:), h4m2L(:,:), &
  h4m2N(:,:), h6m1L(:,:), h6m1N(:,:), h6m2L(:,:), h6m2N(:,:)
  real(kind(0d0)),allocatable:: p1(:),p2(:),p3(:)
  complex(kind(0d0)),allocatable:: g0(:),tabg0(:,:,:,:,:),tabg0der(:,:,:,:,:),tabg0_small_buffer(:,:,:,:)
  complex(kind(0d0)),allocatable:: d0(:)
  complex(kind(0d0)):: g0tmp(2),g0dertmp(2) ! forward
  ! variable for the stack point
  integer,allocatable:: isp(:),issp(:),ilsp(:),jssp(:),jsp(:),ksp(:),lsp(:),lmax_r(:),ismall_r(:),arret(:)
  integer::isdr,jsdr,ildr,cista,cksta
  ! variables for the output stack point
  integer,allocatable:: istazone(:)
  integer,allocatable:: ksta(:)   ! output stack point for g
  integer,allocatable:: jsta(:)   ! output stack point for d
  ! variables for the gridding
  integer,allocatable:: jjdr(:),kkdr(:)
  integer:: jdr,kdr
  real(kind(0d0)),allocatable:: vmin(:),gridpar(:),dzpar(:),maxamp_r(:)
  ! variables for l cut off
  integer:: kc,lsuf,sufzone,ismall,llog,llog0,index_l,lref
  real(kind(0d0)):: maxamp,amp,ampratio
  ! variables for the numerical integration
  complex(kind(0d0)):: anum(4,4,10),bnum(4,4,10)

  ! other variables
  integer:: i,j,nn,ier,itmp,jtmp,mtmp,kkdr0,nn0,ig2
  integer:: ll(12),lli(12),llj(12)
  real(kind(0d0)):: eps,l2,lsq,bazi
  real(kind(0d0)),allocatable:: work(:)
  complex(kind(0d0)), allocatable :: z(:), w(:),cwork(:)

  logical :: we_left_the_loop_through_exit

  !-----------------------------------------------------------------------
  character(len=11) debug

  data eps /-1.d0/

  ! ajout de variables pour MPI --- VM VM

  integer, allocatable :: Ifrq(:),Ifrq2(:,:)
  integer irec
  integer imin_glob, imax_glob
  integer para
  !! VM VM timer
  real time_to_write, time_to_dcsymbdl, time_to_dcsbdlv,time_to_define_matrix,time_total
  real start_time_to_write, start_time_to_dcsymbdl, start_time_to_dcsbdlv, &
       start_time_to_define_matrix,start_time_total
  real finish_time_to_write, finish_time_to_dcsymbdl, finish_time_to_dcsbdlv, &
       finish_time_to_define_matrix,finish_time_total
  real total_global_time
  !--------------------------------------------------------------------------

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nbproc,ierr)
  if (USE_TIMER) then
     call cpu_time(start_time_total)
     time_to_write=0.
     time_to_dcsymbdl=0.
     time_to_dcsbdlv=0.
     time_to_define_matrix=0.
     time_total=0.
  endif
  allocate(Ifrq(nbproc+1))
  allocate(Ifrq2(0:nbproc-1,2))
  if (myrank == 0) call ReadFrqsByProc(Ifrq,Ifrq2,para,nbproc+1)
  call MPI_Bcast(Ifrq,nbproc+1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(Ifrq2,2*nbproc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(para,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  ! si on veut un fichier de debugage
  write(debug,'(a6,i5.5)') 'debug',myrank
  if (myrank == 0) then
  call pinputTra_Part1(outputDir,psvmodel,modelname,stationsinf,tlen,imin_glob,imax_glob,r0min,r0max,r0delta,r0lat, &
  r0lon,itranslat,myrank)
  endif
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
  itranslat=0

  ! ------ VM VM  ! on peut choisir 2 options pour la parallelisation : on lit
  ! les frequences en
  if (para == 1) then
    imin=Ifrq(myrank+1)+ 1
    imax=Ifrq(myrank+2)
    if (myrank == 0) imin=Ifrq(myrank+1) ! on commence au debut du fichier
    imax_global = maxval(Ifrq)
  else
    imin=Ifrq2(myrank,1)
    imax=Ifrq2(myrank,2)
  endif

  imax_global = imax_glob  ! frequence max lue dans le fichier parametres :
                             ! pour calculer la grille de discretisation

  if (myrank == 0) then
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
  if (myrank == 0) then
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
  if (myrank == 0) then
     open (1,file=stationsinf,status='old',action='read',position='rewind')
     open(2,file='recdepth')
     read(2,*) r_n
     read(1,*) nsta_ignored_in_this_code
  endif
  call MPI_Bcast(r_n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  allocate(r_(1:r_n))
  allocate(A0sta(1:r_n))
  allocate(C0sta(1:r_n))
  allocate(F0sta(1:r_n))
  allocate(L0sta(1:r_n))
  allocate(N0sta(1:r_n))
  allocate(theta(1:theta_n))
  allocate(phi(1:theta_n))
  allocate(stla(1:nsta))
  allocate(stlo(1:nsta))
  allocate(updown(1:r_n))
  allocate(idum(1:r_n))
  allocate(lmax_r(r_n),ismall_r(r_n),arret(r_n))
  allocate(maxamp_r(r_n))
!! DK DK moved index 1..maxlmax_g to the first position in order to be able to vectorize a loop
  allocate(tabg0(maxlmax_g,2,6,-2:2,r_n))
  allocate(tabg0der(maxlmax_g,2,6,-2:2,r_n))
  allocate(tabg0_small_buffer(maxlmax_g,2,6,-2:2))

  if (myrank == 0) then
     do i = 1,nsta
        read(1,*) stla(i),stlo(i)
        if (itranslat == 1) call translat(stla(i),stla(i))
        call epitra(r0lat,r0lon,stla(i),stlo(i),theta(i),phi(i),bazi)
     enddo
  endif
  call MPI_Bcast(stla,nsta,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(stlo,nsta,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(theta,nsta,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(phi,nsta,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)


  if (myrank == 0) then
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

  ! source depths
  ! for this moment we don't put any necessary allocation
  r0_n =  int((r0max-r0min)/r0delta)+1
  allocate(r0(1:r0_n))
  do i = 1, r0_n
     r0(i) = r0min + dble(i-1)*r0delta
  enddo
  ir0 = r0_n

  if ( (r0(ir0) < rmin) .or. (r0(ir0) > rmax) ) stop 'Location of the source is improper.'

  ! computation de nombre et la location des points de grid
  call calgrid( nzone,vrmin,vrmax,vpv,vsv,rmin,rmax,imax_global,1,tlen,vmin,gridpar,dzpar )

  call calra_psv(nnlayer,inlayer,jnlayer,jnslay,jnllay,gridpar,dzpar,nzone,vrmin,vrmax,iphase, &
       rmin,rmax,nslay,nllay,nlayer,re )

  allocate(ra(1:nnlayer+nzone+1))
  call calra2_psv(nnlayer,gridpar,dzpar,nzone,vrmin,vrmax,rmin,rmax,nlayer,ra,re,r_n,r_, &
  rrsta,iista,r0(ir0),cista,iphase,istazone,ciista,updown)
  if (myrank == 0) then
     ! ecrire les parametres de la grille
!!$     write(*,*) 'nzone ', nzone
!!$     write(*,*)' gridpar '
!!$     write(*,*) gridpar
!!$     write(*,*) ' dzpar '
!!$     write(*,*) dzpar
     do i=1,nzone
        write(*,*) '             GRID USED FOR MAXIMUM FREQUENCY :',imax_global/tlen
        write(*,'("  ZONE :",i5,"  RMIN:",f12.4,"  RMAX:",f12.4)') i,vrmin(i),vrmax(i)
        write(*,'("   NB PTS :",i15,"   STEP :",f30.15)') nlayer(i),(vrmax(i)-vrmin(i)) / nlayer(i)
        !write(*,*)
     enddo
  endif

  nlay = nnlayer

  allocate(vra(1:nlay+2*nzone+1))
  allocate(rho(1:nlay+2*nzone+1))
  allocate(kappa(1:nlay+2*nzone+1))
  allocate(ecKx(1:nlay+2*nzone+1)) !3*Kx=3A-4N
  allocate(ecKy(1:nlay+2*nzone+1)) !3*Ky=3F+2N
  allocate(ecKz(1:nlay+2*nzone+1)) !3*Kz=2F+C
  allocate(mu(1:nlay+2*nzone+1))
  allocate(ecL(1:nlay+2*nzone+1))
  allocate(ecN(1:nlay+2*nzone+1))
  allocate(rhoinv(1:nlay+2*nzone+1))
  allocate(kappainv(nlay+2*nzone+1))
  allocate(a0(1:4,1:2*(2*(nslay+1)+(nllay+1)+2*nzone)))
  allocate(a1(1:4,1:2*(2*(nslay+1)+(nllay+1)+2*nzone)))
  allocate(a2(1:4,1:2*(2*(nslay+1)+(nllay+1)+2*nzone)))
  !allocate(a(1:4,1:2*(nslay+1)+(nllay+1)))
  !! modif vadim
  allocate(a(1:4,1:2*(2*(nslay+1)+(nllay+1)+2*nzone)))
  !! modif
  !allocate(c(1:2,1:(nslay+1)+(nllay+1)))
  allocate(c(1:2,1:(nslay+1)+(nllay+1) + nzone))
  !allocate(ctmp(1:2,1:(nslay+1)+(nllay+1)))  ! ne sert pas
  allocate(t(1:8*nslay))
  allocate(h1x(1:8*nslay))
  allocate(h1y(1:8*nslay))
  allocate(h1z(1:8*nslay))
  allocate(h2L(1:8*nslay))
  allocate(h2N(1:8*nslay))
  allocate(h3ax(1:8*nslay))
  allocate(h3ay(1:8*nslay))
  allocate(h3az(1:8*nslay))
  allocate(h4aL(1:8*nslay))
  allocate(h4aN(1:8*nslay))
  allocate(h5ax(1:8*nslay))
  allocate(h5ay(1:8*nslay))
  allocate(h5az(1:8*nslay))
  allocate(h6aL(1:8*nslay))
  allocate(h6aN(1:8*nslay))
  allocate(h3x(1:8*nslay))
  allocate(h3y(1:8*nslay))
  allocate(h3z(1:8*nslay))
  allocate(h4L(1:8*nslay))
  allocate(h4N(1:8*nslay))
  allocate(h5x(1:8*nslay))
  allocate(h5y(1:8*nslay))
  allocate(h5z(1:8*nslay))
  allocate(h6L(1:8*nslay))
  allocate(h6N(1:8*nslay))
  allocate(h7x(1:8*nslay))
  allocate(h7y(1:8*nslay))
  allocate(h7z(1:8*nslay))
  allocate(h8L(1:8*nslay))
  allocate(h8N(1:8*nslay))
  allocate(h3mx(-2:1,1:2*(nslay+nzone)))
  allocate(h3my(-2:1,1:2*(nslay+nzone)))
  allocate(h3mz(-2:1,1:2*(nslay+nzone)))
  allocate(h5mx(-1:2,1:2*(nslay+nzone)))
  allocate(h5my(-1:2,1:2*(nslay+nzone)))
  allocate(h5mz(-1:2,1:2*(nslay+nzone)))
  allocate(h4m1L(-1:2,1:2*(nslay+nzone)))
  allocate(h4m1N(-1:2,1:2*(nslay+nzone)))
  allocate(h4m2L(-2:1,1:2*(nslay+nzone)))
  allocate(h4m2N(-2:1,1:2*(nslay+nzone)))
  allocate(h6m1L(-1:2,1:2*(nslay+nzone)))
  allocate(h6m1N(-1:2,1:2*(nslay+nzone)))
  allocate(h6m2L(-2:1,1:2*(nslay+nzone)))
  allocate(h6m2N(-2:1,1:2*(nslay+nzone)))
  allocate(p1(1:8*nllay))
  allocate(p2(1:8*nllay))
  allocate(p3(1:8*nllay))
  allocate(g0(1:2*(nslay+1)+(nllay+1)+nzone))
  allocate(d0(1:(nslay+1)+(nllay+1)+nzone))
  allocate(work(1:8*nslay))
  allocate(z(1:2*(nslay+1)+(nllay+1) + nzone))
  allocate(w(1:2*(nslay+1)+(nllay+1) + nzone))
  allocate(cwork(1:4*(16*nslay+4*nllay)),stat=ier)
  if (ier /= 0) stop 'error: not enough memory to allocate array'

  ! computing the stack points
  call calsp(nzone,ndc,nsl,nll,iphase,nlayer,nslay,nllay,isp,jsp,ksp,issp,ilsp,lsp,jssp,isdr,jsdr,ildr,jdr,kdr )

  ! computing the source location
  call calspo(nlay,nzone,ndc,vrmax,iphase,inlayer,r0(ir0),rmin,rmax,ra,isp,spo,spn )

  ! ******************* Computing the matrix elements *******************
  ! data initialization
  a = 0.d0
  t = 0.d0
  h1x = 0.d0
  h1y = 0.d0
  h1z = 0.d0
  h2L = 0.d0
  h2N = 0.d0
  h3ax = 0.d0
  h3ay = 0.d0
  h3az = 0.d0
  h4aL = 0.d0
  h4aN = 0.d0
  h5ax = 0.d0
  h5ay = 0.d0
  h5az = 0.d0
  h6aL = 0.d0
  h6aN = 0.d0
  h3x = 0.d0
  h3y = 0.d0
  h3z = 0.d0
  h4L = 0.d0
  h4N = 0.d0
  h5x = 0.d0
  h5y = 0.d0
  h5z = 0.d0
  h6L = 0.d0
  h6N = 0.d0
  h7x = 0.d0
  h7y = 0.d0
  h7z = 0.d0
  h8L = 0.d0
  h8N = 0.d0
  h3mx = 0.d0
  h3my = 0.d0
  h3mz = 0.d0
  h5mx = 0.d0
  h5my = 0.d0
  h5mz = 0.d0
  h4m1L = 0.d0
  h4m1N = 0.d0
  h4m2L = 0.d0
  h4m2N = 0.d0
  h6m1L = 0.d0
  h6m1N = 0.d0
  h6m2L = 0.d0
  h6m2N = 0.d0
  p1 = 0.d0
  p2 = 0.d0
  p3 = 0.d0

  ! computing the structure grid points
  call calstg(nlay,nzone,nzone,iphase,rrho,vpv,vph,vsv,vsh,eta,nlayer,ra,rmax,vnp,vra,rho,kappa,ecKx,ecKy, &
  ecKz,mu,ecL,ecN,r0(ir0),spn,ecC0,ecF0,ecL0 )
  call calinv( vnp,rho,kappa,rhoinv,kappainv )
  isl = 0
  ill = 0

  do i=1,ndc+1
     if ( iphase(i) == 1 ) then
        isl = isl + 1
        itmp = isdr+issp(isl)
        call calmatc( nlayer(i),vnp,vra,rho ,2,0,0,ra(isp(i)), t(itmp) )
        call caltl( nlayer(i),vnp,vra,rho,ra(isp(i)),work(itmp) )
        call calt( nlayer(i), t(itmp),work(itmp), t(itmp) )
        call calmatc( nlayer(i),vnp,vra,ecKx,0,0,0,ra(isp(i)),h1x(itmp) )
        call calhl( nlayer(i),vnp,vra,ecKx,ra(isp(i)),work(itmp) )
        call calt( nlayer(i), h1x(itmp),work(itmp),h1x(itmp))
        call calmatc( nlayer(i),vnp,vra,ecKy,0,0,0,ra(isp(i)),h1y(itmp) )
        call calhl( nlayer(i),vnp,vra,ecKy,ra(isp(i)),work(itmp) )
        call calt( nlayer(i), h1y(itmp),work(itmp),h1y(itmp))
        call calmatc( nlayer(i),vnp,vra,ecKz,0,0,0,ra(isp(i)),h1z(itmp) )
        call calhl( nlayer(i),vnp,vra,ecKz,ra(isp(i)),work(itmp) )
        call calt( nlayer(i), h1z(itmp),work(itmp),h1z(itmp))
        call calmatc( nlayer(i),vnp,vra,ecL ,0,0,0,ra(isp(i)),h2L(itmp) )
        call calhl( nlayer(i),vnp,vra,ecL,ra(isp(i)),work(itmp) )
        call calt( nlayer(i), h2L(itmp),work(itmp),h2L(itmp))
        call calmatc( nlayer(i),vnp,vra,ecN ,0,0,0,ra(isp(i)),h2N(itmp) )
        call calhl( nlayer(i),vnp,vra,ecN,ra(isp(i)),work(itmp) )
        call calt( nlayer(i), h2N(itmp),work(itmp),h2N(itmp))
        call calmatc( nlayer(i),vnp,vra,ecKx,1,0,1,ra(isp(i)),h5ax(itmp) )
        call calmatc( nlayer(i),vnp,vra,ecKy,1,0,1,ra(isp(i)),h5ay(itmp) )
        call calmatc( nlayer(i),vnp,vra,ecKz,1,0,1,ra(isp(i)),h5az(itmp) )
        call calmatc( nlayer(i),vnp,vra,ecL,1,0,1,ra(isp(i)),h6aL(itmp) )
        call calmatc( nlayer(i),vnp,vra,ecN,1,0,1,ra(isp(i)),h6aN(itmp) )
        call mtrnp( nlayer(i),h5ax(itmp),h3ax(itmp) )
        call mtrnp( nlayer(i),h5ay(itmp),h3ay(itmp) )
        call mtrnp( nlayer(i),h5az(itmp),h3az(itmp) )
        call mtrnp( nlayer(i),h6aL(itmp),h4aL(itmp) )
        call mtrnp( nlayer(i),h6aN(itmp),h4aN(itmp) )
        call calmatc( nlayer(i),vnp,vra,ecKx,2,1,1,ra(isp(i)),h7x(itmp) )
        call calmatc( nlayer(i),vnp,vra,ecKy,2,1,1,ra(isp(i)), h7y(itmp) )
        call calmatc( nlayer(i),vnp,vra,ecKz,2,1,1,ra(isp(i)), h7z(itmp) )
        call calmatc( nlayer(i),vnp,vra,ecL ,2,1,1,ra(isp(i)), h8L(itmp) )
        call calmatc( nlayer(i),vnp,vra,ecN ,2,1,1,ra(isp(i)), h8N(itmp) )
     else
        ill = ill + 1
        itmp = ildr+ilsp(ill)
        call calmatc( nlayer(i),vnp,vra,rhoinv,2,1,1,ra(isp(i)),p1(itmp) )
        call calmatc( nlayer(i),vnp,vra,rhoinv,0,0,0,ra(isp(i)),p2(itmp) )
        call calhl( nlayer(i),vnp,vra,rhoinv,ra(isp(i)),work(itmp) )
        call calt( nlayer(i),p2(itmp),work(itmp),p2(itmp) )
        call calmatc( nlayer(i),vnp,vra,kappainv,2,0,0,ra(isp(i)),p3(itmp) )
        call caltl( nlayer(i),vnp,vra,kappainv,ra(isp(i)),work(itmp) )
        call calt( nlayer(i),p3(itmp),work(itmp),p3(itmp) )
     endif
  enddo

  ! Computing the modified operator of the 1st derivative
  call caltstg(nlay,nzone,nzone,rrho,vpv,vph,vsv,vsh,eta,nlayer,ra,rmax,vra,kappa,ecKx,ecKy,ecKz,mu,ecL,ecN )
  isl = 0
  do i=1,ndc+1
     if ( iphase(i) == 1 ) then
        isl = isl + 1
        itmp = isdr+issp(isl)
        jtmp = isp(i)+i-1
        call calh5( nlayer(i),vra(jtmp),ecKx(jtmp),ra(isp(i)),work(itmp) )
        call submat( nlayer(i),h5ax(itmp),work(itmp),h5x(itmp) )
        call calh5( nlayer(i),vra(jtmp),ecKy(jtmp),ra(isp(i)),work(itmp) )
        call submat( nlayer(i),h5ay(itmp),work(itmp),h5y(itmp) )
        call calh5( nlayer(i),vra(jtmp),ecKz(jtmp),ra(isp(i)),work(itmp) )
        call submat( nlayer(i),h5az(itmp),work(itmp),h5z(itmp) )
        call calh5( nlayer(i),vra(jtmp),ecL(jtmp),ra(isp(i)),work(itmp) )
        call submat( nlayer(i),h6aL(itmp),work(itmp),h6L(itmp) )
        call calh5( nlayer(i),vra(jtmp),ecN(jtmp),ra(isp(i)),work(itmp) )
        call submat( nlayer(i),h6aN(itmp),work(itmp),h6N(itmp) )
        call mtrnp( nlayer(i),h5x(itmp),h3x(itmp) )
        call mtrnp( nlayer(i),h5y(itmp),h3y(itmp) )
        call mtrnp( nlayer(i),h5z(itmp),h3z(itmp) )
        call mtrnp( nlayer(i),h6L(itmp),h4L(itmp) )
        call mtrnp( nlayer(i),h6N(itmp),h4N(itmp) )
        itmp = jsdr+jssp(isl)
        call calhm1( nlayer(i),vra(jtmp),ecKx(jtmp),ra(isp(i)),h5mx(-1,itmp) )
        call calhm1( nlayer(i),vra(jtmp),ecKy(jtmp),ra(isp(i)),h5my(-1,itmp) )
        call calhm1( nlayer(i),vra(jtmp),ecKz(jtmp),ra(isp(i)),h5mz(-1,itmp) )
        call calhm1( nlayer(i),vra(jtmp),ecL(jtmp),ra(isp(i)),h6m1L(-1,itmp) )
        call calhm1( nlayer(i),vra(jtmp),ecN(jtmp),ra(isp(i)),h6m1N(-1,itmp) )
        call calhm2( nlayer(i),vra(jtmp),ecL(jtmp),ra(isp(i)),h6m2L(-2,itmp) )
        call calhm2( nlayer(i),vra(jtmp),ecN(jtmp),ra(isp(i)),h6m2N(-2,itmp) )
        call mtrnp2( nlayer(i),1,2,h5mx(-1,itmp),h3mx(-2,itmp) )
        call mtrnp2( nlayer(i),1,2,h5my(-1,itmp),h3my(-2,itmp) )
        call mtrnp2( nlayer(i),1,2,h5mz(-1,itmp),h3mz(-2,itmp) )
        call mtrnp2( nlayer(i),1,2,h6m1L(-1,itmp),h4m2L(-2,itmp) )
        call mtrnp2( nlayer(i),1,2,h6m1N(-1,itmp),h4m2N(-2,itmp) )
        call mtrnp2( nlayer(i),2,1,h6m2L(-2,itmp),h4m1L(-1,itmp) )
        call mtrnp2( nlayer(i),2,1,h6m2N(-2,itmp),h4m1N(-1,itmp) )
     endif
  enddo

  !******************** plm reading                 *********************

  ! Record the date and time at the beginning of the job
if (SLOW_DEBUG_MODE) then
  write(list, '(I6.6,".",I6.6)') imin,imax
  list = trim(outputDir)//"/log/calLog"//"."//trim(modelname)//"."//trim(list)

  open(25,file =list, status = 'unknown', form = 'formatted')
  call date_and_time(datex,timex)
  write(25,'(/a,a4,a,a2,a,a2,a,a2,a,a2,a,a4)') '    Starting date and time:                     ', &
       datex(1:4),'-',datex(5:6),'-',datex(7:8),'.  ',timex(1:2),':',timex(3:4),':',timex(5:8)
endif

!******************** Computing the coefficients *********************

  llog = 0

if (SLOW_DEBUG_MODE) then
  write(list1, '(I6.6,".",I6.6)') imin,imax
  list1 = trim(outputDir)//"/log/list"//"."//trim(modelname)//"."//trim(list1)
  open(24, file = list1, status = 'unknown', form = 'formatted')
  write(24,*)

  call date_and_time(datex,timex)
  write(25,'(/a,a4,a,a2,a,a2,a,a2,a,a2,a,a4)') '    PLM calculation done:                     ', &
       datex(1:4),'-',datex(5:6),'-',datex(7:8),'.  ',timex(1:2),':',timex(3:4),':',timex(5:8)
endif
  lref=0
  irec=0
  lmax_r(:)=0
  l=0
!********************************************************************************
  do i=imin,imax            ! f-loop start  loop on the frequencies
!********************************************************************************

if (SLOW_DEBUG_MODE) then
     write(24,*) '------------ FRQ -------------'
     write(24,*)  i, dble(i)/tlen,llog-1
     write(24,*)
endif

!! VM VM add action='write' when open output files
!! DK DK now write the displacement and stress coefficients to the same file
!! DK DK to avoid opening too many files on the file system, which is usually slow
     write(coutfile, '(I5.5,".G0_PSV_",I5.5)') int(r0(ir0)*10.d0),i
     coutfile = trim(modelname)//"."//coutfile
     coutfile = trim(outputDir)//"/Stress/"//coutfile
     open(34,file=coutfile, form='unformatted',action='write')

     ir0 = r0_n

     ! initialisation
     tabg0(:,:,:,:,:) = dcmplx(0.d0)
     tabg0der(:,:,:,:,:) = dcmplx(0.d0)

     omega = 2.d0 * pi * dble(i) / tlen

     if ( i /= 0 ) then
        call callsuf(omega,nzone,vrmax,vsv,lsuf)
        call calcoef( nzone,omega,qmu,qkappa,coef1,coef2,coef)
        mtmp = isp(spn) + int(spo)
        if ( spo == int(spo) ) mtmp = mtmp - 1
        call calabnum( omega,omegai,rmax, rrho(1,spn),vpv(1,spn),vph(1,spn),vsv(1,spn),vsh(1,spn), &
        eta(1,spn),ra(mtmp),r0(ir0),coef1(spn),coef2(spn),anum(1,1,1),bnum(1,1,1) )

        if (USE_TIMER) call cpu_time(start_time_to_define_matrix)
        ! computing the matrix elements independent of l
        isl = 0
        ill = 0
        do j=1,ndc+1
           if ( iphase(j) == 1 ) then
              isl = isl + 1
              itmp = isdr+issp(isl)
              jtmp = jdr+jsp(j)
              mtmp = kdr+ksp(j)
              call cala0( nlayer(j),omega,omegai,t(itmp),h1x(itmp), h1y(itmp),h1z(itmp), h2L(itmp), h2N(itmp), &
              h3ax(itmp),h3ay(itmp),h3az(itmp), h4aL(itmp),h4aN(itmp), h5ax(itmp),h5ay(itmp),h5az(itmp), &
              h6aL(itmp),h6aN(itmp), h7x(itmp),h7y(itmp),h7z(itmp), h8L(itmp), h8N(itmp), coef1(j),coef2(j), &
              cwork(jtmp) )
              call overlapa( nlayer(j),cwork(jtmp),a0(1,mtmp))
              call cala1( nlayer(j), h1x(itmp),h1y(itmp),h1z(itmp), h2L(itmp),h2N(itmp), h3x(itmp), h3y(itmp), &
              h3z(itmp),  h4L(itmp), h4N(itmp), h5x(itmp), h5y(itmp), h5z(itmp), h6L(itmp), h6N(itmp),coef1(j), &
              coef2(j),cwork(jtmp) )
              call overlapa( nlayer(j),cwork(jtmp),a1(1,mtmp))
              call cala2( nlayer(j),h1x(itmp), h1y(itmp),h1z(itmp),h2L(itmp),h2N(itmp),coef1(j),coef2(j),cwork(jtmp) )
              call overlapa( nlayer(j), cwork(jtmp),a2(1,mtmp))
              jtmp = jsdr+jssp(isl)
              call calhml( nlayer(j),coef1(j),coef2(j),h3mx(-2,jtmp),h3my(-2,jtmp),h3mz(-2,jtmp), h5mx(-1,jtmp), &
              h5my(-1,jtmp),h5mz(-1,jtmp),h4m1L(-1,jtmp),h4m1N(-1,jtmp), h4m2L(-2,jtmp),h4m2N(-2,jtmp), &
              h6m1L(-1,jtmp),h6m1N(-1,jtmp), h6m2L(-2,jtmp),h6m2N(-2,jtmp), a1(1,mtmp) )
           else
              ill = ill + 1
              itmp = ildr+ilsp(ill)
              jtmp = jdr+jsp(j)
              mtmp = kdr+ksp(j)
              call calb0( nlayer(j),omega,omegai, p1(itmp),p3(itmp),coef(j),cwork(jtmp) )
              call overlapb( nlayer(j), cwork(jtmp),a0(1,mtmp))
              call calb2( nlayer(j),omega,omegai, p2(itmp),coef(j),cwork(jtmp) )
              call overlapb( nlayer(j), cwork(jtmp),a2(1,mtmp))
           endif
        enddo
        if (USE_TIMER) then
           call cpu_time(finish_time_to_define_matrix)
           time_to_define_matrix = time_to_define_matrix + finish_time_to_define_matrix - start_time_to_define_matrix
        endif

        kc = 1
        ismall = 0
        maxamp = -1.d0
        llog = maxlmax
        maxamp_r(:)=-1.d0
        lmax_r(:)=maxlmax
        ismall_r(:)=0
        arret(:)=1
        index_l=0
        lref=0

!! DK DK added this to make sure we do not write the same last buffer twice
        we_left_the_loop_through_exit = .false.

!********************************************************************************************
        do l=0,maxlmax    ! l-loop start
!********************************************************************************************

           index_l=index_l+1  ! index pour stockage des ceoffs
!! DK DK it is not worth trying to vectorize this loop because it contains only "if" statements
           do ir_=1,r_n
             if (ismall_r(ir_) > 20) then
                if (lmax_r(ir_) > l) then
                   lmax_r(ir_)=l
                   arret(ir_)=0
                endif
             endif
           enddo

           if (SUM(arret(:)) == 0) then
!! DK DK added this to make sure we do not write the same last buffer twice
             we_left_the_loop_through_exit = .true.
             exit
           endif

           l2 = dble(l)*dble(l+1)
           lsq = dsqrt( l2 )

           ! computing the coefficient matrix elements
           ! --- renewing  mdr
           if ( mod(l,50) == 0 ) then
              call calmdr( omega,l,nzone,vrmin,vrmax,vmin,dzpar,rmax,sufzone )
              call calspdr(nzone,nzone,iphase,nlayer,jjdr,kkdr )
!! DK DK r_n is of the order of ~500 and thus it is worth trying to vectorize it
!IBM* ASSERT (MINITERCNT(1000))
!DIR$ loop count min(1000)
              do ir_=1,r_n
                ksta(ir_) = kkdr(istazone(ir_))+2*iista(1,ir_) - 1
              enddo
              cksta = kkdr(istazone(cista))+2*iista(1,cista) - 1
              nn = kkdr(nzone) + 2 * nlayer(nzone) + 1
           endif
           if (USE_TIMER) call cpu_time(start_time_to_define_matrix)
           ! computing the matrix elements
           call cala( nzone,ndc,iphase,nlayer,kkdr,kdr,ksp,l2,lsq,nn,a0,a1,a2,a )
           ! computing the boundary condition elements
           call calbc( nzone,ndc,vrmax,iphase,kkdr,a )
           if (USE_TIMER) then
              call cpu_time(finish_time_to_define_matrix)
              time_to_define_matrix = time_to_define_matrix + finish_time_to_define_matrix - start_time_to_define_matrix
           endif
           jtmp = kkdr(spn) + 2 * int(spo)
           mtmp = isp(spn) + int(spo)
           if ( spo == int(spo) ) then
              jtmp = jtmp - 2
              mtmp = mtmp - 1
           endif

           call calya(anum(1,1,1),bnum(1,1,1),l2,ra(mtmp),r0(ir0),ya,yb,yc,yd)

           do m=-2,2        ! m-loop start

              if ( iabs(m) <= iabs(l) ) then
                 ig2 = 0

!! DK DK this part is negligible in terms of cost because only done for l == 0; thus no need to optimize it
                 if ( l == 0 ) then ! l-branch for calu (l=0)

                    !  rearranging the matrix elements
                    do imt = 1,6

                       call setmt(imt,mt)
                       g0 = cmplx(0.d0)
                       call calg(l,m,coef1(spn),coef2(spn),lsq,ecC0,ecF0,ecL0,ya,yb,yc,yd,ra(mtmp),r0(ir0),mt,g0(jtmp))
                       call rea2(nn,a,g0,c,d0,nzone,iphase,kkdr,spn,kkdr0,nn0,r_n,r_n,istazone,iista,jsta)

                       itmp=1
                       if ( rmin == 0.d0 ) itmp=2

                       ns = kkdr0 + ( nint(spo) - 1 )
                       if (USE_TIMER) call cpu_time(start_time_to_dcsymbdl)
                       call dcsymbdl0_m_equals_1_nn_equals_1(c(1,itmp),nn0-itmp+1,eps,z(itmp),w(itmp),ll,lli,llj,ier)
                       if (USE_TIMER) then
                          call cpu_time( finish_time_to_dcsymbdl)
                          time_to_dcsymbdl = time_to_dcsymbdl +  finish_time_to_dcsymbdl - start_time_to_dcsymbdl
                       endif
                       if ( ((abs(m) == 0) .and. ((imt == 1) .or. (imt == 2) .or. (imt == 3))) .or. &
                           ((abs(m) == 1) .and. ((imt == 4) .or. (imt == 5))) .or. &
                           ((abs(m) == 2) .and. ((imt == 2) .or. (imt == 3) .or. (imt == 6))) ) then
                          if (USE_TIMER) call cpu_time(start_time_to_dcsbdlv)
                          call dcsbdlv0_m_equals_1(c(1,itmp),d0(itmp),nn0-itmp+1,z(itmp))
                          if (USE_TIMER) then
                             call cpu_time( finish_time_to_dcsbdlv)
                             time_to_dcsbdlv = time_to_dcsbdlv +  finish_time_to_dcsbdlv - start_time_to_dcsbdlv
                          endif
                       endif

!! DK DK this loop cannot be vectorized because the call to "interpolate" that it contains cannot be inlined
                       do ir_=1,r_n
                          g0tmp = (0.d0,0.d0)
                          g0dertmp = (0.d0,0.d0)
                          call interpolate(1,0,r_(ir_),rrsta(1,ir_),d0(jsta(ir_)),g0tmp(1))
                          call interpolate(1,1,r_(ir_),rrsta(1,ir_),d0(jsta(ir_)),g0dertmp(1))
                          tabg0(index_l,1,imt,m,ir_)=g0tmp(1)
                          tabg0(index_l,2,imt,m,ir_)=g0tmp(2)
                          tabg0der(index_l,1,imt,m,ir_)=g0dertmp(1)
                          tabg0der(index_l,2,imt,m,ir_)=g0dertmp(2)
                       enddo
                    enddo ! imt-loop


!******************************************************************************************************************
!! DK DK this is the expensive part to optimize, done for all values of l except the first one
!******************************************************************************************************************
                 else ! for l /= 0

                    do imt = 1,6
                       call setmt(imt,mt)
                       g0 = cmplx(0.d0)
                       call calg(l,m,coef1(spn),coef2(spn),lsq,ecC0,ecF0,ecL0,ya,yb,yc,yd,ra(mtmp),r0(ir0),mt,g0(jtmp))

                       ! computing forward propagating component (l != 0)
                       itmp=1
                       if ( rmin == 0.d0 ) itmp=3
                       ns = kkdr(spn) + 2 * ( nint(spo) - 1 )
                       if ( ( m == -2 ) .or. ( m == -l ) ) then
                          if (ig2 == 0) then
!******************************************************************************************************************
                             if (USE_TIMER) call cpu_time(start_time_to_dcsymbdl)
                             call dcsymbdl0_m_equals_3_nn_equals_6(a(1,itmp),nn-itmp+1,eps,z(itmp),w(itmp),ier) ! expensive routine
                             if (USE_TIMER) then
                                call cpu_time( finish_time_to_dcsymbdl)
                                time_to_dcsymbdl = time_to_dcsymbdl +  finish_time_to_dcsymbdl - start_time_to_dcsymbdl
                             endif
!******************************************************************************************************************
                             ig2 = 1
                          endif
                       endif
                       if ( ((abs(m) == 0) .and. ((imt == 1) .or. (imt == 2) .or. (imt == 3))) .or. &
                           ((abs(m) == 1) .and. ((imt == 4) .or. (imt == 5))) .or. &
                           ((abs(m) == 2) .and. ((imt == 2) .or. (imt == 3) .or. (imt == 6))) ) then
!******************************************************************************************************************
                          if (USE_TIMER) call cpu_time(start_time_to_dcsbdlv)
                          call dcsbdlv0_m_equals_3(a(1,itmp),g0(itmp),nn-itmp+1,z(itmp)) ! expensive routine
                          if (USE_TIMER) then
                             call cpu_time( finish_time_to_dcsbdlv)
                             time_to_dcsbdlv = time_to_dcsbdlv +  finish_time_to_dcsbdlv - start_time_to_dcsbdlv
                          endif
!******************************************************************************************************************
                       endif

                       ! test sur l'amplitude
                       if (imt == 1 .and. ir0 == r0_n) then
!! DK DK r_n is of the order of ~500 and thus it is worth trying to vectorize it
!IBM* ASSERT (MINITERCNT(1000))
!DIR$ loop count min(1000)
                          do ir_ = 1,r_n
!                           call calamp(g0(ksta(ir_)-1),l,lsuf,maxamp_r(ir_),ismall_r(ir_),ratl)
!! DK DK inlined the call in order to vectorize the loop
                            ampratio = 0.d0
                            amp = dsqrt(zabs(g0(ksta(ir_)-1))**2 + zabs(g0(ksta(ir_)))**2)
                            if (amp > maxamp_r(ir_)) maxamp_r(ir_) = amp
                            if (amp /= 0.d0 .and. maxamp_r(ir_) /= 0.d0) ampratio = amp / maxamp_r(ir_)
                            if (ampratio < ratl .and. l >= lsuf) then
                              ismall_r(ir_) = ismall_r(ir_) + 1
                            else
                              ismall_r(ir_) = 0
                            endif
                          enddo
                          call calamp(g0(ksta(1)-1),l,lsuf,maxamp,ismall,ratl)
                       endif

!! DK DK this loop cannot be vectorized because the call to "interpolate" that it contains cannot be inlined
                       do ir_=1,r_n ! stack point
                          g0tmp = (0.d0,0.d0)
                          g0dertmp = (0.d0,0.d0)
                          call interpolate(2,0,r_(ir_),rrsta(1,ir_),g0(ksta(ir_)-1),g0tmp(1))
                          call interpolate(2,1,r_(ir_),rrsta(1,ir_),g0(ksta(ir_)-1),g0dertmp(1))
                          tabg0(index_l,1,imt,m,ir_)=g0tmp(1)
                          tabg0(index_l,2,imt,m,ir_)=g0tmp(2)
                          tabg0der(index_l,1,imt,m,ir_)=g0dertmp(1)
                          tabg0der(index_l,2,imt,m,ir_)=g0dertmp(2)
                       enddo   ! stack point
                    enddo ! mt-loop

                 endif   ! l-branch for calu
              endif
           enddo            ! m-loop end

           ! on ecrit ici les fichiers des coefficients
!! From Vadim Monteiller, April 2013:
! je n'ecris pas tous les ir_ en une seule fois car le lmax depend des ir_.
! Au cours de la boucle sur l, au debut, il est vrai que l'on doit ecrire le tableau tabg0 pour tous les ir_,
! mais au fur et a mesure que l grandit, on perd des  ir_ (i.e. les coefs sont consideres comme nuls),
! donc j'evite d'ecrire des zeros, et pour l'instant la seule solution que j'ai trouvee c'est celle que j'ai implementee.
           if (index_l == maxlmax_g) then
if (SLOW_DEBUG_MODE) then
              write(24,*) 'partial writing of the coefficients'
endif
              do ir_=1,r_n
                 llog0=min(maxlmax_g, lmax_r(ir_) - lref)
if (SLOW_DEBUG_MODE) then
                 write(24,*) ir_,llog0,l
endif
                 if (USE_TIMER) call cpu_time(start_time_to_write)
                 if (llog0 > 0) call write_a_block_of_coefs_to_disk(tabg0,tabg0der,tabg0_small_buffer,ir_,llog0,maxlmax_g,r_n)

                 if (USE_TIMER) then
                    call cpu_time(finish_time_to_write)
                    time_to_write = time_to_write + finish_time_to_write - start_time_to_write
                 endif
              enddo

              ! write a separator between blocks and re-initialise arrays to zero
              tabg0=dcmplx(0.d0)
              tabg0der=dcmplx(0.d0)
              write(34) -1,1
              !write(34) tabg0(1,:,:,:,1) !! VM VM this is not useful any more
              !write(34) tabg0(1,:,:,:,1) !! VM VM this is not useful any more

              ! reset the index, for the next block to write later
              index_l=0
              lref=l
           endif

!******************************************************************************
        enddo               ! l-loop end
!******************************************************************************

     endif

if (SLOW_DEBUG_MODE) then
        call date_and_time(datex,timex)
        write(25,'(/a,i5,a,a4,a,a2,a,a2,a,a2,a,a2,a,a4)') '    Frequency-index ', i, ' :', &
             datex(1:4),'-',datex(5:6),'-',datex(7:8),'.  ',timex(1:2),':',timex(3:4),':',timex(5:8)
endif

!! DK DK to Vadim: je crois qu'il y a possibilite de petit bug ici: si jamais maxlmax est un multiple de maxlmax_g
!! (ce qui est le cas actuellement je crois) et que tu n'es pas sorti de la boucle sur l ci-dessus par
!! le test "if (SUM(arret(:)) == 0) exit" mais que la boucle s'est poursuivie jusqu'au bout sans arret,
!! alors je crois que ce que tu ecris en complement de datas ci-dessous aura deja ete ecrit une fois ci-dessous
!! a la fin de la boucle, et apparaitra donc deux fois. Solution: soit imposer que maxlmax ne soit jamais un
!! multiple de maxlmax_g (et s'il l'est, lui ajouter 1 au debut du code), soit tester avec un flag logique
!! si tu es sorti de la boucle par le exit ou non, et si non et si en plus maxlmax est un multiple de maxlmax_g
!! alors il ne faudrait pas faire l'ecriture ci-dessous? Stp verifie.
!********************************************************************************************
!! DK DK added this to make sure we do not write the same last buffer twice.
!! DK DK Basically what we need to test here is that we did not already write that same buffer
!! DK DK to disk at the end of the loop above.
!! DK DK Vadim writes the end of the last buffer here.
!! DK DK but we should do this ONLY if it has not been written before.
   if (mod(maxlmax,maxlmax_g) /= 0 .or. we_left_the_loop_through_exit) then

     if (SLOW_DEBUG_MODE) write(24,*) 'writing last block'

     do ir_=1,r_n
        llog0=min(maxlmax_g,lmax_r(ir_)-lref)
if (SLOW_DEBUG_MODE) then
        write(24,*) ir_,llog0,l
endif
        if (USE_TIMER) call cpu_time(start_time_to_write)
        if (llog0 > 0) call write_a_block_of_coefs_to_disk(tabg0,tabg0der,tabg0_small_buffer,ir_,llog0,maxlmax_g,r_n)
        if (USE_TIMER) then
           call cpu_time(finish_time_to_write)
           time_to_write = time_to_write + finish_time_to_write - start_time_to_write
        endif
     enddo

if (SLOW_DEBUG_MODE) then
     write(24,*) 'end of writing for freq ', i
endif

   endif

     ! reset the array
     tabg0=dcmplx(0.d0)
     tabg0der=dcmplx(0.d0)
     ! close the file (write dummy negative values to indicate the end of the data)
     write(34) -1,1
     write(34) -1,-1 !! VM VM add this line because of bugs occurs sometimes when
                     !! VM VM mod(maxlmax_g,lmax)==0
     !write(34) tabg0(1,:,:,:,1) !! VM VM this is not useful any more
     !write(34) tabg0(1,:,:,:,1) !! VM VM this is not useful any more
     close(34)

!********************************************************************************
  enddo     ! of frequency loop
!********************************************************************************

if (SLOW_DEBUG_MODE) then
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
     time_to_write = time_to_write/3600.
     time_to_dcsymbdl = time_to_dcsymbdl/3600.
     time_to_dcsbdlv = time_to_dcsbdlv/3600.
     time_to_define_matrix = time_to_define_matrix/3600.
     !
     call mpi_allreduce( time_total, total_global_time,1, MPI_REAL ,MPI_SUM, MPI_COMM_WORLD, ierr)
     time_total= total_global_time
     call mpi_allreduce( time_to_write, total_global_time,1, MPI_REAL ,MPI_SUM, MPI_COMM_WORLD, ierr)
     time_to_write= total_global_time
     call mpi_allreduce( time_to_dcsymbdl, total_global_time,1, MPI_REAL ,MPI_SUM, MPI_COMM_WORLD, ierr)
     time_to_dcsymbdl= total_global_time
     call mpi_allreduce( time_to_dcsbdlv, total_global_time,1, MPI_REAL ,MPI_SUM, MPI_COMM_WORLD, ierr)
     time_to_dcsbdlv= total_global_time
     call mpi_allreduce( time_to_define_matrix, total_global_time,1, MPI_REAL ,MPI_SUM, MPI_COMM_WORLD, ierr)
     time_to_define_matrix= total_global_time
     if (myrank == 0) then
        open(25,file='timer_part2.txt')
        write(25,*) 'Total :',time_total
        write(25,*) 'Matrix definition :',time_to_define_matrix,100*time_to_define_matrix/time_total,' %'
        write(25,*) 'Matrix LU decp :', time_to_dcsymbdl,100* time_to_dcsymbdl/time_total
        write(25,*) 'System Solver  :', time_to_dcsbdlv,100*time_to_dcsbdlv/time_total
        write(25,*) ' Write Exp Coef:',time_to_write,100*time_to_write/time_total
        close(25)
     endif
  endif
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call MPI_FINALIZE(ierr)

end program TraPSV

