program TraSH_Write_coefficients


!-----------------------------------------------------------------------
!
!
!
!
!
!                                               2002.10.KAWAI Kenji
!                                               2009.6. FUJI Nobuaki
!                                               2011.9. FUJI Nobuaki
!                                               2013.11 YI WANG (MPI)
!
!
!
!
!-----------------------------------------------------------------------

  implicit none

! MPI -- YW

include 'mpif.h'
!include "./constants.h"

INTEGER myrank,nbproc,ierr

!------------------------- <  < input matrix >>----------------------------------

  character(120) :: inputdir,outputDir,psvmodel,modelname,stationsinf
  character(120) :: list,list1
  character(40) :: datex,timex
  real(kind(0d0)), parameter :: pi=3.1415926535897932d0
  real(kind(0d0)), parameter :: re=1.d-2, ratc=1.d-10, ratl=1.d-4
  integer, parameter :: maxlmax = 25000, maxlmax_g = 1000, nsta = 1, theta_n=nsta
  !change to 30000 for possible high frequency calculation.

  real(kind(0d0)) :: tlen
  real(kind(0d0)) :: r0min, r0max, r0delta  !!! JUST FOR ONE DEPTH FOR THIS MOMENT !!
  real(kind(0d0)) :: r0lat, r0lon,stla_curr,stlo_curr
  real(kind(0d0)), allocatable :: stla(:),stlo(:),stla_g(:),stlo_g(:)
  real(kind(0d0)), allocatable :: r_(:),r0(:),theta(:),phi(:),phi_g(:),theta_g(:)
  real(kind(0d0)), allocatable :: A0sta(:),C0sta(:),F0sta(:),L0sta(:),N0sta(:)
  integer, allocatable :: updown(:),idum(:)
  real(kind(0d0)), allocatable :: rrsta(:,:)
  integer, allocatable :: iista(:,:)
  integer :: r_n,r0_n, ir_,ir0,imt,itheta,nsta_ignored_in_this_code

  character(120) :: coutfile, coutfile1, coutfile2
  integer :: imin,imax
  integer :: itranslat
  integer :: nzone
  integer :: iimax

  integer :: i ,j, ier
  real(kind(0d0)) :: dummy
  real(kind(0d0)), allocatable :: vrmin(:), vrmax(:)
  real(kind(0d0)), allocatable :: rrho(:,:), vsv(:,:), vsh(:,:), qmu(:)
  real(kind(0d0)), allocatable :: vra (:), rho(:), ecL(:), ecN(:)
  real(kind(0d0)), allocatable :: gvra(:,:), grho(:,:), gecL(:,:), gecN(:,:),gra(:,:)
  complex(kind(0d0)), allocatable :: coef(:),cwork(:)
  real(kind(0d0)) :: rmin, rmax
  real(kind(0d0)), allocatable :: vmin(:), gridpar(:), dzpar(:)
  complex(kind(0d0)), allocatable :: tmpc(:)
  real(kind(0d0)) :: maxamp
  real(kind(0d0)) :: omegai
  real(kind(0d0)), allocatable :: ra(:)
  integer :: nnlayer, vnp,nn
  integer, allocatable :: nlayer(:), iphase(:)
  integer :: ioutercore

  ! variables pour des points stackes
  integer, allocatable :: isp(:),jsp(:),ins(:)

  ! variables pour la source

  integer, allocatable :: spn(:),ns(:)
  real(kind(0d0)) :: mt(3,3),lsq
  real(kind(0d0)), allocatable :: mu0(:),spo(:)

!-----------------------------------------------------------------------
  ! variables pour des elements de matrice
  complex(kind(0d0)), allocatable :: a0(:,:), a2(:,:), a(:,:),dr(:),z(:)
  real(kind(0d0)), allocatable :: t(:), h1(:), h2(:), h3(:), h4(:), work(:),w(:)
  real(kind(0d0)), allocatable :: gt(:,:),gh1(:,:),gh2(:,:),gh3(:,:),gh4(:,:)
  complex(kind(0d0)),allocatable :: aa(:,:), ga(:,:),ga2(:,:,:),gdr(:,:)
  complex(kind(0d0)), allocatable :: g0(:)
  complex(kind(0d0)) :: g0tmp, g0dertmp
  ! la frequence
  real(kind(0d0)) :: omega
  integer :: lsuf

  ! des autres
  integer :: lda
  integer :: kc, ismall,ismall0, llog,m, l,ig2,ll(12),lli(12),llj(12)
  real(kind(0d0)) :: eps

!-----------------------------------------------------------------------
  complex(kind(0d0)), allocatable :: bvec(:,:,:),bvecdt(:,:,:),bvecdp(:,:,:)
  real(kind(0d0)), allocatable :: plm(:,:,:),maxamp_r(:)
  complex(kind(0d0)), allocatable :: stress(:,:,:), displacement(:,:,:)
  complex(kind(0e0)), allocatable :: stresssngl(:,:,:), displacementsngl(:,:,:)
  complex(kind(0e0)), allocatable :: stresssngl_global(:,:,:), displacementsngl_global(:,:,:)
  complex(kind(0d0))::u(1:3),udr(1:3),udt(1:3),udp(1:3),uder(1:3,1:3)

  data lda / 2 /
  data eps / -1.d0 /
!-----------------------------------------------------------------------
 ! variables for the MPI communications.
  integer, allocatable :: Ifrq(:),Ifrq2(:,:)
  integer imyrank
  integer imin_glob, imax_glob, imax_global
  integer para
!-----------------------------------------------------------------------
! other variables
complex(kind(0d0)),allocatable:: tabg0(:,:,:,:,:),tabg0der(:,:,:,:,:),tabg0_small_buffer(:,:,:,:)
integer,allocatable:: lmax_r(:),ismall_r(:),arret(:)
integer:: llog0,index_l,lref,irank,nsta_global,ifrequ_min,ifrequ_max
real(kind(0d0)):: amp,ampratio
logical :: we_left_the_loop_through_exit

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nbproc,ierr)

if (myrank == 0)  write(*,*) 'how many processor do i successfully apply:',nbproc

 if (myrank == 0) then
    call pinputTra_Part1(outputDir,psvmodel,modelname,stationsinf,tlen,imin_glob,imax_glob,r0min,r0max,r0delta,r0lat,r0lon,itranslat,myrank)
   if (nbproc == 0) then
        imin=imin_glob
        imax=imax_glob
   endif
   !! VM VM commented that because the file FrqMpi.txt Should be writen before
   !call WriteFrqByproc(nbproc,imin_glob,imax_glob)  !distributing the task into every processors evenly.
   !! VM VM
 endif !only the root can excute this function call.



 allocate(Ifrq(nbproc+1))
 allocate(Ifrq2(0:nbproc-1,2))
 if (myrank == 0)   call ReadFrqsByProc(Ifrq,Ifrq2,para,nbproc+1)  !frequency calculation separation.
     call MPI_Bcast(Ifrq,nbproc+1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(Ifrq2,2*nbproc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(para,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
! write(*,*) 'initial setting is completed'

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
  call MPI_Bcast(itranslat,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  itranslat=0

 call MPI_Barrier(MPI_COMM_WORLD,ierr)

 ! frequency section making now
   if (para == 1) then
     imin=Ifrq(myrank+1)+ 1
     imax=Ifrq(myrank+2)
   if (myrank == 0) imin=Ifrq(myrank+1) ! pay attension for this judgement
     imax_global = maxval(Ifrq)
   else
     imin=Ifrq2(myrank,1)
     imax=Ifrq2(myrank,2)
   endif
!   write(*,*) 'the',myrank,'th processor, frequency section is:',imin,imax
  imax_global = imax_glob


if (myrank == 0) then
  open(20, file = psvmodel, status = 'old', action='read', position='rewind')
  read(20,*) nzone
!  call MPI_Bcast(nzone,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  allocate(vrmin(1:nzone))
  allocate(vrmax(1:nzone))
  allocate(rrho(1:4,1:nzone))
  allocate(vsv(1:4,1:nzone))
  allocate(vsh(1:4,1:nzone))
  allocate(qmu(1:nzone))
  allocate(vmin(1:nzone))
  allocate(gridpar(1:nzone))
  allocate(dzpar(1:nzone))
  allocate(nlayer(1:nzone))
  allocate(iphase(1:nzone))
  allocate(isp(1:nzone))
  allocate(jsp(1:nzone))
  allocate(coef(1:nzone))

   do i = 1, nzone
     read (20, *) vrmin(i), vrmax(i), rrho(1,i), rrho(2,i), rrho(3,i), rrho(4,i), dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, vsv(1,i), vsv(2,i), vsv(3,i), vsv(4,i), vsh(1,i), vsh(2,i), vsh(3,i), vsh(4,i), dummy, dummy, dummy, dummy, qmu(i), dummy
     if ((vsv(1,i) == 0.d0) .and. (vsv(2,i) == 0.d0) .and. (vsv(3,i) == 0.d0) .and. (vsv(4,i) == 0.d0)) then
        iphase(i) = 2
  ioutercore = i
     else
        iphase(i) = 1
     endif
   enddo
   close(20)
! CAUTION: this program can only calculate for solid media (SH) for this moment

   open(20, file = psvmodel, status = 'old', action='read', position='rewind')
   read(20,*) nzone
   nzone = nzone - ioutercore
   deallocate(vrmin,vrmax,rrho,vsv,vsh,qmu,vmin,gridpar,dzpar,nlayer,iphase,isp,jsp,coef)
endif ! only root can read file for inputing the setting parameters.

 call MPI_Bcast(nzone,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) !nzone changed
 call MPI_Bcast(ioutercore,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

 call MPI_Barrier(MPI_COMM_WORLD,ierr)

  allocate(vrmin(1:nzone))
  allocate(vrmax(1:nzone))
  allocate(rrho(1:4,1:nzone))
  allocate(vsv(1:4,1:nzone))
  allocate(vsh(1:4,1:nzone))
  allocate(qmu(1:nzone))
  allocate(vmin(1:nzone))
  allocate(gridpar(1:nzone))
  allocate(dzpar(1:nzone))
  allocate(nlayer(1:nzone))
  allocate(iphase(1:nzone))
  allocate(isp(1:nzone))
  allocate(jsp(1:nzone))
  allocate(coef(1:nzone))

 if (myrank == 0) then
   do i = 1,ioutercore
      read (20, *) dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy
   enddo
   do i = 1, nzone
      read (20, *) vrmin(i), vrmax(i), rrho(1,i), rrho(2,i), rrho(3,i), rrho(4,i), dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, vsv(1,i), vsv(2,i), vsv(3,i), vsv(4,i), vsh(1,i), vsh(2,i), vsh(3,i), vsh(4,i), dummy, dummy, dummy, dummy, qmu(i), dummy
   enddo
   close(20)
 endif

 call MPI_Bcast(vrmin,nzone,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
 call MPI_Bcast(vrmax,nzone,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
 call MPI_Bcast(rrho,4*nzone,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
 call MPI_Bcast(vsv,4*nzone,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
 call MPI_Bcast(vsh,4*nzone,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
 call MPI_Bcast(qmu,nzone,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
 call MPI_Bcast(iphase,nzone,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

 call MPI_Barrier(MPI_COMM_WORLD,ierr)
   rmin = vrmin(1)
   rmax = vrmax(nzone)
   omegai = - dlog(1.d-2) / tlen

 if (myrank == 0) then
   open (1,file=stationsinf,status='old',action='read',position='rewind')
   if (itranslat == 1) call translat (r0lat,r0lat)
   open(2,file='recdepth')
   read(2,*) r_n
   read(1,*) nsta_ignored_in_this_code
 endif

 call MPI_Bcast(r_n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
 call MPI_Bcast(r0lat,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
 call MPI_Bcast(r0lon,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
 ! if itraslat = 1.the r0lat and r0lon will be changed.
 call MPI_Barrier(MPI_COMM_WORLD,ierr)

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
  allocate(idum(1:r_n))
  allocate(updown(1:r_n))
!  allocate(stress(1:6,1:6,1:nsta))
!  allocate(displacement(1:3,1:6,1:nsta))
!  allocate(stresssngl(1:6,1:6,1:nsta))
!  allocate(displacementsngl(1:3,1:6,1:nsta))
  allocate(lmax_r(r_n),ismall_r(r_n),arret(r_n))
  allocate(maxamp_r(r_n))
!! DK DK moved index 1..maxlmax_g to the first position in order to be able to vectorize a loop
  allocate(tabg0(maxlmax_g,1,6,-2:2,r_n))
  allocate(tabg0der(maxlmax_g,1,6,-2:2,r_n))
  allocate(tabg0_small_buffer(maxlmax_g,1,6,-2:2))
!! for now is the computation of SH parts, So the tabg0 is the half of the PSV parts ,because only one
!! kind of coefficients need to be stored.


 if (myrank == 0) then
   do i = 1,nsta
   read(1,*) r_(i),stla(i),stlo(i),updown(i)
     r_(i) = 6371.d0 -r_(i)
     if (itranslat == 1) call translat(stla(i),stla(i))
     call calthetaphi(r0lat,r0lon,stla(i),stlo(i),theta(i),phi(i))
   enddo
   close(1)
 endif
!! in fact ,this part is unused in the writing coefficients routine.
  if (myrank == 0) then
    do i=1,r_n
        read(2,*) r_(i),idum(i),updown(i)
  r_(i) = 6371.d0 -r_(i)
    enddo
    close(1)
    close(2)
  endif


 call MPI_Barrier(MPI_COMM_WORLD,ierr)

      call MPI_Bcast(r_,r_n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      !! please notice the size of the r_(:) array.!!!
      call MPI_Bcast(updown,nsta,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(stla,nsta,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(stlo,nsta,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(theta,nsta,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(phi,nsta,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(updown,r_n,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(idum,r_n,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
! because the 'recdepth' file has 3 colum,so we have to include the index reading for every layer.
! idum is the array to store these index number.

 do i = 1,r_n
!!   call calthetaphi(r0lat,r0lon,stla(i),stlo(i),theta(i),phi(i))
    call calstg4onedepth(nzone,nzone,nzone,vrmin,vrmax,iphase,rrho,0.d0,0.d0,vsv,vsh,1.d0,rmax,r_(i),updown(i),A0sta(i),C0sta(i),F0sta(i),L0sta(i),N0sta(i))
 enddo  !pay attention for the MPI_BCAST.because it's easy to forget to declaire some element by BCAST.

 write(list, '(I6,".",I6)'), imin_glob,imax_glob
 do j = 1,15
    if (list(j:j) == ' ')list(j:j) = '0'
 enddo
 list = trim(outputDir)//"/log/calLogSH"//"."//trim(modelname)//"."//trim(list)

 if (myrank == 0) then
    open(1,file =list, status = 'unknown', form = 'formatted')
    call date_and_time(datex,timex)
    write(1,'(/a,a4,a,a2,a,a2,a,a2,a,a2,a,a4)') &
    '    Starting date and time:                     ', &
    datex(1:4),'-',datex(5:6),'-',datex(7:8),'.  ', &
    timex(1:2),':',timex(3:4),':',timex(5:8)
  close (1)
 endif

  allocate(rrsta(1:3,1:r_n))
  allocate(iista(1:3,1:r_n))

  allocate(bvec(1:3,-2:2,1:theta_n))
  allocate(bvecdt(1:3,-2:2,1:theta_n))
  allocate(bvecdp(1:3,-2:2,1:theta_n))
  allocate(plm(1:3,0:3,1:theta_n))

  ! source depths

    r0_n =  int((r0max-r0min)/r0delta)+1

      allocate(r0(1:r0_n))
      allocate(spo(1:r0_n))
      allocate(spn(1:r0_n))
      allocate(ns(1:r0_n))
      allocate(mu0(1:r0_n))
      allocate(ins(1:r0_n))
      allocate(gra(1:3,1:r0_n))
      allocate(gvra(1:3,1:r0_n))
      allocate(grho(1:3,1:r0_n))
      allocate(gecL(1:3,1:r0_n))
      allocate(gecN(1:3,1:r0_n))
      allocate(gt(1:8,1:r0_n))
      allocate(gh1(1:8,1:r0_n))
      allocate(gh2(1:8,1:r0_n))
      allocate(gh3(1:8,1:r0_n))
      allocate(gh4(1:8,1:r0_n))
      allocate(aa(1:4,1:r0_n))
      allocate(ga(1:8,1:r0_n))
      allocate(ga2(1:2,1:3,1:r0_n))
      allocate(gdr(1:3,r0_n))

  do i = 1, r0_n
     r0(i) = r0min + dble(i-1)*r0delta
  enddo
  ir0 = r0_n
  if ( (r0(ir0) < rmin) .or. (r0(ir0) > rmax) ) stop 'Location of the source is improper.'

 ! computation de nombre et la location des points de grid

 iimax =imax_glob  !!!this is very important for computing multi-station !!!!!!!

   call calgrid( nzone,vrmin,vrmax,vsv,rmin,rmax,iimax,1,tlen,vmin,gridpar,dzpar )
!   if (myrank ==0) write(*,*) 'rmin,rmax,iimax',rmin,rmax,iimax
!   if (myrank ==0) write(*,*) '&&dzpar',dzpar
!   if (myrank ==0) write(*,*) '&&gridpar',gridpar
   call calra(nnlayer,gridpar,dzpar,nzone,vrmin,vrmax,rmin,rmax,nlayer,re )
!   if (myrank ==0) write(*,*) '!!!nlayer',nlayer
   allocate(ra(1:nnlayer+nzone+1))
!   if (myrank ==0) write(*,*) '@@ size of ra s @@',size(ra),'layer numbers',nnlayer
   call calra2(nnlayer,gridpar,dzpar,nzone,vrmin,vrmax,rmin,rmax,nlayer,ra,re,r_n,r_,rrsta, iista,updown)
   ! computation de points stackes et la location de la source
   !!! The every upper boundary should be marked with the updown(:) = -1 value
   call calsp( nzone-1,nlayer,isp,jsp )
   do ir0 = 1,r0_n
      call calspo( nzone-1,vrmax,nnlayer,r0(ir0),rmin,rmax,ra,isp,spo(ir0),spn(ir0) )
      call calgra( isp,ra,r0(ir0),spn(ir0),spo(ir0),gra(1:3,ir0))
   enddo

!!!!test
!   if (myrank == nbproc -1) then
!     do i=1,r_n,1
!      write(*,*) 'index',updown(i),i,'th dh:',(r_(i)-rrsta(1,i)),(r_(i)-rrsta(2,i)),(r_(i)-rrsta(3,i))
!     enddo
!   endif

  ! computation des elements de matrice

  allocate(vra(1:nnlayer+2*nzone+1))
  allocate(rho(1:nnlayer+2*nzone+1))
  allocate(ecL(1:nnlayer+2*nzone+1))
  allocate(ecN(1:nnlayer+2*nzone+1))
  allocate(a0(1:2,1:nnlayer+1))
  allocate(a2(1:2,1:nnlayer+1))
  allocate(a(1:2,1:nnlayer+1))
  allocate(t(1:4*nnlayer))
  allocate(cwork(1:4*nnlayer))
  allocate(h1(1:4*nnlayer))
  allocate(h2(1:4*nnlayer))
  allocate(h3(1:4*nnlayer))
  allocate(h4(1:4*nnlayer))
  allocate(work(1:4*nnlayer))
  allocate(tmpc(nnlayer+1))
  allocate(g0(1:nnlayer+1))
  allocate(dr(1:nnlayer+1))
  allocate(w(1:nnlayer+1))
  allocate(z(1:nnlayer+1),stat=ier)
  if (ier /= 0) stop '<!!!error: not enough memory to allocate array > !!!!'
  call calstg( nzone,rrho,vsv,vsh,nnlayer,nlayer,ra,rmax,vnp,vra,rho,ecL,ecN)
  do ir0 = 1, r0_n
     call calgstg(nzone,nnlayer,spn(ir0),rrho,vsv,vsh,gra(1:3,ir0),gvra(1:3,ir0),rmax,grho(1:3,ir0),gecL(1:3,ir0),gecN(1:3,ir0),r0(ir0),mu0(ir0))
  enddo

  do i= 1, nzone
     call calmatc(nlayer(i),vnp,vra,rho,2,0,0,ra(isp(i)),t(jsp(i)),work(jsp(i)))
     call calmatc(nlayer(i),vnp,vra,ecL,2,1,1,ra(isp(i)),h1(jsp(i)),work(jsp(i)))
     call calmatc(nlayer(i),vnp,vra,ecL,1,1,0,ra(isp(i)),h2(jsp(i)),work(jsp(i)))
     call calmatc(nlayer(i),vnp,vra,ecL,0,0,0,ra(isp(i)),h3(jsp(i)),work(jsp(i)))
     call calmatc(nlayer(i),vnp,vra,ecN,0,0,0,ra(isp(i)),h4(jsp(i)),work(jsp(i)))
     call caltl(nlayer(i),vnp,vra,rho,ra(isp(i)),work(jsp(i)))
     call calt(nlayer(i),t(jsp(i)),work(jsp(i)),t(jsp(i)))
     call calhl(nlayer(i),vnp,vra,ecL,ra(isp(i)),work(jsp(i)))
     call calt(nlayer(i),h3(jsp(i)),work(jsp(i)),h3(jsp(i)))
     call calhl(nlayer(i),vnp,vra,ecN,ra(isp(i)),work(jsp(i)))
     call calt(nlayer(i),h4(jsp(i)),work(jsp(i)),h4(jsp(i)))
  enddo
  do ir0 = 1, r0_n
     call calmatc( 2,3,gvra(1:3,ir0),grho(1:3,ir0),2,0,0,gra(1:3,ir0),gt(1:8,ir0), work )
     call calmatc( 2,3,gvra(1:3,ir0),gecL(1:3,ir0),2,1,1,gra(1:3,ir0),gh1(1:8,ir0),work )
     call calmatc( 2,3,gvra(1:3,ir0),gecL(1:3,ir0),1,1,0,gra(1:3,ir0),gh2(1:8,ir0),work )
     call calmatc( 2,3,gvra(1:3,ir0),gecL(1:3,ir0),0,0,0,gra(1:3,ir0),gh3(1:8,ir0),work )
     call calmatc( 2,3,gvra(1:3,ir0),gecN(1:3,ir0),0,0,0,gra(1:3,ir0),gh4(1:8,ir0),work )
     call caltl(2,3,gvra(1:3,ir0),grho(1:3,ir0),gra(1:3,ir0),work )
     call calt(2,gt(1:8,ir0),work,gt(1:8,ir0))
     call calhl(2,3,gvra(1:3,ir0),gecL(1:3,ir0),gra(1:3,ir0),work)
     call calt( 2,gh3(1:8,ir0),work,gh3(1:8,ir0))
     call calhl(2,3,gvra(1:3,ir0),gecN(1:3,ir0),gra(1:3,ir0),work )
     call calt( 2,gh4(1:8,ir0),work,gh4(1:8,ir0))
  enddo


  !computation de la dislocation
  nn = nnlayer + 1
  do ir0 = 1, r0_n
     ns(ir0) = isp(spn(ir0)) + dint(spo(ir0))
     ins(ir0) = 4 * ns(ir0) - 3
  enddo
  llog = 0


  write(list1, '(I6,".",I6)'), imin_glob,imax_glob  !!unify scanning file for the process
  do j = 1,15
     if (list1(j:j) == ' ')list1(j:j) = '0'
  enddo
  list1 = trim(outputDir)//"/log/listSH"//"."//trim(modelname)//"."//trim(list1)

  open(24, file = list1, status = 'unknown', form = 'formatted')
  write(24,*)
  close(24)



  ! Record the date and time at the beginning of the job
 if (myrank == 0) then
  write(list, '(I6,".",I6)'), imin_glob,imax_glob  !unify scanning file for the process
  do j = 1,15
     if (list(j:j) == ' ')list(j:j) = '0'
  enddo
  list = trim(outputDir)//"/log/calLogSH"//"."//trim(modelname)//"."//trim(list)

  open(1,file =list, status = 'unknown', form = 'formatted')
  call date_and_time(datex,timex)
  write(1,'(/a,a4,a,a2,a,a2,a,a2,a,a2,a,a4)') &
       '    Starting date and time:                     ', &
       datex(1:4),'-',datex(5:6),'-',datex(7:8),'.  ', &
       timex(1:2),':',timex(3:4),':',timex(5:8)
  close (1)
 endif


   ! call clPLM(plm(1:3,0:3,1:theta_n,0:maxlmax),maxlmax,theta(1:theta_n),theta_n)
  !do l = 0,maxlmax
  !   do itheta = 1,theta_n
  !      call caldvec_dejaplm(l,(theta(itheta)/180.d0*pi), (phi(itheta)/180.d0*pi),plm(1:3,0:3,itheta,l),dvec(1:3,-2:2,itheta,l),dvecdt(1:3,-2:2,itheta,l),dvecdp(1:3,-2:2,itheta,l))
  !   enddo
  !enddo
  !deallocate(plm)
! These codes is not used here. they are used in the second part of the program by vectorization.


!******************** Computing the coefficients *********************

  llog = 0
  lref=0
  lmax_r(:)=0
  l=0
!********************************************************************************


!! DK DK now write the displacement and stress coefficients to the same file
!! DK DK to avoid opening too many files on the file system, which is usually slow


    ir0 = r0_n
!write(*,*) myrank,"th !!!!!!!!!!!!!!!start frequency calculation!!!!!!----"
  do i=imin,imax ! each frequency

  write(coutfile, '(I5.5,".G0_SH_",I5.5)') int(r0(r0_n)*10.d0),i
! here we must change the r0(:) index to fix value r0_n,otherwise the error occurs.
  coutfile = trim(modelname)//"."//coutfile
  coutfile = trim(outputDir)//"/Stress/"//coutfile
!  write(*,*) i,'th freq for the filename is:',coutfile,ir0,r0(ir0)
  open(34,file=coutfile, form='unformatted',action='write')

  ! initialisation
  tabg0 = dcmplx(0.d0)
  tabg0der = dcmplx(0.d0)

  omega = 2.d0 * pi * dble(i)/tlen

     if ( i /= 0 ) then
        call callsuf(omega,nzone,vrmax,vsv,lsuf)
        call calcoef( nzone,omega,qmu,coef)
!        plm = 0.d0
        a0 = 0.d0
        a2 = 0.d0
        do j = 1, nzone
           call cala0( nlayer(j),omega,omegai,t(jsp(j)), h1(jsp(j)), &
                & h2(jsp(j)), h3(jsp(j)), h4(jsp(j)),coef(j), cwork(jsp(j)) )
           call overlap( nlayer(j),cwork(jsp(j)),a0( 1,isp(j) ) )
           call cala2( nlayer(j),h4(jsp(j)),coef(j), cwork(jsp(j)) )
           call overlap( nlayer(j),cwork(jsp(j)),a2( 1,isp(j) ) )
        enddo


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

        do l = 0, maxlmax ! l-loop commence

        index_l=index_l+1
        lsq = dsqrt(dble(l)*dble(l+1))
!           do itheta = 1,theta_n
!              call calbvec(l,(theta(itheta)/180.d0*pi),(phi(itheta)/180.d0*pi),plm(1,0,itheta),bvec(1,-2,itheta),bvecdt(1,-2,itheta),bvecdp(1,-2,itheta))
!           enddo   ! this code is not used here which is enable in the second part of the program.
           do ir_=1,r_n
             if (ismall_r(ir_) > 20) then
        if (lmax_r(ir_) > l) then
          lmax_r(ir_)=l
    arret(ir_)=0
        endif
       endif
     enddo

           if (SUM(arret(:)) == 0) then
          llog = l
!! DK DK added this to make sure we do not write the same last buffer twice
                we_left_the_loop_through_exit = .true.
                exit
           endif

!           if (ismall>20) then
!              if (llog>l) then
!                 llog =l
!              endif
!              cycle   !exit
!           endif

           tmpc = 0.d0

           a = 0.d0
           ga2 = 0.d0

           call cala( nn,l,lda,a0,a2,a )
           do ir0 = 1, r0_n
              call calga( 1,omega,omegai,l,t(ins(ir0)),h1(ins(ir0)),h2(ins(ir0)),h3(ins(ir0)),h4(ins(ir0)),coef(spn(ir0)),aa(1:4,ir0))
              call calga( 2,omega,omegai,l,gt(1:8,ir0),gh1(1:8,ir0),gh2(1:8,ir0),gh3(1:8,ir0),gh4(1:8,ir0),coef(spn(ir0)),ga(1:8,ir0))
              call overlap( 2,ga(1:8,ir0),ga2(1:2,1:3,ir0))
           enddo


           do m = -2, 2 ! m-loop commence
     !remembering that the computation for SH has no l=0 components, don't like the PSV part.
              if ( ( m /= 0 ) .and. ( iabs(m) <= iabs(l) ) ) then

                 do ir0 = 1,r0_n
                    ig2 = 0
                    do imt = 2,6
                       g0 = cmplx(0.d0)
                       call setmt(imt,mt)
                       call calg2( l,m,spo(ir0),r0(ir0),mt,mu0(ir0),coef(spn(ir0)),ga(1:8,ir0),aa(1:4,ir0),ga2(1:2,1:3,ir0),gdr(1:3,ir0),g0(isp(spn(ir0))),ig2)

                       if ( (m == -2) .or. (m == -l) ) then
                          if (ig2 == 1) then
!!                             call dclisb0( a,nn,1,lda,g0,eps,dr,z,ier)
!                             call dcsymbdl0( a,1,nn,1,eps,z(1),w(1),ll,lli,llj,ier)
!                             call dcsbdlv0( a,g0,1,nn,eps,z(1),ier )
                              call dcsymbdl0_m_equals_1_nn_equals_1(a,nn,eps,z(1),w(1),ll,lli,llj,ier)
                              call dcsbdlv0_m_equals_1(a,g0,nn,z(1))
                             ig2 = ig2+1
                          else
!!                             call dcsbsub0(a,nn,1,lda,g0,eps,dr,z,ier)
!                             call dcsbdlv0( a,g0,1,nn,eps,z(1),ier )
                              call dcsbdlv0_m_equals_1(a,g0,nn,z(1))
                          endif
                       else
!!                          call dcsbsub0( a,nn,1,lda,g0,eps,dr,z,ier)
!                          call dcsbdlv0( a,g0,1,nn,eps,z(1),ier )
                           call dcsbdlv0_m_equals_1(a,g0,nn,z(1))
                       endif
!dclisb0 is subroutine used for solving simultaneous linear equations with real symmetric positive definite  band matrix by cholesky method.
!  parameters
!    (1) a : 2-dim. array containing the matrix.here the dimension of the 2-dim array is 2*nn, and nn = (nnlayer+1)
!    (2) nn : order of the matrix.
!    (3) 3rd parameter : size of band's half width, excluding diagonal.here is set to 1.
!    (4) lda : row size of the array a in the 'dimension' statement.here lda = 2 which is used for defining the dimension of band matrix a.
!    (5) g0 : 1-dim. array containing the right hand side vector.
!    (6) eps : parameter to check singurarity off the matrix
!              standard value = 1.0d-14
!    (7) dr : 1-dim. working array.
!    (8) z :  1-dim. working array.dimension of the dr and z is nn = (nnlayer+1) too.
!    (9) ier : error code.
!entry dcsbsub0(a, nn, 1, lda, g0, eps, dr, z, ier) is used for forward substitution.

                       if ((imt == 2) .and. (ir0 == r0_n)) then
!! DK DK r_n is of the order of ~500 and thus it is worth trying to vectorize it
!IBM* ASSERT (MINITERCNT(1000))
!!DIR$ loop count min(1000)
        do ir_ = 1,r_n
!! DK DK inlined the call in order to vectorize the loop
          ampratio = 0.d0
!         amp = dsqrt(zabs(g0(nn-1))**2 + zabs(g0(nn))**2)  !this is for the PSV judgement
          amp = abs(g0(nn-1))
        if (amp > maxamp_r(ir_))   maxamp_r(ir_)=amp
        if (amp /= 0.d0 .and. maxamp_r(ir_) /= 0.d0) ampratio = amp / maxamp_r(ir_)
        if (ampratio < ratl .and. l > lsuf) then
          ismall_r(ir_) = ismall_r(ir_) + 1
        else
          ismall_r(ir_) = 0
        endif
        enddo
        call calamp(g0(nn-1),l,lsuf,maxamp,ismall,ratl)
!        call calamp_vector(g0(nn-1),l,lsuf,maxamp_r,ismall_r,ratl,r_n)
                       endif

!! DK DK this loop cannot be vectorized because the call to "interpolate" that it contains cannot be inlined
                       do ir_= 1,r_n
                          g0tmp = (0.d0,0.d0)
                          g0dertmp = (0.d0,0.d0)
                          call interpolate(1,0,r_(ir_),rrsta(1:3,ir_),g0(iista(1:3,ir_)),g0tmp)
                          call interpolate(1,1,r_(ir_),rrsta(1:3,ir_),g0(iista(1:3,ir_)),g0dertmp)
                          tabg0(index_l,1,imt,m,ir_)=g0tmp
                          tabg0der(index_l,1,imt,m,ir_)=g0dertmp
!! this part is different from the PSV part,only SH coefficients is needed.So we don't need so much space to save the g0tmp like PSV routine. in fact the g0tmp and g0dertmp are scalar only.

!                          itheta = ir_
!                             u = cmplx(0.d0)
!                             udr = cmplx(0.d0)
!                             udt = cmplx(0.d0)
!                             udp = cmplx(0.d0)
!                             uder = cmplx(0.d0)
!                             call calu(g0tmp,lsq,bvec(1:3,m,itheta),u(1:3))
!                             call calu(g0dertmp,lsq,bvec(1:3,m,itheta),udr(1:3))
!                             call calu(g0tmp,lsq,bvecdt(1:3,m,itheta),udt(1:3))
!                             call calu(g0tmp,lsq,bvecdp(1:3,m,itheta),udp(1:3))
!                             call locallyCartesianDerivatives(u(1:3),udr(1:3),udt(1:3),udp(1:3),uder(1:3,1:3),r_(ir_),theta(itheta)/180.d0*pi)
!                             call udertoStress(uder(1:3,1:3),stress(1:6,imt,itheta),A0sta(itheta),C0sta(itheta),F0sta(itheta),L0sta(itheta),N0sta(itheta))
!                             displacement(1:3,imt,itheta) = u(1:3)+displacement(1:3,imt,itheta)

                          !enddo
                       enddo ! ir_-loop termine
                    enddo ! imt-loop termine
                 enddo !ir0-loop termine
              endif
           enddo ! m-loop termine

     if (index_l == maxlmax_g) then
              do ir_=1,r_n
          llog0=min(maxlmax_g, lmax_r(ir_) - lref)
                 if (llog0 > 0) call write_a_block_of_coefs_to_disk(tabg0,tabg0der,tabg0_small_buffer,ir_,llog0,maxlmax_g,r_n)

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



        enddo !l-loop termine

        open(24,file =list1, status = 'old',access='append', form = 'formatted')
        write(24,*) i, dble(i)/tlen, llog-1
        close(24)


!        stress(1:6,1:6,1:nsta) = stress(1:6,1:6,1:nsta)/cmplx(0,omega)
!        stresssngl(1:6,1:6,1:nsta) = stress(1:6,1:6,1:nsta)
!        displacementsngl(1:3,1:6,1:nsta) = displacement(1:3,1:6,1:nsta)


     endif


   if (mod(maxlmax,maxlmax_g) /= 0 .or. we_left_the_loop_through_exit) then

    do ir_=1,r_n
       llog0=min(maxlmax_g,lmax_r(ir_)-lref)
       if (llog0 > 0) call write_a_block_of_coefs_to_disk(tabg0,tabg0der,tabg0_small_buffer,ir_,llog0,maxlmax_g,r_n)
    enddo


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


     if (mod(i,8) == 0) then

        open(1,file =list, status = 'old',access='append', form = 'formatted')
        call date_and_time(datex,timex)
        write(1,'(/a,i,a,a4,a,a2,a,a2,a,a2,a,a2,a,a4)') &
             '    Frequency-index ', i, ' :', &
             datex(1:4),'-',datex(5:6),'-',datex(7:8),'.  ', &
             timex(1:2),':',timex(3:4),':',timex(5:8)
        close (1)
     endif


  enddo ! omega-loop termine


  call MPI_Barrier(MPI_COMM_WORLD,ierr)




 if (myrank == nbproc-1) then
   open(1,file =list, status = 'old',access='append', form = 'formatted')
   call date_and_time(datex,timex)
   write(1,'(/a,a4,a,a2,a,a2,a,a2,a,a2,a,a4)') &
        '    Finishing date and time:                     ', &
        datex(1:4),'-',datex(5:6),'-',datex(7:8),'.  ', &
        timex(1:2),':',timex(3:4),':',timex(5:8)
   close (1)
 endif  !using the last processor to write the final statement.


!!call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_FINALIZE(ierr)


  stop






end program TraSH_Write_coefficients
