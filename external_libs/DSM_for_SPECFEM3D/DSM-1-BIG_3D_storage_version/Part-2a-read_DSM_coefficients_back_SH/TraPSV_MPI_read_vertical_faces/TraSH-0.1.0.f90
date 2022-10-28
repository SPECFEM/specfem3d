program TraSH_Zmin


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
include "./constants.h"

INTEGER myrank,nbproc,ierr
integer SubCommunicators,mybigrank,nbbigproc
integer, allocatable :: key(:),color(:)

!------------------------- <  < input matrix >>----------------------------------

  character(120) :: inputdir,outputDir,psvmodel,modelname,stationsinf
  character(120) :: list,list1
  character(40) :: datex,timex
  real(kind(0d0)), parameter :: pi=3.1415926535897932d0
  real(kind(0d0)), parameter :: re=1.d-2, ratc=1.d-10, ratl=1.d-4
  integer, parameter :: maxlmax = 25000, maxlmax_g = 1000
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
  integer :: r_n,r0_n, ir_,ir0,imt,itheta

  character(120) :: coutfile, coutfile1,cinfile,coutfile2,coutfile3
  integer :: imin,imax
  integer :: itranslat
  integer :: nzone


  integer :: i ,j, ier
  real(kind(0d0)) :: dummy
  real(kind(0d0)), allocatable :: vrmin(:), vrmax(:)
  real(kind(0d0)), allocatable :: rrho(:,:), vsv(:,:), vsh(:,:),qmu(:),vpv(:,:),vph(:,:),eta(:,:)


  complex(kind(0d0)), allocatable :: coef(:)
  real(kind(0d0)) :: rmin, rmax
  real(kind(0d0)), allocatable :: vmin(:), gridpar(:), dzpar(:)
  real(kind(0d0)) :: Aval,Cval,Fval,Lval,Nval,rval,thetasin,thetaval
  complex(kind(0d0)) :: inv_lsq,inv_rval,thetacot_complex,inv_thetasin,inverse_of_omega_complex
  real(kind(0d0)) :: omegai
  complex(kind(0d0)) :: u1,u2,u3,udr1,udr2,udr3,udt1,udt2,udt3,udp1,udp2,udp3
  complex(kind(0d0)) :: uder11 , uder12 , uder13 , uder21 , uder22 , uder23 , uder31 , uder32 , uder33

  integer, allocatable :: nlayer(:), iphase(:)
  integer :: ioutercore,ista

  ! variables pour des points stackes
  integer, allocatable :: isp(:),jsp(:)

  ! variables pour la source


  real(kind(0d0)) :: lsq

!-----------------------------------------------------------------------
  ! variables pour des elements de matrice

  ! la frequence
  real(kind(0d0)) :: omega
  complex(kind(0d0)), allocatable :: g0tmp(:,:),g0dertmp(:,:)

  ! des autres

  integer :: llog,m, l,nlay


!-----------------------------------------------------------------------
  complex(kind(0d0)), allocatable :: dvec(:,:,:,:),dvecdt(:,:,:,:),dvecdp(:,:,:,:)
  real(kind(0d0)), allocatable :: plm(:,:,:)
  complex(kind(0d0)), allocatable :: stress(:,:,:,:), displacement(:,:,:,:)
  complex(kind(0e0)), allocatable :: stresssngl(:,:,:,:), displacementsngl(:,:,:,:)
  complex(kind(0e0)), allocatable :: stresssngl_global(:,:,:,:), displacementsngl_global(:,:,:,:)
  complex(kind(0d0))::u(1:3),udr(1:3),udt(1:3),udp(1:3),uder(1:3,1:3)

!-----------------------------------------------------------------------
 ! variables for the MPI communications.
  integer, allocatable :: Ifrq2(:,:)

  integer imin_glob, imax_glob

!-----------------------------------------------------------------------
! other variables
complex(kind(0d0)),allocatable:: tabg0(:,:,:,:,:),tabg0der(:,:,:,:,:)
integer,allocatable:: lmax_r(:)
integer i_color,nb_colors,nb_tot_proc,current_color,imin_current,imax_current,nbproc_current,iread
integer:: index_l,lref,irank,nsta_global,nsta,ifrequ_min,ifrequ_max,lmax_lu
real(kind(0d0)):: bazi,l2
integer, allocatable :: Ind(:),NbInd(:),IndSt(:),NbIndSt(:)
integer :: istamin,istamax,k,ifq,c8,coeff,theta_n
 double precision, parameter :: convert_to_radians = pi/180.d0


!write(*,*) "routine is starting to run!@@@@@@@@@@@"
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mybigrank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nbbigproc,ierr)

if (mybigrank == 0)  write(*,*) 'how many processor do i successfully apply:',nbbigproc
! Now here the variable for the MPI_COMM_WORLD is mybigrank .
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
! the sub-communicator distribution.
call MPI_Bcast(key,nbbigproc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(color,nbbigproc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(Ifrq2,2*nbbigproc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)


 if (mybigrank == 0) then
   call pinputTra_part2(ifrequ_min,ifrequ_max,inputdir,outputDir,psvmodel,modelname,stationsinf,tlen,imin_glob,imax_glob,r0min,r0max,r0delta,r0lat,r0lon,itranslat)
!!here i delete the last myrank para in the vadim's routine,because the myrank para is not used.
   if (nbbigproc == 0) then
    imin=imin_glob
    imax=imax_glob
   endif
!   call WriteFrqByproc(nbproc,imin_glob,imax_glob)  !distributing the task into every processors evenly.
 endif !only the root can excute this function call.


!  allocate(Ifrq(nbproc+1))
!  allocate(Ifrq2(0:nbproc-1,2))
!  if (mybigrank==0)   call ReadFrqsByProc(Ifrq,Ifrq2,para,nbproc+1)  !frequency calculation separation.
!       call MPI_Bcast(Ifrq,nbproc+1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!       call MPI_Bcast(Ifrq2,2*nbproc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!       call MPI_Bcast(para,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!  write(*,*) 'initial setting is completed'

  call MPI_Bcast(inputDir,120,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
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

  call MPI_Barrier(MPI_COMM_WORLD,ierr)

!   write(*,*) 'the',mybigrank,'th processor, frequency section is:',ifrequ_min,ifrequ_max

if (mybigrank == 0) then
  open(20, file = psvmodel, status = 'old', action='read', position='rewind')
  read(20,*) nzone

  !call MPI_Bcast(nzone,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)


  allocate(vrmin(1:nzone))
  allocate(vrmax(1:nzone))
  allocate(rrho(1:4,1:nzone))
  allocate(vsv(1:4,1:nzone))
  allocate(vsh(1:4,1:nzone))
  !allocate(vpv(1:4,1:nzone))
  !allocate(vph(1:4,1:nzone))
  !allocate(eta(1:4,1:nzone))
  allocate(qmu(1:nzone))
  allocate(vmin(1:nzone))
  allocate(gridpar(1:nzone))
  allocate(dzpar(1:nzone))
  allocate(nlayer(1:nzone))
  allocate(iphase(1:nzone))
  allocate(isp(1:nzone))
  allocate(jsp(1:nzone))
  allocate(coef(1:nzone))
  !vpv=0.d0
  !vph=0.d0


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
  allocate(vpv(1:4,1:nzone))
  allocate(vph(1:4,1:nzone))
  allocate(eta(1:4,1:nzone))
  allocate(qmu(1:nzone))
  allocate(vmin(1:nzone))
  allocate(gridpar(1:nzone))
  allocate(dzpar(1:nzone))
  allocate(nlayer(1:nzone))
  allocate(iphase(1:nzone))
  allocate(isp(1:nzone))
  allocate(jsp(1:nzone))
  allocate(coef(1:nzone))
  vpv=0.d0
  vph=0.d0
  eta=1.d0
 if (mybigrank == 0) then
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

 if (mybigrank == 0) then
  open (1,file=stationsinf,status='old',action='read',position='rewind')
  if (itranslat == 1) call translat (r0lat,r0lat)
  open(2,file='recdepth')
  read(2,*) r_n   !now the r_n is not the radius of receivers,it's the index of receriver plane depth index.
  read(1,*) nsta_global
!!now the total nsta_global is stations number in the horizontal direction read from file stationsinf.
!!in fact,the total number for the vertical surface computation is r_n plus nsta_global
!  r_n = nsta
!  theta_n = nsta
 endif

! call MPI_Bcast(nsta,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
 call MPI_Bcast(r_n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
 call MPI_Bcast(nsta_global,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
! call MPI_Bcast(theta_n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
 call MPI_Bcast(r0lat,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
 call MPI_Bcast(r0lon,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!! ! if itraslat = 1.the r0lat and r0lon will be changed.
 call MPI_Barrier(MPI_COMM_WORLD,ierr)

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
! now the size allocated is nsta_global

  allocate(updown(1:r_n))
  allocate(idum(1:r_n))
! it is sufficient for us to save the final result to disk in single precision
  allocate(stresssngl_global(1:6,1:6,1:r_n,1:nsta_global))
  allocate(displacementsngl_global(1:3,1:6,1:r_n,1:nsta_global))
! Here there is more 1 dimension than the zmin routine.
!! Remembering that the order of the dimensions of array stresssngl_global/displacementsngl_global is different
!! from the code of PSV reading vertical part written by Vadim.
  allocate(tabg0(1:maxlmax_g,1,1:6,r_n,-2:2))
  allocate(tabg0der(1:maxlmax_g,1,6,r_n,-2:2))
  allocate(g0tmp(6,r_n))
  allocate(g0dertmp(6,r_n))

! Here there is more 1 dimension than the zmin routine,size of the increased dimension is r_n.
! now the array is different from the PSV routine ,it is half size of the array
! in the PSV routine.

 if (mybigrank == 0) then
 write(*,*) '   PROCESS VERTICAL FACE :',r_n
 write(*,*) '**********************************'
 write(*,*) 'SOURCE :',r0lon,r0lat
 write(*,*) '**********************************'
do i = 1,nsta_global
     read(1,*) stlo_curr,stla_curr
!  be careful, here is the stlo_curr reading first. which is the inverse order compared with the souce position reading
!  so we change the reading order that is compatible with the test station position file.
     stla_g(i)=stla_curr
     stlo_g(i)=stlo_curr
     call epitra(r0lat,r0lon,stla_g(i),stlo_g(i),theta_g(i),phi_g(i),bazi)
 enddo
 endif
!! epitra subroutine is used for transfering the latitude and longitude of the stations
!! to the epicentral distance and azimuth.

!if (mybigrank==0) then
! do i = 1,nsta_global
!  write(*,*) i,theta_g(i),phi_g(i)
! enddo
!endif

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
! receiver plane depth coordinates broadcast and boundary index test.
 do i = 1,r_n
   call calstg4onedepth(nzone,nzone,nzone,vrmin,vrmax,iphase,rrho,0.d0,0.d0,vsv,vsh,1.d0,rmax,r_(i),updown(i),A0sta(i),C0sta(i),F0sta(i),L0sta(i),N0sta(i))
!   call calstg4onedepth(nlay,nzone,nzone,vrmin,vrmax,iphase,rrho,vpv,vph,vsv,vsh,eta,rmax,r_(i), &
!     updown(i),A0sta(i),C0sta(i),F0sta(i),L0sta(i),N0sta(i),idum(i))
!  write(*,*) mybigrank, i,C0sta(i),F0sta(i)
 enddo  !pay attention for the MPI_BCAST.because it's easy to forget to declaire some element by BCAST.
 !call mpi_barrier(MPI_COMM_WORLD,ierr)
 !stop
! depth for stocking the Green function.
 allocate(rrsta(1:3,1:r_n))
 allocate(iista(1:3,1:r_n))

 ! ######################### sous communicateurs #####################"
 ! define the sub-communicators

  call MPI_COMM_SPLIT(MPI_COMM_WORLD,color(mybigrank),key(mybigrank),subcommunicators,ierr)
  call MPI_COMM_RANK(SubCommunicators,myrank,ierr)
  call MPI_COMM_SIZE(SubCommunicators,nbproc,ierr)
!write(*,*) mybigrank,'th old proc.',myrank,'th new proc!!!'
! distribuer les stations sur les procs
  call DistribSta(nsta,nsta_global,myrank,nbproc,istamin,istamax)
!write(*,*) 'myrank  ista :', myrank, istamin,istamax
! nsta is the stations number require to be computed for every proc

  call MPI_Barrier(SubCommunicators,ierr)

  allocate(Ind(nbproc))
  allocate(NbInd(nbproc))
  allocate(IndSt(nbproc))
  allocate(NbIndSt(nbproc))
  call mpi_allgather(nsta,1,MPI_INTEGER,NbInd,1,MPI_INTEGER,SubCommunicators, ierr )

  NbIndSt(:) = NbInd(:)
! NbIndSt and NbInd is the array for storing the computation distribution of every process in sub-group
  ! indices et offsets pour mpi_gatherv
  coeff  = 3*6*r_n
!this is different from the zmin reading routinue. because now there are r_n
!times number of stations in the vertical surface than the zmin surface.
!! in fact,the total number for the vertical surface computation is r_n plus nsta_global
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

! Ind and NbInd array store the start and data-block size for every proc about the displacement.
! IndSt and NbIndSt array store the start and data-block size for every proc about the stress.

  theta_n = nsta
! nsta and theta_n is the stations number require to be computed for every proc
  allocate(theta(1:theta_n))
  allocate(phi(1:theta_n))
  allocate(stlo(theta_n))
  allocate(stla(theta_n))

  allocate(dvec(1:theta_n,1:3,-2:2,0:maxlmax))
  allocate(dvecdt(1:theta_n,1:3,-2:2,0:maxlmax))
  allocate(dvecdp(1:theta_n,1:3,-2:2,0:maxlmax))
!!now the dimension of the matrix dvec et al is increasing for more 1 dimension.
!!the last new dimension is for the l angular order from 0 to maxlmax possible.
!! DK DK moved index 1:theta_n to the first position in order to be able to vectorize some loops
!! DK DK and also got rid of the "l" index because we do not store them for all l any more

  allocate(plm(1:theta_n,1:3,0:3),stat=ier)
  if (ier /= 0) stop 'error: not enough memory to allocate array'

  phi(:)=phi_g(istamin:istamax)
  theta(:)=theta_g(istamin:istamax)
  stla(:)=stla_g(istamin:istamax)
  stlo(:)=stlo_g(istamin:istamax)
!! now the phi,theta,stla,stlo is the relative position array for every proc.

! source depths,currently we do not put any necessary allocation
  r0_n =  int((r0max-r0min)/r0delta)+1
  allocate(r0(1:r0_n))
  do i = 1, r0_n
     r0(i) = r0min + dble(i-1)*r0delta
  enddo
  ir0 = r0_n
!write(*,*) 'the source radiu is :%%%%%%%%',r0

  write(list1, '(I6,".",I6)'), imin_glob,imax_glob  !!unify scanning file for the process
  do j = 1,15
     if (list1(j:j) == ' ')list1(j:j) = '0'
  enddo
  list1 = trim(outputDir)//"/log/listSH"//"."//trim(modelname)//"."//trim(list1)
  open(24, file = list1, status = 'unknown', form = 'formatted')
  write(24,*)
  close(24)

  plm=0.d0
  dvec=(0.d0,0.d0)
  dvecdt=(0.d0,0.d0)
  dvecdp=(0.d0,0.d0)

  do l = 0,maxlmax
!    do itheta = 1,theta_n
!       call caldvec_already_plm(l,(theta(itheta)/180.d0*pi),
!       (phi(itheta)/180.d0*pi),plm(1:3,0:3,itheta,l), &
!               dvec(itheta,1:3,-2:2,l),dvecdt(itheta,1:3,-2:2,l),dvecdp(itheta,1:3,-2:2,l))
!    enddo
!! DK DK created new, vectorized versions of these routines
!! DK DK when l is greater than 5 we can get rid of all the "if" tests in this
!subroutine to make it much faster
call calbvec_vector(l,theta/180.d0*pi,phi/180.d0*pi,plm(1:theta_n,:,:),dvec(1:theta_n,1:3,-2:2,l),dvecdt(1:theta_n,1:3,-2:2,l),dvecdp(1:theta_n,1:3,-2:2,l),theta_n)

  enddo
  deallocate(plm)


  allocate(stress(1:6,1:6,1:r_n,1:nsta))
  allocate(displacement(1:3,1:6,1:r_n,1:nsta))
! it is sufficient for us to save the final result to disk in single precision
!! Remembering that the order of the dimensions of array stress/displacement is different
!! from the code of PSV reading vertical part written by Vadim.
  allocate(stresssngl(1:6,1:6,1:r_n,1:nsta))
  allocate(displacementsngl(1:3,1:6,1:r_n,1:nsta))
! noticing the different order(size) of the dimension.Which is different from
! the routine of PSV part written by vadim.
!! Remembering that the order of the dimensions of array stresssngl/displacementsngl is different
!! from the code of PSV reading vertical part written by Vadim.


  llog = 0
  lref=0
  lmax_r(:)=0

!! DK DK now write the displacement and stress coefficients to the same file
!! DK DK to avoid opening too many files on the file system, which is usually slow


    ir0 = r0_n
!write(*,*) mybigrank,"th !!!!!!!!!!!!!!!start frequency calculation!!!!!!----"

  if (myrank == 0 .and. SLOW_DEBUG_MODE) then
      write(list, '(I6.6,".",I6.6)') ifrequ_min, ifrequ_max
      list =trim(outputDir)//"/log/calLog"//"."//trim(modelname)//"."//trim(list)

     open(25,file =list, status = 'unknown', form = 'formatted')
     call date_and_time(datex,timex)
     write(25,'(/a,a4,a,a2,a,a2,a,a2,a,a2,a,a4)') '    Starting date and time:', &
        datex(1:4),'-',datex(5:6),'-',datex(7:8),'.',timex(1:2),':',timex(3:4),':',timex(5:8)

    write(list1, '(I6.6,".",I6.6)') ifrequ_min, ifrequ_max
    list1 = trim(outputDir)//"/log/list"//"."//trim(modelname)//"."//trim(list1)

     if (SLOW_DEBUG_MODE) then
      open(24, file = list1, status = 'unknown', form = 'formatted')
      write(24,*)
    endif

  endif


  do ifq = ifrequ_min,ifrequ_max
!!frequency loop start now for every proc



   ir0 = r0_n  !Be careful with this declaration
   write(cinfile, '(I5.5,".G0_SH_",I5.5)') int(r0(ir0)*10.d0),ifq
! Noticing that now for the file constrction,i use ifq instead of ifq-1!
   cinfile = trim(modelname)//"."//cinfile
   cinfile = trim(inputdir)//"/Stress/"//cinfile  ! TO DO : il faudrait mettre un inputdir
   if (myrank == 0) then
      open(34,file=cinfile, form='unformatted',action='read')
   endif
!! if (myrank == 0) write(*,*) ir0,'cinfile is:',cinfile

  ! initialisation
     stress = (0.d0,0.d0)
     displacement = (0.d0,0.d0)
     stresssngl = (0.e0,0.e0)
     displacementsngl = (0.e0,0.e0)
     stresssngl_global = (0.e0,0.e0)
     displacementsngl_global = (0.e0,0.e0)

!     plm = 0.d0
!     dvec = (0.d0,0.d0)
!     dvecdp = (0.d0,0.d0)
!     dvecdt = (0.d0,0.d0)


  lmax_r(:)=maxlmax
  omega = 2.d0 * pi * dble(Ifq)/tlen
!for this part of routines.the omega won't be used.
     if ( Ifq /= 0 ) then !!frequency don't equal to 0
        lref=0
        index_l=0
        lmax_lu=0

        do l = 0, maxlmax ! l-loop commence
         call MPI_Barrier(SubCommunicators,ierr)

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
                 write(24,*) 'I already read for l=',l,' and ir_=',ir_,llog
                 write(24,*)
                 write(24,*) mybigrank,'is reading for ir,llog',ir_,llog
              endif

        if (ir_ == -1) then ! this flag indicates the end of all the data to read
                if (SLOW_DEBUG_MODE) then
                   write(24,*) mybigrank,ir_,llog, ' exit '
                endif

          exit
!       else if (ir_ == r_n) then
!! in the case of the Zmin face, we read but ignore all the data from the input file except the last one;
!! in the other code (for vertical faces) we use all the data read, but in this code we use the last one only from the same file
!                read(34) tabg0(1:llog,1,:,:) !here only SH coefficients is read.so second para =1 fixed.
!   read(34) tabg0der(1:llog,1,:,:)
!   lmax_r(ir_) = llog + lref ! recover the computed number of l-order cofficients.
!         lmax_lu=max(lmax_lu,lmax_r(r_n))  ! just for test.

              else !!! VM VM : we need ir_=r_n which is not nessesary the last one thus I test it
          read(34) tabg0(1:llog,1,:,ir_,:) !here only SH coefficients is read.so second para =1 fixed
    read(34) tabg0der(1:llog,1,:,ir_,:)
                lmax_r(ir_) = llog + lref ! lmax courant pour ir_
          lmax_lu=max(lmax_lu,lmax_r(r_n))
                if (SLOW_DEBUG_MODE) then
                  write(24,*) mybigrank,'allready read TABS',ir_,llog,lmax_r(ir_),lmax_lu
                endif

              endif
      enddo
             if (SLOW_DEBUG_MODE) write(24,*) 'end reading block'
     endif ! rank 0 proc reading process is over

      index_l=0
             if (SLOW_DEBUG_MODE .and. myrank == 0)  write(24,*) 'before bcast lmax_r',lmax_r
      call MPI_Bcast(lmax_r,r_n,MPI_INTEGER,0,SubCommunicators,ierr)
      call MPI_Bcast(tabg0,30*r_n*maxlmax_g,MPI_DOUBLE_COMPLEX,0,SubCommunicators,ierr)
      call MPI_Bcast(tabg0der,30*r_n*maxlmax_g,MPI_DOUBLE_COMPLEX,0,SubCommunicators,ierr)
!! here the dize of tabg0/tabg0der array is 30*r_n*maxlmax_g,has been multiplied r_n.
      call MPI_Bcast(ir_,1,MPI_INTEGER,0,SubCommunicators,ierr)
      call MPI_Bcast(llog,1,MPI_INTEGER,0,SubCommunicators,ierr)
      !I have to change the size of tabg0/tabg0der from 60 in PSV to 30 in SH.
     endif ! mod(l,maxlmax_g) == 0 judgement over

     if (l > maxval(lmax_r)) then
                 if (SLOW_DEBUG_MODE .and. myrank == 0) write(24,*) 'I exit because l too big'
                 exit  !! general speaking,l < lmax_r(r_n) when computation is continuing.
           endif
     if (ir_ == -1 .and. llog == -1) then
                if (SLOW_DEBUG_MODE .and. myrank == 0) write(24,*) 'I exit because end of file'
               exit  !! VM VM: I add this test because of bugs
            endif
     !! VM VM: when mod(maxval(lmax_r),maxlmax_g)==0
     !! in fact this condition means that the final data-block of coefficients reading is over now.

           index_l=index_l+1
! the initial value of index_l is set to 0 before the l-loop start. then if the mod(l,maxlmax_g) /= 0
! the index_l will add 1 until the if judgement mod(l,maxlmax_g) =0 valid. if this condition is valid.
! then the index_l is set to 0 again.
! because the loop is start from 0.so at first when index_l = 0,the first llog number of coefficients is
! read and saved. So next the index_l =1.2.3.4.....until l= maxlmax_g. which means the index_l = maxlmax_g
! at the moment. so the coefficients is available in the array tabg0(:,:,:,:).
! and remembering that first para is set from 1 to llog.
     l2 = dble(l)*dble(l+1)
     lsq = dsqrt( l2 )
           inv_lsq = 1.d0 / dcmplx(lsq)

!      if (ifq ==512 .and. l==1800 .and. myrank==0) then
!            !plm = 0.d0
!     !call calbvec_vector_test(l,theta/180.d0*pi,phi/180.d0*pi,plm,dvec,dvecdt,dvecdp,theta_n)
!                  itheta = 1
!      write(*,*) 'before itheta th plm components are:0',plm(itheta,:,0)
!                  write(*,*) 'before itheta th plm components are:1',plm(itheta,:,1)
!                         write(*,*) 'before itheta th plm components are:2',plm(itheta,:,2)
!      write(*,*) 'before itheta th plm components are:3',plm(itheta,:,3)
!      write(*,*) 'before before before'
!      endif

!          call calbvec_vector(l,theta/180.d0*pi,phi/180.d0*pi,plm(1:theta_n,:,:),dvec,dvecdt,dvecdp,theta_n)

!!     do itheta=1,theta_n
!!     call calbvec(l,theta(itheta)/180.d0*pi,phi(itheta)/180.d0*pi,plm(itheta,:,:),dvec(itheta,:,:),dvecdt(itheta,:,:),dvecdp(itheta,:,:))
!!     enddo

!      if (ifq ==512 .and. l==1800 .and. myrank==0) then
!            !plm = 0.d0
!     !call calbvec_vector_test(l,theta/180.d0*pi,phi/180.d0*pi,plm,dvec,dvecdt,dvecdp,theta_n)
!                  itheta = 1
!      write(*,*) 'itheta th plm components are:0',plm(itheta,:,0)
!                  write(*,*) 'itheta th plm components are:1',plm(itheta,:,1)
!                         write(*,*) 'itheta th plm components are:2',plm(itheta,:,2)
!      write(*,*) 'itheta th plm components are:3',plm(itheta,:,3)
!      endif

     !! this is the vectorized calbvec subroutine.but i just make a little change different from PSV case.So i don't need to separate the routine into for_l_less_than_5_ or for_l_more_than_5_ case.
           !! be careful here.i should multiplay the degree to rad factor for the theta and phi array.
           do m = -2, 2 ! m-loop commence
     !remembering that the computation for SH has no l=0 components, don't like the PSV part.
     !and only the nonzero m components are excited by the source.
       g0tmp(1:6,:) = tabg0(index_l,1,1:6,:,m) !noticing again for the second para =1 fixed.
       g0dertmp(1:6,:) = tabg0der(index_l,1,1:6,:,m) !index_l =l+1 actually, is the proper index for the array tabg0/tabg0der.
!! size of dimension for g0tmp/g0dertmp increase 1.
              if ( ( m /= 0 ) .and. ( iabs(m) <= iabs(l) ) ) then

                 do ir0 = 1,r0_n  ! for the source number cycle.in fact is only 1 source here.

!                    do imt = 2,6 !! DK DK this loop cannot be vectorized because the call to "interpolate" that it contains cannot be inlined
                    do ir_= 1,r_n  !it is the vertical surface, so r_n depth is considered.

                      if (l <= lmax_r(ir_)) then ! it is decided by the llog data-block length.

                             Aval = A0sta(ir_)
           Cval = C0sta(ir_)
           Fval = F0sta(ir_)
           Lval = L0sta(ir_)
           Nval = N0sta(ir_)
                                                       rval = r_(ir_)
                              inv_rval = 1.d0 / cmplx(rval)
!                              inv_rval = 1.d0 / cmplx(r_(r_n))


!IBM* ASSERT (MINITERCNT(1000))
!DIR$ loop count min(1000)

                       do itheta = 1,theta_n  ! loop for stations computation by every proc.
!! In fact,above is setting the exit judgement:if (l >= lmax_r(r_n)) exit
!! So l = lmax_r(r_n) can't be valid any more at the moment.

!!#                          g0tmp = (0.d0,0.d0)
!!#                          g0dertmp = (0.d0,0.d0)
!!#                        call interpolate(1,0,r_(ir_),rrsta(1:3,ir_),g0(iista(1:3,ir_)),g0tmp)
!!#                          call interpolate(1,1,r_(ir_),rrsta(1:3,ir_),g0(iista(1:3,ir_)),g0dertmp)
!!#                          tabg0(index_l,1,imt,m,ir_)=g0tmp
!!#                          tabg0der(index_l,1,imt,m,ir_)=g0dertmp
!! this part is different from the PSV part,only SH coefficients is needed.So we don't need so much space to save the g0tmp like PSV routine. in fact the g0tmp and g0dertmp are scalar only.

!        if (ifq ==502 .and. l==1500 .and. m==2 .and. myrank==0 .and. itheta ==1) then
!       write(*,*) 'g0tmp',imt,g0tmp(imt)
!     write(*,*) 'g0dertmp',imt,g0dertmp(imt)
!!   if (imt==6)  write(*,*) 'theta',theta(itheta)
!  if (imt==6)  write(*,*) 'dvec',dvec(itheta,2,:)
!  if (imt==6)  write(*,*) 'dvec',dvec(itheta,3,:)
!!   if (imt==6)  write(*,*) 'dvecdt 2',dvecdt(itheta,3,:)
!!   if (imt==6)  write(*,*) 'dvecdp 2',dvecdp(itheta,3,:)
!  endif

!                          itheta = ir_

              thetaval = theta(itheta) * convert_to_radians
              thetasin = sin(thetaval)
        thetacot_complex = cmplx(cos(thetaval)/thetasin)
        inv_thetasin = 1.d0 / cmplx(thetasin)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! case imt = 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!                               call calup
                                u1 = dcmplx( 0.d0 )
                                u2 = g0tmp(2,ir_)*dvec(itheta,2,m,l)*inv_lsq
                                u3 = g0tmp(2,ir_)*dvec(itheta,3,m,l)*inv_lsq

                                udr1 = dcmplx( 0.d0 )
                                udr2 = g0dertmp(2,ir_)*dvec(itheta,2,m,l)*inv_lsq
                                udr3 = g0dertmp(2,ir_)*dvec(itheta,3,m,l)*inv_lsq

                                udt1 = dcmplx( 0.d0 )
        udt2 = g0tmp(2,ir_)*dvecdt(itheta,2,m,l)*inv_lsq
        udt3 = g0tmp(2,ir_)*dvecdt(itheta,3,m,l)*inv_lsq

                                udp1 = dcmplx( 0.d0 )
        udp2 = g0tmp(2,ir_)*dvecdp(itheta,2,m,l)*inv_lsq
        udp3 = g0tmp(2,ir_)*dvecdp(itheta,3,m,l)*inv_lsq

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
                                stress(1,2,ir_,itheta) = stress(1,2,ir_,itheta)+Cval*uder11+Fval*uder22+Fval*uder33
                                stress(2,2,ir_,itheta) = stress(2,2,ir_,itheta)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder33
                                stress(3,2,ir_,itheta) = stress(3,2,ir_,itheta)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder22

                                stress(4,2,ir_,itheta) = stress(4,2,ir_,itheta)+Lval*(uder12+uder21)
                                stress(5,2,ir_,itheta) = stress(5,2,ir_,itheta)+Lval*(uder13+uder31)
                                stress(6,2,ir_,itheta) = stress(6,2,ir_,itheta)+Nval*(uder23+uder32)

                                displacement(1,2,ir_,itheta) = u1 + displacement(1,2,ir_,itheta)
                                displacement(2,2,ir_,itheta) = u2 + displacement(2,2,ir_,itheta)
                                displacement(3,2,ir_,itheta) = u3 + displacement(3,2,ir_,itheta)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! case imt = 3 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!                               call calup
                                u1 = dcmplx( 0.d0 )
                                u2 = g0tmp(3,ir_)*dvec(itheta,2,m,l)*inv_lsq
                                u3 = g0tmp(3,ir_)*dvec(itheta,3,m,l)*inv_lsq

                                udr1 = dcmplx( 0.d0 )
                                udr2 = g0dertmp(3,ir_)*dvec(itheta,2,m,l)*inv_lsq
                                udr3 = g0dertmp(3,ir_)*dvec(itheta,3,m,l)*inv_lsq

                                udt1 = dcmplx( 0.d0 )
        udt2 = g0tmp(3,ir_)*dvecdt(itheta,2,m,l)*inv_lsq
        udt3 = g0tmp(3,ir_)*dvecdt(itheta,3,m,l)*inv_lsq

                                udp1 = dcmplx( 0.d0 )
        udp2 = g0tmp(3,ir_)*dvecdp(itheta,2,m,l)*inv_lsq
        udp3 = g0tmp(3,ir_)*dvecdp(itheta,3,m,l)*inv_lsq

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
                                stress(1,3,ir_,itheta) = stress(1,3,ir_,itheta)+Cval*uder11+Fval*uder22+Fval*uder33
                                stress(2,3,ir_,itheta) = stress(2,3,ir_,itheta)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder33
                                stress(3,3,ir_,itheta) = stress(3,3,ir_,itheta)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder22

                                stress(4,3,ir_,itheta) = stress(4,3,ir_,itheta)+Lval*(uder12+uder21)
                                stress(5,3,ir_,itheta) = stress(5,3,ir_,itheta)+Lval*(uder13+uder31)
                                stress(6,3,ir_,itheta) = stress(6,3,ir_,itheta)+Nval*(uder23+uder32)

                                displacement(1,3,ir_,itheta) = u1 + displacement(1,3,ir_,itheta)
                                displacement(2,3,ir_,itheta) = u2 + displacement(2,3,ir_,itheta)
                                displacement(3,3,ir_,itheta) = u3 + displacement(3,3,ir_,itheta)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! case imt = 4 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!                               call calup
                                u1 = dcmplx( 0.d0 )
                                u2 = g0tmp(4,ir_)*dvec(itheta,2,m,l)*inv_lsq
                                u3 = g0tmp(4,ir_)*dvec(itheta,3,m,l)*inv_lsq

                                udr1 = dcmplx( 0.d0 )
                                udr2 = g0dertmp(4,ir_)*dvec(itheta,2,m,l)*inv_lsq
                                udr3 = g0dertmp(4,ir_)*dvec(itheta,3,m,l)*inv_lsq

                                udt1 = dcmplx( 0.d0 )
        udt2 = g0tmp(4,ir_)*dvecdt(itheta,2,m,l)*inv_lsq
        udt3 = g0tmp(4,ir_)*dvecdt(itheta,3,m,l)*inv_lsq

                                udp1 = dcmplx( 0.d0 )
        udp2 = g0tmp(4,ir_)*dvecdp(itheta,2,m,l)*inv_lsq
        udp3 = g0tmp(4,ir_)*dvecdp(itheta,3,m,l)*inv_lsq

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
                                stress(1,4,ir_,itheta) = stress(1,4,ir_,itheta)+Cval*uder11+Fval*uder22+Fval*uder33
                                stress(2,4,ir_,itheta) = stress(2,4,ir_,itheta)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder33
                                stress(3,4,ir_,itheta) = stress(3,4,ir_,itheta)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder22

                                stress(4,4,ir_,itheta) = stress(4,4,ir_,itheta)+Lval*(uder12+uder21)
                                stress(5,4,ir_,itheta) = stress(5,4,ir_,itheta)+Lval*(uder13+uder31)
                                stress(6,4,ir_,itheta) = stress(6,4,ir_,itheta)+Nval*(uder23+uder32)

                                displacement(1,4,ir_,itheta) = u1 + displacement(1,4,ir_,itheta)
                                displacement(2,4,ir_,itheta) = u2 + displacement(2,4,ir_,itheta)
                                displacement(3,4,ir_,itheta) = u3 + displacement(3,4,ir_,itheta)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! case imt = 5 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!                               call calup
                                u1 = dcmplx( 0.d0 )
                                u2 = g0tmp(5,ir_)*dvec(itheta,2,m,l)*inv_lsq
                                u3 = g0tmp(5,ir_)*dvec(itheta,3,m,l)*inv_lsq

                                udr1 = dcmplx( 0.d0 )
                                udr2 = g0dertmp(5,ir_)*dvec(itheta,2,m,l)*inv_lsq
                                udr3 = g0dertmp(5,ir_)*dvec(itheta,3,m,l)*inv_lsq

                                udt1 = dcmplx( 0.d0 )
        udt2 = g0tmp(5,ir_)*dvecdt(itheta,2,m,l)*inv_lsq
        udt3 = g0tmp(5,ir_)*dvecdt(itheta,3,m,l)*inv_lsq

                                udp1 = dcmplx( 0.d0 )
        udp2 = g0tmp(5,ir_)*dvecdp(itheta,2,m,l)*inv_lsq
        udp3 = g0tmp(5,ir_)*dvecdp(itheta,3,m,l)*inv_lsq

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
                                stress(1,5,ir_,itheta) = stress(1,5,ir_,itheta)+Cval*uder11+Fval*uder22+Fval*uder33
                                stress(2,5,ir_,itheta) = stress(2,5,ir_,itheta)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder33
                                stress(3,5,ir_,itheta) = stress(3,5,ir_,itheta)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder22

                                stress(4,5,ir_,itheta) = stress(4,5,ir_,itheta)+Lval*(uder12+uder21)
                                stress(5,5,ir_,itheta) = stress(5,5,ir_,itheta)+Lval*(uder13+uder31)
                                stress(6,5,ir_,itheta) = stress(6,5,ir_,itheta)+Nval*(uder23+uder32)

                                displacement(1,5,ir_,itheta) = u1 + displacement(1,5,ir_,itheta)
                                displacement(2,5,ir_,itheta) = u2 + displacement(2,5,ir_,itheta)
                                displacement(3,5,ir_,itheta) = u3 + displacement(3,5,ir_,itheta)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! case imt = 6 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!                               call calup
                                u1 = dcmplx( 0.d0 )
                                u2 = g0tmp(6,ir_)*dvec(itheta,2,m,l)*inv_lsq
                                u3 = g0tmp(6,ir_)*dvec(itheta,3,m,l)*inv_lsq

                                udr1 = dcmplx( 0.d0 )
                                udr2 = g0dertmp(6,ir_)*dvec(itheta,2,m,l)*inv_lsq
                                udr3 = g0dertmp(6,ir_)*dvec(itheta,3,m,l)*inv_lsq

                                udt1 = dcmplx( 0.d0 )
        udt2 = g0tmp(6,ir_)*dvecdt(itheta,2,m,l)*inv_lsq
        udt3 = g0tmp(6,ir_)*dvecdt(itheta,3,m,l)*inv_lsq

                                udp1 = dcmplx( 0.d0 )
        udp2 = g0tmp(6,ir_)*dvecdp(itheta,2,m,l)*inv_lsq
        udp3 = g0tmp(6,ir_)*dvecdp(itheta,3,m,l)*inv_lsq

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
                                stress(1,6,ir_,itheta) = stress(1,6,ir_,itheta)+Cval*uder11+Fval*uder22+Fval*uder33
                                stress(2,6,ir_,itheta) = stress(2,6,ir_,itheta)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder33
                                stress(3,6,ir_,itheta) = stress(3,6,ir_,itheta)+Fval*uder11+Aval*uder22+Aval*uder33-2.d0*Nval*uder22

                                stress(4,6,ir_,itheta) = stress(4,6,ir_,itheta)+Lval*(uder12+uder21)
                                stress(5,6,ir_,itheta) = stress(5,6,ir_,itheta)+Lval*(uder13+uder31)
                                stress(6,6,ir_,itheta) = stress(6,6,ir_,itheta)+Nval*(uder23+uder32)

                                displacement(1,6,ir_,itheta) = u1 + displacement(1,6,ir_,itheta)
                                displacement(2,6,ir_,itheta) = u2 + displacement(2,6,ir_,itheta)
                                displacement(3,6,ir_,itheta) = u3 + displacement(3,6,ir_,itheta)



!                             u = cmplx(0.d0)
!                             udr = cmplx(0.d0)
!                             udt = cmplx(0.d0)
!                             udp = cmplx(0.d0)
!                             uder = cmplx(0.d0)
!                             call calu(g0tmp(imt,ir_),lsq,dvec(itheta,1:3,m,l),u(1:3)) !noticing that the first para of dvec array is the station index now which is different from the original TraSH routine.
!!     if (ifq ==347 .and. l <= 4 .and. myrank==0 .and. itheta ==6 .and. imt==3) then
!!        write(*,*) 'u(2) is',l,m,displacement(2,imt,itheta),u(2)
!!       ! write(*,*) 'u(2) is',l,m,imt,u(2)
!!     endif
!                             call calu(g0dertmp(imt,ir_),lsq,dvec(itheta,1:3,m,l),udr(1:3))
!                             call calu(g0tmp(imt,ir_),lsq,dvecdt(itheta,1:3,m,l),udt(1:3))
!                             call calu(g0tmp(imt,ir_),lsq,dvecdp(itheta,1:3,m,l),udp(1:3))
!                             call locallyCartesianDerivatives(u(1:3),udr(1:3),udt(1:3),udp(1:3),uder(1:3,1:3),r_(ir_),theta(itheta)/180.d0*pi)  !the unit of the theta(:) is degree. Here the r_(:) index is the depth index,
!       !and the index for the theta(:) is using station index.
!!           if (ifq ==47 .and. l==4 .and. myrank==0 .and. itheta ==6 .and. imt==3 .and. ir_==10) then
!!              !write(*,*) l,m,'uder temporary matrix is:',stress(1:6,imt,ir_,itheta)
!!        write(*,*)  m,'COEFFICIENTS ARE:',A0sta(ir_),C0sta(ir_),F0sta(ir_),L0sta(ir_),N0sta(ir_)
!!     endif
!                             call udertoStress(uder(1:3,1:3),stress(1:6,imt,ir_,itheta),A0sta(ir_),C0sta(ir_),F0sta(ir_),L0sta(ir_),N0sta(ir_))
!                             displacement(1:3,imt,ir_,itheta) = u(1:3)+displacement(1:3,imt,ir_,itheta)

                       enddo ! for the  itheta = 1,theta_n  loop
                      endif !for the condition (l <= lmax_r(ir_))
                    enddo ! ir_= 1,r_n  loop termine
!                    enddo ! imt-loop termine
                 enddo !ir0-loop termine
              endif !for the condition  ( m /= 0 ) .and. ( iabs(m) <= iabs(l) )
           enddo ! m-loop termine

        enddo !l-loop termine

        open(24,file =list1, status = 'old',access='append', form = 'formatted')
        write(24,*) ifq, dble(ifq)/tlen, llog-1
        close(24)


   inverse_of_omega_complex = 1.d0 / dcmplx(0,omega)

         do ir_=1,r_n
!IBM* ASSERT (MINITERCNT(1000))
!DIR$ SIMD
!DIR$ loop count min(1000)
           do ista=1,nsta
              stresssngl(1,1,ir_,ista) = stress(1,1,ir_,ista) * inverse_of_omega_complex
        stresssngl(2,1,ir_,ista) = stress(2,1,ir_,ista) * inverse_of_omega_complex
        stresssngl(3,1,ir_,ista) = stress(3,1,ir_,ista) * inverse_of_omega_complex
        stresssngl(4,1,ir_,ista) = stress(4,1,ir_,ista) * inverse_of_omega_complex
        stresssngl(5,1,ir_,ista) = stress(5,1,ir_,ista) * inverse_of_omega_complex
        stresssngl(6,1,ir_,ista) = stress(6,1,ir_,ista) * inverse_of_omega_complex
        stresssngl(1,2,ir_,ista) = stress(1,2,ir_,ista) * inverse_of_omega_complex
        stresssngl(2,2,ir_,ista) = stress(2,2,ir_,ista) * inverse_of_omega_complex
        stresssngl(3,2,ir_,ista) = stress(3,2,ir_,ista) * inverse_of_omega_complex
        stresssngl(4,2,ir_,ista) = stress(4,2,ir_,ista) * inverse_of_omega_complex
        stresssngl(5,2,ir_,ista) = stress(5,2,ir_,ista) * inverse_of_omega_complex
        stresssngl(6,2,ir_,ista) = stress(6,2,ir_,ista) * inverse_of_omega_complex
        stresssngl(1,3,ir_,ista) = stress(1,3,ir_,ista) * inverse_of_omega_complex
        stresssngl(2,3,ir_,ista) = stress(2,3,ir_,ista) * inverse_of_omega_complex
        stresssngl(3,3,ir_,ista) = stress(3,3,ir_,ista) * inverse_of_omega_complex
        stresssngl(4,3,ir_,ista) = stress(4,3,ir_,ista) * inverse_of_omega_complex
        stresssngl(5,3,ir_,ista) = stress(5,3,ir_,ista) * inverse_of_omega_complex
        stresssngl(6,3,ir_,ista) = stress(6,3,ir_,ista) * inverse_of_omega_complex
        stresssngl(1,4,ir_,ista) = stress(1,4,ir_,ista) * inverse_of_omega_complex
        stresssngl(2,4,ir_,ista) = stress(2,4,ir_,ista) * inverse_of_omega_complex
        stresssngl(3,4,ir_,ista) = stress(3,4,ir_,ista) * inverse_of_omega_complex
        stresssngl(4,4,ir_,ista) = stress(4,4,ir_,ista) * inverse_of_omega_complex
        stresssngl(5,4,ir_,ista) = stress(5,4,ir_,ista) * inverse_of_omega_complex
        stresssngl(6,4,ir_,ista) = stress(6,4,ir_,ista) * inverse_of_omega_complex
        stresssngl(1,5,ir_,ista) = stress(1,5,ir_,ista) * inverse_of_omega_complex
        stresssngl(2,5,ir_,ista) = stress(2,5,ir_,ista) * inverse_of_omega_complex
        stresssngl(3,5,ir_,ista) = stress(3,5,ir_,ista) * inverse_of_omega_complex
        stresssngl(4,5,ir_,ista) = stress(4,5,ir_,ista) * inverse_of_omega_complex
        stresssngl(5,5,ir_,ista) = stress(5,5,ir_,ista) * inverse_of_omega_complex
        stresssngl(6,5,ir_,ista) = stress(6,5,ir_,ista) * inverse_of_omega_complex
        stresssngl(1,6,ir_,ista) = stress(1,6,ir_,ista) * inverse_of_omega_complex
        stresssngl(2,6,ir_,ista) = stress(2,6,ir_,ista) * inverse_of_omega_complex
        stresssngl(3,6,ir_,ista) = stress(3,6,ir_,ista) * inverse_of_omega_complex
        stresssngl(4,6,ir_,ista) = stress(4,6,ir_,ista) * inverse_of_omega_complex
        stresssngl(5,6,ir_,ista) = stress(5,6,ir_,ista) * inverse_of_omega_complex
        stresssngl(6,6,ir_,ista) = stress(6,6,ir_,ista) * inverse_of_omega_complex
     enddo
   enddo


         do ir_=1,r_n
!IBM* ASSERT (MINITERCNT(1000))
!DIR$ SIMD
!DIR$ loop count min(1000)
          do ista = 1,nsta
                displacementsngl(1,1,ir_,ista) = displacement(1,1,ir_,ista)
    displacementsngl(2,1,ir_,ista) = displacement(2,1,ir_,ista)
    displacementsngl(3,1,ir_,ista) = displacement(3,1,ir_,ista)
    displacementsngl(1,2,ir_,ista) = displacement(1,2,ir_,ista)
    displacementsngl(2,2,ir_,ista) = displacement(2,2,ir_,ista)
    displacementsngl(3,2,ir_,ista) = displacement(3,2,ir_,ista)
    displacementsngl(1,3,ir_,ista) = displacement(1,3,ir_,ista)
    displacementsngl(2,3,ir_,ista) = displacement(2,3,ir_,ista)
    displacementsngl(3,3,ir_,ista) = displacement(3,3,ir_,ista)
    displacementsngl(1,4,ir_,ista) = displacement(1,4,ir_,ista)
    displacementsngl(2,4,ir_,ista) = displacement(2,4,ir_,ista)
    displacementsngl(3,4,ir_,ista) = displacement(3,4,ir_,ista)
    displacementsngl(1,5,ir_,ista) = displacement(1,5,ir_,ista)
    displacementsngl(2,5,ir_,ista) = displacement(2,5,ir_,ista)
    displacementsngl(3,5,ir_,ista) = displacement(3,5,ir_,ista)
    displacementsngl(1,6,ir_,ista) = displacement(1,6,ir_,ista)
    displacementsngl(2,6,ir_,ista) = displacement(2,6,ir_,ista)
    displacementsngl(3,6,ir_,ista) = displacement(3,6,ir_,ista)
    enddo
   enddo



!        stress(1:6,1:6,1:r_n,1:theta_n) = stress(1:6,1:6,1:r_n,1:theta_n)/cmplx(0,omega)
!        stresssngl(1:6,1:6,1:r_n,1:theta_n) = stress(1:6,1:6,1:r_n,1:theta_n)
!        displacementsngl(1:3,1:6,1:r_n,1:theta_n) = displacement(1:3,1:6,1:r_n,1:theta_n)

        Coeff = 3*6*nsta*r_n  !here the size of data for gethering is multiplying r_n number of depth.
        call mpi_gatherv(displacementsngl,Coeff,MPI_COMPLEX,displacementsngl_global,NbInd,Ind,MPI_COMPLEX,0,SubCommunicators,ierr)
!! here the 0 is the receiving proc in every sub-group.The NbInd is the size of data-block of every proc in sub-group.Ind is similar to the displacement in the mpi_gatherv which represent the display of the receiving data-block in the root to the 1st element of displacementsngl_global array.

        Coeff = 6*6*nsta*r_n  !here the size of data for gethering is multiplying r_n number of depth.
        call mpi_gatherv(stresssngl,Coeff,MPI_COMPLEX,stresssngl_global,NbIndSt,IndSt,MPI_COMPLEX,0,SubCommunicators,ierr)
!! So now the displacementsngl_global and stresssngl_global are the displacement and stress matrix of the every station.but noticing that for different sub-group,the computation frequency is different still.(even though the station_global is the same)

!if (ifq ==347 .and. myrank==0) then
!   ir_ =6
!   write(*,*) 'displacement(2,:,ir_)',displacementsngl_global(2,1:6,ir_)
!endif
     endif !!for the condition ifrq /= 0

   if (myrank == 0) then
    ir0 = r0_n
    write(coutfile, '(I5.5,".Stress_SH_",I5.5)') int(r0(ir0)*10.d0),ifq
    coutfile = trim(modelname)//"."//coutfile
    coutfile = trim(outputDir)//"/Stress/"//coutfile

!! DK DK not a good idea to give indices when writing the whole array: some compilers could create copies or use loops in such a case
!! DK DK unformatted sequential I/Os (because sequential I/O is usually faster than direct-access I/O)
!! DK DK and writing a single big file is also much better
    open(44,file=coutfile,form='unformatted',action='write')
    write(44) stresssngl_global
    write(44) displacementsngl_global
    close(44)

!    write(coutfile1, '(I5.5,".Displacement_SH_",I5.5)') int(r0(ir0)*10.d0),ifq
!    coutfile1 = trim(modelname)//"."//coutfile1
!    coutfile1 = trim(outputDir)//"/Stress/"//coutfile1

!    open(45,file=coutfile1,form='unformatted',action='write')
!    write(45) displacementsngl_global
!    close(45)
   endif

   if (myrank == 0) close(34)
!! Only the myrank==0 proc can start the 34 file reading.


  enddo ! frequency omega-loop termine

!  write(*,*) 'write the radius of the last plane',r_(r_n)  !!test the depth of the computed plane
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call MPI_COMM_FREE(SubCommunicators,ierr)
  call MPI_FINALIZE(ierr)
  stop

if (mybigrank == 0) then
!@@@@ test the displacement compared with the original routine by Kawaii
   ir0 = r0_n
   write(coutfile, '(I5,".Stress_SH")'), int(r0(ir0)*10.d0)
   do j = 1,7
     if (coutfile(j:j) == ' ')coutfile(j:j) = '0'
   enddo

         coutfile = trim(modelname)//"."//coutfile
         coutfile = trim(outputDir)//"Stress/"//coutfile

    write(coutfile1, '(I5,".Displacement_SH")'), int(r0(ir0)*10.d0)
    do j = 1,7
      if (coutfile1(j:j) == ' ')coutfile1(j:j) = '0'
    enddo

         coutfile1 = trim(modelname)//"."//coutfile1
         coutfile1 = trim(outputDir)//"Displacement/"//coutfile1

 write(*,*) 'reading process is begining:the start/end frequency are',imin_glob,imax_glob
 stop  !! VM VM stop because imin_gloab and imax_gloab are not the frequency
       !! VM VM boundaries computed before then bug when doing this loop
 do i=imin_glob,imax_glob
! write(*,*) 'frequency counting',i
     stresssngl_global=cmplx(0.e0)
     displacementsngl_global=cmplx(0.e0)

      write(coutfile2, '(I5.5,".Stress_SH_",I5.5)') int(r0(ir0)*10.d0),i
      coutfile2 = trim(modelname)//"."//coutfile2
      coutfile2 = trim(outputdir)//"/Stress/"//coutfile2  ! TO DO : il faudrait mettre un inputdir
          open(34,file=coutfile2, form='unformatted',action='read')

!      write(coutfile3, '(I5.5,".Displacement_SH_",I5.5)') int(r0(ir0)*10.d0),i
!      coutfile3 = trim(modelname)//"."//coutfile3
!      coutfile3 = trim(outputdir)//"/Stress/"//coutfile3  ! TO DO : il faudrait mettre un inputdir
!          open(35,file=coutfile3, form='unformatted',action='read')

!!coutifle2 and coutfile3 are the file used for storing the Stress and Displacement for SH every frequency computation

       read(34) stresssngl_global(1:6,1:6,1:r_n,1:nsta_global)
       read(34) displacementsngl_global(1:3,1:6,1:r_n,1:nsta_global)
       close(34)

!       read(35) displacementsngl_global(1:3,1:6,1:r_n,1:nsta_global)
!       close(35)

        open(1,file=coutfile,status='unknown',form='unformatted', access = 'direct', recl=2*6*6*kind(0e0)*nsta_global)
!        open(1,file=coutfile,status='unknown',form='unformatted', access = 'direct', recl=2*6*6*kind(0e0)*nsta_global*r_n)
!        write(1,rec=i+1) stresssngl_global(1:6,1:6,1:r_n,1:nsta_global)
        write(1,rec=i+1) stresssngl_global(1:6,1:6,10:10,1:nsta_global)  !just for the single depth test.
        close(1)

        open(1,file=coutfile1,status='unknown',form='unformatted', access = 'direct', recl=2*3*6*kind(0e0)*nsta_global)
!        open(1,file=coutfile1,status='unknown',form='unformatted', access = 'direct', recl=2*3*6*kind(0e0)*nsta_global*r_n)
!        write(1,rec=i+1) displacementsngl_global(1:3,1:6,1:r_n,1:nsta_global)
        write(1,rec=i+1) displacementsngl_global(1:3,1:6,10:10,1:nsta_global) !just for the single depth test.

        close(1)

 enddo  ! end i-frequency cycle

endif



call MPI_COMM_FREE(SubCommunicators,ierr)
call MPI_FINALIZE(ierr)


  stop






end program TraSH_Zmin
