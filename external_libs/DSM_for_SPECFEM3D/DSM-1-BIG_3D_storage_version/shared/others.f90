
subroutine pinputTra_Part1(outputDir,psvmodel,modelname,stationsinf,tlen,imin,imax,r0min, &
r0max,r0delta,r0lat,r0lon,itranslat,myrank)
  implicit none
  character(120) tmpfile
  character(120) :: dummy,outputDir,psvmodel,modelname,stationsinf
  real(kind(0d0)) :: tlen,r0min,r0max,r0delta
  real(kind(0d0)) :: r0lat,r0lon
  integer :: imin,imax,itranslat,myrank

   write(tmpfile,'(a10,i5.5)') 'TmpWrkFile',myrank
  open(20,file='inputIASP.infTra') !! VM
  open(unit=1, file=tmpfile,status='unknown')
100 continue
  read(20,110) dummy  !! VM
110 format(a120)
  if (dummy(1:1) == '#') goto 100
  if (dummy(1:3) == 'end') goto 120
  write(1,110) dummy
  goto 100
120 continue
  close(1)
  close(20)

  open(unit=1,file=tmpfile,status='unknown')
  read(1,110) outputDir
  read(1,110) psvmodel
  read(1,110) modelname
  read(1,110) stationsinf
  outputDir=trim(outputDir)
  psvmodel=trim(psvmodel)
  modelname=trim(modelname)
  stationsinf=trim(stationsinf)
  read(1,*) tlen
  read(1,*) r0min,r0lat,r0lon
  r0min = 6371.d0 -r0min ! because in this version we write the source DEPTH
  r0max=r0min
  r0delta=20.d0
  read(1,*) imin,imax
  read(1,*) itranslat
  close(1)
end subroutine pinputTra_Part1



subroutine pinputTra_Part2(ifrequ_min,ifrequ_max,inputdir,outputDir,psvmodel,modelname,stationsinf,tlen,imin,imax,r0min, &
r0max,r0delta,r0lat,r0lon,itranslat,myrank)
  implicit none
  character(120) tmpfile
  character(120) :: dummy,outputDir,psvmodel,modelname,stationsinf,inputdir,inputIASP
  real(kind(0d0)) :: tlen,r0min,r0max,r0delta
  real(kind(0d0)) :: r0lat,r0lon,m1,m2,m3,m4,m5,m6
  integer :: imin,imax,itranslat,myrank,ifrequ_min,ifrequ_max

  call get_command_argument(1,inputIASP)
  if (trim(inputIASP) == '') then
     write(*,*) 'input file not specified:  cancel '
     stop
  endif
  !write(tmpfile,'(a10,i5.5)')
  tmpfile = trim(inputIASP)//'.read'
  open(20,file=trim(inputIASP)) !! VM
  open(unit=1, file=trim(tmpfile),status='unknown')
100 continue
  read(20,110) dummy  !! VM
110 format(a120)
  if (dummy(1:1) == '#') goto 100
  if (dummy(1:3) == 'end') goto 120
  write(1,110) dummy
  goto 100
120 continue
  close(1)
  close(20)

  open(unit=1,file=tmpfile,status='unknown')
  read(1,110) outputDir
  read(1,110) psvmodel
  read(1,110) modelname
  read(1,110) stationsinf
  outputDir=trim(outputDir)
  psvmodel=trim(psvmodel)
  modelname=trim(modelname)
  stationsinf=trim(stationsinf)
  read(1,*) tlen
  read(1,*) r0min,r0lat,r0lon
  r0min = 6371.d0 -r0min ! because in this version we write the source DEPTH
  r0max=r0min
  r0delta=20.d0
  read(1,*) imin,imax
  read(1,*) itranslat
  read(1,*) m1,m2,m3,m4,m5,m6
  read(1,110) inputdir
  read(1,*) ifrequ_min,ifrequ_max
  close(1)
end subroutine pinputTra_Part2




subroutine pinput(outputDir,psvmodel,modelname,tlen,rmin_,rmax_,rdelta_,r0min,r0max,r0delta, &
thetamin,thetamax,thetadelta,imin,imax,rsgtswitch)
  implicit none
  character(120), parameter :: tmpfile='tmpworkingfile_for_SGTforPinv'
  character(120) :: dummy,outputDir,psvmodel,modelname
  real(kind(0d0)) :: tlen,rmin_,rmax_,rdelta_,r0min,r0max,r0delta
  real(kind(0d0)) :: thetamin,thetamax,thetadelta
  integer :: imin,imax,rsgtswitch

  open(unit=1, file=tmpfile,status='unknown')
100 continue
  read(5,110) dummy
110 format(a120)
  if (dummy(1:1) == '#') goto 100
  if (dummy(1:3) == 'end') goto 120
  write(1,110) dummy
  goto 100
120 continue
  close(1)

  open(unit=1,file=tmpfile,status='unknown')
  read(1,110) outputDir
  read(1,110) psvmodel
  read(1,110) modelname
  outputDir=trim(outputDir)
  psvmodel=trim(psvmodel)
  modelname=trim(modelname)
  read(1,*) tlen
  read(1,*) rmin_,rmax_,rdelta_
  read(1,*) r0min
  r0max=r0min
  r0delta=20.d0
  read(1,*) thetamin,thetamax,thetadelta
  read(1,*) imin,imax
  read(1,*) rsgtswitch
  close(1)
end subroutine pinput






subroutine udertorrsgt(icomp,uder,rsgt)
  implicit none
  integer :: icomp
  complex(kind(0d0)) :: uder(1:3,1:3), rsgt(1:8)

  ! icomp = 1
  ! This is for P inversion

  if (icomp == 1) then ! vertical component
     rsgt(1) = rsgt(1) + uder(1,1)
     rsgt(2) = rsgt(2) + uder(2,1) + uder(1,2)
     rsgt(3) = rsgt(3) + 5.d-1*uder(2,2) - 5.d-1*uder(3,3)
     rsgt(4) = rsgt(4) - 5.d-1*uder(2,2) - 5.d-1*uder(3,3)
  endif

  if (icomp == 2) then ! radial component
     rsgt(5) = rsgt(5) + uder(1,1)
     rsgt(6) = rsgt(6) - 5.d-1*uder(2,2) - 5.d-1*uder(3,3)
     rsgt(7) = rsgt(7) - 5.d-1*uder(2,2) + 5.d-1*uder(3,3)
     rsgt(8) = rsgt(8) - 2.d0*uder(1,2)  - uder(2,1)
  endif


  return
end subroutine udertorrsgt

!
subroutine udertorsgt(icomp,uder,rsgt)
  implicit none
  integer :: icomp
  complex(kind(0d0)) :: uder(1:3,1:3), rsgt(1:2)

  ! icomp = 1
  ! This is for P inversion

  rsgt(icomp) = rsgt(icomp) + uder(1,1) + uder(2,2) + uder(3,3)

  return
end subroutine udertorsgt

!

subroutine udertotsgt(imt,uder,tsgt)
  implicit none
  integer :: imt
  complex(kind(0d0)) :: uder(1:3,1:3),tsgt(1:4)

  ! 1 <= imt <= 4

  ! This is for P inversion

  tsgt(imt) = tsgt(imt) + uder(1,1) + uder(2,2) + uder(3,3)
  return
end subroutine udertotsgt


subroutine udertoStress(uder,stress,A,C,F,L,N)
  implicit none
  complex(kind(0d0)) :: uder(1:3,1:3), stress(1:6)
  real(kind(0d0)) :: A,C,F,L,N
  !write(100,*) 'Stess 1' ,Stress(1)
  stress(1) = stress(1)+C*uder(1,1)+F*uder(2,2)+F*uder(3,3)
  stress(2) = stress(2)+F*uder(1,1)+A*uder(2,2)+A*uder(3,3)-2.d0*N*uder(3,3)
  stress(3) = stress(3)+F*uder(1,1)+A*uder(2,2)+A*uder(3,3)-2.d0*N*uder(2,2)

  stress(4) = stress(4)+L*(uder(1,2)+uder(2,1))
  stress(5) = stress(5)+L*(uder(1,3)+uder(3,1))
  stress(6) = stress(6)+N*(uder(2,3)+uder(3,2))

  return

end subroutine udertoStress


!

subroutine setmt(imt,mt)
  implicit none
  ! We use the  rr, tt, pp, rt, rp, tp order here!!

  real(kind(0d0)) :: mt(3,3)
  integer :: imt
  mt = 0.d0
  if (imt == 1) mt(1,1) = 1.d0
  if (imt == 2) mt(2,2) = 1.d0
  if (imt == 3) mt(3,3) = 1.d0
  if (imt == 4) mt(1,2) = 1.d0
  if (imt == 5) mt(1,3) = 1.d0
  if (imt == 6) mt(2,3) = 1.d0
  return
end subroutine setmt

!

subroutine setmttest(imt,mt)
  implicit none
  ! We use the  rr, tt, pp, rt, rp, tp order here!!

  real(kind(0d0)) :: mt(3,3)
  integer :: imt
  mt = 0.d0
  if (imt == 1) then
     mt(1,1) = 1.d0
     mt(2,2) = 1.d0
     mt(3,3) = 1.d0
  endif

  return
end subroutine setmttest


!

subroutine locallyCartesianDerivatives(u,udr,udt,udp,uder,r,theta)
  implicit none
  complex(kind(0d0)), intent(in) :: u(1:3),udr(1:3),udt(1:3),udp(1:3)
  complex(kind(0d0)), intent(out) :: uder(1:3,1:3)
  real(kind(0d0)), intent(in) :: r,theta
  real(kind(0d0)) :: thetasin,thetacot

  thetasin = sin(theta)
  thetacot = cos(theta)/thetasin

  ! 1,2,3: r,theta,phi; , denotes the partial derivatives

  uder(1,1) = udr(1)
  uder(1,2) = (udt(1)-u(2))/cmplx(r)
  uder(1,3) = (udp(1)/cmplx(thetasin)-u(3))/cmplx(r)

  uder(2,1) = udr(2)
  uder(2,2) = (udt(2)+u(1))/cmplx(r)
  uder(2,3) = (udp(2)/cmplx(thetasin)-u(3)*cmplx(thetacot))/cmplx(r)

  uder(3,1) = udr(3)
  uder(3,2) = udt(3)/cmplx(r)
  uder(3,3) = (udp(3)/cmplx(thetasin)+u(1)+u(2)*cmplx(thetacot))/cmplx(r)

end subroutine locallyCartesianDerivatives

!

subroutine calnl( nzone,vs,iphase,nsl,nll )

  ! counting of nsl and nll.
  implicit none
  integer:: nzone,iphase(*),nsl,nll
  real(kind(0d0)):: vs(4,*)
  integer:: i

  nsl = 0
  nll = 0
  do i=1,nzone
     if ( ( vs(1,i) == 0.d0 ) .and. ( vs(2,i) == 0.d0 ) .and. ( vs(3,i) == 0.d0 ) .and. ( vs(4,i) == 0.d0 ) ) then
        nll = nll + 1
        iphase(i) = 2
     else
        nsl = nsl + 1
        iphase(i) = 1
     endif
  enddo
  return
end subroutine calnl


!


subroutine calgrid( nzone,vrmin,vrmax,vp,vs,rmin,rmax, imax,lmin,tlen,vmin,gridpar,dzpar )
  implicit none
  real(kind(0d0)), parameter:: pi=3.1415926535897932d0
  integer:: nzone,imax,lmin
  real(kind(0d0)):: vrmin(*),vrmax(*),vp(4,*),vs(4,*)
  real(kind(0d0)):: rmin,rmax,tlen,vmin(*),gridpar(*),dzpar(*)
  integer:: izone,i,j
  real(kind(0d0)):: coef1,coef2,v(4),vs1,vs2,rh,omega,amax,gtmp

  do izone=1,nzone
     ! computing the S-velocity at each zone
     if ( vs(1,izone) == 0.d0 ) then
        do i=1,4
           v(i) = vp(i,izone)
        enddo
     else
        do i=1,4
           v(i) = vs(i,izone)
        enddo
     endif
     vs1 = 0.d0
     vs2 = 0.d0
     do j=1,4
        if ( j == 1 ) then
           coef1 = 1.d0
        else
           coef1 = coef1 * ( vrmin(izone) / rmax )
        endif
        if ( j == 1 ) then
           coef2 = 1.d0
        else
           coef2 = coef2 * ( vrmax(izone) / rmax )
        endif
        vs1 = vs1 + v(j) * coef1
        vs2 = vs2 + v(j) * coef2
     enddo
     rh = vrmax(izone) - vrmin(izone)
     ! computing omega,amax
     omega = 2.d0 * pi * dble(imax) / tlen
     if ( vs1 >= vs2 ) then
        vmin(izone) = vs2
     else
        vmin(izone) = vs1
     endif
     amax = vrmax(izone)
     gtmp = ( omega * omega ) / ( vmin(izone) * vmin(izone) )   - ( (dble(lmin)+0.5d0) * (dble(lmin)+0.5d0) )&
     / ( amax * amax )
     if ( gtmp > 0.d0 ) then
        dzpar(izone)   = dsqrt( 1.d0/gtmp )
        gridpar(izone) = rh / dzpar(izone)
     else
        dzpar(izone)   = 0.d0
        gridpar(izone) = 0.d0
     endif
  enddo
  ! rearangement of gridpar
  gtmp = 0.d0
  do izone=1,nzone
     gtmp = gtmp + gridpar(izone)
  enddo
  do izone=1,nzone
     if ( gridpar(izone) > 0.d0 ) then
        gridpar(izone) = gridpar(izone) / gtmp
     else
        rh = vrmax(izone) - vrmin(izone)
        gridpar(izone) = rh / ( rmax - rmin ) * 0.1d0
     endif
  enddo
! re-rearangement of gridpar
  gtmp = 0.d0
  do izone=1,nzone
     gtmp = gtmp + gridpar(izone)
  enddo
  do izone=1,nzone
     gridpar(izone) = gridpar(izone) / gtmp
  enddo
  return
end subroutine calgrid


!

subroutine calra( maxnlay,maxnslay,maxnllay,maxnzone,maxnstack,nlayer,inlayer,jnlayer,jnslay,jnllay,gridpar,dzpar, &
nzone,vrmin,vrmax,iphase,rmin,rmax,nslay,nllay,nnl,ra,re, nsta,rsta,rrsta,istazone,iista,r0,cista)
  ! Computing the number and the location of grid points.

  implicit none
  real(kind(0d0)), parameter:: pi=3.1415926535897932d0
  integer:: maxnlay,maxnslay,maxnllay,maxnzone,maxnstack
  integer:: nlayer,inlayer,jnlayer,jnslay,jnllay
  integer:: nzone,iphase(*),nslay,nllay,nnl(maxnzone)
  real(kind(0d0)):: gridpar(*),dzpar(*),vrmin(*),vrmax(*),rmin,rmax,r0
  real(kind(0d0)):: ra(maxnlay+maxnzone+1)
  integer:: izone,itmp,i,ntmp
  real(kind(0d0)):: rh,re
  integer:: nsta
  real(kind(0d0)):: rsta(maxnstack),rrsta(3,maxnstack)
  real(kind(0d0)):: ctmp               ! distance betwee source and the nearst
  integer:: istazone(maxnstack)
  integer:: iista(3,maxnstack)
  integer:: ista,j,cista

  ctmp = 7000.d0
  ! Initializing the data
  nslay = 0
  nllay = 0
  inlayer = 0
  do i=1,maxnlay+maxnzone+1
     ra(i) = 0.d0
  enddo
  do izone=1,nzone
     nnl(izone) = 0
  enddo
  jnlayer = 0
  jnslay = 0
  jnllay = 0
  do i=1,maxnstack
     do j=1,3
        rrsta(j,i) = 0.d0
        iista(j,i) = 0
     enddo
  enddo

  do i=1,maxnstack
     istazone(i) = 0
  enddo
  ! computing the number and the location of the grid points
  ra(1) = rmin
  itmp = 1
  do izone=1,nzone
     rh = vrmax(izone) - vrmin(izone)
     if (dzpar(izone) == 0.d0) then
        ntmp = 1
     else
        ntmp = int( sqrt(3.3d0 / re ) * rh / dzpar(izone) / 2.d0 / pi  / 7.d-1 + 1 )
     endif
     !                            ! ntmp (see Geller & Takeuchi 1995 6.2)
     nnl(izone) = ntmp
     if ( nnl(izone) < 5 ) nnl(izone)=5
     if ( iphase(izone) == 1 ) nslay = nslay + nnl(izone)
     if ( nslay > maxnslay )  stop  'nslay is too large. (calra)'
     if ( iphase(izone) == 2 ) nllay = nllay + nnl(izone)
     if ( nllay > maxnllay )  stop  'nllay is too large. (calra)'
     do I=1,nnl(izone)
        itmp = itmp + 1
        if ( itmp > maxnlay ) stop  'nlay is too large. (calra)'
        ra(itmp) = vrmin(izone) + rh * dble(i) / dble( nnl(izone) )
     enddo
  enddo

  itmp = 1
  do izone=1,nzone
     do i=1,nnl(izone)
        do ista=1,nsta
           if ( (ra(itmp) < rsta(ista)) .and. (rsta(ista) <= ra(itmp+1)) ) then
              if (i /= nnl(izone)) then
                 istazone(ista) = izone
                 if (iphase(istazone(ista)) == 2) stop 'rsta is in liquid layer. (calra)'
                 rrsta(1,ista) = ra(itmp)
                 rrsta(2,ista) = ra(itmp+1)
                 rrsta(3,ista) = ra(itmp+2)

                 iista(1,ista) = i
                 iista(2,ista) = i + 1
                 iista(3,ista) = i + 2
              else
                 istazone(ista) = izone
                 if (iphase(istazone(ista)) == 2) stop 'rsta is in liquid layer. (calra)'
                 rrsta(1,ista) = ra(itmp-1)
                 rrsta(2,ista) = ra(itmp)
                 rrsta(3,ista) = ra(itmp+1)

                 iista(1,ista) = i - 1
                 iista(2,ista) = i
                 iista(3,ista) = i + 1
              endif
              if (dabs(r0-rsta(ista)) < ctmp) then
                 cista = ista
                 ctmp = dabs(r0-rsta(ista))
              endif
           endif
        enddo
        itmp = itmp + 1
     enddo
  enddo

  ! recouting the total number of grid points
  inlayer = 0
  do izone=1,nzone
     inlayer = inlayer + nnl(izone)
  enddo
  jnlayer = jnlayer + inlayer
  jnslay  = jnslay  + nslay
  jnllay  = jnllay  + nllay

  return
end subroutine calra

!


subroutine calra2( maxnlay,maxnzone,maxnstack,nlayer,inlayer,jnlayer,jnslay,jnllay,gridpar,nzone,vrmin,vrmax, &
iphase, rmin,rmax,r0,nslay,nllay,nnl,ra, nsta,rsta,rrsta,istazone,iista)
  ! Computing the number and the location of grid points.

  implicit none
  integer:: maxnlay,maxnzone,maxnstack
  integer:: nlayer,inlayer,jnlayer,jnslay,jnllay
  integer:: nzone,iphase(*),nslay,nllay,nnl(maxnzone)
  real(kind(0d0)):: gridpar(*),vrmin(*),vrmax(*),rmin,rmax,r0
  real(kind(0d0)):: ra(maxnlay+maxnzone+1)
  integer:: izone,itmp,i
  real(kind(0d0)):: rh

  integer:: nsta
  real(kind(0d0)):: rsta(maxnstack),rrsta(3,maxnstack)
  integer:: istazone(maxnstack)
  integer:: iista(3,maxnstack)
  integer:: ista,j

  ! Initializing the data
  nslay = 0
  nllay = 0
  inlayer = 0
  do i=1,maxnlay+maxnzone+1
     ra(i) = 0.d0
  enddo
  do izone=1,nzone
     nnl(izone) = 0
  enddo
  do i=1,maxnstack
     do j=1,3
        rrsta(j,i) = 0.d0
        iista(j,i) = 0
     enddo
  enddo
  do i=1,maxnstack
     istazone(i) = 0
  enddo
  jnlayer = 0
  jnslay = 0
  jnllay = 0

  !     computing the number and the location of the grid points
  ra(1) = rmin
  itmp = 1
  do izone=1,nzone
     rh = vrmax(izone) - vrmin(izone)
     nnl(izone) = dint( dble(nlayer) * gridpar(izone) )+ 1
     if ( nnl(izone) < 5 ) nnl(izone)=5
     if ( iphase(izone) == 1 ) nslay = nslay + nnl(izone)
     if ( iphase(izone) == 2 ) nllay = nllay + nnl(izone)
     do i=1,nnl(izone)
        itmp = itmp + 1
        ra(itmp) = vrmin(izone) + rh * dble(i) / dble( nnl(izone) )
     enddo
  enddo
  itmp = 1
  do izone=1,nzone
     do i=1,nnl(izone)
        do ista=1,nsta
           if ( (ra(itmp) < rsta(ista)) .and. (rsta(ista) <= ra(itmp+1)) ) then
              if (i /= nnl(izone)) then
                 istazone(ista) = izone
                 if (iphase(istazone(ista)) == 2) stop 'rsta is in liquid layer. (calra)'
                 rrsta(1,ista) = ra(itmp)
                 rrsta(2,ista) = ra(itmp+1)
                 rrsta(3,ista) = ra(itmp+2)

                 iista(1,ista) = i
                 iista(2,ista) = i + 1
                 iista(3,ista) = i + 2
              else
                 istazone(ista) = izone
                 if (iphase(istazone(ista)) == 2) stop 'rsta is in liquid layer. (calra)'
                 rrsta(1,ista) = ra(itmp-1)
                 rrsta(2,ista) = ra(itmp)
                 rrsta(3,ista) = ra(itmp+1)

                 iista(1,ista) = i - 1
                 iista(2,ista) = i
                 iista(3,ista) = i + 1
              endif
           endif
        enddo
        itmp = itmp + 1
     enddo
  enddo
  ! recouting the total number of grid points
  inlayer = 0
  do izone=1,nzone
     inlayer = inlayer + nnl(izone)
  enddo
  jnlayer = jnlayer + inlayer
  jnslay  = jnslay  + nslay
  jnllay  = jnllay  + nllay

  return
end subroutine calra2

!

subroutine calsp( maxnzone,ndc,nsl,nll,iphase,nlayer,nslay,nllay,isp,jsp,ksp,issp,ilsp,lsp,jssp,isdr,jsdr, &
ildr,jdr,kdr )
  ! Computing the stack points.
  implicit none
  integer:: maxnzone
  integer:: ndc,nsl,nll,iphase(*),nlayer(maxnzone)
  integer:: nslay,nllay
  integer:: isp(maxnzone),jsp(maxnzone),ksp(maxnzone)
  integer:: issp(maxnzone),ilsp(maxnzone)
  integer:: lsp(maxnzone),jssp(maxnzone)
  integer:: isdr,jsdr,ildr,jdr,kdr
  integer:: i,isl,ill

  ! Initialization of the data
  do i=1,maxnzone
     isp(i)  = 0
     jsp(i)  = 0
     ksp(i)  = 0
     issp(i) = 0
     ilsp(i) = 0
     lsp(i)  = 0
     jssp(i) = 0
  enddo
  isdr = 0
  jsdr = 0
  ildr = 0
  jdr = 0
  kdr = 0
  ! computation of isp,jsp,ksp,issp,ilsp,lsp
  isp(1)  = 1
  jsp(1)  = 1
  ksp(1)  = 1
  issp(1) = 1
  ilsp(1) = 1
  lsp(1)  = 1
  jssp(1) = 1
  isl = 0
  ill = 0
  do i=1,ndc
     isp(i+1) = isp(i) + nlayer(i)
     if ( iphase(i) == 1 ) then
        jsp(i+1) = jsp(i) + 16 * nlayer(i)
        ksp(i+1) = ksp(i) + 2 * ( nlayer(i) + 1 )
        lsp(i+1) = lsp(i) + 4 * nlayer(i)
        isl = isl + 1
        if ( isl /= nsl ) then
           issp(isl+1) = issp(isl) + 4 * nlayer(i)
           jssp(isl+1) = jssp(isl) + nlayer(i) + 1
        endif
     else
        jsp(i+1) = jsp(i) + 4 * nlayer(i)
        ksp(i+1) = ksp(i) + ( nlayer(i) + 1 )
        lsp(i+1) = lsp(i) + 2 * nlayer(i)
        ill = ill + 1
        if ( ill /= nll ) ilsp(ill+1) = ilsp(ill) + 4 * nlayer(i)
     endif
  enddo
  isdr = 0
  jsdr = 0
  ildr = 0
  jdr = 0
  isdr = isdr + issp(nsl)-1 + 4 * nlayer(ndc+1)
  jsdr = jsdr + jssp(nsl)-1 + nlayer(ndc+1) + 1
  ildr = ildr + 4 * nllay
  jdr =  jdr  + jsp(ndc+1)-1 + 16 * nlayer(ndc+1)
  kdr =  kdr + ksp(ndc+1)-1 + 2 * ( nlayer(ndc+1)+1 )

  return
end subroutine calsp

!


subroutine calspo( maxnlay,maxnzone,ndc,rdc,iphase,inlayer,r0,rmin,rmax,ra,isp,spo,spn )
  ! Computing the source location.
  implicit none
  integer:: maxnlay,maxnzone,ndc,iphase(*)
  integer:: inlayer,isp(maxnzone),spn
  real(kind(0d0)):: rdc(*),r0,rmin,rmax,ra(maxnlay+maxnzone+1),spo
  integer:: itmp

  ! checking the parameter
  if ( (r0 < rmin) .or. (r0 > rmax) ) stop 'The source location is improper.(calspo)'
  spo = 0
  ! computing 'spo'
  if ( r0 == rmax ) then
     spo = dble( inlayer ) - 0.01d0
     r0 = ra(inlayer) + (spo-dble(inlayer-1)) * ( ra(inlayer+1) -ra(inlayer) )
  else
     itmp = 2
110  continue
     if ( r0 < ra(itmp) ) then
        continue
     else
        itmp = itmp + 1
        goto 110
     endif
     spo = dble(itmp-2)  + ( r0-ra(itmp-1) )   / ( ra(itmp)-ra(itmp-1) )
     ! temporal handling
     if ( (spo-dble(itmp-2)) < 0.01d0 ) then
        spo = dble(itmp-2) + 0.01d0
        r0 = ra(itmp-1) + (spo-dble(itmp-2)) * ( ra(itmp)-ra(itmp-1) )
     endif
     if ( (spo-dble(itmp-2)) > 0.99d0 ) then
        spo = dble(itmp-2) + 0.99d0
        r0 = ra(itmp-1) + (spo-dble(itmp-2)) * ( ra(itmp)-ra(itmp-1) )
     endif
  endif
  ! computing 'spn'
  spn = 0
  itmp = 1
 130  continue
  if ( iphase(itmp) == 1 ) then
     spn = spn + 1
     if ( r0 <= rdc(itmp) ) then
        continue
     else
        itmp = itmp + 1
        goto 130
     endif
  else
     spn = spn + 1
     if ( r0 <= rdc(itmp) ) stop 'The source is in the liquid layer.(calspo)'
     itmp = itmp + 1
     goto 130
  endif
  ! changing 'spo'
  spo = spo - dble( isp(spn) - 1 )

  return
end subroutine calspo


!


subroutine calstg( maxnlay,maxnzone,nzone,iphase,rrho,vpv,vph,vsv,vsh,eta,nnl,ra,rmax,vnp,vra,rho, &
kappa,ecKx,ecKy,ecKz, mu,ecL,ecN,r0,spn,ecC0,ecF0,ecL0 )

  ! Computing the structure grid points.
  implicit none
  integer:: maxnlay,maxnzone,nzone,iphase(*),nnl(*),vnp,spn
  real(kind(0d0)):: rrho(4,*),vpv(4,*),vph(4,*),vsv(4,*),vsh(4,*),eta(4,*)
  real(kind(0d0)):: ra(*),rmax
  real(kind(0d0)):: vra(*),rho(*),kappa(*),ecKx(*),ecKy(*),ecKz(*)
  real(kind(0d0)):: mu(*),ecL(*),ecN(*)
  real(kind(0d0)):: ecA,ecC,ecF
  real(kind(0d0)):: r0,ecA0,ecC0,ecF0,ecL0
  real(kind(0d0)):: trho,tvpv,tvph,tvsv,tvsh,teta,coef
  integer:: izone,i,j,itmp,jtmp

  ! initializing the data
  call vecinit( maxnlay+2*maxnzone+1,vra )
  call vecinit( maxnlay+2*maxnzone+1,rho )
  call vecinit( maxnlay+2*maxnzone+1,kappa )
  call vecinit( maxnlay+2*maxnzone+1,ecKx )
  call vecinit( maxnlay+2*maxnzone+1,ecKy )
  call vecinit( maxnlay+2*maxnzone+1,ecKz )
  call vecinit( maxnlay+2*maxnzone+1,mu )
  call vecinit( maxnlay+2*maxnzone+1,ecL )
  call vecinit( maxnlay+2*maxnzone+1,ecN )
  ! computing the structure grid points
  itmp = 0
  jtmp = 0
  do izone=1,nzone
     do i=1,nnl(izone)+1
        itmp = itmp + 1
        jtmp = jtmp + 1
        vra(itmp) = ra(jtmp)
        ! --- evaluating the density and elastic constants at this point
        trho = 0.d0
        tvpv = 0.d0
        tvph = 0.d0
        tvsv = 0.d0
        tvsh = 0.d0
        teta = 0.d0
        do j=1,4
           if ( j == 1 ) then
              coef = 1.d0
           else
              coef = coef * ( vra(itmp) / rmax )
           endif
           trho  = trho  + rrho(j,izone)  * coef
           tvpv  = tvpv  + vpv(j,izone)   * coef
           tvph  = tvph  + vph(j,izone)   * coef
           tvsv  = tvsv  + vsv(j,izone)   * coef
           tvsh  = tvsh  + vsh(j,izone)   * coef
           teta  = teta  + eta(j,izone)   * coef
        enddo
        rho(itmp) = trho
        ecL(itmp)  = rho(itmp) * tvsv * tvsv
        ecN(itmp)  = rho(itmp) * tvsh * tvsh
        ecA = trho * tvph * tvph
        ecC = trho * tvpv * tvpv
        ecF = teta * ( ecA - 2.d0 * ecL(itmp) )
        kappa(itmp) = ( 4.d0 * ecA + ecC  + 4.d0 * ecF - 4.d0 * ecN(itmp) ) / 9.d0
        ecKx(itmp) = ecA - 4.d0 / 3.d0 * ecN(itmp)
        ecKy(itmp) = ecF + 2.d0 / 3.d0 * ecN(itmp)
        ecKz(itmp) = ( ecC + 2.d0 * ecF ) / 3.d0
     enddo
     jtmp = jtmp - 1
  enddo
  vnp = itmp

  trho = 0.d0
  tvpv = 0.d0
  tvph = 0.d0
  tvsv = 0.d0
  tvsh = 0.d0
  teta = 0.d0
  do j=1,4
     if ( j == 1 ) then
        coef = 1.d0
     else
        coef = coef * ( r0 / rmax )
     endif
     trho  = trho  + rrho(j,spn)  * coef
     tvpv  = tvpv  + vpv(j,spn)   * coef
     tvph  = tvph  + vph(j,spn)   * coef
     tvsv  = tvsv  + vsv(j,spn)   * coef
     tvsh  = tvsh  + vsh(j,spn)   * coef
     teta  = teta  + eta(j,spn)   * coef
  enddo
  ecL0 = trho * tvsv * tvsv
  ecA0 = trho * tvph * tvph
  ecC0 = trho * tvpv * tvpv
  ecF0 = teta * ( ecA0 - 2.d0 * ecL0 )

  return
end subroutine calstg

!


subroutine caltstg( maxnlay,maxnzone,nzone,rrho,vpv,vph,vsv,vsh,eta,nnl,ra,rmax,tvra,tkappa,tecKx,tecKy, &
tecKz,tmu,tecL,tecN)

  ! Computing the structure grid points.

  implicit none
  integer:: maxnlay,maxnzone,nzone,nnl(*)
  real(kind(0d0)):: rrho(4,*),vpv(4,*),vph(4,*),vsv(4,*),vsh(4,*),eta(4,*)
  real(kind(0d0)):: ra(*),rmax
  real(kind(0d0)):: tvra(*),tkappa(*),tmu(*)
  real(kind(0d0)):: tecKx(*),tecKy(*),tecKz(*),tecL(*),tecN(*)
  real(kind(0d0)):: trho,tvpv,tvph,tvsv,tvsh,teta,coef
  real(kind(0d0)):: ecA,ecC,ecF
  integer:: izone,i,j,itmp,jtmp

  call vecinit( maxnlay+2*maxnzone+1,tvra )
  call vecinit( maxnlay+2*maxnzone+1,tkappa )
  call vecinit( maxnlay+2*maxnzone+1,tecKx )
  call vecinit( maxnlay+2*maxnzone+1,tecKy )
  call vecinit( maxnlay+2*maxnzone+1,tecKz )
  call vecinit( maxnlay+2*maxnzone+1,tmu )
  call vecinit( maxnlay+2*maxnzone+1,tecL )
  call vecinit( maxnlay+2*maxnzone+1,tecN )
  ! computing the structure grid points
  itmp = 0
  jtmp = 0
  do izone=1,nzone
     do i=1,nnl(izone)+1
        itmp = itmp + 1
        jtmp = jtmp + 1
        tvra(itmp) = ra(jtmp)
        ! --- evaluating the density and elastic constants at this point
        trho = 0.d0
        tvpv = 0.d0
        tvph = 0.d0
        tvsv = 0.d0
        tvsh = 0.d0
        teta = 0.d0
        do j=1,4
           if ( j == 1 ) then
              coef = 1.d0
           else
              coef = coef * ( tvra(itmp) / rmax )
           endif
           trho = trho + rrho(j,izone) * coef
           tvpv  = tvpv  + vpv(j,izone)   * coef
           tvph  = tvph  + vph(j,izone)   * coef
           tvsv  = tvsv  + vsv(j,izone)   * coef
           tvsh  = tvsh  + vsh(j,izone)   * coef
           teta  = teta  + eta(j,izone)   * coef
        enddo
        tecL(itmp)  = trho * tvsv * tvsv
        tecN(itmp)  = trho * tvsh * tvsh
        ecA = trho * tvph * tvph
        ecC = trho * tvpv * tvpv
        ecF = teta * ( ecA - 2.d0 * tecL(itmp) )
        tkappa(itmp) = ( 4.d0 * ecA + ecC + 4.d0 * ecF - 4.d0 * tecN(itmp) )/ 9.d0
        tecKx(itmp) = ecA - 4.d0 / 3.d0 * tecN(itmp)
        tecKy(itmp) = ecF + 2.d0 / 3.d0 * tecN(itmp)
        tecKz(itmp) = ( ecC + 2.d0 * ecF ) / 3.d0
     enddo
     jtmp = jtmp - 1
  enddo

  return
end subroutine caltstg

!

subroutine calinv(vnp,rho,kappa,rhoinv,kappainv)
  ! Computing the inverse of density and elastic constant.
  implicit none
  integer:: vnp,i
  real(kind(0d0)):: rho(*),kappa(*),rhoinv(*),kappainv(*)

  do i=1,vnp
     rhoinv(i)   = 1.d0 / rho(i)
     kappainv(i) = 1.d0 / kappa(i)
  enddo

  return
end subroutine calinv

!


subroutine submat( nlayer,ha,hb,h )

  ! Subtracting matrix `hb' from matrix `ha'.
  implicit none
  integer:: nlayer
  real(kind(0d0)):: ha(*),hb(*),h(*)
  integer:: i

  do i=1,4*nlayer
     h(i) = ha(i) - hb(i)
  enddo
  return
end subroutine submat


!


subroutine calspdr( maxnzone,nzone,iphase,nlayer,jjdr,kkdr )

  implicit none
  integer:: maxnzone,nzone,iphase(*)
  integer:: nlayer(maxnzone),jjdr(*),kkdr(*)
  integer:: izone

  jjdr(1) = 1
  kkdr(1) = 1
  do izone=1,nzone-1
     if ( iphase(izone) == 1 ) then
        jjdr(izone+1) = jjdr(izone) + 16 * nlayer(izone)
        if ( iphase(izone+1) == 1 ) then
           kkdr(izone+1) = kkdr(izone) + 2 * nlayer(izone)
        else
           kkdr(izone+1) = kkdr(izone) + 2 * ( nlayer(izone)+1 )
        endif
     else
        jjdr(izone+1) = jjdr(izone) + 4 * nlayer(izone)
        if ( iphase(izone+1) == 1 ) then
           kkdr(izone+1) = kkdr(izone) + ( nlayer(izone)+1 )
        else
           kkdr(izone+1) = kkdr(izone) + nlayer(izone)
        endif
     endif
  enddo

  return
end subroutine calspdr

!



subroutine calmdr( omega,l,nzone,vrmin,vrmax,vmin,dzpar,rmax,sufzone )
  implicit none
  real(kind(0d0)), parameter:: pi=3.1415926535897932d0
  integer:: l,nzone,sufzone
  real(kind(0d0)):: omega,vrmin(*),vrmax(*),vmin(*),dzpar(*),rmax
  integer:: izone
  real(kind(0d0)):: gtmp,tdzpar
  sufzone = 0
  do izone=1,nzone
     gtmp = ( omega * omega ) / ( vmin(izone) * vmin(izone) )   - ( (dble(l)+0.5d0) * (dble(l)+0.5d0) ) &
     / ( vrmax(izone) * vrmax(izone) )
     if ( gtmp > 0.d0 ) then
        tdzpar = sqrt( 1.d0/gtmp )
     else
        if ( vrmax(izone) > rmax*(1-2.d0*pi/(dble(l)+0.50)) ) then
           tdzpar = 0.d0
        else
           sufzone = izone
           tdzpar = 0.d0
        endif
     endif
  enddo

  return
end subroutine calmdr

!


subroutine calu0( c0,bvec,u )
  implicit none
  complex(kind(0d0)):: c0,bvec,u

  u = u + c0 * bvec

  return
end subroutine calu0

!


subroutine calulcd0( c0,c0der,rsta,theta, bvec,bvecdt,bvecdp,ulcd )

  implicit none
  complex(kind(0d0)):: c0,c0der,bvec(3),bvecdt(3),bvecdp(3),ulcd(9)
  real(kind(0d0)):: rsta,theta
  complex(kind(0d0)):: u1,uder11,uder12,uder13

  u1 = c0 * bvec(1)
  uder11 = c0der * bvec(1)
  uder12 = c0 * bvecdt(1)
  uder13 = c0 * bvecdp(1)

  ulcd(1) = ulcd(1) + uder11
  ulcd(2) = ulcd(2) + uder12 / rsta
  ulcd(3) = ulcd(3) + uder13 / rsta / dsin(theta)
  ulcd(5) = ulcd(5) + u1 / rsta
  ulcd(9) = ulcd(9) + u1 / rsta

  return
end subroutine calulcd0

!


subroutine calu( c0,lsq,bvec,u )
  implicit none
  real(kind(0d0)):: lsq
  complex(kind(0d0)):: c0(2),bvec(3),u(3)

  u(1) = u(1) + c0(1) * bvec(1)
  u(2) = u(2) + c0(2) * bvec(2) / dcmplx(lsq)
  u(3) = u(3) + c0(2) * bvec(3) / dcmplx(lsq)

  return
end subroutine calu


!

subroutine calup(c1,c2,lsq,bvec,u)
  implicit none
  real(kind(0d0)) :: lsq
  complex(kind(0d0)) :: c1,c2, bvec(1:3), u(1:3)

  u(1) = u(1) + c1*bvec(1)
  u(2) = u(2) + c2*bvec(2)/dcmplx(lsq)
  u(3) = u(3) + c2*bvec(3)/dcmplx(lsq)

  return
end subroutine calup

!

subroutine calup0(c1,bvec,u)
  implicit none

  complex(kind(0d0)), intent(inout) :: u(1:3)
  complex(kind(0d0)), intent(in) :: c1,bvec(1:3)

  u(1) = u(1) + c1*bvec(1)

end subroutine calup0


!


subroutine calulcd( c0,c0der,lsq,rsta,theta, bvec,bvecdt,bvecdp,ulcd )

  implicit none
  real(kind(0d0)):: lsq,rsta,theta
  complex(kind(0d0)):: c0(2),c0der(2)
  complex(kind(0d0)):: bvec(3),bvecdt(3),bvecdp(3),ulcd(9)

  complex(kind(0d0)):: u1,u2,u3
  complex(kind(0d0)):: uder11,uder12,uder13
  complex(kind(0d0)):: uder21,uder22,uder23
  complex(kind(0d0)):: uder31,uder32,uder33

  u1 = c0(1) * bvec(1)
  u2 = c0(2) * bvec(2) / dcmplx(lsq)
  u3 = c0(2) * bvec(3) / dcmplx(lsq)
  ! partial derivatives of u
  uder11 = c0der(1) * bvec(1)
  uder12 = c0(1) * bvecdt(1)
  uder13 = c0(1) * bvecdp(1)
  uder21 = c0der(2) * bvec(2) / dcmplx(lsq)
  uder22 = c0(2) * bvecdt(2) / dcmplx(lsq)
  uder23 = c0(2) * bvecdp(2) / dcmplx(lsq)
  uder31 = c0der(2) * bvec(3) / dcmplx(lsq)
  uder32 = c0(2) * bvecdt(3) / dcmplx(lsq)
  uder33 = c0(2) * bvecdp(3) / dcmplx(lsq)
  ! locally Cartesian derivatives of u
  ulcd(1) = ulcd(1) + uder11
  ulcd(2) = ulcd(2) + ( uder12 - u2 ) / rsta
  ulcd(3) = ulcd(3) + ( uder13 / dsin(theta) - u3 ) / rsta
  ulcd(4) = ulcd(4) + uder21
  ulcd(5) = ulcd(5) + ( uder22 + u1 ) / rsta
  ulcd(6) = ulcd(6) + ( uder23 - u3 * dcos(theta) )  / rsta / dsin(theta)
  ulcd(7) = ulcd(7) + uder31
  ulcd(8) = ulcd(8) + uder32 / rsta
  ulcd(9) = ulcd(9)  + ( ( uder33 + u2 * dcos(theta) ) / dsin(theta) + u1 ) / rsta
  return
end subroutine calulcd

!


subroutine matinit( n1,n2,a )
  implicit none
  integer:: n1,n2,i,j
  real(kind(0d0)):: a(n1,*)

  do j=1,n2
     do i=1,n1
        a(i,j) = 0.d0
     enddo
  enddo

  return
end subroutine matinit

!


subroutine cmatinit( n1,n2,a )
  implicit none
  integer:: n1,n2,i,j
  complex(kind(0d0)):: a(n1,*)

  do j=1,n2
     do i=1,n1
        a(i,j) = dcmplx( 0.d0 )
     enddo
  enddo
  return
end subroutine cmatinit

!


subroutine vecinit( nn,b )

  ! Filling zero to the vector 'g'.
  implicit none
  integer:: nn,i
  real(kind(0d0)):: b(*)

  do i=1,nn
     b(i) = 0.d0
  enddo
  return
end subroutine vecinit

!


subroutine cvecinit( nn,b )

  ! Filling zero to the vector 'g'.
  implicit none
  integer:: nn,i
  complex(kind(0d0)):: b(*)

  do i=1,nn
     b(i) = dcmplx( 0.d0 )
  enddo
  return
end subroutine cvecinit

!

subroutine fillinpb( nderiv,b )

  implicit none
  integer:: nderiv
  complex(kind(0d0)):: b(3)

  if ( (nderiv /= 0) .and. (nderiv /= 1) .and. (nderiv /= 2) ) stop 'invalid argument (fillinpb)'
  if (nderiv == 0) then
     b(1) = dcmplx( 1.d0 )
     b(2) = dcmplx( 0.d0 )
     b(3) = dcmplx( 0.d0 )
  else if (nderiv == 1) then
     b(1) = dcmplx( 0.d0 )
     b(2) = dcmplx( 1.d0 )
     b(3) = dcmplx( 0.d0 )
  else if (nderiv == 2) then
     b(1) = dcmplx( 0.d0 )
     b(2) = dcmplx( 0.d0 )
     b(3) = dcmplx( 1.d0 )
  endif

  return
end subroutine fillinpb


!



subroutine calamp( g,l,lsuf,maxamp,ismall,ratl )
  implicit none
  integer:: ismall,l,lsuf
  real(kind(0d0)):: maxamp,ratl
  complex(kind(0d0)):: g(2)
  real(kind(0d0)):: amp,ampratio

  ampratio = 0.d0
  amp = dsqrt( zabs( g(1) )**2 + zabs( g(2) )**2 )
  if ( amp > maxamp ) maxamp = amp
  if ( (amp /= 0.d0) .and. (maxamp /= 0.d0) ) ampratio = amp / maxamp
  if ( ( ampratio < ratl ) .and. ( l >= lsuf ) ) then
     ismall = ismall + 1
  else
     ismall = 0
  endif

end subroutine calamp

!


subroutine callsuf(omega,nzone,vrmax,vsv,lsuf)
  implicit none
  integer:: nzone,lsuf
  real(kind(0d0)):: omega,vrmax(*),vsv(4,*)
  real(kind(0d0)):: tvs,coef
  integer:: i

  tvs = 0.d0
  do i=1,4
     if (i == 1) then
        coef = 1.d0
     else
        coef = coef
     endif
     tvs = tvs + ( vsv(i,nzone) ) * coef
  enddo
  lsuf = int(omega * vrmax(nzone) / tvs - 0.5d0) + 1
  return
end subroutine callsuf


!

subroutine calra_psv(nlayer,inlayer,jnlayer,jnslay,jnllay,gridpar,dzpar,nzone,vrmin,vrmax,iphase, &
rmin,rmax,nslay,nllay,nnl,re )

  implicit none
  real(kind(0d0)), parameter :: pi=3.1415926535897932d0
  integer :: nlayer,inlayer,jnlayer,jnslay,jnllay
  integer :: nzone, iphase(*), nslay, nllay, nnl(nzone)
  real(kind(0d0)) :: gridpar(*),dzpar(*),vrmin(*),vrmax(*),rmin,rmax
  integer :: izone,itmp,ntmp
  real(kind(0d0)) :: rh,re

  nslay=0
  nllay=0
  inlayer=0

  nnl = 0

  jnlayer=0
  jnslay=0
  jnllay=0


  itmp=1
  do izone=1,nzone
     rh = vrmax(izone) - vrmin(izone)

     if (dzpar(izone) == 0.d0) then
        ntmp = 1
     else
        ntmp = int( sqrt(3.3d0 / re ) * rh / dzpar(izone) / 2.d0 / pi  / 7.d-1 + 1 )
     endif
     ! ntmp (see Geller & Takeuchi 1995 6.2)
     nnl(izone) = ntmp
     if ( nnl(izone) < 5 ) nnl(izone)=5
     if (iphase(izone) == 1) nslay = nslay+nnl(izone)
     if (iphase(izone) == 2) nllay=nllay+nnl(izone)
  enddo

  inlayer = 0
  do izone=1,nzone
     inlayer = inlayer + nnl(izone)
  enddo
  nlayer = inlayer
  jnlayer = jnlayer + inlayer
  jnslay = jnslay + nslay
  jnllay = jnllay + nllay

  return

end subroutine calra_psv

!

subroutine calra2_psv_old(nlayer,gridpar,dzpar,nzone,vrmin,vrmax,rmin,rmax,nnl,ra,re,nsta, rsta, &
rrsta, iista, rs, cista,iphase,istazone,ciista)

  implicit none
  real(kind(0d0)), parameter :: pi=3.1415926535897932d0
  integer :: nlayer
  integer :: nzone,nnl(nzone)
  real(kind(0d0)) :: gridpar(*),dzpar(*),vrmin(*),vrmax(*),rmin,rmax,ra(nlayer+nzone+1)
  integer :: izone,itmp,i,nsta,ista, iista(1:3,1:nsta), cista
  real(kind(0d0)) :: rsta(1:nsta),rrsta(1:3,1:nsta)
  real(kind(0d0)) :: rh,re,rs,ctmp
  real(kind(0d0)) :: chikasa
  integer :: iphase(*),istazone(1:nsta),ciista

  chikasa = 0.d0

  ra = 0
  ra(1) = rmin
  ciista = 0
  ctmp = 6371.d0

  itmp = 1
  do izone=1,nzone
     rh = vrmax(izone) - vrmin(izone)
     do i=1,nnl(izone)
        itmp = itmp + 1
        ra(itmp) = vrmin(izone) &
             &  + rh * dble(i) / dble( nnl(izone) )
     enddo
  enddo


  itmp = 1
  do izone=1,nzone
     rh = vrmax(izone) - vrmin(izone)
     do i=1,nnl(izone)
        do ista = 1, nsta
           if ( (ra(itmp) < rsta(ista)) .and. (rsta(ista) <= ra(itmp+1))) then
              if (i /= nnl(izone)) then
                 istazone(ista) = izone
                 if (iphase(istazone(ista)) == 2) stop 'rsta is in liquid layer (calra2_psv)'
                 rrsta(1,ista) = ra(itmp)
                 rrsta(2,ista) = ra(itmp+1)
                 rrsta(3,ista) = ra(itmp+2)

                 iista(1,ista) = i
                 iista(2,ista) = i + 1
                 iista(3,ista) = i + 2
              else
                 istazone(ista) = izone
                 if (iphase(istazone(ista)) == 2) stop 'rsta is in liquid layer (calra2_psv)'
                 rrsta(1,ista) = ra(itmp-1)
                 rrsta(2,ista) = ra(itmp)
                 rrsta(3,ista) = ra(itmp+1)

                 iista(1,ista) = i - 1
                 iista(2,ista) = i
                 iista(3,ista) = i + 1
              endif

              if ((abs(rs-rsta(ista)) < ctmp) .and. (abs(rs-rsta(ista)) >= chikasa)) then
                 cista = ista
                 ciista = itmp
                 ctmp = abs(rs-rsta(ista))

              endif
           endif
        enddo
        itmp = itmp + 1
     enddo
  enddo

  if (cista == 0) cista = 1
end subroutine calra2_psv_old

subroutine calra2_psv(nlayer,gridpar,dzpar,nzone,vrmin,vrmax,rmin,rmax,nnl,ra,re,nsta, &
  rsta, rrsta, iista, rs, cista,iphase,istazone,ciista,updown)

  implicit none

  real(kind(0d0)), parameter :: pi=3.1415926535897932d0

  integer :: izone,itmp,i,nsta,ista,cista
  integer :: nlayer
  integer :: nzone,nnl(nzone),updown(1:nsta)
  real(kind(0d0)) :: gridpar(*),dzpar(*),vrmin(*),vrmax(*),rmin,rmax,ra(nlayer+nzone+1)
  integer :: iista(1:3,1:nsta)
  real(kind(0d0)) :: rsta(1:nsta),rrsta(1:3,1:nsta)
  real(kind(0d0)) :: rh,re,rs,ctmp
  real(kind(0d0)) :: chikasa
  integer :: iphase(*),istazone(1:nsta),ciista

  chikasa = 0.d0

  ra = 0
  ra(1) = rmin
  ciista = 0
  ctmp = 6371.d0

  itmp = 1
  do izone=1,nzone
     rh = vrmax(izone) - vrmin(izone)
     do i=1,nnl(izone)
        itmp = itmp + 1
        ra(itmp) = vrmin(izone) &
             &  + rh * dble(i) / dble( nnl(izone) )
     enddo
  enddo


  itmp = 1
  do izone=1,nzone
     rh = vrmax(izone) - vrmin(izone)
     do i=1,nnl(izone)
        do ista = 1, nsta
           if ( (ra(itmp) < rsta(ista)) .and. (rsta(ista) <= ra(itmp+1))) then
              if (i /= nnl(izone)) then
                 istazone(ista) = izone
                 if (iphase(istazone(ista)) == 2) stop 'rsta is in liquid layer (calra2_psv)'
                 rrsta(1,ista) = ra(itmp)
                 rrsta(2,ista) = ra(itmp+1)
                 rrsta(3,ista) = ra(itmp+2)

                 iista(1,ista) = i
                 iista(2,ista) = i + 1
                 iista(3,ista) = i + 2
              else
                 istazone(ista) = izone
                 if (iphase(istazone(ista)) == 2) stop 'rsta is in liquid layer (calra2_psv)'
                 rrsta(1,ista) = ra(itmp-1)
                 rrsta(2,ista) = ra(itmp)
                 rrsta(3,ista) = ra(itmp+1)

                 iista(1,ista) = i - 1
                 iista(2,ista) = i
                 iista(3,ista) = i + 1
              endif

              if ((abs(rs-rsta(ista)) < ctmp) .and. (abs(rs-rsta(ista)) >= chikasa)) then
                 cista = ista
                 ciista = itmp
                 ctmp = abs(rs-rsta(ista))

              endif
           endif
        enddo
        itmp = itmp + 1
     enddo
  enddo


  ! Nobuaki 2011/12/20 from here
  itmp = 1
  do izone=1,nzone
     rh = vrmax(izone) - vrmin(izone)
     do i=1,nnl(izone)
        do ista = 1, nsta
           if (updown(ista) /= -1) then
              if ( (ra(itmp) <= rsta(ista)) .and. (rsta(ista) < ra(itmp+1))) then

                 istazone(ista) = izone
                 if (iphase(istazone(ista)) == 2) stop 'rsta is in liquid layer (calra2_psv)'
                 rrsta(1,ista) = ra(itmp)
                 rrsta(2,ista) = ra(itmp+1)
                 rrsta(3,ista) = ra(itmp+2)

                 iista(1,ista) = i
                 iista(2,ista) = i + 1
                 iista(3,ista) = i + 2



                 if ((abs(rs-rsta(ista)) < ctmp) .and. (abs(rs-rsta(ista)) >= chikasa)) then
                    cista = ista
                    ciista = itmp
                    ctmp = abs(rs-rsta(ista))
                 endif

              endif
           endif
        enddo
        itmp = itmp + 1
     enddo
  enddo
  ! Nobuaki 2011/12/20 until here

  if (cista == 0) cista = 1
end subroutine calra2_psv

!
subroutine calcutd(nzone,nnlayer,nnl,tmpc,rat,nn,iphase,spo,spn, ra,kkdr,kc)
  implicit none
  integer :: nzone,nn,nnlayer,spn,kkdr(1:nzone),kc,iphase(1:nzone),nnl(1:nzone)
  complex(kind(0d0)) :: tmpc(1:nn)
  real(kind(0d0)) :: rat,spo,ra(1:nnlayer+nzone+1)
  integer :: nc
  real(kind(0d0)) :: cU(nn),cV(nn),rc
  real(kind(0d0)) :: maxamp,amp(nn)
  integer :: iz,jz,jj,i,ml(nzone),tzone

  do jj=1,nn
     cU(jj) = 0.d0
     cV(jj) = 0.d0
  enddo
  iz = 2
  jz = 1
  do jj=1,nn
     if (iz <= nzone) then
        if (jj == kkdr(iz)) then
           if (iphase(iz) /= iphase(iz-1)) jz = jz - 1
           iz = iz + 1
        endif
     endif
     if (iphase(iz-1) == 1) then
        if (mod((jj-kkdr(iz-1)),2) == 1) then ! U
           cU(jz) = cdabs(tmpc(jj))
           jz = jz + 1
        else    ! V
        endif
     else ! U in fluid
        cU(jz) = cdabs(tmpc(jj))
        jz = jz + 1
     endif
  enddo

  maxamp = -1.d0
  do i=1,jz-1
     amp(i) = cU(i)
     if (maxamp < amp(i)) maxamp = amp(i)
  enddo
!
  maxamp = maxamp * rat ! threshold value
!
  nc = 1
  do i=1,jz-1
     if (amp(i) > maxamp) then
        nc = i
        cycle
     endif
  enddo
  i = 1
  do jj=1,nzone
     i = i + nnl(jj)
     ml(jj) = i
  enddo
  do jj=nzone,1,-1
     if (ml(jj) > nc) tzone = jj
  enddo
  rc = ra(nc)

  do i=1,jz-1
     if ( (ra(i) <= rc) .and. (rc < ra(i+1)) ) then
        nc = i
        if (tzone == 1) then ! case(tzone is innermost zone)
           if (iphase(tzone) == 1) kc = 1 + 2 * nc
           if (iphase(tzone) == 2) kc = 1 + nc
        else
           if (iphase(tzone) == 1) then
              kc = kkdr(tzone) + 2 * (nc - ml(tzone-1))
           endif
           if (iphase(tzone) == 2) then
              kc = kkdr(tzone) + nc - ml(tzone-1)
           endif
        endif
     endif
  enddo

  return
end subroutine calcutd




subroutine translat(geodetic,geocentric)

  implicit none
  real(kind(0d0)),parameter ::  flattening = 1.d0 / 298.25d0
  real(kind(0d0)), parameter :: pi = 3.1415926535897932d0
  real(kind(0d0)) :: geocentric, geodetic
  integer :: flag
  flag = 0
  if (geodetic > 90.d0) then
     geodetic = 1.8d2 - geodetic
     flag = 1
  endif

  geodetic = geodetic / 1.8d2 * pi
  geocentric = datan((1.d0-flattening)*(1.d0-flattening)* dtan(geodetic) )
  geocentric = geocentric * 1.8d2 / pi

  if (flag == 1) then
     geocentric = 1.8d2 - geocentric
  endif

  return
end subroutine translat


subroutine calthetaphi(ievla,ievlo,istla,istlo,theta,phi)

  implicit none
  real(kind(0d0)), parameter:: pi = 3.1415926535897932d0

  real(kind(0d0)) ::  ievla,ievlo,istla,istlo
  real(kind(0d0)) ::  evla,evlo,stla,stlo
  real(kind(0d0)) :: theta,phi
  real(kind(0d0)) :: gcarc,az
  real(kind(0d0)) :: tc,ts

  ! transformation to spherical coordinates

  evla = 90.d0 - ievla
  stla = 90.d0 - istla

  evla = evla / 1.8d2 * pi
  evlo = ievlo / 1.8d2 * pi
  stla = stla / 1.8d2 * pi
  stlo = istlo / 1.8d2 * pi

  gcarc = dacos( dcos(evla) * dcos(stla) + dsin(evla) * dsin(stla) * dcos(evlo - stlo) )

  tc = (dcos(stla)*dsin(evla)-dsin(stla)*dcos(evla)*dcos(stlo-evlo))/dsin(gcarc)
  ts = dsin(stla) * dsin(stlo - evlo) / dsin(gcarc)

  az = dacos(tc)
  if ( ts < 0.d0 ) az = -1.d0 * az

  az = az * 1.8d2 / pi

  gcarc = gcarc * 1.8d2 / pi

  theta = gcarc
  phi   = 180.d0 - az
  return
end subroutine calthetaphi




subroutine calstg4onedepth( maxnlay,maxnzone,nzone,vrmin,vrmax,iphase,rrho,vpv,vph,vsv,vsh,eta,rmax,r, &
updown,ecA,ecC,ecF,ecL,ecN,spn)

  ! Computing the structure grid points.
  implicit none
  integer:: maxnlay,maxnzone,nzone,iphase(*),updown
  real(kind(0d0)):: rrho(4,*),vpv(4,*),vph(4,*),vsv(4,*),vsh(4,*),eta(4,*),vrmin(*),vrmax(*)
  real(kind(0d0)):: rmax,r
  real(kind(0d0)):: ecL,ecN
  real(kind(0d0)):: ecA,ecC,ecF
  real(kind(0d0)):: trho,tvpv,tvph,tvsv,tvsh,teta,coef
  real(kind(0d0)):: epsillon
  integer:: j,itmp,jtmp,spn

  epsillon = 0.1d0
  ! computing the structure grid points
  itmp = 0
  jtmp = 0
!  spn = 0

  !write(100,*)
  !write(100,*) 'izone ',r
  !write(100,*) vrmin(1:nzone)
  !write(100,*) vrmax(1:nzone)

  !do izone=1,nzone
  !   write(100,*) izone,vrmin(izone),vrmax(izone)

  !   if ( (r - vrmin(izone) >= 0.d0) .and. ( (vrmax(izone) - epsillon) >r )) then
  !      spn = izone
  !      write(100,*) ' test OK , spn ',spn
  !      write(100,*) vrmax(spn)
  !      write(100,*) r
  !      write(100,*) vrmin(spn)
  !   endif

  !enddo
 ! write(100,*) 'rayon r :', vrmin(spn),r,vrmax(spn)
  !if (vrmax(nzone)==r) spn = nzone
!!  if ((vrmin(spn)==r) .and. (updown==-1)) then
!!  if (updown==-1) then
!!     spn = spn - 1
!!  endif



  trho = 0.d0
  tvpv = 0.d0
  tvph = 0.d0
  tvsv = 0.d0
  tvsh = 0.d0
  teta = 0.d0
  do j=1,4
     if ( j == 1 ) then
        coef = 1.d0
     else
        coef = coef * ( r / rmax )
     endif
     trho  = trho  + rrho(j,spn)  * coef
     tvpv  = tvpv  + vpv(j,spn)   * coef
     tvph  = tvph  + vph(j,spn)   * coef
     tvsv  = tvsv  + vsv(j,spn)   * coef
     tvsh  = tvsh  + vsh(j,spn)   * coef
     teta  = teta  + eta(j,spn)   * coef
  enddo
  ecL = trho * tvsv * tvsv
  ecN = trho * tvsh * tvsh
  ecA = trho * tvph * tvph
  ecC = trho * tvpv * tvpv
  ecF = teta * ( ecA - 2.d0 * ecL )

  return
end subroutine calstg4onedepth


subroutine ReadFrqsByProc(Ifr,Ifr2,para,n)
  integer k,i,n,Ifr(n),Ifr2(0:n-2,2),para

  open(20,file='FrqsMpi.txt')
  read(20,*) para
  if (para == 1) then
   do i=1,n
     read(20,*) Ifr(i)
   enddo
  else
   do i=1,n-1
     read(20,*) k,Ifr2(k,1),Ifr2(k,2)
     write(*,*) k,Ifr2(k,1),Ifr2(k,2)
   enddo
  endif
  close(20)

end subroutine ReadFrqsByProc



subroutine epitra(SourceLat,SourceLong,RecLat,RecLong,delta,azi,bazi)
!.....                                                        ..........
!                                                                      .
! *** Transforms from spherical to source coordinates; angles          .
!      are given in degrees.                                           .
!      tets ...      latitude (theta) of the source                    .
!      phis ...      longitude (phi) of the source                     .
!      tete ...      latitude of the receiver                          .
!      phie ...      longitude of the receiver                         .
!      delta ...      epicentral distance                              .
!      azi ...      azimuth, measured clockwise from north             .
!      bazi ...     back-azimuth, measured clockwise from north        .
!                                                                      .
!.......................................................................
!     imicit none
  real(kind(0d0))  delta,azi,SourceLat,RecLat,SourceLong,RecLong,bazi
  real(kind(0d0)) pi,tete,phie,tets,phis,ctete,ctets,stete,stets,cdel, &
       sdelta,cal,sal,gamma,d2r,r2d,cbe,sbe
  parameter(pi=3.14159265359d0,d2r=Pi/180.d0,r2d=1.d0/d2r)


!-----------------------------------------------------------------
!  Transform from latitude to colatitude                          |
! Transform from degrees to radians and abbreviate sine and       |
! cosine.                                                         |
!-----------------------------------------------------------------

  tets = 90.d0 - dble(SourceLat)
  tete = 90.d0 - dble(RecLat)
  if (tets == 0.d0) then
     delta = sngl(tete)
     azi   = 180.d0-Reclong
     bazi = 0.d0
     return
  else if (tets == 180.d0 ) then
     delta = 180.d0 - (tete)
     azi   = RecLong
     bazi  = 180.d0
     return
  else if (tete == 0.d0 ) then
     delta = (tets)
     azi = 0.d0
     bazi = 180.d0
     return
  endif
  tets = tets*d2r
  phis = dble(SourceLong)*d2r
  tete = tete*d2r
  phie = dble(RecLong)*d2r
  ctets = dcos(tets)
  ctete = dcos(tete)
  stets = dsin(tets)
  stete = dsin(tete)

!-----------------------------------------------------------------
! Use cosine theorem on the sphere and check for antipode.        |
!-----------------------------------------------------------------

  cdel  = ctets*ctete+stets*stete*dcos(phie-phis)
  delta = (dacos(cdel)*r2d)
  sdelta=dsqrt((1.-cdel)*(1.+cdel))
  if (sdelta == 0.d0) then
     azi=0.d0
     bazi=0.d0
     return
  endif

!-----------------------------------------------------------------
! Use cosine and Sine theorem to get azimut and back-azimut.      |
!-----------------------------------------------------------------

  cal=(ctete-ctets*cdel)/(stets*sdelta)
  sal=stete*dsin(phie-phis)/sdelta
  if (cal > 1.d0) cal=1.d0
  if (cal < -1.d0) cal=-1.d0
  gamma=dacos(cal)
  if (sal >= 0.d0) then
     azi=(gamma)
  else
     azi=(2.d0*pi-gamma)
  endif
  azi = azi*(r2d)


  cbe=(ctets-ctete*cdel)/(stete*sdelta)
  sbe=stets*dsin(phie-phis)/sdelta
  if (cbe > 1.d0) cbe=1.d0
  if (cbe < -1.d0) cbe=-1.d0
  gamma=dacos(cbe)
  if (sbe >= 0.d0) then
     bazi=(2.d0*pi-gamma)
  else
     bazi=(gamma)
  endif
  bazi = bazi*(r2d)

  azi = 180.d0 - azi

end subroutine epitra

!-----------------------------------------------------------------

 subroutine DistribSta(n,n_global,myrank,nbproc,irmin,irmax)
   implicit none
   integer coef,n,n_global,myrank,nbproc,irmin,irmax

   n = int(n_global/nbproc)

   coef = mod(n_global,nbproc)

   if (myrank < mod(n_global,nbproc) .and. mod(n_global,nbproc) /= 0) then
      n = n + 1
      coef = 0
   endif

   irmin = (myrank ) * n + 1 + coef
   irmax = (myrank + 1 ) * n + coef

 end subroutine DistribSta

!-----------------------------------------------------------------

  subroutine write_a_block_of_coefs_to_disk(tabg0,tabg0der,tabg0_small_buffer,ir_,llog0,maxlmax_g,r_n)

  implicit none

  integer, intent(in) :: ir_,llog0,maxlmax_g,r_n

  complex(kind(0d0)), dimension(maxlmax_g,2,6,-2:2,r_n), intent(in) :: tabg0,tabg0der

!! DK DK here I purposely declare the buffer with a size llog0 <= maxlmax_g,
!! DK DK this way I can then use a simple write statement to write the whole array to disk
  complex(kind(0d0)), dimension(llog0,2,6,-2:2) :: tabg0_small_buffer

  integer :: illog

! copy the array to the buffer
!IBM* ASSERT (MINITERCNT(1000))
!DIR$ loop count min(1000)
  do illog = 1,llog0

!   tabg0_small_buffer(illog,:,:,:) = tabg0(illog,:,:,:,ir_)

!! DK DK unrolled the loops to allow for vectorization of the outer loop on "illog"
    tabg0_small_buffer(illog,1,1,-2) = tabg0(illog,1,1,-2,ir_)
    tabg0_small_buffer(illog,2,1,-2) = tabg0(illog,2,1,-2,ir_)
    tabg0_small_buffer(illog,1,2,-2) = tabg0(illog,1,2,-2,ir_)
    tabg0_small_buffer(illog,2,2,-2) = tabg0(illog,2,2,-2,ir_)
    tabg0_small_buffer(illog,1,3,-2) = tabg0(illog,1,3,-2,ir_)
    tabg0_small_buffer(illog,2,3,-2) = tabg0(illog,2,3,-2,ir_)
    tabg0_small_buffer(illog,1,4,-2) = tabg0(illog,1,4,-2,ir_)
    tabg0_small_buffer(illog,2,4,-2) = tabg0(illog,2,4,-2,ir_)
    tabg0_small_buffer(illog,1,5,-2) = tabg0(illog,1,5,-2,ir_)
    tabg0_small_buffer(illog,2,5,-2) = tabg0(illog,2,5,-2,ir_)
    tabg0_small_buffer(illog,1,6,-2) = tabg0(illog,1,6,-2,ir_)
    tabg0_small_buffer(illog,2,6,-2) = tabg0(illog,2,6,-2,ir_)
    tabg0_small_buffer(illog,1,1,-1) = tabg0(illog,1,1,-1,ir_)
    tabg0_small_buffer(illog,2,1,-1) = tabg0(illog,2,1,-1,ir_)
    tabg0_small_buffer(illog,1,2,-1) = tabg0(illog,1,2,-1,ir_)
    tabg0_small_buffer(illog,2,2,-1) = tabg0(illog,2,2,-1,ir_)
    tabg0_small_buffer(illog,1,3,-1) = tabg0(illog,1,3,-1,ir_)
    tabg0_small_buffer(illog,2,3,-1) = tabg0(illog,2,3,-1,ir_)
    tabg0_small_buffer(illog,1,4,-1) = tabg0(illog,1,4,-1,ir_)
    tabg0_small_buffer(illog,2,4,-1) = tabg0(illog,2,4,-1,ir_)
    tabg0_small_buffer(illog,1,5,-1) = tabg0(illog,1,5,-1,ir_)
    tabg0_small_buffer(illog,2,5,-1) = tabg0(illog,2,5,-1,ir_)
    tabg0_small_buffer(illog,1,6,-1) = tabg0(illog,1,6,-1,ir_)
    tabg0_small_buffer(illog,2,6,-1) = tabg0(illog,2,6,-1,ir_)
    tabg0_small_buffer(illog,1,1,0) = tabg0(illog,1,1,0,ir_)
    tabg0_small_buffer(illog,2,1,0) = tabg0(illog,2,1,0,ir_)
    tabg0_small_buffer(illog,1,2,0) = tabg0(illog,1,2,0,ir_)
    tabg0_small_buffer(illog,2,2,0) = tabg0(illog,2,2,0,ir_)
    tabg0_small_buffer(illog,1,3,0) = tabg0(illog,1,3,0,ir_)
    tabg0_small_buffer(illog,2,3,0) = tabg0(illog,2,3,0,ir_)
    tabg0_small_buffer(illog,1,4,0) = tabg0(illog,1,4,0,ir_)
    tabg0_small_buffer(illog,2,4,0) = tabg0(illog,2,4,0,ir_)
    tabg0_small_buffer(illog,1,5,0) = tabg0(illog,1,5,0,ir_)
    tabg0_small_buffer(illog,2,5,0) = tabg0(illog,2,5,0,ir_)
    tabg0_small_buffer(illog,1,6,0) = tabg0(illog,1,6,0,ir_)
    tabg0_small_buffer(illog,2,6,0) = tabg0(illog,2,6,0,ir_)
    tabg0_small_buffer(illog,1,1,1) = tabg0(illog,1,1,1,ir_)
    tabg0_small_buffer(illog,2,1,1) = tabg0(illog,2,1,1,ir_)
    tabg0_small_buffer(illog,1,2,1) = tabg0(illog,1,2,1,ir_)
    tabg0_small_buffer(illog,2,2,1) = tabg0(illog,2,2,1,ir_)
    tabg0_small_buffer(illog,1,3,1) = tabg0(illog,1,3,1,ir_)
    tabg0_small_buffer(illog,2,3,1) = tabg0(illog,2,3,1,ir_)
    tabg0_small_buffer(illog,1,4,1) = tabg0(illog,1,4,1,ir_)
    tabg0_small_buffer(illog,2,4,1) = tabg0(illog,2,4,1,ir_)
    tabg0_small_buffer(illog,1,5,1) = tabg0(illog,1,5,1,ir_)
    tabg0_small_buffer(illog,2,5,1) = tabg0(illog,2,5,1,ir_)
    tabg0_small_buffer(illog,1,6,1) = tabg0(illog,1,6,1,ir_)
    tabg0_small_buffer(illog,2,6,1) = tabg0(illog,2,6,1,ir_)
    tabg0_small_buffer(illog,1,1,2) = tabg0(illog,1,1,2,ir_)
    tabg0_small_buffer(illog,2,1,2) = tabg0(illog,2,1,2,ir_)
    tabg0_small_buffer(illog,1,2,2) = tabg0(illog,1,2,2,ir_)
    tabg0_small_buffer(illog,2,2,2) = tabg0(illog,2,2,2,ir_)
    tabg0_small_buffer(illog,1,3,2) = tabg0(illog,1,3,2,ir_)
    tabg0_small_buffer(illog,2,3,2) = tabg0(illog,2,3,2,ir_)
    tabg0_small_buffer(illog,1,4,2) = tabg0(illog,1,4,2,ir_)
    tabg0_small_buffer(illog,2,4,2) = tabg0(illog,2,4,2,ir_)
    tabg0_small_buffer(illog,1,5,2) = tabg0(illog,1,5,2,ir_)
    tabg0_small_buffer(illog,2,5,2) = tabg0(illog,2,5,2,ir_)
    tabg0_small_buffer(illog,1,6,2) = tabg0(illog,1,6,2,ir_)
    tabg0_small_buffer(illog,2,6,2) = tabg0(illog,2,6,2,ir_)

  enddo
  write(34) ir_,llog0
  write(34) tabg0_small_buffer

! copy the array to the buffer
!IBM* ASSERT (MINITERCNT(1000))
!DIR$ loop count min(1000)
  do illog = 1,llog0

!   tabg0_small_buffer(illog,:,:,:) = tabg0der(illog,:,:,:,ir_)

!! DK DK unrolled the loops to allow for vectorization of the outer loop on "illog"
    tabg0_small_buffer(illog,1,1,-2) = tabg0der(illog,1,1,-2,ir_)
    tabg0_small_buffer(illog,2,1,-2) = tabg0der(illog,2,1,-2,ir_)
    tabg0_small_buffer(illog,1,2,-2) = tabg0der(illog,1,2,-2,ir_)
    tabg0_small_buffer(illog,2,2,-2) = tabg0der(illog,2,2,-2,ir_)
    tabg0_small_buffer(illog,1,3,-2) = tabg0der(illog,1,3,-2,ir_)
    tabg0_small_buffer(illog,2,3,-2) = tabg0der(illog,2,3,-2,ir_)
    tabg0_small_buffer(illog,1,4,-2) = tabg0der(illog,1,4,-2,ir_)
    tabg0_small_buffer(illog,2,4,-2) = tabg0der(illog,2,4,-2,ir_)
    tabg0_small_buffer(illog,1,5,-2) = tabg0der(illog,1,5,-2,ir_)
    tabg0_small_buffer(illog,2,5,-2) = tabg0der(illog,2,5,-2,ir_)
    tabg0_small_buffer(illog,1,6,-2) = tabg0der(illog,1,6,-2,ir_)
    tabg0_small_buffer(illog,2,6,-2) = tabg0der(illog,2,6,-2,ir_)
    tabg0_small_buffer(illog,1,1,-1) = tabg0der(illog,1,1,-1,ir_)
    tabg0_small_buffer(illog,2,1,-1) = tabg0der(illog,2,1,-1,ir_)
    tabg0_small_buffer(illog,1,2,-1) = tabg0der(illog,1,2,-1,ir_)
    tabg0_small_buffer(illog,2,2,-1) = tabg0der(illog,2,2,-1,ir_)
    tabg0_small_buffer(illog,1,3,-1) = tabg0der(illog,1,3,-1,ir_)
    tabg0_small_buffer(illog,2,3,-1) = tabg0der(illog,2,3,-1,ir_)
    tabg0_small_buffer(illog,1,4,-1) = tabg0der(illog,1,4,-1,ir_)
    tabg0_small_buffer(illog,2,4,-1) = tabg0der(illog,2,4,-1,ir_)
    tabg0_small_buffer(illog,1,5,-1) = tabg0der(illog,1,5,-1,ir_)
    tabg0_small_buffer(illog,2,5,-1) = tabg0der(illog,2,5,-1,ir_)
    tabg0_small_buffer(illog,1,6,-1) = tabg0der(illog,1,6,-1,ir_)
    tabg0_small_buffer(illog,2,6,-1) = tabg0der(illog,2,6,-1,ir_)
    tabg0_small_buffer(illog,1,1,0) = tabg0der(illog,1,1,0,ir_)
    tabg0_small_buffer(illog,2,1,0) = tabg0der(illog,2,1,0,ir_)
    tabg0_small_buffer(illog,1,2,0) = tabg0der(illog,1,2,0,ir_)
    tabg0_small_buffer(illog,2,2,0) = tabg0der(illog,2,2,0,ir_)
    tabg0_small_buffer(illog,1,3,0) = tabg0der(illog,1,3,0,ir_)
    tabg0_small_buffer(illog,2,3,0) = tabg0der(illog,2,3,0,ir_)
    tabg0_small_buffer(illog,1,4,0) = tabg0der(illog,1,4,0,ir_)
    tabg0_small_buffer(illog,2,4,0) = tabg0der(illog,2,4,0,ir_)
    tabg0_small_buffer(illog,1,5,0) = tabg0der(illog,1,5,0,ir_)
    tabg0_small_buffer(illog,2,5,0) = tabg0der(illog,2,5,0,ir_)
    tabg0_small_buffer(illog,1,6,0) = tabg0der(illog,1,6,0,ir_)
    tabg0_small_buffer(illog,2,6,0) = tabg0der(illog,2,6,0,ir_)
    tabg0_small_buffer(illog,1,1,1) = tabg0der(illog,1,1,1,ir_)
    tabg0_small_buffer(illog,2,1,1) = tabg0der(illog,2,1,1,ir_)
    tabg0_small_buffer(illog,1,2,1) = tabg0der(illog,1,2,1,ir_)
    tabg0_small_buffer(illog,2,2,1) = tabg0der(illog,2,2,1,ir_)
    tabg0_small_buffer(illog,1,3,1) = tabg0der(illog,1,3,1,ir_)
    tabg0_small_buffer(illog,2,3,1) = tabg0der(illog,2,3,1,ir_)
    tabg0_small_buffer(illog,1,4,1) = tabg0der(illog,1,4,1,ir_)
    tabg0_small_buffer(illog,2,4,1) = tabg0der(illog,2,4,1,ir_)
    tabg0_small_buffer(illog,1,5,1) = tabg0der(illog,1,5,1,ir_)
    tabg0_small_buffer(illog,2,5,1) = tabg0der(illog,2,5,1,ir_)
    tabg0_small_buffer(illog,1,6,1) = tabg0der(illog,1,6,1,ir_)
    tabg0_small_buffer(illog,2,6,1) = tabg0der(illog,2,6,1,ir_)
    tabg0_small_buffer(illog,1,1,2) = tabg0der(illog,1,1,2,ir_)
    tabg0_small_buffer(illog,2,1,2) = tabg0der(illog,2,1,2,ir_)
    tabg0_small_buffer(illog,1,2,2) = tabg0der(illog,1,2,2,ir_)
    tabg0_small_buffer(illog,2,2,2) = tabg0der(illog,2,2,2,ir_)
    tabg0_small_buffer(illog,1,3,2) = tabg0der(illog,1,3,2,ir_)
    tabg0_small_buffer(illog,2,3,2) = tabg0der(illog,2,3,2,ir_)
    tabg0_small_buffer(illog,1,4,2) = tabg0der(illog,1,4,2,ir_)
    tabg0_small_buffer(illog,2,4,2) = tabg0der(illog,2,4,2,ir_)
    tabg0_small_buffer(illog,1,5,2) = tabg0der(illog,1,5,2,ir_)
    tabg0_small_buffer(illog,2,5,2) = tabg0der(illog,2,5,2,ir_)
    tabg0_small_buffer(illog,1,6,2) = tabg0der(illog,1,6,2,ir_)
    tabg0_small_buffer(illog,2,6,2) = tabg0der(illog,2,6,2,ir_)

  enddo
  write(34) tabg0_small_buffer

 end subroutine write_a_block_of_coefs_to_disk

