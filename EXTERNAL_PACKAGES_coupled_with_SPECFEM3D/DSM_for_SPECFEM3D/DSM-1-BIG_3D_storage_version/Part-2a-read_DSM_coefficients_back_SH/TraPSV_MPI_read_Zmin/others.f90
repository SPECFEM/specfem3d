subroutine pinputTra(outputDir,psvmodel,modelname,stationsinf,tlen,imin,imax,r0min,r0max,r0delta,r0lat,r0lon,itranslat)
  implicit none
  character(120), parameter :: tmpfile='tmpworkingfile_for_TraPSV'
  character(120) :: dummy,outputDir,psvmodel,modelname,stationsinf
  real(kind(0d0)) :: tlen,r0min,r0max,r0delta
  real(kind(0d0)) :: r0lat,r0lon
  integer :: imin,imax,itranslat

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
end subroutine pinputTra


subroutine pinputTra_Part2(ifrequ_min,ifrequ_max,inputdir,outputDir,psvmodel,modelname,stationsinf,tlen,imin,imax,r0min,r0max,r0delta,r0lat,r0lon,itranslat)
  implicit none
  character(120) tmpfile
  character(120) :: dummy,outputDir,psvmodel,modelname,stationsinf,inputdir,inputIASP
  real(kind(0d0)) :: tlen,r0min,r0max,r0delta
  real(kind(0d0)) :: r0lat,r0lon,m1,m2,m3,m4,m5,m6
  integer :: imin,imax,itranslat,ifrequ_min,ifrequ_max
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


subroutine pinput(outputDir,psvmodel,modelname,tlen,rmin_,rmax_,rdelta_,r0min,r0max,r0delta,thetamin,thetamax,thetadelta,imin,imax,rsgtswitch)
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




subroutine udertoStress(uder,stress,A,C,F,L,N)
  implicit none
  complex(kind(0d0)) :: uder(1:3,1:3), stress(1:6)
  real(kind(0d0)) :: A,C,F,L,N

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

subroutine udertotsgtSH(imt,uder,ttsgt)
  implicit none
  integer :: imt
  complex(kind(0d0)) :: uder(1:3,1:3),tsgt(1:20),ttsgt(1:10)

  tsgt = cmplx(0.d0)

  tsgt(7)  = ttsgt(1)
  tsgt(12) = ttsgt(2)
  tsgt(13) = ttsgt(3)
  tsgt(14) = ttsgt(4)
  tsgt(15) = ttsgt(5)
  tsgt(16) = ttsgt(6)
  tsgt(17) = ttsgt(7)
  tsgt(18) = ttsgt(8)
  tsgt(19) = ttsgt(9)
  tsgt(20) = ttsgt(10)


  if (imt == 1) then ! rr source
     tsgt(1) = tsgt(1) + uder(1,1)
     tsgt(3) = tsgt(3) - 5.d-1*uder(2,2) - 5.d-1*uder(3,3)
     tsgt(6) = tsgt(6) + uder(1,2) + uder(2,1)
     tsgt(10)= tsgt(10)+ 5.d-1*uder(2,2) - 5.d-1*uder(3,3)
  endif
  if (imt == 2) then ! tt source
     tsgt(2) = tsgt(2) - 5.d-1*uder(1,1)
     tsgt(4) = tsgt(4) + 2.5d-1*uder(2,2) + 2.5d-1*uder(3,3)
     tsgt(8) = tsgt(8) - 5.d-1*uder(1,2) - 5.d-1*uder(2,1)
     tsgt(9) = tsgt(9) + 5.d-1*uder(1,1)
     tsgt(11)= tsgt(11)- 2.5d-1*uder(2,2) + 2.5d-1*uder(3,3)
     tsgt(12)= tsgt(12)- 2.5d-1*uder(2,2) - 2.5d-1*uder(3,3)
     tsgt(17)= tsgt(17)+ 2.5d-1*uder(1,2) + 2.5d-1*uder(2,1)
     tsgt(18)= tsgt(18)- 2.5d-1*uder(1,2) - 2.5d-1*uder(2,1)
     tsgt(19)= tsgt(19)+ 1.25d-1*uder(2,2) - 1.25d-1*uder(3,3)
     tsgt(20)= tsgt(20)+ 1.25d-1*uder(2,2) - 1.25d-1*uder(3,3)
  endif
  if (imt == 3) then ! pp source
     tsgt(2) = tsgt(2) - 5.d-1*uder(1,1)
     tsgt(4) = tsgt(4) + 2.5d-1*uder(2,2) + 2.5d-1*uder(3,3)
     tsgt(8) = tsgt(8) - 5.d-1*uder(1,2) - 5.d-1*uder(2,1)
     tsgt(9) = tsgt(9) - 5.d-1*uder(1,1)
     tsgt(11)= tsgt(11) -2.5d-1*uder(2,2) + 2.5d-1*uder(3,3)
     tsgt(12)= tsgt(12) +2.5d-1*uder(2,2) + 2.5d-1*uder(3,3)
     tsgt(17)= tsgt(17) -2.5d-1*uder(1,2) - 2.5d-1*uder(2,1)
     tsgt(18)= tsgt(18) +2.5d-1*uder(1,2) + 2.5d-1*uder(2,1)
     tsgt(19)= tsgt(19) -1.25d-1*uder(2,2) + 1.25d-1*uder(3,3)
     tsgt(20)= tsgt(20) -1.25d-1*uder(2,2) + 1.25d-1*uder(3,3)
  endif
  if (imt == 4) then ! rt source
     tsgt(5) = tsgt(5) - uder(1,1)
     tsgt(7) = tsgt(7) + 5.d-1*uder(2,2) + 5.d-1*uder(3,3)
     tsgt(13)= tsgt(13)- 5.d-1*uder(1,2) - 5.d-1*uder(2,1)
     tsgt(14)= tsgt(14)+ 5.d-1*uder(1,2) + 5.d-1*uder(2,1)
     tsgt(15)= tsgt(15)- 2.5d-1*uder(2,2) + 2.5d-1*uder(3,3)
     tsgt(16)= tsgt(16)+ 2.5d-1*uder(2,2) - 2.5d-1*uder(3,3)
  endif
  if (imt == 5) then ! rp source
     tsgt(13)= tsgt(13)+ 5.d-1*uder(1,3) + 5.d-1*uder(3,1)
     tsgt(14)= tsgt(14)+ 5.d-1*uder(1,3) + 5.d-1*uder(3,1)
     tsgt(15)= tsgt(15)+ 5.d-1*uder(2,3) + 5.d-1*uder(3,2)
     tsgt(16)= tsgt(16)+ 5.d-1*uder(2,3) + 5.d-1*uder(3,2)
  endif
  if (imt == 6) then ! tp source
     tsgt(17)= tsgt(17)- 5.d-1*uder(1,3) - 5.d-1*uder(3,1)
     tsgt(18)= tsgt(18)- 5.d-1*uder(1,3) - 5.d-1*uder(3,1)
     tsgt(19)= tsgt(19)- 5.d-1*uder(2,3) - 5.d-1*uder(3,2)
     tsgt(20)= tsgt(20)+ 5.d-1*uder(2,3) + 5.d-1*uder(3,2)
  endif

  ttsgt(1) = tsgt(7)
  ttsgt(2) = tsgt(12)
  ttsgt(3) = tsgt(13)
  ttsgt(4) = tsgt(14)
  ttsgt(5) = tsgt(15)
  ttsgt(6) = tsgt(16)
  ttsgt(7) = tsgt(17)
  ttsgt(8) = tsgt(18)
  ttsgt(9) = tsgt(19)
  ttsgt(10)= tsgt(20)
end subroutine udertotsgtSH

!

subroutine locallyCartesianDerivatives (u,udr,udt,udp,uder,r,theta)
  implicit none
  complex(kind(0d0)):: u(1:3),udr(1:3),udt(1:3),udp(1:3),uder(1:3,1:3)
  real(kind(0d0)) :: r,theta
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
  return
end subroutine locallyCartesianDerivatives


subroutine normalisetoKM(u,r)
  implicit none
  complex(kind(0d0)) :: u(1:3)
  integer :: i
  real(kind(0d0)) :: r
  do i = 1,3
     u(i) = u(i) / cmplx(r)
  enddo
  return
end subroutine normalisetoKM



subroutine calgrid( nzone,vrmin,vrmax,vs,rmin,rmax, &
     & imax,lmin,tlen,vmin,gridpar,dzpar )
  implicit none
  real(kind(0d0)), parameter :: pi=3.1415926535897932d0
  integer :: nzone,imax,lmin
  real(kind(0d0)) :: vrmin(*),vrmax(*),vs(4,*)
  real(kind(0d0)) :: rmin,rmax,tlen,vmin(*),gridpar(*),dzpar(*)
  integer :: izone,i,j
  real(kind(0d0)) :: coef1,coef2,v(4),vs1,vs2,rh,omega,amax,gtmp
  do izone=1,nzone
!     computing the S-velocity at each zone
     do i=1,4
        v(i) = vs(i,izone)
     enddo
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
     !     computing rh
     rh = vrmax(izone) - vrmin(izone)
     !      computing omega,amax
     omega = 2.d0 * pi * dble(imax) / tlen
     if ( vs1 >= vs2 ) then
        vmin(izone) = vs2
     else
        vmin(izone) = vs1
     endif
     amax = vrmax(izone)
     gtmp = ( omega * omega ) / ( vmin(izone) * vmin(izone) ) &
          &        - ( (dble(lmin)+0.5d0) * (dble(lmin)+0.5d0) )&
          &        / ( amax * amax )
     if ( gtmp > 0.d0 ) then
        dzpar(izone)   = dsqrt( 1.d0/gtmp )
        gridpar(izone) = rh / dzpar(izone)
     else
        dzpar(izone)   = 0.d0
        gridpar(izone) = 0.d0
     endif
  enddo
  !     rearangement of gridpar
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

  !     re-rearangement of gridpar
  gtmp = 0.d0
  do izone=1,nzone
     gtmp = gtmp + gridpar(izone)
  enddo
  do izone=1,nzone
     gridpar(izone) = gridpar(izone) / gtmp
  enddo


end subroutine calgrid

!-----------------------------------------------------------------------------

subroutine calra(nlayer,gridpar,dzpar,nzone,vrmin,vrmax,rmin,rmax,nnl,re )

  implicit none
  real(kind(0d0)), parameter :: pi=3.1415926535897932d0
  integer :: nlayer,inlayer
  integer :: nzone, nnl(nzone)
  real(kind(0d0)) :: gridpar(*),dzpar(*),vrmin(*),vrmax(*),rmin,rmax
  integer :: izone,itmp,ntmp
  real(kind(0d0)) :: rh,re

  inlayer = 0
  itmp = 1
  do izone=1,nzone
     rh = vrmax(izone) - vrmin(izone)

     if (dzpar(izone) == 0.d0) then
        ntmp = 1
     else
        ntmp = int( sqrt(3.3d0 / re ) * rh / dzpar(izone) &
             &                    / 2.d0 / pi  / 7.d-1 + 1 )
     endif
     ! ntmp (see Geller & Takeuchi 1995 6.2)
     nnl(izone) = ntmp
     if ( nnl(izone) < 5 ) nnl(izone)=5
  enddo

  inlayer = 0
  do izone=1,nzone
     inlayer = inlayer + nnl(izone)
  enddo
  nlayer = inlayer

end subroutine calra


!-----------------------------------------------------------------------------

subroutine calra2(nlayer,gridpar,dzpar,nzone,vrmin,vrmax,rmin,rmax,nnl,ra,re,nsta,rsta,rrsta, iista,updown)

  implicit none
  real(kind(0d0)), parameter :: pi=3.1415926535897932d0
  integer :: nlayer
  integer :: nzone,nnl(nzone)
  real(kind(0d0)) :: gridpar(*),dzpar(*),vrmin(*),vrmax(*),rmin,rmax,ra(nlayer+nzone+1)
  integer :: izone,itmp,i,nsta,ista, iista(1:3,1:nsta), ciista,updown(1:nsta)
  real(kind(0d0)) :: rsta(1:nsta),rrsta(1:3,1:nsta)
  real(kind(0d0)) :: rh,re,ctmp
  !real(kind(0d0)) :: chikasa

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

!do ii =1,itmp-1,1
!if (ra(ii) == ra(ii+1) .or. ra(ii) <= 0) write(*,*) ii,'th ra element equal to next ra elment',ra(ii:ii+1)
!enddo

!write(*,*) 'check for the last a few ra elements'
!do ii=itmp,nlayer+nzone+1,1
!write(*,*) ii,'th radiu of the model is:',ra(ii)
!enddo

  itmp = 1
  do izone=1,nzone
     rh = vrmax(izone) - vrmin(izone)
     do i=1,nnl(izone)
        do ista = 1, nsta
           if ( (ra(itmp) < rsta(ista)) .and. (rsta(ista) <= ra(itmp+1))) then
              if (i /= nnl(izone)) then
                 rrsta(1,ista) = ra(itmp)
                 rrsta(2,ista) = ra(itmp+1)
                 rrsta(3,ista) = ra(itmp+2)

                 iista(1,ista) = itmp
                 iista(2,ista) = itmp + 1
                 iista(3,ista) = itmp + 2
              else
                 rrsta(1,ista) = ra(itmp-1)
                 rrsta(2,ista) = ra(itmp)
                 rrsta(3,ista) = ra(itmp+1)

                 iista(1,ista) = itmp - 1
                 iista(2,ista) = itmp
                 iista(3,ista) = itmp + 1
             ! write(*,*) 'i == nnl(izone)',izone,rsta(ista),rrsta(2,ista),ii
              endif

              !if ((abs(rs-rsta(ista)) < ctmp) .and. (abs(rs-rsta(ista)) >= chikasa)) then
              !   ciista = ista
              !   ctmp = abs(rs-rsta(ista))
              !endif
           endif
        enddo
        itmp = itmp + 1
     enddo
  enddo


  itmp = 1
  do izone=1,nzone
     rh = vrmax(izone) - vrmin(izone)
     do i=1,nnl(izone)
        do ista = 1, nsta
         if (updown(ista) /= -1) then
           if ( (ra(itmp) < rsta(ista)) .and. (rsta(ista) <= ra(itmp+1))) then
                 rrsta(1,ista) = ra(itmp)
                 rrsta(2,ista) = ra(itmp+1)
                 rrsta(3,ista) = ra(itmp+2)

                 iista(1,ista) = itmp
                 iista(2,ista) = itmp + 1
                 iista(3,ista) = itmp + 2

           endif
         endif

              !if ((abs(rs-rsta(ista)) < ctmp) .and. (abs(rs-rsta(ista)) >= chikasa)) then
              !   ciista = ista
              !   ctmp = abs(rs-rsta(ista))
              !endif
        enddo
        itmp = itmp + 1
     enddo
  enddo

  if (ciista == 0) ciista = 1

end subroutine calra2


!-----------------------------------------------------------------------------


subroutine calsp (ndc, nlayer, isp, jsp)
  implicit none
  integer :: ndc,nlayer(*)
  integer :: isp(*),jsp(*)
  integer :: i

  ! computation of isp,jsp,ksp,lsp
  isp(1) = 1
  jsp(1) = 1
  do i=1,ndc
     isp(i+1) = isp(i) + nlayer(i)
     jsp(i+1) = jsp(i) + 4 * nlayer(i)
  enddo

  return
end subroutine calsp

!-----------------------------------------------------------------------------

subroutine calspo( ndc,rdc,nlayer,r0,rmin,rmax,ra, isp,spo,spn )
  ! computation de la source
  implicit none
  integer :: ndc,nlayer,isp(*),spn
  real(kind(0d0)) :: rdc(*),r0,rmin,rmax,ra(*),spo
  integer :: itmp


  ! checquer des parameters
  if ( (r0 < rmin) .or. (r0 > rmax) ) &
       & pause 'The source location is improper.(calspo)'
  ! computing 'spo'
  if ( r0 == rmax ) then
     spo = dble(nlayer) - 0.01d0
     r0 = ra(nlayer) &
          & + (spo-dble(nlayer-1)) * ( ra(nlayer+1)-ra(nlayer) )
  else
     do itmp = 2, ndc+nlayer+2
        if ( r0 < ra(itmp) ) exit
     enddo
     spo = dble(itmp-2) + ( r0-ra(itmp-1) ) / ( ra(itmp)-ra(itmp-1) )
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
  do itmp = 1, ndc+1
     spn = spn + 1
     if ( r0 <= rdc(itmp) ) exit
  enddo
  ! changing 'spo'
  spo = spo - dble( isp(spn) - 1 )



end subroutine calspo


!-----------------------------------------------------------------------------


subroutine calgra( isp,ra,r0,spn,spo,gra )

  integer :: isp(*),spn,itmp
  real(kind(0d0)) :: ra(*),r0,spo,gra(*)

  itmp = isp(spn) + dint( spo )
  gra(1) = ra(itmp)
  gra(2) = r0
  gra(3) = ra(itmp+1)

end subroutine calgra

!-----------------------------------------------------------------------------

subroutine calstg( nzone,rrho,vsv,vsh,nlayer,nnl,ra,rmax,vnp,vra,rho,ecL,ecN )
  implicit none
  integer:: nzone,nlayer,nnl(nzone),vnp
  real(kind(0d0)) :: rrho(4,nzone),vsv(4,nzone),vsh(4,nzone),ra(nlayer+nzone+1),rmax
  real(kind(0d0)) :: vra(nlayer+2*nzone+1),rho(nlayer+2*nzone+1)
  real(kind(0d0)) :: ecL(nlayer+2*nzone+1),ecN(nlayer+2*nzone+1)
  real(kind(0d0)) :: trho,tvsv,tvsh,coef
  integer :: izone,i,j,itmp,jtmp

  vra = 0
  rho = 0
  ecL = 0
  ecN = 0

  itmp = 0
  jtmp = 0
  do izone=1,nzone
     do i=1,nnl(izone)+1
        itmp = itmp + 1
        jtmp = jtmp + 1
        vra(itmp) = ra(jtmp)
        ! --- evaluating the density and elastic constants at this point
        trho = 0.d0
        tvsv = 0.d0
        tvsh = 0.d0
        do j=1,4
           if ( j == 1 ) then
              coef = 1.d0
           else
              coef = coef * ( vra(itmp) / rmax )
           endif
           trho = trho + rrho(j,izone) * coef
           tvsv  = tvsv  + vsv(j,izone)   * coef
           tvsh  = tvsh  + vsh(j,izone)   * coef
        enddo
        rho(itmp) = trho
        ecL(itmp)  = rho(itmp) * tvsv * tvsv
        ecN(itmp)  = rho(itmp) * tvsh * tvsh

     enddo
     jtmp = jtmp - 1
  enddo
  vnp = itmp
end subroutine calstg


!-----------------------------------------------------------------------------

subroutine calgstg(nzone,nlayer,spn,rrho,vsv,vsh, ra,vra,rmax,rho,ecL,ecN,r0,mu0 )

  implicit none
  integer :: nlayer, nzone
  integer :: spn
  real(kind(0d0)) :: rrho(4,nzone),vsv(4,nzone),vsh(4,nzone)
  real(kind(0d0)) :: ra(3),rmax
  real(kind(0d0)) :: vra(3),rho(3)
  real(kind(0d0)) :: ecL(3),ecN(3)
  real(kind(0d0)) :: r0,mu0
  real(kind(0d0)) :: trho,tvsv,tvsh,coef
  integer :: i,j

  vra = 0.d0
  rho = 0.d0
  ecL = 0.d0
  ecN = 0.d0

  do i=1,3
     vra(i) = ra(i)
     trho = 0.d0
     tvsv = 0.d0
     tvsh = 0.d0
     do j=1,4
        if ( j == 1 ) then
           coef = 1.d0
        else
           coef = coef * ( vra(i) / rmax )
        endif
        trho = trho + rrho(j,spn) * coef
        tvsv  = tvsv  + vsv(j,spn)   * coef
        tvsh  = tvsh  + vsh(j,spn)   * coef
     enddo
     rho(i) = trho
     ecL(i)  = rho(i) * tvsv * tvsv
     ecN(i)  = rho(i) * tvsh * tvsh

  enddo

  mu0 = ecL(2)

  return

end subroutine calgstg



!-----------------------------------------------------------------------------

subroutine callsuf(omega,nzone,vrmax,vsv,lsuf)

  implicit none
  integer :: nzone, lsuf
  real(kind(0d0)) :: omega, vrmax(*), vsv(4,*)
  real(kind(0d0)) :: tvs, coef
  integer :: i

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
end subroutine callsuf


!-----------------------------------------------------------------------------


subroutine calcoef( nzone,omega,q,coef )
  implicit none
  real(kind(0d0)), parameter ::  pi = 3.1415926535897932d0
  integer :: izone,nzone
  real(kind(0d0)) :: omega,q(*)
  complex(kind(0d0)) :: coef(*)
  real(kind(0d0)) :: aa,bb

  do izone=1,nzone
     if (q(izone) <= 0.d0) then
        coef(izone) = dcmplx(1.d0)
     else
        if ( omega == 0.d0 ) then
           aa = 1.d0
        else
           aa = 1.d0 + dlog( omega / ( 2.d0 * pi ) ) / ( pi * q(izone) )
        endif
        bb = 1.d0 / ( 2.d0 * Q(izone) )
        coef(izone) = dcmplx( aa, bb ) * dcmplx( aa, bb )
     endif
  enddo
end subroutine calcoef


!-----------------------------------------------------------------------------

subroutine calamp(g,l,lsuf,maxamp,ismall,ratl)

  implicit none
  integer :: l,lsuf,ismall
  real(kind(0d0)) :: maxamp,ratl
  complex(kind(0d0)) :: g
  real(kind(0d0)) :: amp,ampratio
  ampratio = 0.d0
  amp = abs(g)
  if ( amp > maxamp ) maxamp = amp
  if ( (amp /= 0.d0) .and. (maxamp /= 0.d0) ) then
     ampratio = amp / maxamp
  endif
  if ( (ampratio < ratl) .and. (l > lsuf) ) then
     ismall = ismall + 1
  else
     ismall = 0
  endif

end subroutine calamp


!-----------------------------------------------------------------------------

subroutine calu(c0,lsq,bvec,u)
  implicit none
  real(kind(0d0)) :: lsq
  complex(kind(0d0)) :: c0, bvec(3), u(1:3)

  u(1) = dcmplx( 0.d0 )
  u(2) = u(2) + c0 * bvec(2) / dcmplx(lsq)
  u(3) = u(3) + c0 * bvec(3) / dcmplx(lsq)

end subroutine calu


!-----------------------------------------------------------------------------

subroutine interpolate( ncomp,nderiv,rsta,rrsta,g,u )

  implicit none
  integer :: ncomp,nderiv
  real(kind(0d0)) :: rsta,rrsta(3)
  complex(kind(0d0)) :: g(3*ncomp),u(ncomp)
  real(kind(0d0)):: dh(3)
  integer :: ip(3),ier,i,itmp,icomp
  complex(kind(0d0)) :: a(3,3),b(3),wk(3)
  real(kind(0d0)) :: eps
  data eps / -1.d0 /

  do icomp=1,ncomp
     u(icomp) = dcmplx(0.d0)
  enddo

  do i=1,3
     dh(i) = rrsta(i) - rsta
  enddo


  if ( (dh(2) == 0.d0) .and. (nderiv == 0)) then
     itmp = ncomp + 1
     do icomp=1,ncomp
        u(icomp) = g(itmp)
        itmp = itmp + 1
     enddo
     return
  endif

  do i=1,3
     a(1,i) = dcmplx( 1.d0 )
     a(2,i) = dcmplx( dh(i) )
     a(3,i) = dcmplx( dh(i) * dh(i) / 2.d0 )
  enddo

  call fillinpb(nderiv,b)

  call glu(a,3,3,b,eps,wk,ip,ier)


  do icomp=1,ncomp
     do i=1,3
        u(icomp) = u(icomp) + b(i) * g( ncomp * (i-1) + icomp )
     enddo
  enddo

  return
end subroutine interpolate


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




subroutine fillinpb( nderiv,b )

  integer :: nderiv
  complex(kind(0d0)) :: b(3)

  if ( (nderiv /= 0) .and. (nderiv /= 1) .and. (nderiv /= 2) ) &
       &     pause 'invalid argument (fillinpb)'
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



!------------------------------------------------------------------------------


subroutine calcutd(nzone,nnl,tmpr,rat,nn,ra,kc)

  implicit none
  integer :: nzone,nn,kc,nnl(*)
  complex(kind(0d0)) :: tmpr(*)
  real(kind(0d0)) :: rat,ra(*)
  integer  :: nc
  real(kind(0d0)):: cU(nn),rc
  real(kind(0d0)) :: maxamp,amp(nn)
  integer :: iz,jz,jj,i,ml(nzone),tzone

  do jj=1,nn
     cU(jj) = 0.d0
  enddo

  iz = 2
  jz = 1
  do jj=1,nn
     cU(jj) = tmpr(jj)
  enddo

  maxamp = -1.d0
  do i=1,nn
     amp(i) = cU(i)
     if (maxamp < amp(i)) maxamp = amp(i)
  enddo
  maxamp = maxamp * rat ! threshold value
  if (maxamp == 0.d0) then
     kc = 1
     return
  endif

  do i=1,nn
     if (amp(i) > maxamp) then
        nc = i
        exit
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
  kc = nc


  return
end subroutine calcutd




subroutine calstg4onedepth( maxnlay,maxnzone,nzone,vrmin,vrmax,iphase,rrho,vpv,vph,vsv,vsh,eta,rmax,r,updown,ecA,ecC,ecF,ecL,ecN)

  ! Computing the structure grid points.
  implicit none
  integer:: maxnlay,maxnzone,nzone,iphase(*),updown
  real(kind(0d0)):: rrho(4,*),vpv(4,*),vph(4,*),vsv(4,*),vsh(4,*),eta(4,*),vrmin(*),vrmax(*)
  real(kind(0d0)):: rmax,r
!  real(kind(0d0)):: rho,kappa,ecKx,ecKy,ecKz
  real(kind(0d0)):: ecL,ecN
  real(kind(0d0)):: ecA,ecC,ecF
  real(kind(0d0)):: trho,tvpv,tvph,tvsv,tvsh,teta,coef
  integer:: izone,j,itmp,jtmp,spn

  ! computing the structure grid points
  itmp = 0
  jtmp = 0
  spn = 0

  do izone=1,nzone
     if ((vrmin(izone) <= r) .and. (vrmax(izone) > r)) then
        spn = izone
     endif
  enddo

  if (vrmax(nzone) == r) spn = nzone
  if ((vrmin(spn) == r) .and. (updown == -1)) then
     spn = spn - 1
  endif



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


subroutine calstg4onedepth_psv(maxnlay,maxnzone,nzone,vrmin,vrmax,iphase,rrho,vpv,vph,vsv,vsh,eta,rmax,r, &
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
end subroutine calstg4onedepth_psv




subroutine translat(geodetic,geocentric)

  implicit none
  real(kind(0d0)),parameter ::  flattening = 1.d0 / 298.25d0
  real(kind(0d0)), parameter :: pi = 3.1415926535897932d0
  real(kind(0d0)) :: geocentric, geodetic
!  real(kind(0d0)) :: tmp
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


!!!!!!!NEW function FOR THE MPI

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

subroutine WriteFrqByproc(proc,imin,imax)
 integer k,proc,imin,imax,dis,redun,summ,f1,f2
open(20,file='FrqsMpi.txt',status = 'replace',form = 'formatted')
dis = int(abs(imax-imin)/proc)
dis =dis+1
write(*,*) 'dis, imin and imax is :',dis,imin,imax
write(20,"(i2)") 2


if ( abs(imax-imin)-proc*dis > 0 ) then
  write(*,*) 'wrong task distribution.please check the code WriteFrqByproc!'
  stop
else
  redun = abs(imax-imin)-proc*(dis-1)
   summ = 0
   f1 = 0
   f2 = dis
  do k=0,redun-1,1
   summ = f2
  write(20,"(i4,i6,i6)") k,f1,f2
  write(*,*) 'summ is',summ
   f1 = f2+1
   f2 = f2+dis
  enddo

  if (redun == 0) then
   f2 = f1+dis-1
  else
   f2 = f1+dis-2
  endif

  do k=redun,proc-1,1
     summ = f2
     write(*,*) 'summ is',summ
     write(20,"(i4,i6,i6)") k,f1,f2
     f1 = f2+1
     f2 = f2+dis-1
  enddo

  if ( summ /= abs(imax-imin) ) then
  write(*,*) 'wrong task distribution.please check the code WriteFrqByproc!'
  stop
  endif

endif
close(20)
return

end subroutine WriteFrqByproc

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

          if ( myrank < mod(n_global,nbproc) .and. mod(n_global,nbproc) /= 0) then
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
 complex(kind(0d0)), dimension(maxlmax_g,1,6,-2:2,r_n), intent(in) :: tabg0,tabg0der
!! DK DK here I purposely declare the buffer with a size llog0 <= maxlmax_g,
!! DK DK this way I can then use a simple write statement to write the whole array to disk
  complex(kind(0d0)), dimension(llog0,1,6,-2:2) :: tabg0_small_buffer
  integer :: illog

    ! copy the array to the buffer
    !IBM* ASSERT (MINITERCNT(1000))
    !DIR$ loop count min(1000)
      do illog = 1,llog0

      !   tabg0_small_buffer(illog,:,:,:) = tabg0(illog,:,:,:,ir_)
      !! DK DK unrolled the loops to allow for vectorization of the outer loop on "illog"
        tabg0_small_buffer(illog,1,1,-2) = tabg0(illog,1,1,-2,ir_)
  tabg0_small_buffer(illog,1,2,-2) = tabg0(illog,1,2,-2,ir_)
  tabg0_small_buffer(illog,1,3,-2) = tabg0(illog,1,3,-2,ir_)
  tabg0_small_buffer(illog,1,4,-2) = tabg0(illog,1,4,-2,ir_)
  tabg0_small_buffer(illog,1,5,-2) = tabg0(illog,1,5,-2,ir_)
  tabg0_small_buffer(illog,1,6,-2) = tabg0(illog,1,6,-2,ir_)
  tabg0_small_buffer(illog,1,1,-1) = tabg0(illog,1,1,-1,ir_)
  tabg0_small_buffer(illog,1,2,-1) = tabg0(illog,1,2,-1,ir_)
  tabg0_small_buffer(illog,1,3,-1) = tabg0(illog,1,3,-1,ir_)
  tabg0_small_buffer(illog,1,4,-1) = tabg0(illog,1,4,-1,ir_)
  tabg0_small_buffer(illog,1,5,-1) = tabg0(illog,1,5,-1,ir_)
  tabg0_small_buffer(illog,1,6,-1) = tabg0(illog,1,6,-1,ir_)
  tabg0_small_buffer(illog,1,1,0) = tabg0(illog,1,1,0,ir_)
  tabg0_small_buffer(illog,1,2,0) = tabg0(illog,1,2,0,ir_)
  tabg0_small_buffer(illog,1,3,0) = tabg0(illog,1,3,0,ir_)
  tabg0_small_buffer(illog,1,4,0) = tabg0(illog,1,4,0,ir_)
  tabg0_small_buffer(illog,1,5,0) = tabg0(illog,1,5,0,ir_)
  tabg0_small_buffer(illog,1,6,0) = tabg0(illog,1,6,0,ir_)
  tabg0_small_buffer(illog,1,1,1) = tabg0(illog,1,1,1,ir_)
  tabg0_small_buffer(illog,1,2,1) = tabg0(illog,1,2,1,ir_)
  tabg0_small_buffer(illog,1,3,1) = tabg0(illog,1,3,1,ir_)
  tabg0_small_buffer(illog,1,4,1) = tabg0(illog,1,4,1,ir_)
  tabg0_small_buffer(illog,1,5,1) = tabg0(illog,1,5,1,ir_)
  tabg0_small_buffer(illog,1,6,1) = tabg0(illog,1,6,1,ir_)
  tabg0_small_buffer(illog,1,1,2) = tabg0(illog,1,1,2,ir_)
  tabg0_small_buffer(illog,1,2,2) = tabg0(illog,1,2,2,ir_)
  tabg0_small_buffer(illog,1,3,2) = tabg0(illog,1,3,2,ir_)
  tabg0_small_buffer(illog,1,4,2) = tabg0(illog,1,4,2,ir_)
  tabg0_small_buffer(illog,1,5,2) = tabg0(illog,1,5,2,ir_)
  tabg0_small_buffer(illog,1,6,2) = tabg0(illog,1,6,2,ir_)

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
        tabg0_small_buffer(illog,1,2,-2) = tabg0der(illog,1,2,-2,ir_)
  tabg0_small_buffer(illog,1,3,-2) = tabg0der(illog,1,3,-2,ir_)
  tabg0_small_buffer(illog,1,4,-2) = tabg0der(illog,1,4,-2,ir_)
  tabg0_small_buffer(illog,1,5,-2) = tabg0der(illog,1,5,-2,ir_)
  tabg0_small_buffer(illog,1,6,-2) = tabg0der(illog,1,6,-2,ir_)
  tabg0_small_buffer(illog,1,1,-1) = tabg0der(illog,1,1,-1,ir_)
  tabg0_small_buffer(illog,1,2,-1) = tabg0der(illog,1,2,-1,ir_)
  tabg0_small_buffer(illog,1,3,-1) = tabg0der(illog,1,3,-1,ir_)
  tabg0_small_buffer(illog,1,4,-1) = tabg0der(illog,1,4,-1,ir_)
  tabg0_small_buffer(illog,1,5,-1) = tabg0der(illog,1,5,-1,ir_)
  tabg0_small_buffer(illog,1,6,-1) = tabg0der(illog,1,6,-1,ir_)
  tabg0_small_buffer(illog,1,1,0) = tabg0der(illog,1,1,0,ir_)
  tabg0_small_buffer(illog,1,2,0) = tabg0der(illog,1,2,0,ir_)
  tabg0_small_buffer(illog,1,3,0) = tabg0der(illog,1,3,0,ir_)
        tabg0_small_buffer(illog,1,4,0) = tabg0der(illog,1,4,0,ir_)
  tabg0_small_buffer(illog,1,5,0) = tabg0der(illog,1,5,0,ir_)
  tabg0_small_buffer(illog,1,6,0) = tabg0der(illog,1,6,0,ir_)
  tabg0_small_buffer(illog,1,1,1) = tabg0der(illog,1,1,1,ir_)
  tabg0_small_buffer(illog,1,2,1) = tabg0der(illog,1,2,1,ir_)
  tabg0_small_buffer(illog,1,3,1) = tabg0der(illog,1,3,1,ir_)
  tabg0_small_buffer(illog,1,4,1) = tabg0der(illog,1,4,1,ir_)
  tabg0_small_buffer(illog,1,5,1) = tabg0der(illog,1,5,1,ir_)
  tabg0_small_buffer(illog,1,6,1) = tabg0der(illog,1,6,1,ir_)
  tabg0_small_buffer(illog,1,1,2) = tabg0der(illog,1,1,2,ir_)
  tabg0_small_buffer(illog,1,2,2) = tabg0der(illog,1,2,2,ir_)
  tabg0_small_buffer(illog,1,3,2) = tabg0der(illog,1,3,2,ir_)
  tabg0_small_buffer(illog,1,4,2) = tabg0der(illog,1,4,2,ir_)
  tabg0_small_buffer(illog,1,5,2) = tabg0der(illog,1,5,2,ir_)
  tabg0_small_buffer(illog,1,6,2) = tabg0der(illog,1,6,2,ir_)
  enddo
    write(34) tabg0_small_buffer

     end subroutine write_a_block_of_coefs_to_disk
