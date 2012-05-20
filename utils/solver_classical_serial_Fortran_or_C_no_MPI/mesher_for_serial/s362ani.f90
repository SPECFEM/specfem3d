
  subroutine evradker(depth,string,nker,vercof,dvercof,ierror)

  implicit none

  integer :: nker,ierror

  real(kind=4) :: chebyshev(100)
  real(kind=4) :: chebyshev2(100)
  real(kind=4) :: vercof(nker)
  real(kind=4) :: dvercof(nker)
  real(kind=4) :: splpts(100)

  character(len=80) string

  logical upper,upper_650
  logical lower,lower_650

  real(kind=4), parameter :: r0=6371.
  real(kind=4), parameter :: rmoho=6371.0-24.4
  real(kind=4), parameter :: r670=6371.-670.
  real(kind=4), parameter :: r650=6371.-650.
  real(kind=4), parameter :: rcmb=3480.0

  integer :: i,nspl,nskip,nlower,nupper,iker,lstr

  real(kind=4) :: u,u2,ddep,radius2,radius,depth

  ierror=0
  lstr=len_trim(string)

  radius=r0-depth
  ddep=0.1
  radius2=r0-depth+ddep
  upper=.false.
  lower=.false.
  if(radius > rcmb.and.radius < r670) then
  lower=.true.
  else if(radius >= r670.and.radius < rmoho) then
  upper=.true.
  endif
  upper_650=.false.
  lower_650=.false.
  if(radius > rcmb.and.radius < r650) then
  lower_650=.true.
  else if(radius >= r650.and.radius < rmoho) then
  upper_650=.true.
  endif
  do iker=1,nker
  vercof(iker)=0.
  dvercof(iker)=0.
  enddo

  if(string(1:16) == 'WDC+SPC_U4L8CHEB') then
  nupper=5
  nlower=9
  nskip=2
  if(upper) then
    u=(radius+radius-rmoho-r670)/(rmoho-r670)
    u2=(radius2+radius2-rmoho-r670)/(rmoho-r670)
!          write(6,"('upper mantle:',2f10.3)") u,u2
    call chebyfun(u,13,chebyshev)
    do i=1+nskip,nskip+nupper
      vercof(i)=chebyshev(i-nskip)
    enddo
    call chebyfun(u2,13,chebyshev2)
    do i=1+nskip,nskip+nupper
      dvercof(i)=(chebyshev2(i-nskip)-chebyshev(i-nskip))/ddep
    enddo
  else if(lower) then
    u=(radius+radius-r670-rcmb)/(r670-rcmb)
    u2=(radius2+radius2-r670-rcmb)/(r670-rcmb)
!          write(6,"('lower mantle:',2f10.3)") u,u2
    call chebyfun(u,13,chebyshev)
    do i=1+nskip+nupper,nskip+nupper+nlower
      vercof(i)=chebyshev(i-nskip-nupper)
    enddo
    call chebyfun(u2,13,chebyshev2)
    do i=1+nskip+nupper,nskip+nupper+nlower
      dvercof(i)=(chebyshev2(i-nskip-nupper)- &
                    chebyshev(i-nskip-nupper))/ddep
    enddo
  endif
  else if(string(1:13) == 'WDC+SHSVWM20A') then
  nspl=20
  splpts(1)=0.
  splpts(2)=50.
  splpts(3)=100.
  splpts(4)=150.
  splpts(5)=200.
  splpts(6)=250.
  splpts(7)=300.
  splpts(8)=400.
  splpts(9)=500.
  splpts(10)=600.
  splpts(11)=700.
  splpts(12)=850.
  splpts(13)=1050.
  splpts(14)=1300.
  splpts(15)=1600.
  splpts(16)=1900.
  splpts(17)=2200.
  splpts(18)=2500.
  splpts(19)=2700.
  splpts(20)=2891.
  call vbspl(depth,nspl,splpts,vercof(2),dvercof(2))
  do i=22,27
    vercof(i)=vercof(i-20)
    dvercof(i)=dvercof(i-20)
  enddo
  vercof(1)=1.
  else if(string(1:16) == 'WDC+XBS_362_U6L8') then
  if(upper) then
   nspl=6
   splpts(1)=24.4
   splpts(2)=100.
   splpts(3)=225.
   splpts(4)=350.
   splpts(5)=500.
   splpts(6)=670.
   call vbspl(depth,nspl,splpts,vercof(2),dvercof(2))
  else if(lower) then
 nspl=8
   splpts(1)=670.
   splpts(2)=820.
   splpts(3)=1320.
   splpts(4)=1820.
   splpts(5)=2320.
   splpts(6)=2550.
   splpts(7)=2791.
   splpts(8)=2891.
   call vbspl(depth,nspl,splpts,vercof(8),dvercof(8))
  endif
  vercof(1)=1.
!        vercof(16)=1.
!        vercof(17)=1.
!      else if(string(1:21) == 'WDC+ANI_362_U6L8_TOPO') then
!        if(upper) then
!         nspl=6
!         splpts(1)=24.4
!         splpts(2)=100.
!         splpts(3)=225.
!         splpts(4)=350.
!         splpts(5)=500.
!         splpts(6)=670.
!         call vbspl(depth,nspl,splpts,vercof(2),dvercof(2))
!         do i=16,21
!          vercof(i)=vercof(i-14)
!          dvercof(i)=dvercof(i-14)
!         enddo
!     else if(lower) then
!      nspl=8
!         splpts(1)=670.
!         splpts(2)=820.
!         splpts(3)=1320.
!         splpts(4)=1820.
!         splpts(5)=2320.
!         splpts(6)=2550.
!         splpts(7)=2791.
!         splpts(8)=2891.
!         call vbspl(depth,nspl,splpts,vercof(8),dvercof(8))
!     endif
!        vercof(1)=1.
!        vercof(22)=1.
!        vercof(23)=1.
!        vercof(24)=1.
!        vercof(25)=1.
  else if( &
       (string(1:lstr) == 'WDC+ANI_362_U6L8'.and.lstr == 16) &
       .or. &
           (string(1:lstr) == 'WDC+ANI_362_U6L8_TOPO'.and.lstr == 21) &
       ) then
  if(upper) then
   nspl=6
   splpts(1)=24.4
   splpts(2)=100.
   splpts(3)=225.
   splpts(4)=350.
   splpts(5)=500.
   splpts(6)=670.
   call vbspl(depth,nspl,splpts,vercof(2),dvercof(2))
   do i=16,21
    vercof(i)=vercof(i-14)
    dvercof(i)=dvercof(i-14)
   enddo
  else if(lower) then
 nspl=8
   splpts(1)=670.
   splpts(2)=820.
   splpts(3)=1320.
   splpts(4)=1820.
   splpts(5)=2320.
   splpts(6)=2550.
   splpts(7)=2791.
   splpts(8)=2891.
   call vbspl(depth,nspl,splpts,vercof(8),dvercof(8))
  endif
  vercof(1)=1.
  vercof(22)=1.
  vercof(23)=1.
  else if(string(1:lstr) == 'WDC+WM_362_U6L8'.and.lstr == 15) then
  if(upper) then
   nspl=6
   splpts(1)=24.4
   splpts(2)=100.
   splpts(3)=225.
   splpts(4)=350.
   splpts(5)=500.
   splpts(6)=670.
   call vbspl(depth,nspl,splpts,vercof(2),dvercof(2))
   do i=16,21
    vercof(i)=vercof(i-14)
    dvercof(i)=dvercof(i-14)
   enddo
  else if(lower) then
 nspl=8
   splpts(1)=670.
   splpts(2)=820.
   splpts(3)=1320.
   splpts(4)=1820.
   splpts(5)=2320.
   splpts(6)=2550.
   splpts(7)=2791.
   splpts(8)=2891.
   call vbspl(depth,nspl,splpts,vercof(8),dvercof(8))
   do i=22,29
    vercof(i)=vercof(i-14)
    dvercof(i)=dvercof(i-14)
   enddo
  endif
  vercof(1)=1.
  vercof(30)=1.
  vercof(31)=1.
  vercof(32)=1.
  else if( &
     (string(1:lstr) == 'WDC+ANI_362_U6L8_650'.and.lstr == 20) &
     .or. &
         (string(1:lstr) == 'WDC+ANI_362_U6L8_TOPO_650'.and.lstr == 25) &
     ) then
  if(upper_650) then
   nspl=6
   splpts(1)=24.4
   splpts(2)=100.
   splpts(3)=225.
   splpts(4)=350.
   splpts(5)=500.
   splpts(6)=650.
   call vbspl(depth,nspl,splpts,vercof(2),dvercof(2))
   do i=16,21
    vercof(i)=vercof(i-14)
    dvercof(i)=dvercof(i-14)
   enddo
  else if(lower_650) then
 nspl=8
   splpts(1)=650.
   splpts(2)=820.
   splpts(3)=1320.
   splpts(4)=1820.
   splpts(5)=2320.
   splpts(6)=2550.
   splpts(7)=2791.
   splpts(8)=2891.
   call vbspl(depth,nspl,splpts,vercof(8),dvercof(8))
  endif
  vercof(1)=1.
  vercof(22)=1.
  vercof(23)=1.
  else if(string(1:lstr) == 'WDC+WM_362_U6L8_650' &
       .and.lstr == 19) then
  if(upper_650) then
   nspl=6
   splpts(1)=24.4
   splpts(2)=100.
   splpts(3)=225.
   splpts(4)=350.
   splpts(5)=500.
   splpts(6)=650.
   call vbspl(depth,nspl,splpts,vercof(2),dvercof(2))
   do i=16,21
    vercof(i)=vercof(i-14)
    dvercof(i)=dvercof(i-14)
   enddo
  else if(lower_650) then
 nspl=8
   splpts(1)=650.
   splpts(2)=820.
   splpts(3)=1320.
   splpts(4)=1820.
   splpts(5)=2320.
   splpts(6)=2550.
   splpts(7)=2791.
   splpts(8)=2891.
   call vbspl(depth,nspl,splpts,vercof(8),dvercof(8))
   do i=22,29
    vercof(i)=vercof(i-14)
    dvercof(i)=dvercof(i-14)
   enddo
  endif
  vercof(1)=1.
  vercof(30)=1.
  vercof(31)=1.
  vercof(32)=1.
  else if(string(1:lstr) == 'WDC+U8L8_650'.and.lstr == 12) then
  if(upper_650) then
   nspl=8
   splpts(1)=24.4
   splpts(2)=75.
   splpts(3)=150.
   splpts(4)=225.
   splpts(5)=300.
   splpts(6)=410.
   splpts(7)=530.
   splpts(8)=650.
   call vbspl(depth,nspl,splpts,vercof(2),dvercof(2))
   do i=18,25
    vercof(i)=vercof(i-16)
    dvercof(i)=dvercof(i-16)
   enddo
  else if(lower_650) then
 nspl=8
   splpts(1)=650.
   splpts(2)=820.
   splpts(3)=1320.
   splpts(4)=1820.
   splpts(5)=2320.
   splpts(6)=2550.
   splpts(7)=2791.
   splpts(8)=2891.
   call vbspl(depth,nspl,splpts,vercof(10),dvercof(10))
   do i=26,33
    vercof(i)=vercof(i-16)
    dvercof(i)=dvercof(i-16)
   enddo
  endif
  vercof(1)=1.
  vercof(34)=1.
  vercof(35)=1.
  vercof(36)=1.
  else if(string(1:lstr) == 'WDC+U8L8_670'.and.lstr == 12) then
  if(upper) then
   nspl=8
   splpts(1)=24.4
   splpts(2)=75.
   splpts(3)=150.
   splpts(4)=225.
   splpts(5)=300.
   splpts(6)=410.
   splpts(7)=530.
   splpts(8)=670.
   call vbspl(depth,nspl,splpts,vercof(2),dvercof(2))
   do i=18,25
    vercof(i)=vercof(i-16)
    dvercof(i)=dvercof(i-16)
   enddo
  else if(lower) then
 nspl=8
   splpts(1)=670.
   splpts(2)=820.
   splpts(3)=1320.
   splpts(4)=1820.
   splpts(5)=2320.
   splpts(6)=2550.
   splpts(7)=2791.
   splpts(8)=2891.
   call vbspl(depth,nspl,splpts,vercof(10),dvercof(10))
   do i=26,33
    vercof(i)=vercof(i-16)
    dvercof(i)=dvercof(i-16)
   enddo
  endif
  vercof(1)=1.
  vercof(34)=1.
  vercof(35)=1.
  vercof(36)=1.
  else if( &
      (string(1:lstr) == 'WDC+U8L8_I1D_650'.and.lstr == 16) &
      .or. &
      (string(1:lstr) == 'WDC+U8L8_I3D_650'.and.lstr == 16) &
      ) then
  if(upper_650) then
   nspl=8
   splpts(1)=24.4
   splpts(2)=75.
   splpts(3)=150.
   splpts(4)=225.
   splpts(5)=300.
   splpts(6)=410.
   splpts(7)=530.
   splpts(8)=650.
   call vbspl(depth,nspl,splpts,vercof(2),dvercof(2))
   do i=18,25
    vercof(i)=vercof(i-16)
    dvercof(i)=dvercof(i-16)
   enddo
   do i=37,40
    vercof(i)=vercof(i-35)
    dvercof(i)=dvercof(i-35)
   enddo
   do i=41,44
    vercof(i)=vercof(i-39)
    dvercof(i)=dvercof(i-39)
   enddo
   do i=45,48
    vercof(i)=vercof(i-43)
    dvercof(i)=dvercof(i-43)
   enddo
   do i=49,52
    vercof(i)=vercof(i-47)
    dvercof(i)=dvercof(i-47)
   enddo
  else if(lower_650) then
 nspl=8
   splpts(1)=650.
   splpts(2)=820.
   splpts(3)=1320.
   splpts(4)=1820.
   splpts(5)=2320.
   splpts(6)=2550.
   splpts(7)=2791.
   splpts(8)=2891.
   call vbspl(depth,nspl,splpts,vercof(10),dvercof(10))
   do i=26,33
    vercof(i)=vercof(i-16)
    dvercof(i)=dvercof(i-16)
   enddo
  endif
  vercof(1)=1.
  vercof(34)=1.
  vercof(35)=1.
  vercof(36)=1.
  else if((string(1:lstr) == 'WDC+I1D_650'.and.lstr == 11).or. &
          (string(1:lstr) == 'WDC+I3D_650'.and.lstr == 11)) then
  if(upper_650) then
   nspl=8
   splpts(1)=24.4
   splpts(2)=75.
   splpts(3)=150.
   splpts(4)=225.
   splpts(5)=300.
   splpts(6)=410.
   splpts(7)=530.
   splpts(8)=650.
   call vbspl(depth,nspl,splpts,vercof(2),dvercof(2))
   do i=18,25
    vercof(i)=vercof(i-16)
    dvercof(i)=dvercof(i-16)
   enddo
   do i=37,44
    vercof(i)=vercof(i-35)
    dvercof(i)=dvercof(i-35)
   enddo
   do i=53,60
    vercof(i)=vercof(i-51)
    dvercof(i)=dvercof(i-51)
   enddo
   do i=69,76
    vercof(i)=vercof(i-67)
    dvercof(i)=dvercof(i-67)
   enddo
   do i=85,92
    vercof(i)=vercof(i-83)
    dvercof(i)=dvercof(i-83)
   enddo
  else if(lower_650) then
 nspl=8
   splpts(1)=650.
   splpts(2)=820.
   splpts(3)=1320.
   splpts(4)=1820.
   splpts(5)=2320.
   splpts(6)=2550.
   splpts(7)=2791.
   splpts(8)=2891.
   call vbspl(depth,nspl,splpts,vercof(10),dvercof(10))
   do i=26,33
    vercof(i)=vercof(i-16)
    dvercof(i)=dvercof(i-16)
   enddo
   do i=45,52
    vercof(i)=vercof(i-35)
    dvercof(i)=dvercof(i-35)
   enddo
   do i=61,68
    vercof(i)=vercof(i-51)
    dvercof(i)=dvercof(i-51)
   enddo
   do i=77,84
    vercof(i)=vercof(i-67)
    dvercof(i)=dvercof(i-67)
   enddo
   do i=93,100
    vercof(i)=vercof(i-83)
    dvercof(i)=dvercof(i-83)
   enddo
  endif
  vercof(1)=1.
  vercof(34)=1.
  vercof(35)=1.
  vercof(36)=1.
  else if(string(1:lstr) == 'V16A4_V7A4'.and.lstr == 10) then
  if(upper_650) then
   nspl=8
   splpts(1)=24.4
   splpts(2)=75.
   splpts(3)=150.
   splpts(4)=225.
   splpts(5)=300.
   splpts(6)=410.
   splpts(7)=530.
   splpts(8)=650.
   call vbspl(depth,nspl,splpts,vercof(1),dvercof(1))
   do i=17,20
    vercof(i)=vercof(i-16)
    dvercof(i)=dvercof(i-16)
   enddo
   do i=23,29
    vercof(i)=vercof(i-22)
    dvercof(i)=dvercof(i-22)
   enddo
   do i=30,33
    vercof(i)=vercof(i-29)
    dvercof(i)=dvercof(i-29)
   enddo
  else if(lower_650) then
 nspl=8
   splpts(1)=650.
   splpts(2)=820.
   splpts(3)=1320.
   splpts(4)=1820.
   splpts(5)=2320.
   splpts(6)=2550.
   splpts(7)=2791.
   splpts(8)=2891.
   call vbspl(depth,nspl,splpts,vercof(9),dvercof(9))
  endif
  vercof(21)=1.
  vercof(22)=1.
  else
  write(6,"('problem 4')")
  write(6,"(a)")string(1:len_trim(string))
  stop
  endif

  end subroutine evradker

! ---

  subroutine chebyfun(u,kmax,f)

  implicit none

  integer :: kmax

  real(kind=4) :: chebycoeff(0:13),f(0:kmax),u

  integer :: k

  real(kind=4) :: twou

  data chebycoeff / &
   0.70710678118655,1.2247448713916,1.0350983390135,1.0145993123918, &
   1.00803225754840,1.0050890913907,1.0035149493262,1.0025740068320, &
   1.00196657023780,1.0015515913133,1.0012554932754,1.0010368069141, &
   1.00087070107920,1.0007415648034 /

  if(kmax > 13)then
   write(*,"(' kmax exceeds the limit in chebyfun')")
   stop
  endif

  f(0)=1.0
  f(1)=u
  twou=2.0*u

  do k=2,kmax
   f(k) = twou*f(k-1)-f(k-2)
  enddo

  do k=0,kmax
   f(k)=f(k)*chebycoeff(k)
  enddo

  end subroutine chebyfun


  subroutine gt3dmodl(lu,targetfile, &
      maxhpa,maxker,maxcoe, &
      numhpa,numker,numcoe,lmxhpa, &
      ihpakern,itypehpa,coe, &
      itpspl,xlatspl,xlonspl,radispl, &
      numvar,ivarkern,varstr, &
      refmdl,kerstr,hsplfl,dskker,ierror)

  implicit none

  integer, parameter :: mxhpar=2
  integer, parameter :: mxkern=200
  integer, parameter :: mxcoef=2000

  character(len=80) refmodel
  character(len=80) kernstri
  character(len=40) desckern(mxkern)
  character(len=80) hsplfile(mxhpar)

  integer ihorpar(mxkern)
  integer ityphpar(mxhpar)
  integer ixlspl(mxcoef,mxhpar)
  integer lmaxhor(mxhpar)
  integer ncoefhor(mxhpar)

  real(kind=4) coef(mxcoef,mxkern)
  real(kind=4) xlaspl(mxcoef,mxhpar)
  real(kind=4) xlospl(mxcoef,mxhpar)
  real(kind=4) xraspl(mxcoef,mxhpar)

  character(len=128) targetfile

  integer numhpa,numker,maxhpa,maxker,maxcoe

  integer numcoe(maxhpa)
  integer lmxhpa(maxhpa)
  integer ihpakern(maxker)
  integer itypehpa(maxhpa)
  integer itpspl(maxcoe,maxhpa)
  integer ivarkern(maxker)

  real(kind=4) coe(maxcoe,maxker)
  real(kind=4) xlatspl(maxcoe,maxhpa)
  real(kind=4) xlonspl(maxcoe,maxhpa)
  real(kind=4) radispl(maxcoe,maxhpa)

  character(len=80) refmdl
  character(len=80) kerstr
  character(len=80) hsplfl(maxhpa)
  character(len=40) dskker(maxker)
  character(len=40) string
  character(len=40) varstr(maxker)

  integer numvar,ierror,lu,nhorpar,nmodkern,i,j,lstr,k

  ierror=0
  call rd3dmodl(lu,targetfile,ierror, &
    nmodkern,nhorpar,ityphpar, &
    ihorpar,lmaxhor,ncoefhor, &
    xlaspl,xlospl,xraspl,ixlspl,coef, &
    hsplfile,refmodel,kernstri,desckern)

  if(nhorpar <= maxhpa) then
  numhpa=nhorpar
  else
  ierror=ierror+1
  endif

  if(nmodkern <= maxker) then
  numker=nmodkern
  else
  ierror=ierror+1
  endif

  do i=1,nmodkern
  ihpakern(i)=ihorpar(i)
  dskker(i)=desckern(i)
  do j=1,ncoefhor(ihpakern(i))
    coe(j,i)=coef(j,i)
!          if(j == 1) then
!            write(6,"(e12.4)") coe(j,i)
!          endif
  enddo
  enddo

  do i=1,nhorpar
  numcoe(i)=ncoefhor(i)
  lmxhpa(i)=lmaxhor(i)
  itypehpa(i)=ityphpar(i)
  if(itypehpa(i) == 2) then
    do j=1,ncoefhor(i)
      itpspl(j,i)=ixlspl(j,i)
      xlatspl(j,i)=xlaspl(j,i)
      xlonspl(j,i)=xlospl(j,i)
      radispl(j,i)=xraspl(j,i)
    enddo
  endif
  hsplfl(i)=hsplfile(i)
  enddo

  numvar=0
  do i=1,nmodkern
  string=dskker(i)
  lstr=len_trim(string)
  j=1
  do while(string(j:j) /= ','.and.j < lstr)
    j=j+1
  enddo
  ivarkern(i)=0
  do k=1,numvar
    if(string(1:j) == varstr(k)(1:j)) then
      ivarkern(i)=k
    endif
  enddo
  if(ivarkern(i) == 0) then
    numvar=numvar+1
    varstr(numvar)=string(1:j)
    ivarkern(i)=numvar
  endif
  enddo

  refmdl=refmodel
  kerstr=kernstri

  end subroutine gt3dmodl


  subroutine rd3dmodl(lu,filename,ierror, &
    nmodkern,nhorpar,ityphpar, &
    ihorpar,lmaxhor,ncoefhor, &
    xlaspl,xlospl,xraspl,ixlspl,coef, &
    hsplfile,refmodel,kernstri,desckern)

  implicit none

  integer, parameter :: mxhpar=2
  integer, parameter :: mxkern=200
  integer, parameter :: mxcoef=2000

  character(len=80) refmodel
  character(len=80) kernstri
  character(len=40) desckern(mxkern)
  character(len=80) hsplfile(mxhpar)

  integer ihorpar(mxkern)
  integer ityphpar(mxhpar)
  integer ixlspl(mxcoef,mxhpar)
  integer lmaxhor(mxhpar)
  integer ncoefhor(mxhpar)

  real(kind=4) coef(mxcoef,mxkern)
  real(kind=4) xlaspl(mxcoef,mxhpar)
  real(kind=4) xlospl(mxcoef,mxhpar)
  real(kind=4) xraspl(mxcoef,mxhpar)

  character(len=128) filename

  character(len=128) string
  character(len=128) substr

  integer :: lu,ierror

  integer :: ncoef,i,ihor,ifst,ilst,ifst1,ios,lstr,nmodkern,idummy,nhorpar,lmax

  open(lu,file=filename,iostat=ios)
  if(ios /= 0) then
  stop 'error opening 3-d model'
  endif
  do while (ios == 0)
  read(lu,"(a)",iostat=ios) string
  lstr=len_trim(string)
  if(ios == 0) then
    if(string(1:16) == 'REFERENCE MODEL:') then
      substr=string(17:lstr)
      ifst=1
      ilst=len_trim(substr)
      do while (substr(ifst:ifst) == ' '.and.ifst < ilst)
        ifst=ifst+1
      enddo
      if(ilst-ifst <= 0) then
        stop 'error reading model 1'
      else
        refmodel=substr(ifst:ilst)
      endif
    else if(string(1:11) == 'KERNEL SET:') then
      substr=string(12:len_trim(string))
      ifst=1
      ilst=len_trim(substr)
      do while (substr(ifst:ifst) == ' '.and.ifst < ilst)
        ifst=ifst+1
      enddo
      if(ilst-ifst <= 0) then
        stop 'error reading model 2'
      else
        kernstri=substr(ifst:ilst)
      endif
    else if(string(1:25) == 'RADIAL STRUCTURE KERNELS:') then
      substr=string(26:len_trim(string))
      read(substr,*,iostat=ierror) nmodkern
      if(ierror /= 0) then
        stop 'error reading model 3'
      endif
    else if(string(1:4) == 'DESC'.and.string(9:9) == ':') then
      read(string(5:8),"(i4)") idummy
      substr=string(10:len_trim(string))
      ifst=1
      ilst=len_trim(substr)
      do while (substr(ifst:ifst) == ' '.and.ifst < ilst)
        ifst=ifst+1
      enddo
      if(ilst-ifst <= 0) then
        stop 'error reading model 4'
      else
        desckern(idummy)=substr(ifst:ilst)
      endif
    else if(string(1:29) == 'HORIZONTAL PARAMETERIZATIONS:') then
      substr=string(30:len_trim(string))
      read(substr,*,iostat=ierror) nhorpar
      if(ierror /= 0) then
        stop 'error reading model 5'
      endif
    else if(string(1:4) == 'HPAR'.and.string(9:9) == ':') then
      read(string(5:8),"(i4)") idummy
      ifst=10
      ilst=len_trim(string)
      do while (string(ifst:ifst) == ' '.and.ifst < ilst)
        ifst=ifst+1
      enddo
      if(ilst-ifst <= 0) then
        stop 'error reading model 6'
      else if(string(ifst:ifst+19) == 'SPHERICAL HARMONICS,') then
        substr=string(20+ifst:len_trim(string))
        read(substr,*) lmax
        ityphpar(idummy)=1
        lmaxhor(idummy)=lmax
        ncoefhor(idummy)=(lmax+1)**2
      else if(string(ifst:ifst+17) == 'SPHERICAL SPLINES,') then
        ifst1=ifst+18
        ifst=len_trim(string)
        ilst=len_trim(string)
        do while(string(ifst:ifst) /= ',')
          ifst=ifst-1
        enddo
        read(string(ifst+1:ilst),*) ncoef
        substr=string(ifst1:ifst-1)
        do while (string(ifst1:ifst1) == ' '.and.ifst1 < ifst)
          ifst1=ifst1+1
        enddo
        hsplfile(idummy)=string(ifst1:ifst-1)
        ityphpar(idummy)=2
        lmaxhor(idummy)=0
        ncoefhor(idummy)=ncoef
        do i=1,ncoef
          read(lu,*) ixlspl(i,idummy),xlaspl(i,idummy), &
             xlospl(i,idummy),xraspl(i,idummy)
        enddo
      endif
    else if(string(1:4) == 'STRU'.and.string(9:9) == ':') then
      read(string(5:8),"(i4)") idummy
      substr=string(10:len_trim(string))
      read(substr,*) ihor
      ihorpar(idummy)=ihor
      ncoef=ncoefhor(ihor)
      read(lu,"(6e12.4)") (coef(i,idummy),i=1,ncoef)
    endif
  endif
  enddo
  close(lu)

  end subroutine rd3dmodl


   subroutine read_model_s362ani(THREE_D_MODEL, &
              THREE_D_MODEL_S362ANI,THREE_D_MODEL_S362WMANI, &
              THREE_D_MODEL_S362ANI_PREM,THREE_D_MODEL_S29EA, &
              numker,numhpa,ihpa,lmxhpa,itypehpa,ihpakern,numcoe,ivarkern,itpspl, &
              xlaspl,xlospl,radspl,coe,hsplfl,dskker,kerstr,varstr,refmdl)

  implicit none

  integer THREE_D_MODEL,THREE_D_MODEL_S362ANI
  integer THREE_D_MODEL_S362WMANI
  integer THREE_D_MODEL_S362ANI_PREM,THREE_D_MODEL_S29EA

  integer lu
  character(len=128) modeldef
  logical exists
  integer numvar
  integer ierror

  integer, parameter :: maxker=200
  integer, parameter :: maxl=72
  integer, parameter :: maxcoe=2000
  integer, parameter :: maxver=1000
  integer, parameter :: maxhpa=2

  integer numker
  integer numhpa
  integer ihpa
  integer lmxhpa(maxhpa)
  integer itypehpa(maxhpa)
  integer ihpakern(maxker)
  integer numcoe(maxhpa)
  integer ivarkern(maxker)
  integer itpspl(maxcoe,maxhpa)

  real(kind=4) xlaspl(maxcoe,maxhpa)
  real(kind=4) xlospl(maxcoe,maxhpa)
  real(kind=4) radspl(maxcoe,maxhpa)
  real(kind=4) coe(maxcoe,maxker)
  character(len=80) hsplfl(maxhpa)
  character(len=40) dskker(maxker)

  character(len=80) kerstr
  character(len=80) refmdl
  character(len=40) varstr(maxker)

! -------------------------------------

  lu=1                    ! --- log unit: input 3-D model
  if(THREE_D_MODEL  ==  THREE_D_MODEL_S362ANI) then
    modeldef='DATA/s362ani/S362ANI'
  elseif(THREE_D_MODEL  ==  THREE_D_MODEL_S362WMANI) then
    modeldef='DATA/s362ani/S362WMANI'
  elseif(THREE_D_MODEL  ==  THREE_D_MODEL_S362ANI_PREM) then
    modeldef='DATA/s362ani/S362ANI_PREM'
  elseif(THREE_D_MODEL  ==  THREE_D_MODEL_S29EA) then
    modeldef='DATA/s362ani/S2.9EA'
  else
    stop 'unknown 3D model in read_model_s362ani'
  endif
  inquire(file=modeldef,exist=exists)
  if(exists) then
    call gt3dmodl(lu,modeldef, &
        maxhpa,maxker,maxcoe, &
        numhpa,numker,numcoe,lmxhpa, &
        ihpakern,itypehpa,coe, &
        itpspl,xlaspl,xlospl,radspl, &
        numvar,ivarkern,varstr, &
        refmdl,kerstr,hsplfl,dskker,ierror)
  else
    write(6,"('the model ',a,' does not exits')") modeldef(1:len_trim(modeldef))
  endif

!         --- check arrays

  if(numker > maxker) stop 'numker > maxker'
  do ihpa=1,numhpa
    if(itypehpa(ihpa) == 1) then
      if(lmxhpa(ihpa) > maxl) stop 'lmxhpa(ihpa) > maxl'
    else if(itypehpa(ihpa) == 2) then
      if(numcoe(ihpa) > maxcoe) stop 'numcoe(ihpa) > maxcoe'
    else
      stop 'problem with itypehpa'
    endif
  enddo

  end subroutine read_model_s362ani


  subroutine splcon(xlat,xlon,nver,verlat,verlon,verrad,ncon,icon,con)

  implicit none

  integer icon(1)

  real(kind=4) verlat(1)
  real(kind=4) verlon(1)
  real(kind=4) verrad(1)
  real(kind=4) con(1)

  double precision dd
  double precision rn
  double precision dr
  double precision xrad
  double precision ver8
  double precision xla8

  integer :: ncon,iver,nver

  real(kind=4) :: xlat,xlon

  xrad=3.14159265358979/180.d0

  ncon=0

  do iver=1,nver
  if(xlat > verlat(iver)-2.*verrad(iver)) then
    if(xlat < verlat(iver)+2.*verrad(iver)) then
      ver8=xrad*(verlat(iver))
      xla8=xrad*(xlat)
      dd=sin(ver8)*sin(xla8)
      dd=dd+cos(ver8)*cos(xla8)* cos(xrad*(xlon-verlon(iver)))
      dd=acos(dd)/xrad
      if(dd > (verrad(iver))*2.d0) then
      else
        ncon=ncon+1
        icon(ncon)=iver
        rn=dd/(verrad(iver))
        dr=rn-1.d0
        if(rn <= 1.d0) then
          con(ncon)=(0.75d0*rn-1.5d0)*(rn**2)+1.d0
        else if(rn > 1.d0) then
          con(ncon)=((-0.25d0*dr+0.75d0)*dr-0.75d0)*dr+0.25d0
        else
          con(ncon)=0.
        endif
      endif
    endif
  endif
  enddo

  end subroutine splcon


! --- evaluate perturbations in per cent

  subroutine subshsv(xcolat,xlon,xrad,dvsh,dvsv,dvph,dvpv, &
    numker,numhpa,numcof,ihpa,lmax,nylm, &
    lmxhpa,itypehpa,ihpakern,numcoe,ivarkern, &
    nconpt,iver,iconpt,conpt,xlaspl,xlospl,radspl, &
    coe,vercof,vercofd,ylmcof,wk1,wk2,wk3,kerstr,varstr)

  implicit none

  integer, parameter :: maxker=200
  integer, parameter :: maxl=72
  integer, parameter :: maxcoe=2000
  integer, parameter :: maxver=1000
  integer, parameter :: maxhpa=2

  integer numker
  integer numhpa,numcof
  integer ihpa,lmax,nylm
  integer lmxhpa(maxhpa)
  integer itypehpa(maxhpa)
  integer ihpakern(maxker)
  integer numcoe(maxhpa)
  integer ivarkern(maxker)

  integer nconpt(maxhpa),iver
  integer iconpt(maxver,maxhpa)
  real(kind=4) conpt(maxver,maxhpa)

  real(kind=4) xlaspl(maxcoe,maxhpa)
  real(kind=4) xlospl(maxcoe,maxhpa)
  real(kind=4) radspl(maxcoe,maxhpa)
  real(kind=4) coe(maxcoe,maxker)
  real(kind=4) vercof(maxker)
  real(kind=4) vercofd(maxker)

  real(kind=4) ylmcof((maxl+1)**2,maxhpa)
  real(kind=4) wk1(maxl+1)
  real(kind=4) wk2(maxl+1)
  real(kind=4) wk3(maxl+1)

  character(len=80) kerstr
  character(len=40) varstr(maxker)

  real(kind=4) :: xcolat,xlon,xrad
  real(kind=4) :: dvsh,dvsv,dvph,dvpv

! --- model evaluation

  integer ish ! --- 0 if SV, 1 if SH
  integer ieval     ! --- 1 for velocity, 2 for anisotropy
  real(kind=4) :: valu(2)    ! --- valu(1) if S; valu(1)=velo, valu(2)=aniso
  real(kind=4) :: value      ! --- used in single evaluation of perturbation
  integer isel      ! --- if variable should be included
  real(kind=4) :: depth      ! --- depth
  real(kind=4) :: x,y  ! --- lat lon
  real(kind=4) :: vsh3drel   ! --- relative perturbation
  real(kind=4) :: vsv3drel   ! --- relative perturbation

! ---

  integer iker,i
  character(len=40) vstr
  integer lstr
  integer ierror

! -------------------------------------

  depth=6371.0-xrad
  call evradker (depth,kerstr,numker,vercof,vercofd,ierror)
  if(ierror /= 0) stop 'ierror evradker'

! --- loop over sv and sh (sv=0,sh=1)

  do ish=0,1

!       --- contributing horizontal basis functions at xlat,xlon

  y=90.0-xcolat
  x=xlon
  do ihpa=1,numhpa
      if(itypehpa(ihpa) == 1) then
        lmax=lmxhpa(ihpa)
        call ylm(y,x,lmax,ylmcof(1,ihpa),wk1,wk2,wk3)
      else if(itypehpa(ihpa) == 2) then
        numcof=numcoe(ihpa)
        call splcon(y,x,numcof,xlaspl(1,ihpa), &
              xlospl(1,ihpa),radspl(1,ihpa), &
              nconpt(ihpa),iconpt(1,ihpa),conpt(1,ihpa))
      else
        write(6,"('problem 1')")
      endif
  enddo

!         --- evaluate 3-D perturbations in velocity and anisotropy

  valu(1)=0. ! --- velocity
  valu(2)=0. ! --- anisotropy

  do ieval=1,2
    value=0.
    do iker=1,numker
      isel=0
      lstr=len_trim(varstr(ivarkern(iker)))
      vstr=(varstr(ivarkern(iker)))
      if(ieval == 1) then
        if(vstr(1:lstr) == 'UM (SH+SV)*0.5,'.or. &
                 vstr(1:lstr) == 'LM (SH+SV)*0.5,'.or. &
                 vstr(1:lstr) == 'EA (SH+SV)*0.5,') then
          isel=1
      endif
      else if(ieval == 2) then
        if(vstr(1:lstr) == 'UM SH-SV,'.or. &
                       vstr(1:lstr) == 'LM SH-SV,'.or. &
                       vstr(1:lstr) == 'EA SH-SV,') then
          isel=1
        endif
      endif

      if(isel == 1) then
        if(vercof(iker) /= 0.) then
            if(itypehpa(ihpakern(iker)) == 1) then
          ihpa=ihpakern(iker)
              nylm=(lmxhpa(ihpakern(iker))+1)**2
              do i=1,nylm
                value=value+vercof(iker)*ylmcof(i,ihpa) &
                          *coe(i,iker)
              enddo
            else if(itypehpa(ihpakern(iker)) == 2) then
          ihpa=ihpakern(iker)
              do i=1,nconpt(ihpa)
                iver=iconpt(i,ihpa)
                value=value+vercof(iker)*conpt(i,ihpa) &
                          *coe(iver,iker)
              enddo
            else
              write(6,"('problem 2')")
              stop
            endif ! --- itypehpa
        endif ! --- vercof(iker) /= 0.
      endif ! --- isel == 1
    enddo ! --- end of do iker=1,numker

    valu(ieval)=value
  enddo ! --- ieval

!       --- evaluate perturbations in vsh and vsv

  if(ish == 1) then
    vsh3drel=valu(1)+0.5*valu(2)
  else if(ish == 0) then
    vsv3drel=valu(1)-0.5*valu(2)
  else
    stop 'something wrong'
  endif

  enddo ! --- by ish

! --- evaluate perturbations in per cent

  dvsh=vsh3drel
  dvsv=vsv3drel
  dvph=0.55*dvsh    ! --- scaling used in the inversion
  dvpv=0.55*dvsv    ! --- scaling used in the inversion

  end subroutine subshsv


! --- evaluate depressions of the 410- and 650-km discontinuities in km

  subroutine subtopo(xcolat,xlon,topo410,topo650, &
                     numker,numhpa,numcof,ihpa,lmax,nylm, &
                     lmxhpa,itypehpa,ihpakern,numcoe,ivarkern, &
                     nconpt,iver,iconpt,conpt,xlaspl,xlospl,radspl, &
                     coe,ylmcof,wk1,wk2,wk3,varstr)

  implicit none

  integer, parameter :: maxker=200
  integer, parameter :: maxl=72
  integer, parameter :: maxcoe=2000
  integer, parameter :: maxver=1000
  integer, parameter :: maxhpa=2

  integer numker
  integer numhpa,numcof
  integer ihpa,lmax,nylm
  integer lmxhpa(maxhpa)
  integer itypehpa(maxhpa)
  integer ihpakern(maxker)
  integer numcoe(maxhpa)
  integer ivarkern(maxker)

  integer nconpt(maxhpa),iver
  integer iconpt(maxver,maxhpa)
  real(kind=4) conpt(maxver,maxhpa)

  real(kind=4) xlaspl(maxcoe,maxhpa)
  real(kind=4) xlospl(maxcoe,maxhpa)
  real(kind=4) radspl(maxcoe,maxhpa)
  real(kind=4) coe(maxcoe,maxker)

  real(kind=4) ylmcof((maxl+1)**2,maxhpa)
  real(kind=4) wk1(maxl+1)
  real(kind=4) wk2(maxl+1)
  real(kind=4) wk3(maxl+1)

  character(len=40) varstr(maxker)

  real(kind=4) :: xcolat,xlon
  real(kind=4) :: topo410,topo650

! --- model evaluation

  integer ieval     ! --- 1 for velocity, 2 for anisotropy
  real(kind=4) :: valu(2)    ! --- valu(1) if S; valu(1)=velo, valu(2)=aniso
  real(kind=4) :: value      ! --- used in single evaluation of perturbation
  integer isel      ! --- if variable should be included
  real(kind=4) :: x,y  ! --- lat lon

! ---
  integer iker,i
  character(len=40) vstr
  integer lstr

! -------------------------------------

!       --- contributing horizontal basis functions at xlat,xlon

  y=90.0-xcolat
  x=xlon
  do ihpa=1,numhpa
      if(itypehpa(ihpa) == 1) then
        lmax=lmxhpa(ihpa)
        call ylm(y,x,lmax,ylmcof(1,ihpa),wk1,wk2,wk3)
      else if(itypehpa(ihpa) == 2) then
        numcof=numcoe(ihpa)
        call splcon(y,x,numcof,xlaspl(1,ihpa), &
              xlospl(1,ihpa),radspl(1,ihpa), &
              nconpt(ihpa),iconpt(1,ihpa),conpt(1,ihpa))
      else
        write(6,"('problem 1')")
      endif
  enddo

!         --- evaluate topography (depression) in km

  valu(1)=0. ! --- 410
  valu(2)=0. ! --- 650

  do ieval=1,2
    value=0.
    do iker=1,numker
      isel=0
      lstr=len_trim(varstr(ivarkern(iker)))
      vstr=(varstr(ivarkern(iker)))
      if(ieval == 1) then
        if(vstr(1:lstr) == 'Topo 400,') then
          isel=1
      endif
      else if(ieval == 2) then
        if(vstr(1:lstr) == 'Topo 670,') then
          isel=1
        endif
      endif

      if(isel == 1) then
            if(itypehpa(ihpakern(iker)) == 1) then
          ihpa=ihpakern(iker)
              nylm=(lmxhpa(ihpakern(iker))+1)**2
              do i=1,nylm
                value=value+ylmcof(i,ihpa)*coe(i,iker)
              enddo
            else if(itypehpa(ihpakern(iker)) == 2) then
          ihpa=ihpakern(iker)
              do i=1,nconpt(ihpa)
                iver=iconpt(i,ihpa)
                value=value+conpt(i,ihpa)*coe(iver,iker)
              enddo
            else
              write(6,"('problem 2')")
              stop
            endif ! --- itypehpa
      endif ! --- isel == 1
    enddo ! --- end of do iker=1,numker

    valu(ieval)=value
  enddo ! --- ieval

  topo410=valu(1)
  topo650=valu(2)

  end subroutine subtopo

  subroutine vbspl(x,np,xarr,splcon,splcond)
!
!---- this subroutine returns the spline contributions at a particular value of x
!
  implicit none

  integer :: np

  real(kind=4) :: xarr(np),x
  real(kind=4) :: splcon(np)
  real(kind=4) :: splcond(np)

  real(kind=4) :: r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13
  real(kind=4) :: r1d,r2d,r3d,r4d,r5d,r6d,r7d,r8d,r9d,r10d,r11d,r12d,r13d,val,vald

  real(kind=4) :: rr1,rr2,rr3,rr4,rr5,rr6,rr7,rr8,rr9,rr10,rr11,rr12
  real(kind=4) :: rr1d,rr2d,rr3d,rr4d,rr5d,rr6d,rr7d,rr8d,rr9d,rr10d,rr11d,rr12d

  integer :: iflag,interval,ik,ib

!
!---- iflag=1 ==>> second derivative is 0 at end points
!---- iflag=0 ==>> first derivative is 0 at end points
!
  iflag=1
!
!---- first, find out within which interval x falls
!
  interval=0
  ik=1
  do while(interval == 0.and.ik < np)
  ik=ik+1
  if(x >= xarr(ik-1).and.x <= xarr(ik)) interval=ik-1
  enddo
  if(x > xarr(np)) then
  interval=np
  endif

  if(interval == 0) then
!        write(6,"('low value:',2f10.3)") x,xarr(1)
  else if(interval > 0.and.interval < np) then
!        write(6,"('bracket:',i5,3f10.3)") interval,xarr(interval),x,xarr(interval+1)
  else
!        write(6,"('high value:',2f10.3)") xarr(np),x
  endif

  do ib=1,np
  val=0.
  vald=0.
  if(ib == 1) then

    r1=(x-xarr(1))/(xarr(2)-xarr(1))
    r2=(xarr(3)-x)/(xarr(3)-xarr(1))
    r4=(xarr(2)-x)/(xarr(2)-xarr(1))
    r5=(x-xarr(1))/(xarr(2)-xarr(1))
    r6=(xarr(3)-x)/(xarr(3)-xarr(1))
   r10=(xarr(2)-x)/(xarr(2)-xarr(1))
   r11=(x-xarr(1))  /(xarr(2)-xarr(1))
   r12=(xarr(3)-x)/(xarr(3)-xarr(2))
   r13=(xarr(2)-x)/(xarr(2)-xarr(1))

    r1d=1./(xarr(2)-xarr(1))
    r2d=-1./(xarr(3)-xarr(1))
    r4d=-1./(xarr(2)-xarr(1))
    r5d=1./(xarr(2)-xarr(1))
    r6d=-1./(xarr(3)-xarr(1))
   r10d=-1./(xarr(2)-xarr(1))
   r11d=1./(xarr(2)-xarr(1))
   r12d=-1./(xarr(3)-xarr(2))
   r13d=-1./(xarr(2)-xarr(1))

    if(interval == ib.or.interval == 0) then
         if(iflag == 0) then
           val=r1*r4*r10 + r2*r5*r10 + r2*r6*r11 +r13**3
           vald=r1d*r4*r10+r1*r4d*r10+r1*r4*r10d
           vald=vald+r2d*r5*r10+r2*r5d*r10+r2*r5*r10d
           vald=vald+r2d*r6*r11+r2*r6d*r11+r2*r6*r11d
           vald=vald+3.*r13d*r13**2
         else if(iflag == 1) then
           val=0.6667*(r1*r4*r10 + r2*r5*r10 + r2*r6*r11 &
                    + 1.5*r13**3)
           vald=r1d*r4*r10+r1*r4d*r10+r1*r4*r10d
           vald=vald+r2d*r5*r10+r2*r5d*r10+r2*r5*r10d
           vald=vald+r2d*r6*r11+r2*r6d*r11+r2*r6*r11d
           vald=vald+4.5*r13d*r13**2
           vald=0.6667*vald
         endif
    else if(interval == ib+1) then
         if(iflag == 0) then
           val=r2*r6*r12
           vald=r2d*r6*r12+r2*r6d*r12+r2*r6*r12d
         else if(iflag == 1) then
           val=0.6667*r2*r6*r12
           vald=0.6667*(r2d*r6*r12+r2*r6d*r12+r2*r6*r12d)
         endif
    else
      val=0.
    endif

  else if(ib == 2) then

    rr1=(x-xarr(1))/(xarr(2)-xarr(1))
    rr2=(xarr(3)-x)/(xarr(3)-xarr(1))
    rr4=(xarr(2)-x)/(xarr(2)-xarr(1))
    rr5=(x-xarr(1))/(xarr(2)-xarr(1))
    rr6=(xarr(3)-x)/(xarr(3)-xarr(1))
   rr10=(xarr(2)-x)/(xarr(2)-xarr(1))
   rr11=(x-xarr(1))  /(xarr(2)-xarr(1))
   rr12=(xarr(3)-x)/(xarr(3)-xarr(2))

    rr1d=1./(xarr(2)-xarr(1))
    rr2d=-1./(xarr(3)-xarr(1))
    rr4d=-1./(xarr(2)-xarr(1))
    rr5d=1./(xarr(2)-xarr(1))
    rr6d=-1./(xarr(3)-xarr(1))
   rr10d=-1./(xarr(2)-xarr(1))
   rr11d=1./(xarr(2)-xarr(1))
   rr12d=-1./(xarr(3)-xarr(2))

    r1=(x-xarr(ib-1))/(xarr(ib+1)-xarr(ib-1))
    r2=(xarr(ib+2)-x)/(xarr(ib+2)-xarr(ib-1))
    r3=(x-xarr(ib-1))/(xarr(ib)-xarr(ib-1))
    r4=(xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib-1))
    r5=(x-xarr(ib-1))/(xarr(ib+1)-xarr(ib-1))
    r6=(xarr(ib+2)-x)/(xarr(ib+2)-xarr(ib))
    r8=(xarr(ib)-x)/  (xarr(ib)-xarr(ib-1))
    r9=(x-xarr(ib-1))/(xarr(ib)-xarr(ib-1))
   r10=(xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib))
   r11=(x-xarr(ib))  /(xarr(ib+1)-xarr(ib))
   r12=(xarr(ib+2)-x)/(xarr(ib+2)-xarr(ib+1))

    r1d=1./(xarr(ib+1)-xarr(ib-1))
    r2d=-1./(xarr(ib+2)-xarr(ib-1))
    r3d=1./(xarr(ib)-xarr(ib-1))
    r4d=-1./(xarr(ib+1)-xarr(ib-1))
    r5d=1./(xarr(ib+1)-xarr(ib-1))
    r6d=-1./(xarr(ib+2)-xarr(ib))
    r8d=-1./  (xarr(ib)-xarr(ib-1))
    r9d=1./(xarr(ib)-xarr(ib-1))
   r10d=-1./(xarr(ib+1)-xarr(ib))
   r11d=1./(xarr(ib+1)-xarr(ib))
   r12d=-1./(xarr(ib+2)-xarr(ib+1))

    if(interval == ib-1.or.interval == 0) then
         val=r1*r3*r8 + r1*r4*r9 + r2*r5*r9
         vald=r1d*r3*r8+r1*r3d*r8+r1*r3*r8d
         vald=vald+r1d*r4*r9+r1*r4d*r9+r1*r4*r9d
         vald=vald+r2d*r5*r9+r2*r5d*r9+r2*r5*r9d
         if(iflag == 1) then
           val=val+0.3333*(rr1*rr4*rr10 + rr2*rr5*rr10 + &
                     rr2*rr6*rr11)
           vald=vald+0.3333*(rr1d*rr4*rr10+rr1*rr4d*rr10+ &
                    rr1*rr4*rr10d)
           vald=vald+0.3333*(rr2d*rr5*rr10+rr2*rr5d*rr10+ &
                    rr2*rr5*rr10d)
           vald=vald+0.3333*(rr2d*rr6*rr11+rr2*rr6d*rr11+ &
                    rr2*rr6*rr11d)
         endif
    else if(interval == ib) then
         val=r1*r4*r10 + r2*r5*r10 + r2*r6*r11
         vald=r1d*r4*r10+r1*r4d*r10+r1*r4*r10d
         vald=vald+r2d*r5*r10+r2*r5d*r10+r2*r5*r10d
         vald=vald+r2d*r6*r11+r2*r6d*r11+r2*r6*r11d
         if(iflag == 1) then
           val=val+0.3333*rr2*rr6*rr12
           vald=vald+0.3333*(rr2d*rr6*rr12+rr2*rr6d*rr12+ &
                    rr2*rr6*rr12d)
         endif
    else if(interval == ib+1) then
         val=r2*r6*r12
         vald=r2d*r6*r12+r2*r6d*r12+r2*r6*r12d
    else
         val=0.
    endif
  else if(ib == np-1) then

    rr1=(x-xarr(np-2))/(xarr(np)-xarr(np-2))
    rr2=(xarr(np)-x)/(xarr(np)-xarr(np-1))
    rr3=(x-xarr(np-2))/(xarr(np)-xarr(np-2))
    rr4=(xarr(np)-x)/(xarr(np)-xarr(np-1))
    rr5=(x-xarr(np-1))/(xarr(np)-xarr(np-1))
    rr7=(x-xarr(np-2))/(xarr(np-1)-xarr(np-2))
    rr8=(xarr(np)-x)/  (xarr(np)-xarr(np-1))
    rr9=(x-xarr(np-1))/(xarr(np)-xarr(np-1))

    rr1d=1./(xarr(np)-xarr(np-2))
    rr2d=-1./(xarr(np)-xarr(np-1))
    rr3d=1./(xarr(np)-xarr(np-2))
    rr4d=-1./(xarr(np)-xarr(np-1))
    rr5d=1./(xarr(np)-xarr(np-1))
    rr7d=1./(xarr(np-1)-xarr(np-2))
    rr8d=-1./  (xarr(np)-xarr(np-1))
    rr9d=1./(xarr(np)-xarr(np-1))

    r1=(x-xarr(ib-2))/(xarr(ib+1)-xarr(ib-2))
    r2=(xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib-1))
    r3=(x-xarr(ib-2))/(xarr(ib)-xarr(ib-2))
    r4=(xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib-1))
    r5=(x-xarr(ib-1))/(xarr(ib+1)-xarr(ib-1))
    r6=(xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib))
    r7=(x-xarr(ib-2))/(xarr(ib-1)-xarr(ib-2))
    r8=(xarr(ib)-x)/  (xarr(ib)-xarr(ib-1))
    r9=(x-xarr(ib-1))/(xarr(ib)-xarr(ib-1))
   r10=(xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib))
   r11=(x-xarr(ib))  /(xarr(ib+1)-xarr(ib))

    r1d=1./(xarr(ib+1)-xarr(ib-2))
    r2d=-1./(xarr(ib+1)-xarr(ib-1))
    r3d=1./(xarr(ib)-xarr(ib-2))
    r4d=-1./(xarr(ib+1)-xarr(ib-1))
    r5d=1./(xarr(ib+1)-xarr(ib-1))
    r6d=-1./(xarr(ib+1)-xarr(ib))
    r7d=1./(xarr(ib-1)-xarr(ib-2))
    r8d=-1./(xarr(ib)-xarr(ib-1))
    r9d=1./(xarr(ib)-xarr(ib-1))
   r10d=-1./(xarr(ib+1)-xarr(ib))
   r11d=1./(xarr(ib+1)-xarr(ib))

    if(interval == ib-2) then
         val=r1*r3*r7
         vald=r1d*r3*r7+r1*r3d*r7+r1*r3*r7d
    else if(interval == ib-1) then
         val=r1*r3*r8 + r1*r4*r9 + r2*r5*r9
         vald=r1d*r3*r8+r1*r3d*r8+r1*r3*r8d
         vald=vald+r1d*r4*r9+r1*r4d*r9+r1*r4*r9d
         vald=vald+r2d*r5*r9+r2*r5d*r9+r2*r5*r9d
         if(iflag == 1) then
           val=val+0.3333*rr1*rr3*rr7
           vald=vald+0.3333*(rr1d*rr3*rr7+rr1*rr3d*rr7+ &
                    rr1*rr3*rr7d)
         endif
    else if(interval == ib.or.interval == np) then
         val=r1*r4*r10 + r2*r5*r10 + r2*r6*r11
         vald=r1d*r4*r10+r1*r4d*r10+r1*r4*r10d
         vald=vald+r2d*r5*r10+r2*r5d*r10+r2*r5*r10d
         vald=vald+r2d*r6*r11+r2*r6d*r11+r2*r6*r11d
         if(iflag == 1) then
           val=val+0.3333*(rr1*rr3*rr8 + rr1*rr4*rr9 + &
                     rr2*rr5*rr9)
           vald=vald+0.3333*(rr1d*rr3*rr8+rr1*rr3d*rr8+ &
                    rr1*rr3*rr8d)
           vald=vald+0.3333*(rr1d*rr4*rr9+rr1*rr4d*rr9+ &
                    rr1*rr4*rr9d)
           vald=vald+0.3333*(rr2d*rr5*rr9+rr2*rr5d*rr9+ &
                    rr2*rr5*rr9d)
         endif
    else
      val=0.
    endif
  else if(ib == np) then

    r1=(x-xarr(np-2))/(xarr(np)-xarr(np-2))
    r2=(xarr(np)-x)/(xarr(np)-xarr(np-1))
    r3=(x-xarr(np-2))/(xarr(np)-xarr(np-2))
    r4=(xarr(np)-x)/(xarr(np)-xarr(np-1))
    r5=(x-xarr(np-1))/(xarr(np)-xarr(np-1))
    r7=(x-xarr(np-2))/(xarr(np-1)-xarr(np-2))
    r8=(xarr(np)-x)/  (xarr(np)-xarr(np-1))
    r9=(x-xarr(np-1))/(xarr(np)-xarr(np-1))
    r13=(x-xarr(np-1))/(xarr(np)-xarr(np-1))

    r1d=1./(xarr(np)-xarr(np-2))
    r2d=-1./(xarr(np)-xarr(np-1))
    r3d=1./(xarr(np)-xarr(np-2))
    r4d=-1./(xarr(np)-xarr(np-1))
    r5d=1./(xarr(np)-xarr(np-1))
    r7d=1./(xarr(np-1)-xarr(np-2))
    r8d=-1./  (xarr(np)-xarr(np-1))
    r9d=1./(xarr(np)-xarr(np-1))
    r13d=1./(xarr(np)-xarr(np-1))

    if(interval == np-2) then
         if(iflag == 0) then
           val=r1*r3*r7
           vald=r1d*r3*r7+r1*r3d*r7+r1*r3*r7d
         else if(iflag == 1) then
           val=0.6667*r1*r3*r7
           vald=0.6667*(r1d*r3*r7+r1*r3d*r7+r1*r3*r7d)
         endif
    else if(interval == np-1.or.interval == np) then
         if(iflag == 0) then
           val=r1*r3*r8 + r1*r4*r9 + r2*r5*r9 + r13**3
           vald=r1d*r3*r8+r1*r3d*r8+r1*r3*r8d
           vald=vald+r1d*r4*r9+r1*r4d*r9+r1*r4*r9d
           vald=vald+r2d*r5*r9+r2*r5d*r9+r2*r5*r9d
           vald=vald+3.*r13d*r13**2
         else if(iflag == 1) then
           val=0.6667*(r1*r3*r8 + r1*r4*r9 + r2*r5*r9 + &
                     1.5*r13**3)
           vald=r1d*r3*r8+r1*r3d*r8+r1*r3*r8d
           vald=vald+r1d*r4*r9+r1*r4d*r9+r1*r4*r9d
           vald=vald+r2d*r5*r9+r2*r5d*r9+r2*r5*r9d
           vald=vald+4.5*r13d*r13**2
           vald=0.6667*vald
         endif
    else
      val=0.
    endif
  else

    r1=(x-xarr(ib-2))/(xarr(ib+1)-xarr(ib-2))
    r2=(xarr(ib+2)-x)/(xarr(ib+2)-xarr(ib-1))
    r3=(x-xarr(ib-2))/(xarr(ib)-xarr(ib-2))
    r4=(xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib-1))
    r5=(x-xarr(ib-1))/(xarr(ib+1)-xarr(ib-1))
    r6=(xarr(ib+2)-x)/(xarr(ib+2)-xarr(ib))
    r7=(x-xarr(ib-2))/(xarr(ib-1)-xarr(ib-2))
    r8=(xarr(ib)-x)/  (xarr(ib)-xarr(ib-1))
    r9=(x-xarr(ib-1))/(xarr(ib)-xarr(ib-1))
   r10=(xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib))
   r11=(x-xarr(ib))  /(xarr(ib+1)-xarr(ib))
   r12=(xarr(ib+2)-x)/(xarr(ib+2)-xarr(ib+1))

    r1d=1./(xarr(ib+1)-xarr(ib-2))
    r2d=-1./(xarr(ib+2)-xarr(ib-1))
    r3d=1./(xarr(ib)-xarr(ib-2))
    r4d=-1./(xarr(ib+1)-xarr(ib-1))
    r5d=1./(xarr(ib+1)-xarr(ib-1))
    r6d=-1./(xarr(ib+2)-xarr(ib))
    r7d=1./(xarr(ib-1)-xarr(ib-2))
    r8d=-1./  (xarr(ib)-xarr(ib-1))
    r9d=1./(xarr(ib)-xarr(ib-1))
   r10d=-1./(xarr(ib+1)-xarr(ib))
   r11d=1./(xarr(ib+1)-xarr(ib))
   r12d=-1./(xarr(ib+2)-xarr(ib+1))

    if(interval == ib-2) then
         val=r1*r3*r7
         vald=r1d*r3*r7+r1*r3d*r7+r1*r3*r7d
    else if(interval == ib-1) then
         val=r1*r3*r8 + r1*r4*r9 + r2*r5*r9
         vald=r1d*r3*r8+r1*r3d*r8+r1*r3*r8d
         vald=vald+r1d*r4*r9+r1*r4d*r9+r1*r4*r9d
         vald=vald+r2d*r5*r9+r2*r5d*r9+r2*r5*r9d
    else if(interval == ib) then
         val=r1*r4*r10 + r2*r5*r10 + r2*r6*r11
         vald=r1d*r4*r10+r1*r4d*r10+r1*r4*r10d
         vald=vald+r2d*r5*r10+r2*r5d*r10+r2*r5*r10d
         vald=vald+r2d*r6*r11+r2*r6d*r11+r2*r6*r11d
    else if(interval == ib+1) then
         val=r2*r6*r12
         vald=r2d*r6*r12+r2*r6d*r12+r2*r6*r12d
    else
      val=0.
    endif
  endif
  splcon(ib)=val
  splcond(ib)=vald
  enddo

  end subroutine vbspl


  subroutine ylm(XLAT,XLON,LMAX,Y,WK1,WK2,WK3)

  implicit none

  complex TEMP,FAC,DFAC

  real(kind=4) WK1(1),WK2(1),WK3(1),Y(1),XLAT,XLON

  integer :: LMAX

!
!     WK1,WK2,WK3 SHOULD BE DIMENSIONED AT LEAST (LMAX+1)*4
!
  real(kind=4), parameter :: RADIAN = 57.2957795

  integer :: IM,IL1,IND,LM1,L

  real(kind=4) :: THETA,PHI

  THETA=(90.-XLAT)/RADIAN
  PHI=XLON/RADIAN

  IND=0
  LM1=LMAX+1

  DO IL1=1,LM1

  L=IL1-1
  CALL legndr(THETA,L,L,WK1,WK2,WK3)

  FAC=(1.,0.)
  DFAC=CEXP(CMPLX(0.,PHI))

  do IM=1,IL1
    TEMP=FAC*CMPLX(WK1(IM),0.)
    IND=IND+1
    Y(IND)=REAL(TEMP)
    IF(IM == 1) GOTO 20
    IND=IND+1
    Y(IND)=AIMAG(TEMP)
 20 FAC=FAC*DFAC
  enddo

  enddo

  end subroutine ylm

!------------------------------------

  subroutine legndr(THETA,L,M,X,XP,XCOSEC)

  implicit none

  real(kind=4) :: X(2),XP(2),XCOSEC(2)

  double precision :: SMALL,SUM,COMPAR,CT,ST,FCT,COT,X1,X2,X3,F1,F2,XM,TH

  double precision, parameter :: FPI = 12.56637062D0

  integer :: i,M,MP1,k,l,LP1

  real(kind=4) :: THETA,DSFL3,COSEC,SFL3

!!!!!! illegal statement, removed by Dimitri Komatitsch   DFLOAT(I)=FLOAT(I)

  SUM=0.D0
  LP1=L+1
  TH=THETA
  CT=DCOS(TH)
  ST=DSIN(TH)
  MP1=M+1
  FCT=DSQRT(dble(2*L+1)/FPI)
  SFL3=SQRT(FLOAT(L*(L+1)))
  COMPAR=dble(2*L+1)/FPI
  DSFL3=SFL3
  SMALL=1.D-16*COMPAR

  do I=1,MP1
    X(I)=0.
    XCOSEC(I)=0.
    XP(I)=0.
  enddo

  IF(L > 1.AND.ABS(THETA) > 1.E-5) GO TO 3
  X(1)=FCT
  IF(L == 0) RETURN
  X(1)=CT*FCT
  X(2)=-ST*FCT/DSFL3
  XP(1)=-ST*FCT
  XP(2)=-.5D0*CT*FCT*DSFL3
  IF(ABS(THETA) < 1.E-5) XCOSEC(2)=XP(2)
  IF(ABS(THETA) >= 1.E-5) XCOSEC(2)=X(2)/ST
  RETURN

 3 X1=1.D0
  X2=CT

  do I=2,L
    X3=(dble(2*I-1)*CT*X2-dble(I-1)*X1)/dble(I)
    X1=X2
    X2=X3
  enddo

  COT=CT/ST
  COSEC=1./ST
  X3=X2*FCT
  X2=dble(L)*(X1-CT*X2)*FCT/ST
  X(1)=X3
  X(2)=X2
  SUM=X3*X3
  XP(1)=-X2
  XP(2)=dble(L*(L+1))*X3-COT*X2
  X(2)=-X(2)/SFL3
  XCOSEC(2)=X(2)*COSEC
  XP(2)=-XP(2)/SFL3
  SUM=SUM+2.D0*X(2)*X(2)
  IF(SUM-COMPAR > SMALL) RETURN
  X1=X3
  X2=-X2/DSQRT(dble(L*(L+1)))

  do I=3,MP1
    K=I-1
    F1=DSQRT(dble((L+I-1)*(L-I+2)))
    F2=DSQRT(dble((L+I-2)*(L-I+3)))
    XM=K
    X3=-(2.D0*COT*(XM-1.D0)*X2+F2*X1)/F1
    SUM=SUM+2.D0*X3*X3
    IF(SUM-COMPAR > SMALL.AND.I /= LP1) RETURN
    X(I)=X3
    XCOSEC(I)=X(I)*COSEC
    X1=X2
    XP(I)=-(F1*X2+XM*COT*X3)
    X2=X3
  enddo

  end subroutine legndr

