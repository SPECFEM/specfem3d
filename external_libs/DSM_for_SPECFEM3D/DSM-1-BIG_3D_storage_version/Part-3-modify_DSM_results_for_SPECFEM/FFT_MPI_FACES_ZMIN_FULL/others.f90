subroutine pinputTra(outputDir,psvmodel,modelname,stationsinf,tlen,imin,imax,r0min, &
  r0max,r0delta,r0lat,r0lon,itranslat,mt,myrank,samplingHz,start,ending,fgauss)
  implicit none
  character(120) tmpfile
  character(120) :: dummy,outputDir,psvmodel,modelname,stationsinf,inputdir,inputIASP
  real(kind(0d0)) :: tlen,r0min,r0max,r0delta
  real(kind(0d0)) :: r0lat,r0lon,samplingHz,start,ending,fgauss
  real(kind(0d0)) :: mt(6)
  integer :: imin,imax,ifrqmin,ifrqmax,itranslat,myrank


  call get_command_argument(1,inputIASP)
  if (trim(inputIASP) == '') then
     write(*,*) 'input file not specified:  cancel '
     stop
  endif
  tmpfile = trim(inputIASP)//'.read'
  !write(tmpfile,'(a10,i5.5)') 'TmpWrkFile',myrank
  open(unit=1, file=trim(tmpfile),status='unknown')
  open(10,file=trim(inputIASP))
100 continue
  read(10,110) dummy
110 format(a120)
  if (dummy(1:1) == '#') goto 100
  if (dummy(1:3) == 'end') goto 120
  write(1,110) dummy
  goto 100
120 continue
  close(1)
  close(10)

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
  ! MT reading for TraFFT
  read(1,*) mt(1),mt(2),mt(3),mt(4),mt(5),mt(6)
  read(1,110) inputdir
  read(1,*) ifrqmin,ifrqmax
  read(1,*) samplingHz,start,ending
  read(1,*) fgauss!,fh

  close(1)
end subroutine pinputTra

subroutine tensorFFT_real(n,imin,np1,cvec,rvec,omegai,tlen)
  ! this subroutine particularly calculates the FFT of the given tensor and make a float tensor
  implicit none
  integer :: i,j,n,imin,np1,n1,m1
  complex(kind(0d0)) :: cvec(1:n,0:2*np1-1)
  real(kind(0d0)) :: rvec(1:n,0:2*np1-1)
  real(kind(0d0)), parameter :: pi = 3.141592653589793d0
  real(kind(0d0)) :: omegai, tlen,samplingHz

  samplingHz = dble(2*np1)/tlen
  do j = 1,n
     do i = imin, np1-1
        n1 = np1 +i
        m1 = np1 -i
        cvec(j,n1) = conjg(cvec(j,m1))
     enddo
  enddo




  do j = 1,n
     call cdft(4*np1,dcos(pi/(2*np1)),dsin(pi/(2*np1)), cvec(j,0:2*np1-1))
     do i = 0, 2*np1-1
        cvec(j,i) = dble(cvec(j,i))*dble(dexp(omegai*dble(i)/samplingHz))/tlen*1.d3
        rvec(j,i) = dreal(cvec(j,i))
     enddo
  enddo


  return
end subroutine tensorFFT_real



!---------------------------------------------------------------

subroutine lsmoothfinder (tlen, np0, freq, lsmooth)

  implicit none
  real(kind(0d0)) :: tlen, freq
  integer :: np0, np, lsmooth, i
  np = 1
  do while (np < np0)
     np = np*2
  enddo
  lsmooth = int(0.5*tlen*freq/dble(np))
  i = 1

  do while (i < lsmooth)
     i = i*2
  enddo

  lsmooth = i

  return

end subroutine lsmoothfinder


!---------------------------------------------------------------

!---------------------------------------------------------------

subroutine fft(data, nn, isign)

  implicit none
  integer :: nn, isign
  real(kind(0d0)) :: data(nn)
  real(kind(0d0)) :: wr, wi, wpr, wpi, wtemp, theta
  real(kind(0d0)) :: tempr, tempi
  integer :: n,i,j,m,mmax,istep

  n = 2*nn
  j = 1
  do i = 1,n,2
     if (j > i) then
        tempr = data(j)
        tempi = data(j+1)
        data(j) = data(i)
        data(j+1) = data(i+1)
        data(i) = tempr
        data(i+1) = tempi
     endif
     m = n/2

     do while ((m >= 2) .and. (j > m))
        j = j-m
        m = m/2
     enddo

     j = j+m

  enddo

  mmax = 2

  do while (n > mmax)
     istep = 2*mmax
     theta = 6.28318530717959d0/(isign*mmax)
     wpr = -2.d0*dsin(0.5d0*theta)**2
     wpi = dsin(theta)
     wr = 1.d0
     wi = 0.d0
     do m = 1, mmax, 2
        do i = m, n, istep
           j = i+mmax
           tempr = wr*data(j) -wi*data(j+1)
           tempi = wr*data(j+1)+wi*data(j)
           data(j) = data(i) - tempr
           data(j+1) = data(i+1) -tempi
           data(i) = data(i) + tempr
           data(i+1) = data(i+1) + tempi
        enddo
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
     enddo
     mmax = istep
  enddo
  return
end subroutine fft


!---------------------------------------------------------------


subroutine bwfilt (x, y, dt, n, irek, norder, f1, f2)

  ! recursive filtering of data with butterworth filter
  ! x: input array
  ! y: output array
  ! dt: time increment
  ! n: number of data points

  ! irek=0: forward filtering only
  ! irek=1: forward and backward filtering

  ! norder: order of butterworth filter
  ! norder=0: only filtering, no determination of coefficients
  ! norder < 0: no starplots of transfer function and impulse response

  ! f1: low cutoff frequency (Hz)
  ! f1=0: low pass filter

  ! f2: high cutoff frequency (Hz)
  ! f2>0.5/dt: high pass filter

  implicit none

  real(kind(0d0)), dimension(1)::x,y
  real(kind(0d0)), dimension (10) ::  a, b1, b2
  real(kind(0d0)) :: dt,f1,f2
  integer :: iunit, npoles,norder,irek,n,lx
  !real(kind(0d0)) :: x(n),y(n)

   iunit = 3

   if (norder /= 0) then
      npoles=iabs(norder)
      !determination of filter coefficients
      call bpcoeff(f1,f2,npoles, dt, a,b1, b2)
      if (norder >= 0) then
         !plot of transfer function and impuulse response
         lx = 100
         !filtering
      endif
   endif


   if (n /= 0) then
      call rekurs(x,y,n,a,b1,b2,npoles,irek)
   endif
   return
 end subroutine bwfilt

!---------------------------------------------------------------


subroutine rekurs(x,y,ndat,a,b1,b2,npoles,iflag)
  ! performs recursive filtering of data in array x of length ndat
  ! filtered output in y
  ! a, b1, b2 are the filtercoefficients previously determined in bwcoef
  ! npoles is the number of poles
  ! iflag=0: forward filtering only
  ! iflag /= 0: forward and backward filtering

  implicit none

  real(kind(0d0)), dimension(10) :: z,z1,z2 ,a,b1,b2
  real(kind(0d0)) ::  x1,x2
  integer :: ndat, npoles, iflag, n,i
  real(kind(0d0)) :: x(ndat), y(ndat)

  !forward

  x1 = 0.d0
  x2 = 0.d0

  do i = 1, npoles
     z1(i) = 0.d0
     z2(i) = 0.d0
  enddo

  do n = 1, ndat
     z(1) = a(1)*(x(n)-x2) -b1(1)*z1(1) -b2(1)*z2(1)
     do i = 2, npoles
        z(i) = a(i)*(z(i-1)-z2(i-1))-b1(i)*z1(i)-b2(i)*z2(i)
     enddo
     x2=x1
     x1=x(n)
     do i = 1, npoles
        z2(i) =z1(i)
        z1(i) =z(i)
     enddo
     y(n) = z(npoles)
  enddo

  if (iflag == 0) then
     return
  endif

  !backward

  x1 =0.d0
  x2 =0.d0

  do i = 1, npoles
     z1(i) = 0.d0
     z2(i) = 0.d0
  enddo

  do n = ndat, 1, -1
     z(1) = a(1)*(y(n)-x2)-b1(1)*z1(1)-b2(1)*z2(1)
     do i =2, npoles
        z(i) = a(i)*(z(i-1)-z2(i-1))-b1(i)*z1(i)-b2(i)*z2(i)
     enddo
     x2=x1
     x1=y(n)
     do i = 1,npoles
        z2(i)=z1(i)
        z1(i)=z(i)
     enddo
     y(n) = z(npoles)
  enddo
  return
end subroutine rekurs



!---------------------------------------------------------------


subroutine bpcoeff(f1,f2,npoles,dt,a,b1,b2)
  !determines filtercoefficients for recursive bandpassfilter

  real(kind(0d0)),dimension(10) :: a,b1,b2
  complex(kind(0d0)) :: s(20), t1,t2,p
  real(kind(0d0)), parameter :: pi = 3.141592653589793d0
  real(kind(0d0)) :: f1,f2,dt,d2,w0,w1,w2,ssum, sprod,fact1,fact2,fact3
  integer :: i,npol2,n,npoles


  if (npoles > 10) then
     stop ' npoles greater than 10: STOP '
  endif

  d2= 2.d0/dt
  w1=d2*tan(2.d0*pi*f1/d2)
  w2=d2*tan(2.d0*pi*f2/d2)
  w0=0.5*(w2-w1)

  i=1
  npol2=npoles/2+1
  do n =1,npoles
     p = cexp(cmplx(0.d0,dble(2*n-1+npoles)*pi/dble(2*npoles)))
     t1 = p*cmplx(w0,0.d0)
     t2 = sqrt(t1*t1-cmplx(w1*w2,0.d0))
     s(i)=t1+t2
     s(i+1)=t1-t2
     i=i+2
  enddo

  do n=1,npoles
     ssum=2*real(s(n))
     sprod=dble(s(n)*conjg(s(n)))
     fact1=d2*d2-d2*ssum+sprod
     fact2=2.d0*(sprod-d2*d2)
     fact3=d2*d2+d2*ssum+sprod
     a(n)=2.d0*d2*w0/fact1
     b1(n)=fact2/fact1
     b2(n)=fact3/fact1
  enddo
  return
end subroutine bpcoeff




!---------------------------------------------------------------

subroutine cdft(n, wr, wi, c)

  integer :: n, i, j, k, l, m
  real(kind(0d0)) :: wr, wi, a(0 : n - 1), wmr, wmi, wkr, wki
  real(kind(0d0)) ::wdr, wdi, ss, xr, xi
  complex(kind(0d0)) :: c(0:n/2-1)

  do i = 0, n/2-1
     a(2*i) = dble(c(i))
     a(2*i+1) = imag(c(i))
  enddo


  wmr = wr
  wmi = wi
  m = n
  do while (m > 4)
     l = m / 2
     wkr = 1
     wki = 0
     wdr = 1 - 2 * wmi * wmi
     wdi = 2 * wmi * wmr
     ss = 2 * wdi
     wmr = wdr
     wmi = wdi
     do j = 0, n - m, m
        i = j + l
        xr = a(j) - a(i)
        xi = a(j + 1) - a(i + 1)
        a(j) = a(j) + a(i)
        a(j + 1) = a(j + 1) + a(i + 1)
        a(i) = xr
        a(i + 1) = xi
        xr = a(j + 2) - a(i + 2)
        xi = a(j + 3) - a(i + 3)
        a(j + 2) = a(j + 2) + a(i + 2)
        a(j + 3) = a(j + 3) + a(i + 3)
        a(i + 2) = wdr * xr - wdi * xi
        a(i + 3) = wdr * xi + wdi * xr
     enddo
     do k = 4, l - 4, 4
        wkr = wkr - ss * wdi
        wki = wki + ss * wdr
        wdr = wdr - ss * wki
        wdi = wdi + ss * wkr
        do j = k, n - m + k, m
           i = j + l
           xr = a(j) - a(i)
           xi = a(j + 1) - a(i + 1)
           a(j) = a(j) + a(i)
           a(j + 1) = a(j + 1) + a(i + 1)
           a(i) = wkr * xr - wki * xi
           a(i + 1) = wkr * xi + wki * xr
           xr = a(j + 2) - a(i + 2)
           xi = a(j + 3) - a(i + 3)
           a(j + 2) = a(j + 2) + a(i + 2)
           a(j + 3) = a(j + 3) + a(i + 3)
           a(i + 2) = wdr * xr - wdi * xi
           a(i + 3) = wdr * xi + wdi * xr
        enddo
     enddo
     m = l
  enddo
  if (m > 2) then
     do j = 0, n - 4, 4
        xr = a(j) - a(j + 2)
        xi = a(j + 1) - a(j + 3)
        a(j) = a(j) + a(j + 2)
        a(j + 1) = a(j + 1) + a(j + 3)
        a(j + 2) = xr
        a(j + 3) = xi
     enddo
  endif
  if (n > 4) call bitrv2_old(n, a)


  do i = 0, n/2-1
     c(i) = dcmplx(a(2*i), a(2*i+1))
  enddo


end subroutine cdft



subroutine bitrv2_old(n, a)
  integer :: n, j, j1, k, k1, l, m, m2, n2
  real(kind(0d0)) :: a(0 : n - 1), xr, xi

  m = n / 4
  m2 = 2 * m
  n2 = n - 2
  k = 0
  do j = 0, m2 - 4, 4
     if (j < k) then
        xr = a(j)
        xi = a(j + 1)
        a(j) = a(k)
        a(j + 1) = a(k + 1)
        a(k) = xr
        a(k + 1) = xi
     else if (j > k) then
        j1 = n2 - j
        k1 = n2 - k
        xr = a(j1)
        xi = a(j1 + 1)
        a(j1) = a(k1)
        a(j1 + 1) = a(k1 + 1)
        a(k1) = xr
        a(k1 + 1) = xi
     endif
     k1 = m2 + k
     xr = a(j + 2)
     xi = a(j + 3)
     a(j + 2) = a(k1)
     a(j + 3) = a(k1 + 1)
     a(k1) = xr
     a(k1 + 1) = xi
     l = m
     do while (k >= l)
        k = k - l
        l = l / 2
     enddo
     k = k + l
  enddo
end subroutine bitrv2_old



!

subroutine calcijkl(omega,coef1,coef2,rho,ecKx,ecKy,ecKz,ecL,ecN,A,C,F,L,N,H)
  implicit none
  real(kind(0d0)):: rho,ecKx,ecKy,ecKz
  real(kind(0d0)):: ecL,ecN,qmu,qkappa
  complex(kind(0d0)):: A,C,F,L,N,H,Kx,Ky,Kz
  complex(kind(0d0))::omega,cijkl_acf(3,3),cijkl_lln(3),coef1,coef2
  real(kind(0d0)):: aa,bb
  real(kind(0d0)), parameter:: pi =3.1415926535897932d0

  ! l'effet d l'anelasticite

  qmu = dble(coef1)
  qkappa = dble(coef2)
  coef1 = 0.d0
  coef2 = 0.d0


  if ( qmu <= 0.d0 ) then
     coef1 = cmplx( 1.d0 )
  else
     if ( omega == 0.d0 ) then
        aa = 1.d0
     else
        aa = 1.d0 + log( omega / ( 2.d0 * pi ) ) / ( pi * qmu )
     endif
     bb = 1.d0 / ( 2.d0 * qmu)
     !coef1 = cmplx( aa, bb ) * cmplx( aa, bb )
     coef1 = cmplx(1.d0,1.d0/qmu)
  endif
  if ( qkappa <= 0.d0 ) then
     coef2 = cmplx( 1.d0 )
     !coef = cmplx( 1.d0 ) / coef2(i)
  else
     if ( omega == 0.d0 ) then
        aa = 1.d0
     else
        aa = 1.d0 + log( omega / ( 2.d0 * pi ) ) / ( pi * qkappa)
     endif
     bb = 1.d0 / ( 2.d0 * qkappa)
     !coef2 = cmplx( aa, bb ) * cmplx( aa, bb )
     !coef(i) = dcmplx( 1.d0 ) / coef2(i)
     coef2 = cmplx(1.d0,1.d0/qkappa)
  endif


  coef1 = 1.d0
  coef2 = 1.d0

  Kx = coef2 * ecKx
  Ky = coef2 * ecKy
  Kz = coef2 * ecKz
  N = coef1 * ecN
  L = coef1 * ecN

  A = Kx + 4.d0 / 3.d0 * N
  F = Ky - 2.d0 / 3.d0 * N
  C = 3.d0 * Kz - 2.d0 * F

  H = A - 2.d0 * N

  return

end subroutine calcijkl

!


subroutine findInterpolationCoefficients(r00,n,r0,ipCoef,ir00)
  implicit none
  integer :: i,n,ir00(1:3)
  real(kind(0d0)) :: r00, r0(1:n), rrsta(1:3), ipCoef(1:3)
  complex(kind(0d0)):: g(1:3), u(1:3)

  ! determination of rrsta

  u = cmplx(0.d0)

  if (r00 == r0(1)) then
     rrsta(1:3) = r0(1:3)
     ir00(1) = 1
     ir00(2) = 2
     ir00(3) = 3
  endif


  do i = 1, n-1
     if ((r0(i) < r00) .and. (r00 <= r0(i+1))) then
        if (i /= n-1) then
           rrsta(1:3) = r0(i:i+2)
           ir00(1) = i
           ir00(2) = i+1
           ir00(3) = i+2
        else
           rrsta(1:3) = r0(i-1:i+1)
           ir00(1) = i-1
           ir00(2) = i
           ir00(3) = i+1
        endif
     endif
  enddo


  do i = 1,3
     g = cmplx(0.d0)
     g(i) = 1.d0
     call interpolate(1,0,r00,rrsta,g,u(i))
  enddo
  ipCoef = dble(u)
  return
end subroutine findInterpolationCoefficients

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


!




subroutine fillinpb( nderiv,b )

  integer :: nderiv
  complex(kind(0d0)) :: b(3)

  if ( (nderiv /= 0) .and. (nderiv /= 1) .and. (nderiv /= 2) ) &
       &     stop 'invalid argument (fillinpb)'
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
subroutine glu( a, n, n1, b, eps, wk, ip, ier )
!
!        glu, glusub
!               copyright : h.hasegawa, aug. 26 1989 v.1
!
!               solves simultaneous linear equations
!               by Gaussian elimination method.
!
!        input - -
!             a(n1,n)  r *8  : 2-dim. array containing the coefficients.
!             n        i *4  : order of matrix.
!             n1       i *4  : size of array a.
!             b(n)     r *8  : 1-dim. array containing the right-hand
!                              side vector.
!             eps      r *8  : parameter to check singularity of the
!                              matrix. ( standard value 3.52d-15 )
!        output - -
!             a(n1,n)        : result of Gaussian elimination.
!             b(n)           : solution.
!             ip(n)    i *4  : pivot number.
!             ier      i *4  : = 0,  for normal execution.
!                              = 1,  for singular matrix.
!                              = 2,  for singular original matrix.
!                              = 3,  for invalid arguement.
!        working  -
!             wk(n)    r *8  : 1-dim. array.
!
  implicit none
  integer :: n,n1
  integer :: ip(n),ier
  real(kind(0d0)) :: eps
  complex(kind(0d0)) :: a(n1,n),b(n),wk(n)
  integer ::i,j,k,ipk
  complex(kind(0d0)) :: amax,aik,w,t
  !             left-hand side

  if ( eps < 0.d0 )  eps = 3.52d-15

  if ( ( n1 < n ) .or. ( n <= 0 ) ) then
     ier = 3
     print *, '  (subr. glu)  invalid argument.  n1, n =', n1, n
     return
  endif
  !            check original matrix.

  do i = 1, n
     wk(i) = dcmplx(abs(a(i,1)),0.d0)
  enddo


  do j = 2, n
     do i = 1, n
        wk(i) = max(   (abs( wk(i)) ),( abs(a(i,j))) )
     enddo
  enddo

  do i = 1, n
     if ( abs( wk(i) ) < eps ) then
        ier = 2
        print *, '  (subr. glu)  original matrix is singular.'
        return
     endif
  enddo


  ier = 0
  do k = 1, n
     !             find maximum element in the k-th column.
     amax = abs(a(k,k))
     ipk = k
     do i = k+1, n
        aik = abs(a(i,k))
        if ( abs( aik ) > abs( amax ) ) then
           ipk = i
           amax = aik
        endif
     enddo

     ip(k) = ipk

     if ( abs( amax ) > eps ) then
        if ( ipk /= k ) then
           w = a(ipk,k)
           a(ipk,k) = a(k,k)
           a(k,k) = w
        endif
        !             compute alfa
        do i = k+1, n
           a(i,k) = -a(i,k)/a(k,k)
           wk(i) = a(i,k)
        enddo
        do j = k+1, n
           if ( ipk /= k ) then
              w = a(ipk,j)
              a(ipk,j) = a(k,j)
              a(k,j) = w
           endif
           !             Gaussian elimination
           t = a(k,j)
           do i = k+1, n
              a(i,j) = a(i,j) + wk(i)*t
           enddo
        enddo

        !             matrix is singular.
     else
        ier = 1
        ip(k) = k
        do i = k+1, n
           a(i,k) = 0.0d0
        enddo
        print *, '  (subr. glu)  matrix is singular at k =', k
        return
     endif
  enddo

!             right-hand side
  !  entry glusub( a, b )
  entry glusub( a, n, n1, b, eps, wk, ip, ier )
  !forward elimination process
  do k = 1, n
     if ( ip(k) /= k ) then
        w = b(ip(k))
        b(ip(k)) = b(k)
        b(k) = w
     endif

     t = b(k)
     do i = k+1, n
        b(i) = b(i) + a(i,k)*t
     enddo
  enddo
  !          backward substitution process
  b(n) = b(n)/a(n,n)
  do k = n-1, 1, -1
     t = b(k+1)
     do i = 1, k
        b(i) = b(i) - a(i,k+1)*t
     enddo
     b(k) = b(k)/a(k,k)
  enddo
end subroutine glu


!


subroutine mtMult(g0,imt,mt)
  implicit none
  complex(kind(0d0)) :: g0
  real(kind(0d0)):: mt(3,3)
  integer :: imt

  if (imt == 1) g0 = g0 * cmplx(mt(1,1))
  if (imt == 2) g0 = g0 * cmplx(mt(1,2))
  if (imt == 3) g0 = g0 * cmplx(mt(1,3))
  if (imt == 4) g0 = g0 * cmplx(mt(2,2))
  if (imt == 5) g0 = g0 * cmplx(mt(2,3))
  if (imt == 6) g0 = g0 * cmplx(mt(3,3))

end subroutine mtMult


!
subroutine calstg_for_card(r,nzone,vrmin,vrmax,rrho,vpv,vph,vsv,vsh,eta,qmu,qkappa,array)

  ! Computing the structure grid points.
  implicit none
  integer:: nzone
  real(kind(0d0)):: r,rrho(4,nzone),vpv(4,nzone),vph(4,nzone),vsv(4,nzone),vsh(4,nzone),eta(4,nzone)
  real(kind(0d0)):: qmu(nzone), qkappa(nzone),vrmin(nzone),vrmax(nzone)
  real(kind(0d0)), parameter:: rmax  = 6371.d0
  real(kind(0d0)):: rho,ecKx,ecKy,ecKz
  real(kind(0d0)):: ecL,ecN
  real(kind(0d0)):: ecA,ecC,ecF
  real(kind(0d0)):: trho,tvpv,tvph,tvsv,tvsh,teta,coef
  complex(kind(0d0)):: coef1,coef2
  integer:: izone,i,j,itmp,jtmp
  real(kind(0d0)):: array(1:9)

  array = 0.d0
  do izone = 1, nzone
     if ((r > vrmin(izone)) .and. (r <= vrmax(izone))) then

        coef1 = cmplx(qmu(izone))
        coef2 = cmplx(qkappa(izone))
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
              coef = coef * (r / rmax )
           endif
           trho  = trho  + rrho(j,izone)  * coef
           tvpv  = tvpv  + vpv(j,izone)   * coef
           tvph  = tvph  + vph(j,izone)   * coef
           tvsv  = tvsv  + vsv(j,izone)   * coef
           tvsh  = tvsh  + vsh(j,izone)   * coef
           teta  = teta  + eta(j,izone)   * coef
        enddo
        rho = trho
        ecL  = rho * tvsv * tvsv
        ecN  = rho * tvsh * tvsh
        ecA = trho * tvph * tvph
        ecC = trho * tvpv * tvpv
        ecF = teta * ( ecA - 2.d0 * ecL )
        !kappa(itmp) = ( 4.d0 * ecA + ecC  + 4.d0 * ecF - 4.d0 * ecN(itmp) ) / 9.d0
        ecKx = ecA - 4.d0 / 3.d0 * ecN
        ecKy = ecF + 2.d0 / 3.d0 * ecN
        ecKz = ( ecC + 2.d0 * ecF ) / 3.d0


        array(1) = 1.d3 * r
        array(2) = 1.d3 * rho
        array(3) = 1.d3 * tvpv
        array(4) = 1.d3 * tvsv
        array(5) = qkappa(izone)
        array(6) = qmu(izone)
        array(7) = 1.d3 * tvph
        array(8) = 1.d3 * tvsh
        array(9) = teta
     endif
  enddo



  return
end subroutine calstg_for_card





subroutine translat(geodetic,geocentric)

  implicit none
  real(kind(0d0)),parameter ::  flattening = 1.d0 / 298.25d0
  real(kind(0d0)), parameter :: pi = 3.1415926535897932d0
  real(kind(0d0)) :: geocentric, geodetic
  real(kind(0d0)) :: tmp
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




subroutine calstg4onedepth( maxnlay,maxnzone,nzone,vrmin,vrmax,iphase,rrho,vpv,vph,vsv,vsh,eta,rmax,r,updown,ecA,ecC,ecF,ecL,ecN)

  ! Computing the structure grid points.
  implicit none
  integer:: maxnlay,maxnzone,nzone,iphase(*),updown
  real(kind(0d0)):: rrho(4,*),vpv(4,*),vph(4,*),vsv(4,*),vsh(4,*),eta(4,*),vrmin(*),vrmax(*)
  real(kind(0d0)):: rmax,r
  real(kind(0d0)):: rho,kappa,ecKx,ecKy,ecKz
  real(kind(0d0)):: mu,ecL,ecN
  real(kind(0d0)):: ecA,ecC,ecF
  real(kind(0d0)):: trho,tvpv,tvph,tvsv,tvsh,teta,coef
  integer:: izone,i,j,itmp,jtmp,spn

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
!!
!!  VM --- VM 2011.10
!!
!!  Coordonnees du stress
!!
!!  1 : rr : zz (bas - haut)
!!  2 : tt : XX (radial = source->recepteur)
!!  3 : pp : YY (?signe + vers le Nord?) (c'est apres la rotation que l'on choisi la sens)
!!  4 : rt : ZX
!!  5 : rp : ZY
!!  6 : tp : XY
!!
!!
!!
!!
!!
!!
!!
!!
!!
!!
!!
!!
!!
subroutine Local2chunk(vec,n,np1,Orign_Lon_Chunk,Origin_Lat_Chunk,azi_chunk,Lat_current,Lon_current,Bazi,imin,imax)
  implicit none
  integer n,np1
  complex(kind(0d0)) vec(n,0:2*np1-1)
  real(kind(0d0)) Bazi_chunk,Bazi,Orign_Lon_Chunk,Origin_Lat_Chunk,Lat_current,Lon_current,azi_chunk!,normal(3)
  complex(kind(0d0)) normal_RT(3),normal_local(3),Pass_local_chunk(3,3),Pass_chunk_local(3,3),normal(3)
  complex(kind(0d0)) Rot_RT_NE(3,3),Pass_Chunk_global(3,3),Pass_global_chunk(3,3),Rot_azi_chunk(3,3)
  complex(kind(0d0)) Pass_local_global(3,3),Pass_global_local(3,3),Stress(3,3),Disp(3)
  complex(kind(0d0)) Trac_local(3),Trac_RT(3),Trac_chunk(3),Disp_local(3),Disp_chunk(3),Trac_global(3),Disp_global(3)
  integer i,j,k
  integer icomp,itime,imin,imax

  call  RotationMatrix_RT2Local(Rot_RT_NE,Bazi)
  Bazi_chunk=270.d0-azi_chunk ! car on fait cette operation dans la subroutine ci dessous ca annule le truc, reflechir
  ! a un truc pour mieux faire !! !! VM VM
  call RotationMatrix_RT2local(Rot_azi_chunk,Bazi_chunk)  !! c'est l'azimuth du chunk


  call  RotationMatrix(Origin_Lat_Chunk,Orign_Lon_Chunk,Pass_Chunk_global)
  call  RotationMatrix(Lat_current,lon_current,Pass_local_global)

  call InverseMatrix(Pass_local_global,Pass_global_local)
  call InverseMatrix(Pass_chunk_global,Pass_global_chunk)


  call  MatrixMultInverse(Pass_local_global,Pass_global_chunk,Pass_local_chunk)
  call  MatrixMultInverse(Pass_chunk_global,Pass_global_local,Pass_chunk_local)
   !write(*,*) '-----------------------------------------'
   !write(*,*) Rot_local_global
   !write(*,*) Rot_Chunk_global
   !write(*,*) Rot_local_chunk
   !write(*,*) '-----------------------------------------'

  ! calcul de la normale dans repere RT
  !normal_local = inv(Rot_local_chunk) * normal

 !*! call MultDir(Pass_chunk_local,normal,normal_local)

   normal_local(1) = dcmplx(-1.d0)
   normal_local(2) = dcmplx(0.d0)
   normal_local(3) = dcmplx(0.d0)

  !call MultDir(Pass_chunk_local,normal,normal_local)


  !if (itime>1000 .and. itime < 1002) then
  !write(100,*) '----------NORMALES---------------------'
  !write(100,*) 'normal chunk ' ,normal_local
  !write(100,*) 'normal chunk rot  ',  Disp_chunk
  !endif
  ! normal_RT = inv(Rot_RT_NE) * normal_local
  call MultInv(Rot_RT_NE,normal_local,normal_RT)
  !if (itime>10 .and. itime < 12) then
  !write(100,*) 'normal RT ', normal_RT
  !write(*,*)
  !endif
  ! write(100,*) Pass_chunk_local(1,1),Pass_chunk_local(1,2),Pass_chunk_local(1,3)
  ! write(100,*) Pass_chunk_local(2,2),Pass_chunk_local(2,3),Pass_chunk_local(3,3)
  !write(100,*)
  !write(100,*) Pass_local_global
  !write(100,*)
  !write(100,*) Pass_global_chunk
  !write(100,*)
  !write(100,*) Pass_chunk_local
     do itime=imin,imax !0,2*np1-1

       ! repere ZRT
      !write(100,*) '------------- ',itime
      !write(100,*) vec(:,itime)
       Stress(1,1)=(vec(1,itime))
       Stress(1,2)=(vec(4,itime))
       Stress(1,3)=(vec(5,itime))
       Stress(2,1)=(vec(4,itime))
       Stress(2,2)=(vec(2,itime))
       Stress(2,3)=(vec(6,itime))
       Stress(3,1)=(vec(5,itime))
       Stress(3,2)=(vec(6,itime))
       Stress(3,3)=(vec(3,itime))

       ! repere ZRT
       Disp(1)=(vec(7,itime))
       Disp(2)=(vec(8,itime))
       Disp(3)=(vec(9,itime))

     !if (itime>400 .and. itime < 402) then
     !*write(*,*) 'vec ',vec(:,itime)
     !*write(*,*) 'Stress ',Stress
     !write(*,*) vec(1,itime),vec(4,itime),vec(5,itime)
     !write(*,*) '------------------------'
    !endif
      !!!!!!!!!!!!!!!!!   A FAIRE !!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!! PB de signe avec la composante  NS ou SN

      ! Trac_RT=Stress*normal_RT
      !if (itime>1000 .and. itime < 1002)
      !write(100,*) '------------normal RT',normal_RT
      call MultDir(Stress,normal_RT,Trac_RT)
      !if (itime>100 .and. itime < 102)
      !write(100,*) 'Trac RT ',Trac_RT
      !! ROTATION DU REPRE RT -> CHUNK
      ! Trac_local = Rot_RT_NE * Trac_RT
      call MultDir(Rot_RT_NE,Trac_RT,Trac_local)
      !if (itime>100 .and. itime < 102)
      !write(100,*) 'Trac local ',Trac_local
      ! Trac_chunk = Rot_local_chunk * Trac_local
      call MultDir(Pass_local_chunk,Trac_local,Trac_global)
      !write(100,*) 'Trac global ',Trac_global
      ! rotation azimuthale du chunk
       call MultDir(Rot_azi_chunk,Trac_global,Trac_chunk)
      !if (itime>400+imin .and. itime < 402+imin)
      !write(100,*) 'Trac chunk ',Trac_chunk
      ! Disp_local = Rot_RT_NE * Disp
      ! Disp_chunk = Rot_local_chunk * Disp_local

      ! traction (1,2,3) (X,Y,Z chunk)
       vec(1,itime)=(Trac_chunk(1))
       vec(2,itime)=(Trac_chunk(2))
       vec(3,itime)=(Trac_chunk(3))
       !vec(4,itime)=
       !vec(5,itime)=
       !vec(6,itime)=
       !if (itime==imax) then
       !  write(100,*)
       !  write(100,*) Disp_local
       !endif

       !!  ROTATION DU REPRE RT -> CHUNK
       call MultDir(Rot_RT_NE,Disp,Disp_local)
       call MultDir(Pass_local_chunk,Disp_local,Disp_global)
       call MultDir(Rot_azi_chunk,Disp_global,Disp_chunk)

       ! save for debug
        vec(4,itime)=vec(7,itime)
        vec(5,itime)=vec(8,itime)
        vec(6,itime)=vec(9,itime)
       ! velocity (X,Y,Z chunk)
       vec(7,itime)=(Disp_chunk(1))
       vec(8,itime)=(Disp_chunk(2))
       vec(9,itime)=(Disp_chunk(3))
       !if (itime==imax) then
       !  write(100,*) Disp_chunk
       !endif


     enddo

   !*write(*,*)
   !*write(*,*)
   !*write(*,*)
   !*write(*,*) '-------------------------------------------'
   !*


end subroutine Local2chunk


subroutine InverseMatrix(R,Ri)
  implicit none
  complex(kind(0d0)) R(3,3), Ri(3,3),Rtmp(3,3)

!  permutation col 2 et 3
  Rtmp(:,1) = R(:,1)
  Rtmp(:,2) = R(:,3)
  Rtmp(:,3) = R(:,2)
! transposition
  Ri(:,:) = Rtmp(:,:)
  Ri(2,3) = Rtmp(3,2)
  Ri(3,2) = Rtmp(2,3)
  Ri(1,2) = Rtmp(2,1)
  Ri(2,1) = Rtmp(1,2)
  Ri(1,3) = Rtmp(3,1)
  Ri(3,1) = Rtmp(1,3)
! permutaion col 2 et 3
  Rtmp(:,:) = Ri(:,:)
  Ri(2,:) = Rtmp(3,:)
  Ri(3,:) = Rtmp(2,:)

end subroutine InverseMatrix

subroutine MultDir(R,x,y)
 implicit none
 complex(kind(0d0)) R(3,3),x(3),y(3)
 integer i,j

  y(:) = 0.d0
  do i=1,3
     do j=1,3
       y(i) = y(i) + R(i,j) * x(j)
     enddo
  enddo
end subroutine MultDir



subroutine MultInv(R,x,y)
 implicit none
 complex(kind(0d0)) R(3,3),x(3),y(3)
 integer i,j

  y(:) = 0.d0
  do i=1,3
     do j=1,3
       y(i) = y(i) + R(j,i) * x(j)
     enddo
  enddo
end subroutine MultInv



subroutine RotationMatrix_RT2Local(R,bazi)
!!
!!
!!  matrice de rotation (R,T) -> (WE,SN)
!!
!!
!!
!!
 implicit none
 complex(kind(0d0)) R(3,3)
 real(kind(0d0)) bazi,deg2rad,t

 deg2rad= 3.141592653589793d0/180.d0


 !t = 180.d0 - bazi

 t =  270.d0 - bazi


 R(1,1)=dcmplx(1.d0)
 R(1,2)=dcmplx(0.d0)
 R(1,3)=dcmplx(0.d0)

 R(2,1)=dcmplx(0.d0)
 R(2,2)=dcmplx(dcos(t*deg2rad))
 R(2,3)=dcmplx(-dsin(t*deg2rad))

 R(3,1)=dcmplx(0.d0)
 R(3,2)=dcmplx(dsin(t*deg2rad))
 R(3,3)=dcmplx(dcos(t*deg2rad))



end subroutine  RotationMatrix_RT2Local





subroutine MatrixMultInverse(R1,R2,R)

  implicit none
  complex(kind(0d0)) R(3,3),R1(3,3),R2(3,3)
  integer i,j,k

  R(:,:) = dcmplx(0.d0)

  do i=1,3
    do j=1,3
      do k=1,3
        R(i,j) = R(i,j)  + R2(i,k) * R1(k,j)
       enddo
     enddo
   enddo

end subroutine MatrixMultInverse



subroutine RotationMatrix(lat,lon,R)

 implicit none
 real(kind(0d0)) lat,lon
 real(kind(0d0))  t,p,deg2rad
 complex(kind(0d0))  cp,sp,st,ct,R(3,3)
 integer i,j

 deg2rad = 3.141592653589793d0/180.d0

 t =  ( 90.d0 - lat ) * deg2rad
 p =  lon * deg2rad

 ct = dcmplx(dcos(t))
 st = dcmplx(dsin(t))

 cp = dcmplx(dcos(p))
 sp = dcmplx(dsin(p))

 R(1,1)=st*cp
 R(2,1)=st*sp
 R(3,1)=ct

 !! vm changement de signe
 R(1,3)= - ct*cp
 R(2,3)= - ct*sp
 R(3,3)=  st
 !!
 R(1,2)=-sp
 R(2,2)=cp
 R(3,2)=dcmplx(0.d0)


end subroutine  RotationMatrix


subroutine InputForTraction(Orign_Lon_Chunk,Origin_Lat_Chunk,azi_chunk,delta_lat,delta_lon)
 implicit none
 real(kind(0.d0)) Orign_Lon_Chunk,Origin_Lat_Chunk,normal_real(3),azi_chunk,delta_lat,delta_lon
 !complex(kind(0.d0)) normal(3)

 open(10,file='OrigRepSpecfm')
 read(10,*) Orign_Lon_Chunk,Origin_Lat_Chunk
 read(10,*) azi_chunk,delta_lon,delta_lat
 close(10)
 !open(10,file='normal.txt')
 !read(10,*) normal(1),normal(2),normal(3)
 !read(10,*) normal_real(2),normal_real(3),normal_real(1)
 !normal(1)=dcmplx(normal_real(1))
 !normal(2)=dcmplx(normal_real(2))
 !normal(3)=dcmplx(normal_real(3))
 !!            theta   phi           r
 !!      Y: comp S-N/*/
 !!      Z :comp B-H/*/
 !!      X :comp W-E
 !close(10)

end subroutine InputForTraction




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
!
 double precision delta,azi,SourceLat,RecLat,SourceLong,RecLong,bazi
 double precision pi,tete,phie,tets,phis,ctete,ctets,stete,stets,cdel,sdelta,cal,sal,gamma,d2r,r2d,cbe,sbe
 parameter(pi=3.141592653589793d0,d2r=Pi/180.d0,r2d=1.d0/d2r)


!-----------------------------------------------------------------
!  Transform from latitude to colatitude                          |
! Transform from degrees to radians and abbreviate sine and       |
! cosine.                                                         |
!-----------------------------------------------------------------

      tets = 90.d0 - (SourceLat)
      tete = 90.d0 - (RecLat)
      if (tets == 0.d0) then
         delta = tete
         azi   = 180.-Reclong
         bazi = 0.
         return
      else if (tets == 180.d0 ) then
         delta = 180. - (tete)
         azi   = RecLong
         bazi  = 180.
         return
      else if (tete == 0.d0 ) then
         delta = (tets)
         azi = 0.
         bazi = 180.
         return
      endif
      tets = tets*d2r
      phis = (SourceLong)*d2r
      tete = tete*d2r
      phie = (RecLong)*d2r
      ctets = dcos(tets)
      ctete = dcos(tete)
      stets = dsin(tets)
      stete = dsin(tete)

!-----------------------------------------------------------------
! Use cosine theorem on the sphere and check for antipode.        |
!-----------------------------------------------------------------

      cdel  = ctets*ctete+stets*stete*dcos(phie-phis)
      delta = (dacos(cdel)*r2d)
      sdelta=dsqrt((1.d0-cdel)*(1.d0+cdel))
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


 end subroutine epitra


 subroutine DistribDepth(n,n_global,myrank,nbproc,irmin,irmax)
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



   !write(100,*) myrank, irmin, irmax
   !if (myrank <= mod(n_global,nbproc) .and. mod(n_global,nbproc) /= 0) then
   !   n = n + 1
   !   irmin = irmin + myrank
   !   irmax = irmax + myrank
   !endif



 end subroutine DistribDepth


!----------------------------------------------------------------------
      subroutine filtbutter(x,npts,dt,fl,fh,norder)


        implicit none

        integer npts,norder
        real(kind(0d0)) dt,fl,fh,fnyq,x(npts)
        real(kind(0d0)) b(2,10),bfact
!-----
!                       set up filter coefficients which
!                       depend upon dt
!-----
        fnyq = 0.5/dt
        if (fh < 0.0 .or. fh > fnyq)fh = fnyq
        call filcof(fl,fh,dt,norder,b,bfact)
!-----
!                   convolve

!-----
!-----
!                   initialize the time series
!-----
        call filt(x,npts,b,bfact,norder)
!-----
!                   output the filtered trace
!-----


        return
      end subroutine filtbutter
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
      subroutine filt(dat,n,b,bfact,npole)

        implicit none

!-----
!       data    R*4 - array to be filtered, same array used for
!                   input and output
!       n   I*4 - number of points in time series to be filtered
!-----


        integer i,k,n,npole
        real(kind(0d0)) y(3,10),x(3)
        real(kind(0d0)) b(2,10),bfact
        real(kind(0d0)) dum
        real(kind(0d0)) dat(n)

        x = 0.d0
        y = 0.d0
        do k = 1,n

           x(1) = dat(k)
           y(1,1) = x(1) - x(3) - b(1,1)*y(2,1) - b(2,1)*y(3,1)
           do i = 2,npole
              y(1,i) = y(1,i-1)-y(3,i-1)-b(1,i)*y(2,i)-b(2,i)*y(3,i)
           enddo
           !-----
           !           update history
           !-----
           x(3) = x(2)
           x(2) = x(1)
           do i = 1,npole
              y(3,i) = y(2,i)
              y(2,i) = y(1,i)
           enddo
           dum = bfact * y(1,npole)
           if (dabs(dum) < 1.d-36)dum = 0.d0
           dat(k) = dum
        !       print *,k,x(1),data(k)
        enddo

        return
      end subroutine filt
!------------------------------------------------------------------------
!-------------------------------------------------------------------------
      subroutine filcof(fl,fh,dt,norder,b,bfact)

        implicit none

        integer i,npole,npole2,norder
        real(kind(0d0)) fl,fh,dt,pi
        real(kind(0d0))  b(2,10),bfact
        real(kind(0d0)) wdl,wdh,x,wa,s1,sr,si,xn,xd,ang
        complex(kind(0d0)) s2,az,cz,z1,z2,zc,xc,bfac
        !DOUBLE COMPLEX CSQRT

        pi=3.1415926535d0
        bfac = dcmplx(1.0d+00,0.0d+00)
        npole = norder
        wdl=pi*fl*dt
        wdh=pi*fh*dt
        xn=dcos(wdl+wdh)
        xd=dcos(wdh-wdl)
        x=xn/xd
        wa = dsin(wdh-wdl)/dcos(wdh-wdl)
        s1=dabs(wa)
        npole2 = (npole + 1 ) /2
        do i=1,npole2,1
           ang = 0.5d0*pi *(1.d0 + (2.d0*i-1.d0)/npole )
           sr = -s1 * dcos(ang)
           si = -s1 * dsin(ang)
           s2 = dcmplx(sr,si)
           xc=dcmplx(x,0.0d+00)
           az=(1.0d+00,0.0d+00)+s2
           cz=(1.0d+00,0.0d+00)-s2
           zc=xc*xc-az*cz
           zc=CDSQRT(zc)
           z1=(xc+zc)/az
           z2=(xc-zc)/az
           if (i /= (npole+1-i)) then
              b(1,i) = - 2.d0*dreal(z1)
              b(2,i) = CDABS(z1)**2
              b(1,npole+1-i) = -2.d0*dreal(z2)
              b(2,npole+1-i) = CDABS(z2)**2
              bfac = bfac * s2/(1.d0+s2)
              bfac = bfac * dconjg(s2)/(1.d0+dconjg(s2))
           else
              b(1,i) = -2.d0*x/(1.d0 + s1)
              b(2,i) = ( 1.d0 - s1)/(1.d0 + s1)
              bfac = bfac * s1/(1.d0 + s1)
           endif
        enddo

        bfact = dreal(bfac)

        return
      end subroutine filcof


!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
      subroutine cdft_other(n, wr, wi, a)
      integer n, i, j, k, l, m
      real(kind(0d0)) wr, wi, a(0 : n - 1), wmr, wmi, wkr, wki, wdr, wdi, &
         ss, xr, xi
      wmr = wr
      wmi = wi
      m = n
      do while (m > 4)
          l = m / 2
          wkr = 1
          wki = 0
          wdr = 1 - 2 * wmi * wmi
          wdi = 2 * wmi * wmr
          ss = 2 * wdi
          wmr = wdr
          wmi = wdi
          do j = 0, n - m, m
              i = j + l
              xr = a(j) - a(i)
              xi = a(j + 1) - a(i + 1)
              a(j) = a(j) + a(i)
              a(j + 1) = a(j + 1) + a(i + 1)
              a(i) = xr
              a(i + 1) = xi
              xr = a(j + 2) - a(i + 2)
              xi = a(j + 3) - a(i + 3)
              a(j + 2) = a(j + 2) + a(i + 2)
              a(j + 3) = a(j + 3) + a(i + 3)
              a(i + 2) = wdr * xr - wdi * xi
              a(i + 3) = wdr * xi + wdi * xr
          enddo
          do k = 4, l - 4, 4
              wkr = wkr - ss * wdi
              wki = wki + ss * wdr
              wdr = wdr - ss * wki
              wdi = wdi + ss * wkr
              do j = k, n - m + k, m
                  i = j + l
                  xr = a(j) - a(i)
                  xi = a(j + 1) - a(i + 1)
                  a(j) = a(j) + a(i)
                  a(j + 1) = a(j + 1) + a(i + 1)
                  a(i) = wkr * xr - wki * xi
                  a(i + 1) = wkr * xi + wki * xr
                  xr = a(j + 2) - a(i + 2)
                  xi = a(j + 3) - a(i + 3)
                  a(j + 2) = a(j + 2) + a(i + 2)
                  a(j + 3) = a(j + 3) + a(i + 3)
                  a(i + 2) = wdr * xr - wdi * xi
                  a(i + 3) = wdr * xi + wdi * xr
              enddo
          enddo
          m = l
      enddo
      if (m > 2) then
          do j = 0, n - 4, 4
              xr = a(j) - a(j + 2)
              xi = a(j + 1) - a(j + 3)
              a(j) = a(j) + a(j + 2)
              a(j + 1) = a(j + 1) + a(j + 3)
              a(j + 2) = xr
              a(j + 3) = xi
          enddo
      endif
      if (n > 4) call bitrv2(n, a)
      end subroutine cdft_other
!--------------------------
      subroutine bitrv2(n, a)
      integer n, j, j1, k, k1, l, m, m2, n2
      real(kind(0d0)) a(0 : n - 1), xr, xi
      m = n / 4
      m2 = 2 * m
      n2 = n - 2
      k = 0
      do j = 0, m2 - 4, 4
          if (j < k) then
              xr = a(j)
              xi = a(j + 1)
              a(j) = a(k)
              a(j + 1) = a(k + 1)
              a(k) = xr
              a(k + 1) = xi
          else if (j > k) then
              j1 = n2 - j
              k1 = n2 - k
              xr = a(j1)
              xi = a(j1 + 1)
              a(j1) = a(k1)
              a(j1 + 1) = a(k1 + 1)
              a(k1) = xr
              a(k1 + 1) = xi
          endif
          k1 = m2 + k
          xr = a(j + 2)
          xi = a(j + 3)
          a(j + 2) = a(k1)
          a(j + 3) = a(k1 + 1)
          a(k1) = xr
          a(k1 + 1) = xi
          l = m
          do while (k >= l)
              k = k - l
              l = l / 2
          enddo
          k = k + l
      enddo
      end subroutine bitrv2



subroutine epitra1(SourceLat,SourceLong,RecLat,RecLong,delta,azi,bazi)
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
      !if (azi>180d.0) azi =


end subroutine epitra1

