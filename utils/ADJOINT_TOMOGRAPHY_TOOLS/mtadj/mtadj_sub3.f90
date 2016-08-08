! this module contains the following subroutine/functions:
! subroutine fft(n,xi,zign,dt)
! subroutine fftinv(npow,s,zign,dt,r)

! subroutine staper(nt, fw, nev, v, ndim, a, w)
!   subroutine tsturm(nt,n,a,b,w,nev,r,ndim,ev,ipar)
!   subroutine root(u,el,elam,a,bb,w,n,ik)
! subroutine costaper(ipoint, ndata, tas)
! subroutine boxcar(ipoint, ndata, tas)

module mtadj_sub3

  implicit none

! TOLERRANCE CONTROL
  double precision, parameter ::  TOL=1e-7

contains

!------------------------------------------------------------------
  subroutine fft(n,xi,zzign,dt)
! Fourier transform
! This inputs AND outputs a complex function.
! The convention is FFT --> e^(-iwt)
! numerical factor for Plancherel theorem: planch_fac = dble(NPT * dt * dt)
!------------------------------------------------------------------
      complex*16, dimension(*) :: xi
      integer :: n
      double precision :: dt

      double precision, parameter :: PI = 3.141592653589793d+00
      complex*16 :: wk, hold, q
      double precision :: m(25)
      double precision :: zzign,zign,flx,v
      integer :: lblock,k,fk,jh,ii,istart
      integer :: l,iblock,nblock,i,lbhalf,j,lx

      ! sign must be +1. or -1.
      if (zzign >= 0.) then
        zign = 1.
      else
        zign = -1.
      endif

      lx = 2**n
      do 1 i=1,n
    1 m(i) = 2**(n-i)
      do 4 l=1,n
      nblock = 2**(l-1)
      lblock = lx/nblock
      lbhalf = lblock/2
      k = 0
      do 4 iblock=1,nblock
      fk = k
      flx = lx

      v = zign*2.*PI*fk/flx         ! Fourier convention

      wk = cmplx(cos(v),-sin(v))    ! sign change to -sin(v) 17-Nov-2006
      istart = lblock*(iblock-1)

      do 2 i=1,lbhalf
      j  = istart+i
      jh = j+lbhalf
      q = xi(jh)*wk
      xi(jh) = xi(j)-q
      xi(j)  = xi(j)+q
    2 continue

      do 3 i=2,n
      ii = i
      if (k < m(i)) goto 4
    3 k = k-m(i)
    4 k = k+m(ii)
      k = 0
      do 7 j=1,lx
      if (k < j) goto 5
      hold = xi(j)
      xi(j) = xi(k+1)
      xi(k+1) = hold
    5 do 6 i=1,n
      ii = i
      if (k < m(i)) goto 7
    6 k = k-m(i)
    7 k = k+m(ii)

      ! final steps deal with dt factors
      if (zign > 0.) then       ! FORWARD FFT
         do i = 1,lx
            xi(i) = xi(i)*dt   ! multiplication by dt
         enddo

      else                     ! REVERSE FFT
         flx = flx*dt
         do i = 1,lx
            xi(i) = xi(i)/flx  ! division by dt
         enddo
      endif

  end subroutine fft

!------------------------------------------------------------------
  subroutine fftinv(npow,s,zzign,dt,r)
! inverse Fourier transform -- calls fft
!------------------------------------------------------------------

      !implicit real*8(a-h,o-z)
      !dimension r(4096*4)
      !complex s(4096*4)

      complex*16, intent(in) :: s(*)
      double precision, intent(out) :: r(*)   ! note this is REAL

      double precision :: dt,zzign,zign
      integer :: npow, nsmp, nhalf, i

      nsmp = 2**npow
      nhalf = nsmp/2
      call rspec(s,nhalf)   ! re-structuring

      zign=zzign
      call fft(npow,s,zign,dt)    ! Fourier transform

      do i = 1,nsmp
        r(i) = real(s(i))     ! REAL part
      enddo

  end subroutine fftinv

!------------------------------------------------------------------
  subroutine rspec(s,np2)
!------------------------------------------------------------------

      !implicit real*8(a-h,o-z)
      !complex s(4096*4)

      complex*16 :: s(*)
      integer :: np2,n,n1,i

      n = 2*np2
      n1 = np2+1

      s(n1) = 0.
!     s(1)  = 0.
      s(1)  = cmplx( real(s(1)),0.)

      do i = 1,np2
         s(np2+i) = conjg(s(np2+2-i))
      enddo

  end subroutine rspec

!------------------------------------------------------------------
  subroutine staper(nt, fw, nev, v, ndim, a, w)
!------------------------------------------------------------------
!$$$$ calls tsturm, root
!  Slepian - Thomson multi-taper procedure
!  Slepian, D.     1978  Bell Sys Tech J v57 n5 1371-1430
!  Thomson, D. J.  1982  Proc IEEE v70 n9 1055-1096
!    nt    the number of points in the series
!    fw    the time-bandwidth product (number of Rayleigh bins)
!    nev   the desired number of tapers
!    v     the eigenvectors (tapers) are returned in v(.,nev)
!    a, w  work arrays dimensioned at least nt long (nt+1, nt odd)
!    a(1..nev) contains bandwidth retention factors on output.
!  The tapers are the eigenvectors of the tridiagonal matrix sigma(i,j)
!  [see Slepian(1978) eq 14 and 25.] They are also the eigenvectors of
!  the Toeplitz matrix eq. 18. We solve the tridiagonal system in
!  tsturm for the tapers and use them in Slepians eq 18 to get the
!  bandwidth retention factors (i.e. the eigenvalues) Thomson's
!  normalisation is used with no attention to sign.
      !implicit real*8(a-h,o-z)
      !dimension a(*),w(*),v(ndim,*)
      !parameter (pi=3.14159265358979d0,r2=1.414213562373095d0)

      integer :: nt, nev, ndim
      double precision :: fw
      double precision :: v(ndim,*), a(*), w(*)

      double precision, parameter :: PI = 3.141592653589793d+00
      integer :: i,j,k,m
      integer :: nxi, lh, lp1, neven, nodd, ntot, kk, kmax, nlow, nup
      double precision :: r2,om,com,hn,asav,rbd,dc,sm,s,sn,vmax

      !-------------------------

      r2 = sqrt(2.)

      if (nt < 2) return
      nxi=mod(nt,2)
      lh=(nt/2)+nxi
      lp1=nt+1
      om=2.*PI*fw/nt
      com=cos(om)
      hn=0.5*dble(lp1)
      do 10 i=1,lh
        a(i)=com*(i-hn)**2
   10   w(i)=0.5*dble(i*(nt-i))
      if (nxi == 0) then
        asav=a(lh)-w(lh)
        a(lh)=a(lh)+w(lh)
        rbd=1./(a(lh)+w(lh-1))
      else
        asav=w(lh-1)
        rbd=1./(w(lh)+w(lh-1))
        w(lh-1)=r2*w(lh-1)
      endif
      do 15 i=1,lh
        a(i+lh)=w(i)*rbd
        w(i)=a(i+lh)**2
   15   a(i)=a(i)*rbd
      neven=max0((nev+1)/2,1)
      nodd=nev-neven
!  Do the even tapers
      call tsturm(nt,lh,a,a(lh+1),w,neven,v,ndim,w(lh+1),0)
      do 20 i=1,neven
        k=2*i-1
        if (nxi == 1) v(lh,k)=r2*v(lh,k)
          do 20 j=1,lh
   20     v(lp1-j,k)=v(j,k)
      if (nodd <= 0) goto 34
!  Do the odd tapers
      if (nxi == 0) then
        a(lh)=asav*rbd
      else
        a(nt)=asav*rbd
        w(lh-1)=asav*asav
      endif
      call tsturm(nt,lh-nxi,a,a(lh+1),w,nodd,v,ndim,w(lh+1),1)
      do 30 i=1,nodd
        k=2*i
        if (nxi == 1) v(lh,k)=0.
          do 30 j=1,lh
   30     v(lp1-j,k)=-v(j,k)
   34 ntot=neven+nodd
!  Calculate bandwidth retention parameters
      dc=2.*com
      sm=0.
      s=sin(om)
      w(1)=om/PI
      w(2)=s/PI
      do 35 j=3,nt
        sn=dc*s-sm
        sm=s
        s=sn
   35   w(j)=s/(PI*(j-1))
      do 55 m=1,ntot
        vmax=abs(v(1,m))
        kmax=1
        do 40 kk=2,lh
          if (abs(v(kk,m)) <= vmax) goto 40
          kmax=kk
          vmax=abs(v(kk,m))
   40     continue
        a(m)=0.
        nlow=kmax-1
          do 45 j=1,nlow
   45     a(m)=a(m)+w(j+1)*v(nlow+1-j,m)
        nup=nt-nlow
          do 50 j=1,nup
   50     a(m)=a(m)+w(j)*v(nlow+j,m)
   55 a(m)=a(m)/v(kmax,m)
      return

  end subroutine staper

!------------------------------------------------------------------
  subroutine tsturm(nt,n,a,b,w,nev,r,ndim,ev,ipar)
!------------------------------------------------------------------
!$$$$ calls root
!  Uses bisection and Sturm counting to isolate the eigenvalues of the
!  symmetric tridiagonal matrix with main diagonal a(.) and sub/super
!  diagonal b(.).  Newton's method is used to refine the eigenvalue in
!  subroutine root then direct recursion is used to get the eigenvector
!  as this is always stable.  Note  ipar=0 for even tapers   =1 for odd
!  tapers
      !implicit real*8(a-h,o-z)
      !parameter (epsi=1.d-15,epsi1=5.d-15)
      !dimension a(*),b(*),ev(*),w(*),r(ndim,*)

      double precision, parameter :: epsi = 1.d-15, epsi1 = 5.d-15

      double precision, dimension(ndim) :: a, b, w, ev
      double precision, dimension(ndim,*) :: r
      integer :: nt,n,ndim,nev,ipar

      double precision, dimension(ndim) :: bb
      double precision :: q,el,elam,u,umeps,x,ddot,rnorm
      integer :: i,j,ik,iag,m,jk,jm1

      !-------------------------

      if (n <= 0 .or. nev <= 0) return
      umeps=1.-epsi
      do 5 i=1,nev
    5 ev(i)=-1.
      u=1.
      do 1000 ik=1,nev
      if (ik > 1) u=ev(ik-1)*umeps
      el=min(ev(ik),u)
   10 elam=0.5*(u+el)
      if (abs(u-el) <= epsi1) goto 35
      iag=0
      q=a(1)-elam
      if (q >= 0.) iag=iag+1
      do 15 i=2,n
      if (q == 0.) x=abs(b(i-1))/epsi
      if (q /= 0.) x=w(i-1)/q
      q=a(i)-elam-x
      if (q >= 0.) iag=iag+1
      if (iag > nev) goto 20
   15 continue
      if (iag >= ik) goto 20
      u=elam
      goto 10
   20 if (iag == ik) goto 30
      m=ik+1
      do 25 i=m,iag
   25 ev(i)=elam
      el=elam
      goto 10
   30 el=elam
      call root(u,el,elam,a,b,w,n,ik)
   35 ev(ik)=elam
      jk=2*ik+ipar-1
      r(1,jk)=1.
      r(2,jk)=-(a(1)-ev(ik))/b(1)
      ddot=1.+r(2,jk)*r(2,jk)
      jm1=2
      do 45 j=3,n
      r(j,jk)=-((a(jm1)-ev(ik))*r(jm1,jk)+b(j-2)*r(j-2,jk))/b(jm1)
      ddot=ddot+r(j,jk)*r(j,jk)
   45 jm1=j
      rnorm=sqrt(nt/(2.*ddot))
      do 50 j=1,n
   50 r(j,jk)=r(j,jk)*rnorm
 1000 continue
      return

  end subroutine tsturm

!------------------------------------------------------------------
  subroutine root(u,el,elam,a,bb,w,n,ik)
!------------------------------------------------------------------

      !implicit real*8(a-h,o-z)
      !parameter (epsi = 1.d-15, epsi1 = 5.d-15)
      !dimension a(*),bb(*),w(*)

      double precision, parameter :: epsi = 1.d-15, epsi1 = 5.d-15
      double precision :: u,el,elam
      double precision, dimension(*) :: a,bb,w
      integer :: n,ik

      double precision :: an,b,bm,bn,del,x
      integer :: i,iag

      !----------------------

    5 elam=0.5*(u+el)
   10 if (abs(u-el) <= 1.5*epsi1) return
      an=a(1)-elam
      b=0.
      bn=-1./an
      iag=0
      if (an >= 0.) iag=iag+1
      do 20 i=2,n
      if (an == 0.) x=abs(bb(i-1))/epsi
      if (an /= 0.) x=w(i-1)/an
      an=a(i)-elam-x
      if (an == 0.) an=epsi
      bm=b
      b=bn
      bn=((a(i)-elam)*b-bm*x-1.)/an
      if (an >= 0.) iag=iag+1
   20 continue
      if (iag == ik) goto 25
      u=elam
      goto 30
   25 el=elam
   30 del=1./bn
      if (abs(del) <= epsi1) del=sign(epsi1,del)
      elam=elam-del
      if (elam >= u .or. elam <= el) goto 5
      goto 10

  end subroutine root

!     ------------------------------------------------------------------
!     subroutine costaper(ipoint, ndata, tas)
!     ------------------------------------------------------------------
      subroutine costaper(ipoint, ndata, tas)
      implicit none

      integer ipoint, ndata
      double precision tas(ndata,*)
      double precision sum, pi
      integer i

      pi = asin(1.0d0)*2
      sum = 0.
      do i =1,ipoint
      tas(i,1) = 1 -  cos( 2*pi*i/ipoint)
      tas(i,1) = tas(i,1) / sqrt(1.5)
      enddo
      return
      end subroutine costaper

!     ------------------------------------------------------------------
!     subroutine boxcar(ipoint, ndata, tas)
!     ------------------------------------------------------------------
      subroutine boxcar(ipoint, ndata, tas)

      integer ipoint, ndata
      double precision tas(ndata,*)
      integer i

      do i =1,ipoint
      tas(i,1) = 1.0
      enddo
      return
      end subroutine boxcar

!     ----------------------------
      subroutine rmean(dat,npts)
!     ---------------------------

!  remove mean from data

      implicit none

      integer npts
      real dat(npts),sum
      integer i


      sum=0
      do i = 1, npts
         sum = sum + dat(i)
      enddo
      do i = 1, npts
         dat(i) = dat(i) - sum/npts
      enddo

      return

      end subroutine rmean

!     ----------------------------
      subroutine rtrend(y, n)
!     ---------------------------

!      remove trend a*i + b
!      revised from lupei's c code

      implicit none

      integer n
      real y(n)
      integer i
      real*8  a, b, a11, a12, a22, y1, y2

      y1 = 0.
      y2 = 0.
      do i = 1, n
        y1 = y1 + (i-1)*y(i)
        y2 = y2 + y(i);
      enddo
      a12 = 0.5*n*(n-1)
      a11 = a12*(2*n-1)/3.
      a22 = n;
      b = a11*a22-a12*a12
      a = (a22*y1-a12*y2)/b
      b = (a11*y2-a12*y1)/b
      do i = 1, n
       y(i) = y(i) - a * (i-1) - b
      enddo

      return
      end subroutine rtrend

!     ---------------------------------------
      subroutine taper (a, n, start, end, b)
!     ---------------------------------------
        integer n
        real start, end
        real a(*), b(*)

! taper a(1:n) to b(1:n)
! cosine taper the beginning portion 'start' and end portion 'end' between [0,1]

        real :: an,ang,xi,cs
        integer :: i, m1,m2,m3,m4,m5

      an = n
      m1 = an*start+0.5
      m2 = m1 + 1
      if (m1 > 0) then
        ang = 3.1415926 / float(m1)
        do i=1, m1
          xi = i
          cs = (1.0-cos(xi*ang))/2.0
          b(i) = a(i)*cs
        enddo
      endif

      m3 = an*end+0.5
      m5 = n-m3
      m4 = m5 + 1

      if (m3 > 0) then
        ang = 3.1415926 / float (m3)
        do i=m4,n
          xi = i-n-1
          cs = (1.0-cos(xi*ang))/2.0
          b(i) = a(i)*cs
        enddo
      endif
      do i=m2, m5
        b(i) = a(i)
      enddo

      return
    end subroutine taper

    end module mtadj_sub3
