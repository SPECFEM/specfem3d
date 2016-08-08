!=====================================================================
!
!  Routines from "Numerical Recipes: the Art of Scientific Computing"
!  W. H. Press et al., Cambridge University Press
!
!=====================================================================

! double precision routines

  double precision function erf(x)
  double precision x

! this routine uses routine gammp
  double precision gammp

  if (x < 0.) then
    erf=-gammp(0.5d0,x**2)
  else
    erf=gammp(0.5d0,x**2)
  endif

  end function erf

! ---------------------------------

  double precision function gammp(a,x)
  double precision a,x

! this routine uses routines gcf and gser
  double precision gammcf,gamser,gln

  if (x < 0.d0 .or. a <= 0.d0) stop 'bad arguments in gammp'

  if (x < a+1.d0) then
    call gser(gamser,a,x,gln)
    gammp=gamser
  else
    call gcf(gammcf,a,x,gln)
    gammp=1.d0-gammcf
  endif

  end function gammp

! ---------------------------------

  subroutine gcf(gammcf,a,x,gln)

  double precision a,gammcf,gln,x

  double precision, parameter :: EPS=3.d-7,FPMIN=1.d-30
  integer, parameter :: ITMAX=100

! this routine uses routine gammln

  integer i
  double precision an,b,c,d,del,h

  double precision, external :: gammln

  gln=gammln(a)
  b=x+1.d0-a
  c=1.d0/FPMIN
  d=1.d0/b
  h=d
  do i=1,ITMAX
    an=-i*(i-a)
    b=b+2.d0
    d=an*d+b
    if (dabs(d) < FPMIN)d=FPMIN
    c=b+an/c
    if (dabs(c) < FPMIN)c=FPMIN
    d=1.d0/d
    del=d*c
    h=h*del
    if (dabs(del-1.d0) < EPS) then
      gammcf=exp(-x+a*log(x)-gln)*h
      return
    endif
  enddo

  stop 'a too large, ITMAX too small in gcf'

  end subroutine gcf

! ---------------------------------

  subroutine gser(gamser,a,x,gln)

  double precision a,gamser,gln,x

  integer, parameter :: ITMAX=100
  double precision, parameter :: EPS=3.d-7

! this routine uses routine gammln

  integer n
  double precision ap,del,sumval

  double precision, external :: gammln

  gln=gammln(a)

  if (x <= 0.d0) then
    if (x < 0.d0) stop 'x < 0 in gser'
    gamser=0.d0
    return
  endif

  ap=a
  sumval=1.d0/a
  del=sumval

  do n=1,ITMAX
    ap=ap+1.d0
    del=del*x/ap
    sumval=sumval+del
    if (dabs(del) < dabs(sumval)*EPS) then
      gamser=sumval*exp(-x+a*log(x)-gln)
      return
    endif
  enddo

  stop 'a too large, ITMAX too small in gser'

  end subroutine gser

! ---------------------------------

  double precision function gammln(xx)

  double precision xx

  integer j
  double precision ser,stp,tmp,x,y,cof(6)

  cof(1) = 76.18009172947146d0
  cof(2) = -86.50532032941677d0
  cof(3) = 24.01409824083091d0
  cof(4) = -1.231739572450155d0
  cof(5) = 0.1208650973866179d-2
  cof(6) = -0.5395239384953d-5

  stp = 2.5066282746310005d0

  x=xx
  y=x
  tmp=x+5.5d0
  tmp=(x+0.5d0)*log(tmp)-tmp
  ser=1.000000000190015d0
  do j=1,6
    y=y+1.d0
    ser=ser+cof(j)/y
  enddo
  gammln=tmp+log(stp*ser/x)

  end function gammln

! ---------------------------------

! compute spline coefficients (Numerical Recipes)
! modified to use dynamic allocation

  subroutine spline(x,y,n,yp1,ypn,y2)

  implicit none

  integer n
  double precision x(n),y(n),y2(n)
  double precision yp1,ypn
  double precision, dimension(:), allocatable :: u

  integer i,k
  double precision sig,p,qn,un

  allocate(u(n))

  y2(1)=-0.5d0
  u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)

  do i=2,n-1
    sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
    p=sig*y2(i-1)+2.d0
    y2(i)=(sig-1.d0)/p
    u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
                /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
  enddo

  qn=0.5d0
  un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
  y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)

  do k=n-1,1,-1
    y2(k)=y2(k)*y2(k+1)+u(k)
  enddo

  deallocate(u)

  end subroutine spline

! --------------

! evaluate spline (Numerical Recipes)

  subroutine splint(xa,ya,y2a,n,x,y)

  implicit none

  integer n
  double precision XA(N),YA(N),Y2A(N)
  double precision x,y

  integer k,klo,khi
  double precision h,a,b

  KLO=1
  KHI=N
 1 if (KHI-KLO > 1) then
    K=(KHI+KLO)/2
    if (XA(K) > X) then
      KHI=K
    ELSE
      KLO=K
    endif
  goto 1
  endif
  H=XA(KHI)-XA(KLO)
  if (H == 0.d0) stop 'Bad input in spline evaluation'
  A=(XA(KHI)-X)/H
  B=(X-XA(KLO))/H

  Y=A*YA(KLO)+B*YA(KHI)+((A**3-A)*Y2A(KLO) + (B**3-B)*Y2A(KHI))*(H**2)/6.d0

  end subroutine splint

