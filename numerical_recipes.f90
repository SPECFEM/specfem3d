!=====================================================================
!
!  Routines from "Numerical Recipes: the Art of Scientific Computing"
!  W. H. Press et al., Cambridge University Press
!
!=====================================================================

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
    u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
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

! evaluate spline (adapted from Numerical Recipes)

  subroutine splint(xa,ya,y2a,n,x,y)

  implicit none

  integer n
  double precision, dimension(n) :: XA,YA,Y2A
  double precision x,y

  integer k,klo,khi
  double precision h,a,b

  KLO = 1
  KHI = N

  do while (KHI-KLO > 1)
    K=(KHI+KLO)/2
    if(XA(K) > X) then
      KHI=K
    else
      KLO=K
    endif
  enddo

  H = XA(KHI) - XA(KLO)
  IF (H == 0.d0) stop 'bad input in spline evaluation'

  A = (XA(KHI)-X) / H
  B = (X-XA(KLO)) / H

  Y = A*YA(KLO) + B*YA(KHI) + ((A**3-A)*Y2A(KLO) + (B**3-B)*Y2A(KHI))*(H**2)/6.d0

  end subroutine splint

