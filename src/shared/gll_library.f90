!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

!=======================================================================
!
!  Library to compute the Gauss-Lobatto-Legendre points and weights
!  Based on Gauss-Lobatto routines from M.I.T.
!  Department of Mechanical Engineering
!
!=======================================================================

  double precision function endw1(n,alpha,beta)

  implicit none

  integer n
  double precision alpha,beta

  double precision, parameter :: zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0
  double precision apb,f1,fint1,fint2,f2,di,abn,abnn,a1,a2,a3,f3
  double precision, external :: gammaf
  integer i

  f3 = zero
  apb   = alpha+beta
  if (n == 0) then
   endw1 = zero
   return
  endif
  f1   = gammaf(alpha+two)*gammaf(beta+one)/gammaf(apb+three)
  f1   = f1*(apb+two)*two**(apb+two)/two
  if (n == 1) then
   endw1 = f1
   return
  endif
  fint1 = gammaf(alpha+two)*gammaf(beta+one)/gammaf(apb+three)
  fint1 = fint1*two**(apb+two)
  fint2 = gammaf(alpha+two)*gammaf(beta+two)/gammaf(apb+four)
  fint2 = fint2*two**(apb+three)
  f2    = (-two*(beta+two)*fint1 + (apb+four)*fint2) * (apb+three)/four
  if (n == 2) then
   endw1 = f2
   return
  endif
  do i=3,n
   di   = dble(i-1)
   abn  = alpha+beta+di
   abnn = abn+di
   a1   = -(two*(di+alpha)*(di+beta))/(abn*abnn*(abnn+one))
   a2   =  (two*(alpha-beta))/(abnn*(abnn+two))
   a3   =  (two*(abn+one))/((abnn+two)*(abnn+one))
   f3   =  -(a2*f2+a1*f1)/a3
   f1   = f2
   f2   = f3
  enddo
  endw1  = f3

  end function endw1

!
!=======================================================================
!

  double precision function endw2(n,alpha,beta)

  implicit none

  integer n
  double precision alpha,beta

  double precision, parameter :: zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0
  double precision apb,f1,fint1,fint2,f2,di,abn,abnn,a1,a2,a3,f3
  double precision, external :: gammaf
  integer i

  apb   = alpha+beta
  f3 = zero
  if (n == 0) then
   endw2 = zero
   return
  endif
  f1   = gammaf(alpha+one)*gammaf(beta+two)/gammaf(apb+three)
  f1   = f1*(apb+two)*two**(apb+two)/two
  if (n == 1) then
   endw2 = f1
   return
  endif
  fint1 = gammaf(alpha+one)*gammaf(beta+two)/gammaf(apb+three)
  fint1 = fint1*two**(apb+two)
  fint2 = gammaf(alpha+two)*gammaf(beta+two)/gammaf(apb+four)
  fint2 = fint2*two**(apb+three)
  f2    = (two*(alpha+two)*fint1 - (apb+four)*fint2) * (apb+three)/four
  if (n == 2) then
   endw2 = f2
   return
  endif
  do i=3,n
   di   = dble(i-1)
   abn  = alpha+beta+di
   abnn = abn+di
   a1   =  -(two*(di+alpha)*(di+beta))/(abn*abnn*(abnn+one))
   a2   =  (two*(alpha-beta))/(abnn*(abnn+two))
   a3   =  (two*(abn+one))/((abnn+two)*(abnn+one))
   f3   =  -(a2*f2+a1*f1)/a3
   f1   = f2
   f2   = f3
  enddo
  endw2  = f3

  end function endw2

!
!=======================================================================
!

  double precision function gammaf (x)

  implicit none

  double precision, parameter :: pi = 3.141592653589793d0

  double precision x

  double precision, parameter :: half=0.5d0,one=1.d0,two=2.d0

  gammaf = one

  if (x == -half) gammaf = -two*dsqrt(pi)
  if (x ==  half) gammaf =  dsqrt(pi)
  if (x ==  one ) gammaf =  one
  if (x ==  two ) gammaf =  one
  if (x ==  1.5d0) gammaf =  dsqrt(pi)/2.d0
  if (x ==  2.5d0) gammaf =  1.5d0*dsqrt(pi)/2.d0
  if (x ==  3.5d0) gammaf =  2.5d0*1.5d0*dsqrt(pi)/2.d0
  if (x ==  3.d0 ) gammaf =  2.d0
  if (x ==  4.d0 ) gammaf = 6.d0
  if (x ==  5.d0 ) gammaf = 24.d0
  if (x ==  6.d0 ) gammaf = 120.d0

  end function gammaf

!
!=====================================================================
!

  subroutine jacg (xjac,np,alpha,beta)

!=======================================================================
!
! computes np Gauss points, which are the zeros of the
! Jacobi polynomial with parameters alpha and beta
!
!                  .alpha = beta =  0.0  ->  Legendre points
!                  .alpha = beta = -0.5  ->  Chebyshev points
!
!=======================================================================

  implicit none

  integer np
  double precision alpha,beta
  double precision xjac(np)

  ! local parameters
  integer k,j,i,jmin,jm,n
  double precision xlast,dth,x,x1,x2,recsum,delx,xmin,swap
  double precision p,pd,pm1,pdm1,pm2,pdm2

  integer, parameter :: K_MAX_ITER = 10
  double precision, parameter :: zero = 0.d0, eps = 1.0d-12

  pm1 = zero
  pm2 = zero
  pdm1 = zero
  pdm2 = zero

  xlast = 0.d0
  n   = np-1
  dth = 4.d0*datan(1.d0)/(2.d0*dble(n)+2.d0)
  p = 0.d0
  pd = 0.d0

  do j=1,np
   if(j == 1) then
      x = dcos((2.d0*(dble(j)-1.d0)+1.d0)*dth)
   else
      x1 = dcos((2.d0*(dble(j)-1.d0)+1.d0)*dth)
      x2 = xlast
      x  = (x1+x2)/2.d0
   endif

   do k=1,K_MAX_ITER
      call jacobf (p,pd,pm1,pdm1,pm2,pdm2,np,alpha,beta,x)
      recsum = 0.d0
      jm = j-1
      do i=1,jm
         recsum = recsum+1.d0/(x-xjac(np-i+1))
      enddo
      delx = -p/(pd-recsum*p)
      x    = x+delx

      ! exits loop if increment too small
      if(abs(delx) < eps) exit

   enddo

   ! checks bounds
   if( np-j+1 < 1 .or. np-j+1 > np ) stop 'error np-j+1-index in jacg'

   xjac(np-j+1) = x
   xlast        = x
  enddo

  jmin = 0

  ! orders xjac array in increasing values
  do i=1,np
   xmin = 2.d0
   jmin = i

   ! looks for index with minimum value
   do j=i,np
      ! note: some compilers (cray) might be too aggressive in optimizing this loop,
      !       thus we need this temporary array value x to store and compare values
      x = xjac(j)

      if( x < xmin) then
         xmin = x
         jmin = j
      endif
   enddo

   ! checks bounds
   if(jmin < 1 .or. jmin > np ) stop 'error j-index in jacg'

   if(jmin /= i) then
      swap = xjac(i)
      xjac(i) = xjac(jmin)
      xjac(jmin) = swap
   endif

  enddo

  end subroutine jacg

!
!=====================================================================
!

  subroutine jacobf (poly,pder,polym1,pderm1,polym2,pderm2,n,alp,bet,x)

!=======================================================================
!
! Computes the Jacobi polynomial of degree n and its derivative at x
!
!=======================================================================

  implicit none

  double precision poly,pder,polym1,pderm1,polym2,pderm2,alp,bet,x
  integer n

  double precision apb,polyl,pderl,dk,a1,a2,b3,a3,a4,polyn,pdern,psave,pdsave
  integer k

  apb  = alp+bet
  poly = 1.d0
  pder = 0.d0
  psave = 0.d0
  pdsave = 0.d0

  if (n == 0) return

  polyl = poly
  pderl = pder
  poly  = (alp-bet+(apb+2.d0)*x)/2.d0
  pder  = (apb+2.d0)/2.d0
  if (n == 1) return

  do k=2,n
    dk = dble(k)
    a1 = 2.d0*dk*(dk+apb)*(2.d0*dk+apb-2.d0)
    a2 = (2.d0*dk+apb-1.d0)*(alp**2-bet**2)
    b3 = (2.d0*dk+apb-2.d0)
    a3 = b3*(b3+1.d0)*(b3+2.d0)
    a4 = 2.d0*(dk+alp-1.d0)*(dk+bet-1.d0)*(2.d0*dk+apb)
    polyn  = ((a2+a3*x)*poly-a4*polyl)/a1
    pdern  = ((a2+a3*x)*pder-a4*pderl+a3*poly)/a1
    psave  = polyl
    pdsave = pderl
    polyl  = poly
    poly   = polyn
    pderl  = pder
    pder   = pdern
  enddo

  polym1 = polyl
  pderm1 = pderl
  polym2 = psave
  pderm2 = pdsave

  end subroutine jacobf

!
!------------------------------------------------------------------------
!

  double precision FUNCTION PNDLEG (Z,N)

!------------------------------------------------------------------------
!
!     Compute the derivative of the Nth order Legendre polynomial at Z.
!     Based on the recursion formula for the Legendre polynomials.
!
!------------------------------------------------------------------------
  implicit none

  double precision z
  integer n

  double precision P1,P2,P1D,P2D,P3D,FK,P3
  integer k

  P1   = 1.d0
  P2   = Z
  P1D  = 0.d0
  P2D  = 1.d0
  P3D  = 1.d0

  do K = 1, N-1
    FK  = dble(K)
    P3  = ((2.d0*FK+1.d0)*Z*P2 - FK*P1)/(FK+1.d0)
    P3D = ((2.d0*FK+1.d0)*P2 + (2.d0*FK+1.d0)*Z*P2D - FK*P1D) / (FK+1.d0)
    P1  = P2
    P2  = P3
    P1D = P2D
    P2D = P3D
  enddo

  PNDLEG = P3D

  end function pndleg

!
!------------------------------------------------------------------------
!

  double precision FUNCTION PNLEG (Z,N)

!------------------------------------------------------------------------
!
!     Compute the value of the Nth order Legendre polynomial at Z.
!     Based on the recursion formula for the Legendre polynomials.
!
!------------------------------------------------------------------------
  implicit none

  double precision z
  integer n

  double precision P1,P2,P3,FK
  integer k

  P1   = 1.d0
  P2   = Z
  P3   = P2

  do K = 1, N-1
    FK  = dble(K)
    P3  = ((2.d0*FK+1.d0)*Z*P2 - FK*P1)/(FK+1.d0)
    P1  = P2
    P2  = P3
  enddo

  PNLEG = P3

  end function pnleg

!
!------------------------------------------------------------------------
!

  double precision function pnormj (n,alpha,beta)

  implicit none

  double precision alpha,beta
  integer n

  double precision one,two,dn,const,prod,dindx,frac
  double precision, external :: gammaf
  integer i

  one   = 1.d0
  two   = 2.d0
  dn    = dble(n)
  const = alpha+beta+one

  if (n <= 1) then
    prod   = gammaf(dn+alpha)*gammaf(dn+beta)
    prod   = prod/(gammaf(dn)*gammaf(dn+alpha+beta))
    pnormj = prod * two**const/(two*dn+const)
    return
  endif

  prod  = gammaf(alpha+one)*gammaf(beta+one)
  prod  = prod/(two*(one+const)*gammaf(const+one))
  prod  = prod*(one+alpha)*(two+alpha)
  prod  = prod*(one+beta)*(two+beta)

  do i=3,n
    dindx = dble(i)
    frac  = (dindx+alpha)*(dindx+beta)/(dindx*(dindx+alpha+beta))
    prod  = prod*frac
  enddo

  pnormj = prod * two**const/(two*dn+const)

  end function pnormj

!
!------------------------------------------------------------------------
!

  subroutine zwgjd(z,w,np,alpha,beta)

!=======================================================================
!
!     Z w g j d : Generate np Gauss-Jacobi points and weights
!                 associated with Jacobi polynomial of degree n = np-1
!
!     Note : Coefficients alpha and beta must be greater than -1.
!     ----
!=======================================================================

  implicit none

  double precision, parameter :: zero=0.d0,one=1.d0,two=2.d0

  integer np
  double precision z(np),w(np)
  double precision alpha,beta

  ! local paraeters
  integer n,np1,np2,i
  double precision p,pd,pm1,pdm1,pm2,pdm2
  double precision apb,dnp1,dnp2,fac1,fac2,fac3,fnorm,rcoef
  double precision, external :: gammaf,pnormj

  pd = zero
  pm1 = zero
  pm2 = zero
  pdm1 = zero
  pdm2 = zero

  n    = np-1
  apb  = alpha+beta
  p    = zero
  pdm1 = zero

  if (np <= 0) stop 'minimum number of Gauss points is 1'

  if ((alpha <= -one) .or. (beta <= -one)) stop 'alpha and beta must be greater than -1'

  if (np == 1) then
   z(1) = (beta-alpha)/(apb+two)
   w(1) = gammaf(alpha+one)*gammaf(beta+one)/gammaf(apb+two) * two**(apb+one)
   return
  endif

  call jacg(z,np,alpha,beta)

  np1   = n+1
  np2   = n+2
  dnp1  = dble(np1)
  dnp2  = dble(np2)
  fac1  = dnp1+alpha+beta+one
  fac2  = fac1+dnp1
  fac3  = fac2+one
  fnorm = pnormj(np1,alpha,beta)
  rcoef = (fnorm*fac2*fac3)/(two*fac1*dnp2)
  do i=1,np
    call jacobf(p,pd,pm1,pdm1,pm2,pdm2,np2,alpha,beta,z(i))
    w(i) = -rcoef/(p*pdm1)
  enddo

  end subroutine zwgjd

!
!------------------------------------------------------------------------
!

  subroutine zwgljd(z,w,np,alpha,beta)

!=======================================================================
!
!     Z w g l j d : Generate np Gauss-Lobatto-Jacobi points and the
!     -----------   weights associated with Jacobi polynomials of degree
!                   n = np-1.
!
!     Note : alpha and beta coefficients must be greater than -1.
!            Legendre polynomials are special case of Jacobi polynomials
!            just by setting alpha and beta to 0.
!
!=======================================================================

  implicit none

  double precision, parameter :: zero=0.d0,one=1.d0,two=2.d0

  integer np
  double precision alpha,beta
  double precision z(np), w(np)

  ! local parameters
  integer n,nm1,i
  double precision p,pd,pm1,pdm1,pm2,pdm2
  double precision alpg,betg
  double precision, external :: endw1,endw2

  p = zero
  pm1 = zero
  pm2 = zero
  pdm1 = zero
  pdm2 = zero

  n   = np-1
  nm1 = n-1
  pd  = zero

  if (np <= 1) stop 'minimum number of Gauss-Lobatto points is 2'

! with spectral elements, use at least 3 points
  if (np <= 2) stop 'minimum number of Gauss-Lobatto points for the SEM is 3'

  if ((alpha <= -one) .or. (beta <= -one)) stop 'alpha and beta must be greater than -1'

  if (nm1 > 0) then
    alpg  = alpha+one
    betg  = beta+one
    call zwgjd(z(2:n),w(2:n),nm1,alpg,betg)
  endif

  z(1)  = - one
  z(np) =  one

  do i=2,np-1
   w(i) = w(i)/(one-z(i)**2)
  enddo

  call jacobf(p,pd,pm1,pdm1,pm2,pdm2,n,alpha,beta,z(1))
  w(1)  = endw1(n,alpha,beta)/(two*pd)

  call jacobf(p,pd,pm1,pdm1,pm2,pdm2,n,alpha,beta,z(np))
  w(np) = endw2(n,alpha,beta)/(two*pd)

  end subroutine zwgljd

