!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2008
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

  subroutine lgndr(l,c,s,x,dx)

! computes Legendre function x(l,m,theta)
! theta=colatitude,c=cos(theta),s=sin(theta),l=angular order,
! sin(theta) restricted so that sin(theta) > 1.e-7
! x(1) contains m=0, x(2) contains m=1, x(k+1) contains m=k
! m=azimuthal (longitudinal) order 0 <= m <= l
! dx=dx/dtheta
!
! subroutine originally came from Physics Dept. Princeton through
! Peter Davis, modified by Jeffrey Park

  implicit none

! argument variables
  integer l
  double precision x(2*l+1),dx(2*l+1)
  double precision c,s

! local variables
  integer i,lp1,lpsafe,lsave
  integer m,maxsin,mmm,mp1

  double precision sqroot2over2,c1,c2,cot
  double precision ct,d,f1,f2
  double precision f3,fac,g1,g2
  double precision g3,rfpi,sqroot3,sos
  double precision ss,stom,t,tol
  double precision v,y

  tol = 1.d-05
  rfpi = 0.282094791773880d0
  sqroot3 = 1.73205080756890d0
  sqroot2over2 = 0.707106781186550d0

  if(s >= 1.0d0-tol) s=1.0d0-tol
  lsave=l
  if(l<0) l=-1-l
  if(l>0) goto 1
  x(1)=rfpi
  dx(1)=0.0d0
  l=lsave
  return
 1 if(l /= 1) goto 2
  c1=sqroot3*rfpi
  c2=sqroot2over2*c1
  x(1)=c1*c
  x(2)=-c2*s
  dx(1)=-c1*s
  dx(2)=-c2*c
  l=lsave
  return
    2 sos=s
  if(s<tol) s=tol
  cot=c/s
  ct=2.0d0*c
  ss=s*s
  lp1=l+1
  g3=0.0d0
  g2=1.0d0
  f3=0.0d0

! evaluate m=l value, sans (sin(theta))**l
  do i=1,l
    g2=g2*(1.0d0-1.0d0/(2.0d0*i))
  enddo
  g2=rfpi*dsqrt((2*l+1)*g2)
  f2=l*cot*g2
  x(lp1)=g2
  dx(lp1)=f2
  v=1.0d0
  y=2.0d0*l
  d=dsqrt(v*y)
  t=0.0d0
  mp1=l
  m=l-1

! these recursions are similar to ordinary m-recursions, but since we
! have taken the s**m factor out of the xlm's, the recursion has the powers
! of sin(theta) instead
    3 g1=-(ct*mp1*g2+ss*t*g3)/d
  f1=(mp1*(2.0d0*s*g2-ct*f2)-t*ss*(f3+cot*g3))/d-cot*g1
  x(mp1)=g1
  dx(mp1)=f1
  if(m == 0) goto 4
  mp1=m
  m=m-1
  v=v+1.0d0
  y=y-1.0d0
  t=d
  d=dsqrt(v*y)
  g3=g2
  g2=g1
  f3=f2
  f2=f1
 goto 3
! explicit conversion to integer added
    4 maxsin=int(-72.0d0/log10(s))

! maxsin is the max exponent of sin(theta) without underflow
  lpsafe=min0(lp1,maxsin)
  stom=1.0d0
  fac=sign(1.0d0,dble((l/2)*2-l) + 0.50d0)

! multiply xlm by sin**m
  do m=1,lpsafe
    x(m)=fac*x(m)*stom
    dx(m)=fac*dx(m)*stom
    stom=stom*s
  enddo

! set any remaining xlm to zero
  if(maxsin <= l) then
    mmm=maxsin+1
    do m=mmm,lp1
      x(m)=0.0d0
      dx(m)=0.0d0
    enddo
  endif

  s=sos
  l=lsave

  end subroutine lgndr

