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

  subroutine intgrl(sum,r,nir,ner,f,s1,s2,s3)

! Computes the integral of f[i]*r[i]*r[i] from i=nir to i=ner for
! radii values as in model PREM_an640

  implicit none

! Argument variables
  integer ner,nir
  double precision f(640),r(640),s1(640),s2(640)
  double precision s3(640),sum

! Local variables
  integer i,j,n,kdis(28)
  integer ndis,nir1
  double precision rji,yprime(640)

  double precision, parameter :: third = 1.0d0/3.0d0
  double precision, parameter :: fifth = 1.0d0/5.0d0
  double precision, parameter :: sixth = 1.0d0/6.0d0

  data kdis/163,323,336,517,530,540,565,590,609,619,626,633,16*0/

  ndis = 12
  n = 640

  call deriv(f,yprime,n,r,ndis,kdis,s1,s2,s3)
  nir1 = nir + 1
  sum = 0.0d0
  do i=nir1,ner
    j = i-1
    rji = r(i) - r(j)
    sum=sum+r(j)*r(j)*rji*(f(j)+rji*(.50d0*s1(j)+rji*(third*s2(j)+rji* &
      .250d0*s3(j))))+2.0d0*r(j)*rji*rji*(.50d0*f(j)+rji*(third*s1(j)+rji* &
     (.250d0*s2(j)+rji*fifth*s3(j))))+rji*rji*rji*(third*f(j)+rji* &
     (.250d0*s1(j)+rji*(fifth*s2(j)+rji*sixth*s3(j))))
  enddo

  end subroutine intgrl

! -------------------------------

  subroutine deriv(y,yprime,n,r,ndis,kdis,s1,s2,s3)

  implicit none

! Argument variables
  integer kdis(28),n,ndis
  double precision r(n),s1(n),s2(n),s3(n)
  double precision y(n),yprime(n)

! Local variables
  integer i,j,j1,j2
  integer k,nd,ndp
  double precision a0,b0,b1
  double precision f(3,1000),h,h2,h2a
  double precision h2b,h3a,ha,s13
  double precision s21,s32,yy(3)

  yy(1) = 0.d0
  yy(2) = 0.d0
  yy(3) = 0.d0

  ndp=ndis+1
  do 3 nd=1,ndp
  if(nd == 1) goto 4
  if(nd == ndp) goto 5
  j1=kdis(nd-1)+1
  j2=kdis(nd)-2
  goto 6
    4 j1=1
  j2=kdis(1)-2
  goto 6
    5 j1=kdis(ndis)+1
  j2=n-2
    6 if((j2+1-j1)>0) goto 11
  j2=j2+2
  yy(1)=(y(j2)-y(j1))/(r(j2)-r(j1))
  s1(j1)=yy(1)
  s1(j2)=yy(1)
  s2(j1)=yy(2)
  s2(j2)=yy(2)
  s3(j1)=yy(3)
  s3(j2)=yy(3)
  goto 3
   11 a0=0.0d0
  if(j1 == 1) goto 7
  h=r(j1+1)-r(j1)
  h2=r(j1+2)-r(j1)
  yy(1)=h*h2*(h2-h)
  h=h*h
  h2=h2*h2
  b0=(y(j1)*(h-h2)+y(j1+1)*h2-y(j1+2)*h)/yy(1)
  goto 8
 7 b0=0.0d0
 8 b1=b0

  if(j2 > 1000) stop 'error in subroutine deriv for j2'

  do i=j1,j2
    h=r(i+1)-r(i)
    yy(1)=y(i+1)-y(i)
    h2=h*h
    ha=h-a0
    h2a=h-2.0d0*a0
    h3a=2.0d0*h-3.0d0*a0
    h2b=h2*b0
    s1(i)=h2/ha
    s2(i)=-ha/(h2a*h2)
    s3(i)=-h*h2a/h3a
    f(1,i)=(yy(1)-h*b0)/(h*ha)
    f(2,i)=(h2b-yy(1)*(2.0d0*h-a0))/(h*h2*h2a)
    f(3,i)=-(h2b-3.0d0*yy(1)*ha)/(h*h3a)
    a0=s3(i)
    b0=f(3,i)
  enddo

  i=j2+1
  h=r(i+1)-r(i)
  yy(1)=y(i+1)-y(i)
  h2=h*h
  ha=h-a0
  h2a=h*ha
  h2b=h2*b0-yy(1)*(2.d0*h-a0)
  s1(i)=h2/ha
  f(1,i)=(yy(1)-h*b0)/h2a
  ha=r(j2)-r(i+1)
  yy(1)=-h*ha*(ha+h)
  ha=ha*ha
  yy(1)=(y(i+1)*(h2-ha)+y(i)*ha-y(j2)*h2)/yy(1)
  s3(i)=(yy(1)*h2a+h2b)/(h*h2*(h-2.0d0*a0))
  s13=s1(i)*s3(i)
  s2(i)=f(1,i)-s13

  do j=j1,j2
    k=i-1
    s32=s3(k)*s2(i)
    s1(i)=f(3,k)-s32
    s21=s2(k)*s1(i)
    s3(k)=f(2,k)-s21
    s13=s1(k)*s3(k)
    s2(k)=f(1,k)-s13
    i=k
  enddo

  s1(i)=b1
  j2=j2+2
  s1(j2)=yy(1)
  s2(j2)=yy(2)
  s3(j2)=yy(3)
 3 continue

  do i=1,n
    yprime(i)=s1(i)
  enddo

  end subroutine deriv

