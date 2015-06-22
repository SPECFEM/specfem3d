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

  subroutine read_mantle_model(D3MM_V)

  implicit none

  include "constants.h"

! three_d_mantle_model_variables
  type three_d_mantle_model_variables
    sequence
    double precision dvs_a(0:NK,0:NS,0:NS)
    double precision dvs_b(0:NK,0:NS,0:NS)
    double precision dvp_a(0:NK,0:NS,0:NS)
    double precision dvp_b(0:NK,0:NS,0:NS)
    double precision spknt(NK+1)
    double precision qq0(NK+1,NK+1)
    double precision qq(3,NK+1,NK+1)
  end type three_d_mantle_model_variables

  type (three_d_mantle_model_variables) D3MM_V
! three_d_mantle_model_variables

  integer k,l,m

  character(len=150) S20RTS, P12

  call get_value_string(S20RTS, 'model.S20RTS', 'DATA/s20rts/S20RTS.dat')
  call get_value_string(P12, 'model.P12', 'DATA/s20rts/P12.dat')

! S20RTS degree 20 S model from Ritsema
  open(unit=10,file=S20RTS,status='old',action='read')
  do k=0,NK
    do l=0,NS
      read(10,*) D3MM_V%dvs_a(k,l,0),(D3MM_V%dvs_a(k,l,m),D3MM_V%dvs_b(k,l,m),m=1,l)
    enddo
  enddo
  close(10)

! P12 degree 12 P model from Ritsema
  open(unit=10,file=P12,status='old',action='read')
  do k=0,NK
    do l=0,12
      read(10,*) D3MM_V%dvp_a(k,l,0),(D3MM_V%dvp_a(k,l,m),D3MM_V%dvp_b(k,l,m),m=1,l)
    enddo
    do l=13,NS
      D3MM_V%dvp_a(k,l,0) = 0.0d0
      do m=1,l
        D3MM_V%dvp_a(k,l,m) = 0.0d0
        D3MM_V%dvp_b(k,l,m) = 0.0d0
      enddo
    enddo
  enddo
  close(10)

! set up the splines used as radial basis functions by Ritsema
  call splhsetup(D3MM_V)

  end subroutine read_mantle_model

!---------------------------

  subroutine mantle_model(radius,theta,phi,dvs,dvp,drho,D3MM_V)

  implicit none

  include "constants.h"

! three_d_mantle_model_variables
  type three_d_mantle_model_variables
    sequence
    double precision dvs_a(0:NK,0:NS,0:NS)
    double precision dvs_b(0:NK,0:NS,0:NS)
    double precision dvp_a(0:NK,0:NS,0:NS)
    double precision dvp_b(0:NK,0:NS,0:NS)
    double precision spknt(NK+1)
    double precision qq0(NK+1,NK+1)
    double precision qq(3,NK+1,NK+1)
  end type three_d_mantle_model_variables

  type (three_d_mantle_model_variables) D3MM_V
! three_d_mantle_model_variables

! factor to convert perturbations in shear speed to perturbations in density
  double precision, parameter :: SCALE_RHO = 0.40d0

  double precision radius,theta,phi,dvs,dvp,drho

  double precision, parameter :: RMOHO_ = 6346600.d0
  double precision, parameter :: RCMB_ = 3480000.d0
  double precision, parameter :: R_EARTH_ = 6371000.d0
  double precision, parameter :: ZERO_ = 0.d0

  integer l,m,k
  double precision r_moho,r_cmb,xr
  double precision dvs_alm,dvs_blm
  double precision dvp_alm,dvp_blm
  double precision rsple,radial_basis(0:NK)
  double precision sint,cost,x(2*NS+1),dx(2*NS+1)

  dvs = ZERO_
  dvp = ZERO_
  drho = ZERO_

  r_moho = RMOHO_ / R_EARTH_
  r_cmb = RCMB_ / R_EARTH_
  if(radius>=r_moho .or. radius <= r_cmb) return

  xr=-1.0d0+2.0d0*(radius-r_cmb)/(r_moho-r_cmb)
  do k=0,NK
    radial_basis(k)=rsple(1,NK+1,D3MM_V%spknt(1),D3MM_V%qq0(1,NK+1-k),D3MM_V%qq(1,1,NK+1-k),xr)
  enddo

  do l=0,NS
    sint=dsin(theta)
    cost=dcos(theta)
    call lgndr(l,cost,sint,x,dx)
    dvs_alm=0.0d0
    dvp_alm=0.0d0
    do k=0,NK
      dvs_alm=dvs_alm+radial_basis(k)*D3MM_V%dvs_a(k,l,0)
      dvp_alm=dvp_alm+radial_basis(k)*D3MM_V%dvp_a(k,l,0)
    enddo
    dvs=dvs+dvs_alm*x(1)
    dvp=dvp+dvp_alm*x(1)
    do m=1,l
      dvs_alm=0.0d0
      dvp_alm=0.0d0
      dvs_blm=0.0d0
      dvp_blm=0.0d0
      do k=0,NK
        dvs_alm=dvs_alm+radial_basis(k)*D3MM_V%dvs_a(k,l,m)
        dvp_alm=dvp_alm+radial_basis(k)*D3MM_V%dvp_a(k,l,m)
        dvs_blm=dvs_blm+radial_basis(k)*D3MM_V%dvs_b(k,l,m)
        dvp_blm=dvp_blm+radial_basis(k)*D3MM_V%dvp_b(k,l,m)
      enddo
      dvs=dvs+(dvs_alm*dcos(dble(m)*phi)+dvs_blm*dsin(dble(m)*phi))*x(m+1)
      dvp=dvp+(dvp_alm*dcos(dble(m)*phi)+dvp_blm*dsin(dble(m)*phi))*x(m+1)
    enddo
  enddo

  drho = SCALE_RHO*dvs

  end subroutine mantle_model

!----------------------------------

  subroutine splhsetup(D3MM_V)!!!!!!!!!!!!!!(spknt,qq0,qq)

  implicit none
  include "constants.h"

!!!!!!!!!!!!!!!!!!!  double precision spknt(NK+1),qq0(NK+1,NK+1),qq(3,NK+1,NK+1)

! three_d_mantle_model_variables
  type three_d_mantle_model_variables
    sequence
    double precision dvs_a(0:NK,0:NS,0:NS)
    double precision dvs_b(0:NK,0:NS,0:NS)
    double precision dvp_a(0:NK,0:NS,0:NS)
    double precision dvp_b(0:NK,0:NS,0:NS)
    double precision spknt(NK+1)
    double precision qq0(NK+1,NK+1)
    double precision qq(3,NK+1,NK+1)
  end type three_d_mantle_model_variables

  type (three_d_mantle_model_variables) D3MM_V
! three_d_mantle_model_variables


  integer i,j
  double precision qqwk(3,NK+1)

  D3MM_V%spknt(1) = -1.00000d0
  D3MM_V%spknt(2) = -0.78631d0
  D3MM_V%spknt(3) = -0.59207d0
  D3MM_V%spknt(4) = -0.41550d0
  D3MM_V%spknt(5) = -0.25499d0
  D3MM_V%spknt(6) = -0.10909d0
  D3MM_V%spknt(7) = 0.02353d0
  D3MM_V%spknt(8) = 0.14409d0
  D3MM_V%spknt(9) = 0.25367d0
  D3MM_V%spknt(10) = 0.35329d0
  D3MM_V%spknt(11) = 0.44384d0
  D3MM_V%spknt(12) = 0.52615d0
  D3MM_V%spknt(13) = 0.60097d0
  D3MM_V%spknt(14) = 0.66899d0
  D3MM_V%spknt(15) = 0.73081d0
  D3MM_V%spknt(16) = 0.78701d0
  D3MM_V%spknt(17) = 0.83810d0
  D3MM_V%spknt(18) = 0.88454d0
  D3MM_V%spknt(19) = 0.92675d0
  D3MM_V%spknt(20) = 0.96512d0
  D3MM_V%spknt(21) = 1.00000d0

  do i=1,NK+1
    do j=1,NK+1
      if(i == j) then
        D3MM_V%qq0(j,i)=1.0d0
      else
        D3MM_V%qq0(j,i)=0.0d0
      endif
    enddo
  enddo
  do i=1,NK+1
    call rspln(1,NK+1,D3MM_V%spknt(1),D3MM_V%qq0(1,i),D3MM_V%qq(1,1,i),qqwk(1,1))
  enddo

  end subroutine splhsetup

!----------------------------------

! changed the obsolecent f77 features in the two routines below
! now still awful Fortran, but at least conforms to f90 standard

  double precision function rsple(I1,I2,X,Y,Q,S)

  implicit none

! rsple returns the value of the function y(x) evaluated at point S
! using the cubic spline coefficients computed by rspln and saved in Q.
! If S is outside the interval (x(i1),x(i2)) rsple extrapolates
! using the first or last interpolation polynomial. The arrays must
! be dimensioned at least - x(i2), y(i2), and q(3,i2).

      integer i1,i2
      double precision  X(*),Y(*),Q(3,*),s

      integer i,ii
      double precision h

      i = 1
      II=I2-1

!   GUARANTEE I WITHIN BOUNDS.
      I=MAX0(I,I1)
      I=MIN0(I,II)

!   SEE IF X IS INCREASING OR DECREASING.
      IF(X(I2)-X(I1) <  0) goto 1
      IF(X(I2)-X(I1) >= 0) goto 2

!   X IS DECREASING.  CHANGE I AS NECESSARY.
 1    IF(S-X(I) <= 0) goto 3
      IF(S-X(I) >  0) goto 4

 4    I=I-1

      IF(I-I1 <  0) goto 11
      IF(I-I1 == 0) goto 6
      IF(I-I1 >  0) goto 1

 3    IF(S-X(I+1) <  0) goto 5
      IF(S-X(I+1) >= 0) goto 6

 5    I=I+1

      IF(I-II <  0) goto 3
      IF(I-II == 0) goto 6
      IF(I-II >  0) goto 7

!   X IS INCREASING.  CHANGE I AS NECESSARY.
 2    IF(S-X(I+1) <= 0) goto 8
      IF(S-X(I+1) >  0) goto 9

 9    I=I+1

      IF(I-II <  0) goto 2
      IF(I-II == 0) goto 6
      IF(I-II >  0) goto 7

 8    IF(S-X(I) <  0) goto 10
      IF(S-X(I) >= 0) goto 6

 10   I=I-1
      IF(I-I1 <  0) goto 11
      IF(I-I1 == 0) goto 6
      IF(I-I1 >  0) goto 8

 7    I=II
      GOTO 6
 11   I=I1

!   CALCULATE RSPLE USING SPLINE COEFFICIENTS IN Y AND Q.
 6    H=S-X(I)
      RSPLE=Y(I)+H*(Q(1,I)+H*(Q(2,I)+H*Q(3,I)))

      end function rsple

!----------------------------------

  subroutine rspln(I1,I2,X,Y,Q,F)

  implicit none

! Subroutine rspln computes cubic spline interpolation coefficients
! for y(x) between grid points i1 and i2 saving them in q.The
! interpolation is continuous with continuous first and second
! derivatives. It agrees exactly with y at grid points and with the
! three point first derivatives at both end points (i1 and i2).
! X must be monotonic but if two successive values of x are equal
! a discontinuity is assumed and separate interpolation is done on
! each strictly monotonic segment. The arrays must be dimensioned at
! least - x(i2), y(i2), q(3,i2), and f(3,i2).
! F is working storage for rspln.

      integer i1,i2
      double precision X(*),Y(*),Q(3,*),F(3,*)

      integer i,j,k,j1,j2
      double precision y0,a0,b0,b1,h,h2,ha,h2a,h3a,h2b
      double precision YY(3),small
      equivalence (YY(1),Y0)
      data SMALL/1.0d-08/,YY/0.0d0,0.0d0,0.0d0/

      J1=I1+1
      Y0=0.0d0

!   BAIL OUT IF THERE ARE LESS THAN TWO POINTS TOTAL
      IF(I2-I1  < 0) return
      IF(I2-I1 == 0) goto 17
      IF(I2-I1  > 0) goto 8

 8    A0=X(J1-1)
!   SEARCH FOR DISCONTINUITIES.
      DO 3 I=J1,I2
      B0=A0
      A0=X(I)
      IF(DABS((A0-B0)/DMAX1(A0,B0)) < SMALL) GOTO 4
 3    CONTINUE
 17   J1=J1-1
      J2=I2-2
      GOTO 5
 4    J1=J1-1
      J2=I-3
!   SEE IF THERE ARE ENOUGH POINTS TO INTERPOLATE (AT LEAST THREE).
 5    IF(J2+1-J1 <  0) goto 9
      IF(J2+1-J1 == 0) goto 10
      IF(J2+1-J1 >  0) goto 11

!   ONLY TWO POINTS.  USE LINEAR INTERPOLATION.
 10   J2=J2+2
      Y0=(Y(J2)-Y(J1))/(X(J2)-X(J1))
      DO J=1,3
        Q(J,J1)=YY(J)
        Q(J,J2)=YY(J)
      enddo
      GOTO 12

!   MORE THAN TWO POINTS.  DO SPLINE INTERPOLATION.
 11   A0=0.
      H=X(J1+1)-X(J1)
      H2=X(J1+2)-X(J1)
      Y0=H*H2*(H2-H)
      H=H*H
      H2=H2*H2
!   CALCULATE DERIVITIVE AT NEAR END.
      B0=(Y(J1)*(H-H2)+Y(J1+1)*H2-Y(J1+2)*H)/Y0
      B1=B0

!   EXPLICITLY REDUCE BANDED MATRIX TO AN UPPER BANDED MATRIX.
      DO I=J1,J2
        H=X(I+1)-X(I)
        Y0=Y(I+1)-Y(I)
        H2=H*H
        HA=H-A0
        H2A=H-2.0d0*A0
        H3A=2.0d0*H-3.0d0*A0
        H2B=H2*B0
        Q(1,I)=H2/HA
        Q(2,I)=-HA/(H2A*H2)
        Q(3,I)=-H*H2A/H3A
        F(1,I)=(Y0-H*B0)/(H*HA)
        F(2,I)=(H2B-Y0*(2.0d0*H-A0))/(H*H2*H2A)
        F(3,I)=-(H2B-3.0d0*Y0*HA)/(H*H3A)
        A0=Q(3,I)
        B0=F(3,I)
      enddo

!   TAKE CARE OF LAST TWO ROWS.
      I=J2+1
      H=X(I+1)-X(I)
      Y0=Y(I+1)-Y(I)
      H2=H*H
      HA=H-A0
      H2A=H*HA
      H2B=H2*B0-Y0*(2.0d0*H-A0)
      Q(1,I)=H2/HA
      F(1,I)=(Y0-H*B0)/H2A
      HA=X(J2)-X(I+1)
      Y0=-H*HA*(HA+H)
      HA=HA*HA

!   CALCULATE DERIVATIVE AT FAR END.
      Y0=(Y(I+1)*(H2-HA)+Y(I)*HA-Y(J2)*H2)/Y0
      Q(3,I)=(Y0*H2A+H2B)/(H*H2*(H-2.0d0*A0))
      Q(2,I)=F(1,I)-Q(1,I)*Q(3,I)

!   SOLVE UPPER BANDED MATRIX BY REVERSE ITERATION.
      DO J=J1,J2
        K=I-1
        Q(1,I)=F(3,K)-Q(3,K)*Q(2,I)
        Q(3,K)=F(2,K)-Q(2,K)*Q(1,I)
        Q(2,K)=F(1,K)-Q(1,K)*Q(3,K)
        I=K
      enddo
      Q(1,I)=B1
!   FILL IN THE LAST POINT WITH A LINEAR EXTRAPOLATION.
 9    J2=J2+2
      DO J=1,3
        Q(J,J2)=YY(J)
      enddo

!   SEE IF THIS DISCONTINUITY IS THE LAST.
 12   IF(J2-I2 < 0) then
        goto 6
      else
        return
      endif

!   NO.  GO BACK FOR MORE.
 6    J1=J2+2
      IF(J1-I2 <= 0) goto 8
      IF(J1-I2 >  0) goto 7

!   THERE IS ONLY ONE POINT LEFT AFTER THE LATEST DISCONTINUITY.
 7    DO J=1,3
        Q(J,I2)=YY(J)
      enddo

      end subroutine rspln

