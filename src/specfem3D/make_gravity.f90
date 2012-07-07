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
!
! United States and French Government Sponsorship Acknowledged.


  subroutine make_gravity(nspl,rspl,gspl,gspl2, &
                          ROCEAN,RMIDDLE_CRUST,RMOHO,R80,R220,R400,R600,R670, &
                          R771,RTOPDDOUBLEPRIME,RCMB,RICB)

! creates a spline for the gravity profile in PREM
! radius and density are non-dimensional

  implicit none

  include "constants.h"

  integer:: nspl
  double precision:: rspl(NR),gspl(NR),gspl2(NR)
  double precision ROCEAN,RMIDDLE_CRUST,RMOHO,R80,R220,R400,R600,R670, &
                   R771,RTOPDDOUBLEPRIME,RCMB,RICB

  ! local parameters
  double precision r_icb,r_cmb,r_topddoubleprime,r_771,r_670,r_600
  double precision r_400,r_220,r_80,r_moho,r_middle_crust,r_ocean,r_0
  double precision r(NR),rho(NR),g(NR),i_rho
  double precision s1(NR),s2(NR),s3(NR)
  double precision yp1,ypn
  integer i

! PREM
  ! values in (m)
  ROCEAN = 6368000.d0
  RMIDDLE_CRUST = 6356000.d0
  RMOHO = 6346600.d0 ! PREM moho depth at 24.4 km
  R80  = 6291000.d0
  R220 = 6151000.d0
  R400 = 5971000.d0
  R600 = 5771000.d0
  R670 = 5701000.d0
  R771 = 5600000.d0
  RTOPDDOUBLEPRIME = 3630000.d0
  RCMB = 3480000.d0
  RICB = 1221000.d0

! non-dimensionalize
  r_icb = RICB/R_EARTH_GRAVITY
  r_cmb = RCMB/R_EARTH_GRAVITY
  r_topddoubleprime = RTOPDDOUBLEPRIME/R_EARTH_GRAVITY
  r_771 = R771/R_EARTH_GRAVITY
  r_670 = R670/R_EARTH_GRAVITY
  r_600 = R600/R_EARTH_GRAVITY
  r_400 = R400/R_EARTH_GRAVITY
  r_220 = R220/R_EARTH_GRAVITY
  r_80 = R80/R_EARTH_GRAVITY
  r_moho = RMOHO/R_EARTH_GRAVITY
  r_middle_crust = RMIDDLE_CRUST/R_EARTH_GRAVITY
  r_ocean = ROCEAN_GRAVITY/R_EARTH_GRAVITY
  r_0 = 1.d0

  do i=1,163
    r(i) = r_icb*dble(i-1)/dble(162)
  enddo
  do i=164,323
    r(i) = r_icb+(r_cmb-r_icb)*dble(i-164)/dble(159)
  enddo
  do i=324,336
    r(i) = r_cmb+(r_topddoubleprime-r_cmb)*dble(i-324)/dble(12)
  enddo
  do i=337,517
    r(i) = r_topddoubleprime+(r_771-r_topddoubleprime)*dble(i-337)/dble(180)
  enddo
  do i=518,530
    r(i) = r_771+(r_670-r_771)*dble(i-518)/dble(12)
  enddo
  do i=531,540
    r(i) = r_670+(r_600-r_670)*dble(i-531)/dble(9)
  enddo
  do i=541,565
    r(i) = r_600+(r_400-r_600)*dble(i-541)/dble(24)
  enddo
  do i=566,590
    r(i) = r_400+(r_220-r_400)*dble(i-566)/dble(24)
  enddo
  do i=591,609
    r(i) = r_220+(r_80-r_220)*dble(i-591)/dble(18)
  enddo
  do i=610,619
    r(i) = r_80+(r_moho-r_80)*dble(i-610)/dble(9)
  enddo
  do i=620,626
    r(i) = r_moho+(r_middle_crust-r_moho)*dble(i-620)/dble(6)
  enddo
  do i=627,633
    r(i) = r_middle_crust+(r_ocean-r_middle_crust)*dble(i-627)/dble(6)
  enddo
  do i=634,NR
    r(i) = r_ocean+(r_0-r_ocean)*dble(i-634)/dble(6)
  enddo

! use PREM to get the density profile for ellipticity (fine for other 1D reference models)
  do i=1,NR
    call prem_density(r(i),rho(i), &
                    RICB,RCMB,RTOPDDOUBLEPRIME, &
                    R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN_GRAVITY)
  enddo

  g(1)=0.0d0
  do i=2,NR
    call intgrl(i_rho,r,1,i,rho,s1,s2,s3)
    g(i)=4.0d0*i_rho/(r(i)*r(i))
  enddo

!
! get ready to spline g
!
  nspl=1
  rspl(1)=r(1)
  gspl(1)=g(1)
  do i=2,NR
    if(r(i)/=r(i-1)) then
      nspl=nspl+1
      rspl(nspl)=r(i)
      gspl(nspl)=g(i)
    endif
  enddo
  yp1=(4.0d0/3.0d0)*rho(1)
  ypn=4.0d0*rho(NR)-2.0d0*g(NR)/r(NR)
  call spline_construction(rspl,gspl,nspl,yp1,ypn,gspl2)

  end subroutine make_gravity


!--------------------------------------------------------------------------------------------------
!
! PREM [Dziewonski and Anderson, 1981].
!
! A. M. Dziewonski and D. L. Anderson.
! Preliminary reference Earth model.
! Phys. Earth Planet. Inter., 25:297–356, 1981.
!
! Isotropic (iso) and transversely isotropic (aniso) version of the
! spherically symmetric Preliminary Reference Earth Model
!
!--------------------------------------------------------------------------------------------------


  subroutine model_prem_iso(x,rho,drhodr,vp,vs,Qkappa,Qmu,&
                          RICB,RCMB,RTOPDDOUBLEPRIME, &
                          R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN)

  implicit none

  include "constants.h"

! given a normalized radius x, gives the non-dimensionalized density rho,
! speeds vp and vs, and the quality factors Qkappa and Qmu

  double precision x,rho,drhodr,vp,vs,Qkappa,Qmu,RICB,RCMB,RTOPDDOUBLEPRIME, &
      R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN

  double precision r,scaleval

  logical,parameter :: ONE_CRUST = .false.
  logical,parameter :: CRUSTAL = .false.

! compute real physical radius in meters
  r = x * R_EARTH

!
!--- inner core
!
  if(r >= 0.d0 .and. r <= RICB) then
    drhodr=-2.0d0*8.8381d0*x
    rho=13.0885d0-8.8381d0*x*x
    vp=11.2622d0-6.3640d0*x*x
    vs=3.6678d0-4.4475d0*x*x
    Qmu=84.6d0
    Qkappa=1327.7d0
!
!--- outer core
!
  else if(r > RICB .and. r <= RCMB) then
    drhodr=-1.2638d0-2.0d0*3.6426d0*x-3.0d0*5.5281d0*x*x
    rho=12.5815d0-1.2638d0*x-3.6426d0*x*x-5.5281d0*x*x*x
    vp=11.0487d0-4.0362d0*x+4.8023d0*x*x-13.5732d0*x*x*x
    vs=0.0d0
    Qmu=0.0d0
    Qkappa=57827.0d0
!
!--- D" at the base of the mantle
!
  else if(r > RCMB .and. r <= RTOPDDOUBLEPRIME) then
    drhodr=-6.4761d0+2.0d0*5.5283d0*x-3.0d0*3.0807d0*x*x
    rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
    vp=15.3891d0-5.3181d0*x+5.5242d0*x*x-2.5514d0*x*x*x
    vs=6.9254d0+1.4672d0*x-2.0834d0*x*x+0.9783d0*x*x*x
    Qmu=312.0d0
    Qkappa=57827.0d0
!
!--- mantle: from top of D" to d670
!
  else if(r > RTOPDDOUBLEPRIME .and. r <= R771) then
    drhodr=-6.4761d0+2.0d0*5.5283d0*x-3.0d0*3.0807d0*x*x
    rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
    vp=24.9520d0-40.4673d0*x+51.4832d0*x*x-26.6419d0*x*x*x
    vs=11.1671d0-13.7818d0*x+17.4575d0*x*x-9.2777d0*x*x*x
    Qmu=312.0d0
    Qkappa=57827.0d0
  else if(r > R771 .and. r <= R670) then
    drhodr=-6.4761d0+2.0d0*5.5283d0*x-3.0d0*3.0807d0*x*x
    rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
    vp=29.2766d0-23.6027d0*x+5.5242d0*x*x-2.5514d0*x*x*x
    vs=22.3459d0-17.2473d0*x-2.0834d0*x*x+0.9783d0*x*x*x
    Qmu=312.0d0
    Qkappa=57827.0d0
!
!--- mantle: above d670
!
  else if(r > R670 .and. r <= R600) then
    drhodr=-1.4836d0
    rho=5.3197d0-1.4836d0*x
    vp=19.0957d0-9.8672d0*x
    vs=9.9839d0-4.9324d0*x
    Qmu=143.0d0
    Qkappa=57827.0d0
  else if(r > R600 .and. r <= R400) then
    drhodr=-8.0298d0
    rho=11.2494d0-8.0298d0*x
    vp=39.7027d0-32.6166d0*x
    vs=22.3512d0-18.5856d0*x
    Qmu=143.0d0
    Qkappa=57827.0d0
  else if(r > R400 .and. r <= R220) then
    drhodr=-3.8045d0
    rho=7.1089d0-3.8045d0*x
    vp=20.3926d0-12.2569d0*x
    vs=8.9496d0-4.4597d0*x
    Qmu=143.0d0
    Qkappa=57827.0d0
  else if(r > R220 .and. r <= R80) then
    drhodr=0.6924d0
    rho=2.6910d0+0.6924d0*x
    vp=4.1875d0+3.9382d0*x
    vs=2.1519d0+2.3481d0*x
    Qmu=80.0d0
    Qkappa=57827.0d0
  else
    if(CRUSTAL) then
    ! fill with PREM mantle and later add CRUST2.0
      if(r > R80) then
        ! density/velocity from mantle just below moho
        drhodr=0.6924d0
        rho=2.6910d0+0.6924d0*x
        vp=4.1875d0+3.9382d0*x
        vs=2.1519d0+2.3481d0*x
        ! shear attenuation for R80 to surface
        Qmu=600.0d0
        Qkappa=57827.0d0
      endif
    else
    ! use PREM crust
      if(r > R80 .and. r <= RMOHO) then
        drhodr=0.6924d0
        rho=2.6910d0+0.6924d0*x
        vp=4.1875d0+3.9382d0*x
        vs=2.1519d0+2.3481d0*x
        Qmu=600.0d0
        Qkappa=57827.0d0

      else if(r > RMOHO .and. r <= RMIDDLE_CRUST) then
        drhodr=0.0d0
        rho=2.9d0
        vp=6.8d0
        vs=3.9d0
        Qmu=600.0d0
        Qkappa=57827.0d0

    ! same properties everywhere in PREM crust if we decide to define only one layer in the crust
        if(ONE_CRUST) then
          drhodr=0.0d0
          rho=2.6d0
          vp=5.8d0
          vs=3.2d0
          Qmu=600.0d0
          Qkappa=57827.0d0
        endif

      else if(r > RMIDDLE_CRUST .and. r <= ROCEAN) then
        drhodr=0.0d0
        rho=2.6d0
        vp=5.8d0
        vs=3.2d0
        Qmu=600.0d0
        Qkappa=57827.0d0
    ! for density profile for gravity, we do not check that r <= R_EARTH
      else if(r > ROCEAN) then
        drhodr=0.0d0
        rho=2.6d0
        vp=5.8d0
        vs=3.2d0
        Qmu=600.0d0
        Qkappa=57827.0d0

      endif
    endif
  endif

! non-dimensionalize
! time scaling (s^{-1}) is done with scaleval
  scaleval=dsqrt(PI*GRAV*RHOAV)
  drhodr=drhodr*1000.0d0/RHOAV
  rho=rho*1000.0d0/RHOAV
  vp=vp*1000.0d0/(R_EARTH*scaleval)
  vs=vs*1000.0d0/(R_EARTH*scaleval)

  end subroutine model_prem_iso

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prem_density(x,rho, &
                         RICB,RCMB,RTOPDDOUBLEPRIME, &
                         R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN)

  implicit none

  include "constants.h"

  double precision x,rho,RICB,RCMB,RTOPDDOUBLEPRIME, &
      R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN
  double precision r

  logical,parameter :: ONE_CRUST = .false.

  ! compute real physical radius in meters
  r = x * R_EARTH

  ! calculates density according to radius
  if(r <= RICB) then
    rho=13.0885d0-8.8381d0*x*x
  else if(r > RICB .and. r <= RCMB) then
    rho=12.5815d0-1.2638d0*x-3.6426d0*x*x-5.5281d0*x*x*x
  else if(r > RCMB .and. r <= RTOPDDOUBLEPRIME) then
    rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
  else if(r > RTOPDDOUBLEPRIME .and. r <= R771) then
    rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
  else if(r > R771 .and. r <= R670) then
    rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
  else if(r > R670 .and. r <= R600) then
    rho=5.3197d0-1.4836d0*x
  else if(r > R600 .and. r <= R400) then
    rho=11.2494d0-8.0298d0*x
  else if(r > R400 .and. r <= R220) then
    rho=7.1089d0-3.8045d0*x
  else if(r > R220 .and. r <= R80) then
    rho=2.6910d0+0.6924d0*x
  else
    if(r > R80 .and. r <= RMOHO) then
      rho=2.6910d0+0.6924d0*x
    else if(r > RMOHO .and. r <= RMIDDLE_CRUST) then
      if(ONE_CRUST) then
        rho=2.6d0
      else
        rho=2.9d0
      endif
    else if(r > RMIDDLE_CRUST .and. r <= ROCEAN) then
      rho=2.6d0
    else if(r > ROCEAN) then
      rho=2.6d0
    endif
  endif

  rho=rho*1000.0d0/RHOAV

  end subroutine prem_density

!
!-------------------------------------------------------------------------------------------------
!

  subroutine intgrl(sum,r,nir,ner,f,s1,s2,s3)

! Computes the integral of f[i]*r[i]*r[i] from i=nir to i=ner for
! radii values as in model PREM_an640

  implicit none

! Argument variables
  integer ner,nir
  double precision f(640),r(640),s1(640),s2(640)
  double precision s3(640),sum

! Local variables
  double precision, parameter :: third = 1.0d0/3.0d0
  double precision, parameter :: fifth = 1.0d0/5.0d0
  double precision, parameter :: sixth = 1.0d0/6.0d0

  double precision rji,yprime(640)
  double precision s1l,s2l,s3l

  integer i,j,n,kdis(28)
  integer ndis,nir1

  data kdis/163,323,336,517,530,540,565,590,609,619,626,633,16*0/

  ndis = 12
  n = 640

  call deriv(f,yprime,n,r,ndis,kdis,s1,s2,s3)
  nir1 = nir + 1
  sum = 0.0d0
  do i=nir1,ner
    j = i-1
    rji = r(i) - r(j)
    s1l = s1(j)
    s2l = s2(j)
    s3l = s3(j)
    sum = sum + r(j)*r(j)*rji*(f(j) &
              + rji*(0.5d0*s1l + rji*(third*s2l + rji*0.25d0*s3l))) &
              + 2.0d0*r(j)*rji*rji*(0.5d0*f(j) + rji*(third*s1l + rji*(0.25d0*s2l + rji*fifth*s3l))) &
              + rji*rji*rji*(third*f(j) + rji*(0.25d0*s1l + rji*(fifth*s2l + rji*sixth*s3l)))
  enddo

  end subroutine intgrl

!
!-------------------------------------------------------------------------------------------------
!

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

!
!-------------------------------------------------------------------------------------------------
!

! compute spline coefficients

  subroutine spline_construction(xpoint,ypoint,npoint,tangent_first_point,tangent_last_point,spline_coefficients)

  implicit none

! tangent to the spline imposed at the first and last points
  double precision, intent(in) :: tangent_first_point,tangent_last_point

! number of input points and coordinates of the input points
  integer, intent(in) :: npoint
  double precision, dimension(npoint), intent(in) :: xpoint,ypoint

! spline coefficients output by the routine
  double precision, dimension(npoint), intent(out) :: spline_coefficients

  integer :: i

  double precision, dimension(:), allocatable :: temporary_array

  allocate(temporary_array(npoint))

  spline_coefficients(1) = - 1.d0 / 2.d0

  temporary_array(1) = (3.d0/(xpoint(2)-xpoint(1)))*((ypoint(2)-ypoint(1))/(xpoint(2)-xpoint(1))-tangent_first_point)

  do i = 2,npoint-1

    spline_coefficients(i) = ((xpoint(i)-xpoint(i-1))/(xpoint(i+1)-xpoint(i-1))-1.d0) &
       / ((xpoint(i)-xpoint(i-1))/(xpoint(i+1)-xpoint(i-1))*spline_coefficients(i-1)+2.d0)

    temporary_array(i) = (6.d0*((ypoint(i+1)-ypoint(i))/(xpoint(i+1)-xpoint(i)) &
       - (ypoint(i)-ypoint(i-1))/(xpoint(i)-xpoint(i-1)))/(xpoint(i+1)-xpoint(i-1)) &
       - (xpoint(i)-xpoint(i-1))/(xpoint(i+1)-xpoint(i-1))*temporary_array(i-1)) &
       / ((xpoint(i)-xpoint(i-1))/(xpoint(i+1)-xpoint(i-1))*spline_coefficients(i-1)+2.d0)

  enddo

  spline_coefficients(npoint) = ((3.d0/(xpoint(npoint)-xpoint(npoint-1))) &
      * (tangent_last_point-(ypoint(npoint)-ypoint(npoint-1))/(xpoint(npoint)-xpoint(npoint-1))) &
      - 1.d0/2.d0*temporary_array(npoint-1))/(1.d0/2.d0*spline_coefficients(npoint-1)+1.d0)

  do i = npoint-1,1,-1
    spline_coefficients(i) = spline_coefficients(i)*spline_coefficients(i+1) + temporary_array(i)
  enddo

  deallocate(temporary_array)

  end subroutine spline_construction

!
!-------------------------------------------------------------------------------------------------
!

! evaluate a spline

  subroutine spline_evaluation(xpoint,ypoint,spline_coefficients,npoint,x_evaluate_spline,y_spline_obtained)

  implicit none

! number of input points and coordinates of the input points
  integer, intent(in) :: npoint
  double precision, dimension(npoint), intent(in) :: xpoint,ypoint

! spline coefficients to use
  double precision, dimension(npoint), intent(in) :: spline_coefficients

! abscissa at which we need to evaluate the value of the spline
  double precision, intent(in):: x_evaluate_spline

! ordinate evaluated by the routine for the spline at this abscissa
  double precision, intent(out):: y_spline_obtained

  integer :: index_loop,index_lower,index_higher

  double precision :: coef1,coef2

! initialize to the whole interval
  index_lower = 1
  index_higher = npoint

! determine the right interval to use, by dichotomy
  do while (index_higher - index_lower > 1)
! compute the middle of the interval
    index_loop = (index_higher + index_lower) / 2
    if(xpoint(index_loop) > x_evaluate_spline) then
      index_higher = index_loop
    else
      index_lower = index_loop
    endif
  enddo

! test that the interval obtained does not have a size of zero
! (this could happen for instance in the case of duplicates in the input list of points)
  if(xpoint(index_higher) == xpoint(index_lower)) stop 'incorrect interval found in spline evaluation'

  coef1 = (xpoint(index_higher) - x_evaluate_spline) / (xpoint(index_higher) - xpoint(index_lower))
  coef2 = (x_evaluate_spline - xpoint(index_lower)) / (xpoint(index_higher) - xpoint(index_lower))

  y_spline_obtained = coef1*ypoint(index_lower) + coef2*ypoint(index_higher) + &
        ((coef1**3 - coef1)*spline_coefficients(index_lower) + &
         (coef2**3 - coef2)*spline_coefficients(index_higher))*((xpoint(index_higher) - xpoint(index_lower))**2)/6.d0

  end subroutine spline_evaluation


