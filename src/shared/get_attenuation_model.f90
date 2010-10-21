!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  1 . 4
!               ---------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
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

  subroutine get_attenuation_model(myrank,iattenuation, &
         tau_mu,tau_sigma,beta,one_minus_sum_beta,factor_scale)

! return attenuation mechanisms Q_mu using standard linear solids
! frequency range: 20.000000 -- 1000.000000 mHz
! period range: 1.000000 -- 50.000000 s
! central logarithmic frequency: 0.141421 Hz
! the Tau values computed by Jeroen's code are used
! number of relaxation mechanisms: 3

! in the future when more memory is available on computers
! it would be more accurate to use four mechanisms instead of three

  implicit none

  include "constants.h"

! define central frequency of source in seconds using values from Jeroen's code
! logarithmic mean of frequency interval
  double precision, parameter :: f_c_source = 0.141421d0

! reference frequency for target velocity values in velocity model
! arbitrarily set to typical resolution of model (3 sec)
  double precision, parameter :: f0_REFERENCE = 0.3d0

  integer iattenuation,myrank

  double precision, dimension(N_SLS) :: tau_mu,tau_sigma,beta
  double precision one_minus_sum_beta

  integer i

  double precision Q_mu,w_c_source
  double precision factor_scale_mu0,factor_scale_mu,factor_scale
  double precision a_val,b_val,big_omega

! check number of SLS is okay
  if(N_SLS /= 3) call exit_MPI(myrank,'wrong number of SLS for attenuation, must be 3')

! clear arrays
  tau_mu(:) = 0.d0
  tau_sigma(:) = 0.d0

! tau sigma evenly spaced in log frequency, does not depend on value of Q
  tau_sigma( 1) =  7.957747154594766669788441504352d0
  tau_sigma( 2) =  1.125395395196382652969191440206d0
  tau_sigma( 3) =  0.159154943091895345608222100964d0

! determine in which region we are based upon doubling flag

  select case(iattenuation)

!--- sediments

! select value needed here, from Q_mu = 40 to Q_mu = 150

  case(IATTENUATION_SEDIMENTS_40)

  Q_mu = 40.000000d0
  tau_mu( 1) = 8.207413221956890936326090013608d0
  tau_mu( 2) = 1.161729745747647424281012717984d0
  tau_mu( 3) = 0.165834182312059152941685624683d0

  case(IATTENUATION_SEDIMENTS_50)

  Q_mu = 50.000000d0
  tau_mu( 1) = 8.169307711419302009403509146068d0
  tau_mu( 2) = 1.153839195800796080249028818798d0
  tau_mu( 3) = 0.164437605011117371489604011003d0

  case(IATTENUATION_SEDIMENTS_60)

  Q_mu = 60.000000d0
  tau_mu( 1) = 8.140254475505114939437589782756d0
  tau_mu( 2) = 1.148759228190431747052002720011d0
  tau_mu( 3) = 0.163522774234807849458306350243d0

  case(IATTENUATION_SEDIMENTS_70)

  Q_mu = 70.000000d0
  tau_mu( 1) = 8.117833196570874321196242817678d0
  tau_mu( 2) = 1.145216760190841176481058028003d0
  tau_mu( 3) = 0.162877472647593862786763452277d0

  case(IATTENUATION_SEDIMENTS_80)

  Q_mu = 80.000000d0
  tau_mu( 1) = 8.100148465407393416626291582361d0
  tau_mu( 2) = 1.142606124533341649396334105404d0
  tau_mu( 3) = 0.162398031255151509277823151933d0

  case(IATTENUATION_SEDIMENTS_90)

  Q_mu = 90.000000d0
  tau_mu( 1) = 8.085897732468197318667080253363d0
  tau_mu( 2) = 1.140602642076625095057806902332d0
  tau_mu( 3) = 0.162027854074084459723437134926d0

  case(IATTENUATION_SEDIMENTS_100)

  Q_mu = 100.000000d0
  tau_mu( 1) = 8.074193745349216300155603676103d0
  tau_mu( 2) = 1.139016691991711960341149278975d0
  tau_mu( 3) = 0.161733443689579814428469717313d0

  case(IATTENUATION_SEDIMENTS_110)

  Q_mu = 110.000000d0
  tau_mu( 1) = 8.064421863800781409281626110896d0
  tau_mu( 2) = 1.137730132230029722606445830024d0
  tau_mu( 3) = 0.161493715940844051459635011270d0

  case(IATTENUATION_SEDIMENTS_120)

  Q_mu = 120.000000d0
  tau_mu( 1) = 8.056146565814696458573962445371d0
  tau_mu( 2) = 1.136665532765689157201904890826d0
  tau_mu( 3) = 0.161294740739552050490246415393d0

  case(IATTENUATION_SEDIMENTS_130)

  Q_mu = 130.000000d0
  tau_mu( 1) = 8.049052148467024991873586259317d0
  tau_mu( 2) = 1.135770035674695810357093250786d0
  tau_mu( 3) = 0.161126946571733903335044146843d0

  case(IATTENUATION_SEDIMENTS_140)

  Q_mu = 140.000000d0
  tau_mu( 1) = 8.042904857756342451580167107750d0
  tau_mu( 2) = 1.135006327178704310654211440124d0
  tau_mu( 3) = 0.160983540254336005004276444197d0

  case(IATTENUATION_SEDIMENTS_150)

  Q_mu = 150.000000d0
  tau_mu( 1) = 8.037528252037535736462814384140d0
  tau_mu( 2) = 1.134347316535732730358176922891d0
  tau_mu( 3) = 0.160859567464536307168643247678d0

!--- bedrock

  case(IATTENUATION_BEDROCK)

    tau_mu( 1) = 7.959142154402283786396310460987d0
    tau_mu( 2) = 1.125540477911388892451327592426d0
    tau_mu( 3) = 0.159182872336587483141912002793d0

    Q_mu = 9000.d0

  case default

    call exit_MPI(myrank,'wrong attenuation flag in mesh')

  end select

!--- compute beta
    beta(:) = 1.d0 - tau_mu(:) / tau_sigma(:)

!--- compute central angular frequency of source
    w_c_source = TWO_PI * f_c_source

!--- quantity by which to scale mu_0 to get mu
    factor_scale_mu0 = ONE + TWO * log(f_c_source / f0_REFERENCE) / (PI * Q_mu)

!--- compute a, b and Omega parameters, also compute one minus sum of betas
  a_val = ONE
  b_val = ZERO
  one_minus_sum_beta = ONE

  do i = 1,N_SLS
    a_val = a_val - w_c_source * w_c_source * tau_mu(i) * &
      (tau_mu(i) - tau_sigma(i)) / (1.d0 + w_c_source * w_c_source * tau_mu(i) * tau_mu(i))
    b_val = b_val + w_c_source * (tau_mu(i) - tau_sigma(i)) / &
      (1.d0 + w_c_source * w_c_source * tau_mu(i) * tau_mu(i))
    one_minus_sum_beta = one_minus_sum_beta - beta(i)
  enddo

  big_omega = a_val*(sqrt(1.d0 + b_val*b_val/(a_val*a_val))-1.d0)

!--- quantity by which to scale mu to get mu_relaxed
  factor_scale_mu = b_val * b_val / (TWO * big_omega)

!--- total factor by which to scale mu0
  factor_scale = factor_scale_mu * factor_scale_mu0

!--- check that the correction factor is close to one
  if(factor_scale < 0.9 .or. factor_scale > 1.1) &
    call exit_MPI(myrank,'incorrect correction factor in attenuation model')

  end subroutine get_attenuation_model


!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_attenuation_model_olsen( vs_val, iselected )

! uses scaling rule similar to Olsen et al. (2003) to determine attenuation medium
!
! returns: selected sediment iselected
  
  implicit none
  
  include "constants.h"
  
  real(kind=CUSTOM_REAL) :: vs_val  
  integer :: iselected

!local parameters
  real(kind=CUSTOM_REAL) :: Q_mu
  integer :: int_Q_mu,iattenuation_sediments
  
  ! use rule Q_mu = constant * v_s
  Q_mu = OLSEN_ATTENUATION_RATIO * vs_val
  int_Q_mu = 10 * nint(Q_mu / 10.)
  
  if(int_Q_mu < 40) int_Q_mu = 40
  if(int_Q_mu > 150) int_Q_mu = 150

  if(int_Q_mu == 40) then
    iattenuation_sediments = IATTENUATION_SEDIMENTS_40
  else if(int_Q_mu == 50) then
    iattenuation_sediments = IATTENUATION_SEDIMENTS_50
  else if(int_Q_mu == 60) then
    iattenuation_sediments = IATTENUATION_SEDIMENTS_60
  else if(int_Q_mu == 70) then
    iattenuation_sediments = IATTENUATION_SEDIMENTS_70
  else if(int_Q_mu == 80) then
    iattenuation_sediments = IATTENUATION_SEDIMENTS_80
  else if(int_Q_mu == 90) then
    iattenuation_sediments = IATTENUATION_SEDIMENTS_90
  else if(int_Q_mu == 100) then
    iattenuation_sediments = IATTENUATION_SEDIMENTS_100
  else if(int_Q_mu == 110) then
    iattenuation_sediments = IATTENUATION_SEDIMENTS_110
  else if(int_Q_mu == 120) then
    iattenuation_sediments = IATTENUATION_SEDIMENTS_120
  else if(int_Q_mu == 130) then
    iattenuation_sediments = IATTENUATION_SEDIMENTS_130
  else if(int_Q_mu == 140) then
    iattenuation_sediments = IATTENUATION_SEDIMENTS_140
  else if(int_Q_mu == 150) then
    iattenuation_sediments = IATTENUATION_SEDIMENTS_150
  else
    stop 'incorrect attenuation coefficient'
  endif
  
  ! return sediment number
  iselected = iattenuation_sediments  
  
  end subroutine get_attenuation_model_olsen
