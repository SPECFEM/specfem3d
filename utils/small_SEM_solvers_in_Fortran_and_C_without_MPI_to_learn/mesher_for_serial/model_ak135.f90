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

  subroutine model_ak135(x,rho,vp,vs,Qkappa,Qmu,iregion_code,Mak135_V)

  implicit none

  include "constants.h"

! model_ak135_variables
  type model_ak135_variables
    sequence
    double precision, dimension(NR_AK135) :: radius_ak135
    double precision, dimension(NR_AK135) :: density_ak135
    double precision, dimension(NR_AK135) :: vp_ak135
    double precision, dimension(NR_AK135) :: vs_ak135
    double precision, dimension(NR_AK135) :: Qkappa_ak135
    double precision, dimension(NR_AK135) :: Qmu_ak135
  end type model_ak135_variables

 type (model_ak135_variables) Mak135_V
! model_ak135_variables

! input:
! radius r: meters

! output:
! density rho: kg/m^3
! compressional wave speed vp: km/s
! shear wave speed vs: km/s

  integer iregion_code

  double precision x,rho,vp,vs,Qmu,Qkappa

  integer i

  double precision r,frac,scaleval

!! DK DK UGLY implementation of model ak135 below and its radii in
!! DK DK UGLY subroutine read_parameter_file.f90 has not been thoroughly
!! DK DK UGLY checked yet

! compute real physical radius in meters
  r = x * R_EARTH

  i = 1
  do while(r >= Mak135_V%radius_ak135(i) .and. i /= NR_AK135)
    i = i + 1
  enddo

! make sure we stay in the right region and never take a point above
! and a point below the ICB or the CMB and interpolate between them,
! which would lead to a wrong value (keeping in mind that we interpolate
! between points i-1 and i below)
  if(iregion_code == IREGION_INNER_CORE .and. i > 25) i = 25

  if(iregion_code == IREGION_OUTER_CORE .and. i < 27) i = 27
  if(iregion_code == IREGION_OUTER_CORE .and. i > 71) i = 71

  if(iregion_code == IREGION_CRUST_MANTLE .and. i < 73) i = 73

  if(i == 1) then
    rho = Mak135_V%density_ak135(i)
    vp = Mak135_V%vp_ak135(i)
    vs = Mak135_V%vs_ak135(i)
    Qmu = Mak135_V%Qmu_ak135(i)
    Qkappa = Mak135_V%Qkappa_ak135(i)
  else

! interpolate from radius_ak135(i-1) to r using the values at i-1 and i
    frac = (r-Mak135_V%radius_ak135(i-1))/(Mak135_V%radius_ak135(i)-Mak135_V%radius_ak135(i-1))

    rho = Mak135_V%density_ak135(i-1) + frac * (Mak135_V%density_ak135(i)-Mak135_V%density_ak135(i-1))
    vp = Mak135_V%vp_ak135(i-1) + frac * (Mak135_V%vp_ak135(i)-Mak135_V%vp_ak135(i-1))
    vs = Mak135_V%vs_ak135(i-1) + frac * (Mak135_V%vs_ak135(i)-Mak135_V%vs_ak135(i-1))
    Qmu = Mak135_V%Qmu_ak135(i-1) + frac * (Mak135_V%Qmu_ak135(i)-Mak135_V%Qmu_ak135(i-1))
    Qkappa = Mak135_V%Qkappa_ak135(i-1) + frac * (Mak135_V%Qkappa_ak135(i)-Mak135_V%Qkappa_ak135(i-1))

  endif

! make sure Vs is zero in the outer core even if roundoff errors on depth
! also set fictitious attenuation to a very high value (attenuation is not used in the fluid)
  if(iregion_code == IREGION_OUTER_CORE) then
    vs = 0.d0
    Qkappa = 3000.d0
    Qmu = 3000.d0
  endif

! non-dimensionalize
! time scaling (s^{-1}) is done with scaleval
  scaleval=dsqrt(PI*GRAV*RHOAV)
  rho=rho*1000.0d0/RHOAV
  vp=vp*1000.0d0/(R_EARTH*scaleval)
  vs=vs*1000.0d0/(R_EARTH*scaleval)

  end subroutine model_ak135

!-------------------

  subroutine define_model_ak135(USE_EXTERNAL_CRUSTAL_MODEL,Mak135_V)

  implicit none
  include "constants.h"

! model_ak135_variables
  type model_ak135_variables
    sequence
    double precision, dimension(NR_AK135) :: radius_ak135
    double precision, dimension(NR_AK135) :: density_ak135
    double precision, dimension(NR_AK135) :: vp_ak135
    double precision, dimension(NR_AK135) :: vs_ak135
    double precision, dimension(NR_AK135) :: Qkappa_ak135
    double precision, dimension(NR_AK135) :: Qmu_ak135
  end type model_ak135_variables

 type (model_ak135_variables) Mak135_V
! model_ak135_variables

  logical USE_EXTERNAL_CRUSTAL_MODEL

  integer i

! define all the values in the model

  Mak135_V%radius_ak135(  1) =  0.000000000000000
  Mak135_V%radius_ak135(  2) =   50710.0000000000
  Mak135_V%radius_ak135(  3) =   101430.000000000
  Mak135_V%radius_ak135(  4) =   152140.000000000
  Mak135_V%radius_ak135(  5) =   202850.000000000
  Mak135_V%radius_ak135(  6) =   253560.000000000
  Mak135_V%radius_ak135(  7) =   304280.000000000
  Mak135_V%radius_ak135(  8) =   354990.000000000
  Mak135_V%radius_ak135(  9) =   405700.000000000
  Mak135_V%radius_ak135( 10) =   456410.000000000
  Mak135_V%radius_ak135( 11) =   507130.000000000
  Mak135_V%radius_ak135( 12) =   557840.000000000
  Mak135_V%radius_ak135( 13) =   608550.000000000
  Mak135_V%radius_ak135( 14) =   659260.000000000
  Mak135_V%radius_ak135( 15) =   709980.000000000
  Mak135_V%radius_ak135( 16) =   760690.000000000
  Mak135_V%radius_ak135( 17) =   811400.000000000
  Mak135_V%radius_ak135( 18) =   862110.000000000
  Mak135_V%radius_ak135( 19) =   912830.000000000
  Mak135_V%radius_ak135( 20) =   963540.000000000
  Mak135_V%radius_ak135( 21) =   1014250.00000000
  Mak135_V%radius_ak135( 22) =   1064960.00000000
  Mak135_V%radius_ak135( 23) =   1115680.00000000
  Mak135_V%radius_ak135( 24) =   1166390.00000000
  Mak135_V%radius_ak135( 25) =   1217500.00000000
  Mak135_V%radius_ak135( 26) =   1217500.00000000
  Mak135_V%radius_ak135( 27) =   1267430.00000000
  Mak135_V%radius_ak135( 28) =   1317760.00000000
  Mak135_V%radius_ak135( 29) =   1368090.00000000
  Mak135_V%radius_ak135( 30) =   1418420.00000000
  Mak135_V%radius_ak135( 31) =   1468760.00000000
  Mak135_V%radius_ak135( 32) =   1519090.00000000
  Mak135_V%radius_ak135( 33) =   1569420.00000000
  Mak135_V%radius_ak135( 34) =   1619750.00000000
  Mak135_V%radius_ak135( 35) =   1670080.00000000
  Mak135_V%radius_ak135( 36) =   1720410.00000000
  Mak135_V%radius_ak135( 37) =   1770740.00000000
  Mak135_V%radius_ak135( 38) =   1821070.00000000
  Mak135_V%radius_ak135( 39) =   1871400.00000000
  Mak135_V%radius_ak135( 40) =   1921740.00000000
  Mak135_V%radius_ak135( 41) =   1972070.00000000
  Mak135_V%radius_ak135( 42) =   2022400.00000000
  Mak135_V%radius_ak135( 43) =   2072730.00000000
  Mak135_V%radius_ak135( 44) =   2123060.00000000
  Mak135_V%radius_ak135( 45) =   2173390.00000000
  Mak135_V%radius_ak135( 46) =   2223720.00000000
  Mak135_V%radius_ak135( 47) =   2274050.00000000
  Mak135_V%radius_ak135( 48) =   2324380.00000000
  Mak135_V%radius_ak135( 49) =   2374720.00000000
  Mak135_V%radius_ak135( 50) =   2425050.00000000
  Mak135_V%radius_ak135( 51) =   2475380.00000000
  Mak135_V%radius_ak135( 52) =   2525710.00000000
  Mak135_V%radius_ak135( 53) =   2576040.00000000
  Mak135_V%radius_ak135( 54) =   2626370.00000000
  Mak135_V%radius_ak135( 55) =   2676700.00000000
  Mak135_V%radius_ak135( 56) =   2727030.00000000
  Mak135_V%radius_ak135( 57) =   2777360.00000000
  Mak135_V%radius_ak135( 58) =   2827700.00000000
  Mak135_V%radius_ak135( 59) =   2878030.00000000
  Mak135_V%radius_ak135( 60) =   2928360.00000000
  Mak135_V%radius_ak135( 61) =   2978690.00000000
  Mak135_V%radius_ak135( 62) =   3029020.00000000
  Mak135_V%radius_ak135( 63) =   3079350.00000000
  Mak135_V%radius_ak135( 64) =   3129680.00000000
  Mak135_V%radius_ak135( 65) =   3180010.00000000
  Mak135_V%radius_ak135( 66) =   3230340.00000000
  Mak135_V%radius_ak135( 67) =   3280680.00000000
  Mak135_V%radius_ak135( 68) =   3331010.00000000
  Mak135_V%radius_ak135( 69) =   3381340.00000000
  Mak135_V%radius_ak135( 70) =   3431670.00000000
  Mak135_V%radius_ak135( 71) =   3479500.00000000
  Mak135_V%radius_ak135( 72) =   3479500.00000000
  Mak135_V%radius_ak135( 73) =   3531670.00000000
  Mak135_V%radius_ak135( 74) =   3581330.00000000
  Mak135_V%radius_ak135( 75) =   3631000.00000000
  Mak135_V%radius_ak135( 76) =   3631000.00000000
  Mak135_V%radius_ak135( 77) =   3681000.00000000
  Mak135_V%radius_ak135( 78) =   3731000.00000000
  Mak135_V%radius_ak135( 79) =   3779500.00000000
  Mak135_V%radius_ak135( 80) =   3829000.00000000
  Mak135_V%radius_ak135( 81) =   3878500.00000000
  Mak135_V%radius_ak135( 82) =   3928000.00000000
  Mak135_V%radius_ak135( 83) =   3977500.00000000
  Mak135_V%radius_ak135( 84) =   4027000.00000000
  Mak135_V%radius_ak135( 85) =   4076500.00000000
  Mak135_V%radius_ak135( 86) =   4126000.00000000
  Mak135_V%radius_ak135( 87) =   4175500.00000000
  Mak135_V%radius_ak135( 88) =   4225000.00000000
  Mak135_V%radius_ak135( 89) =   4274500.00000000
  Mak135_V%radius_ak135( 90) =   4324000.00000000
  Mak135_V%radius_ak135( 91) =   4373500.00000000
  Mak135_V%radius_ak135( 92) =   4423000.00000000
  Mak135_V%radius_ak135( 93) =   4472500.00000000
  Mak135_V%radius_ak135( 94) =   4522000.00000000
  Mak135_V%radius_ak135( 95) =   4571500.00000000
  Mak135_V%radius_ak135( 96) =   4621000.00000000
  Mak135_V%radius_ak135( 97) =   4670500.00000000
  Mak135_V%radius_ak135( 98) =   4720000.00000000
  Mak135_V%radius_ak135( 99) =   4769500.00000000
  Mak135_V%radius_ak135(100) =   4819000.00000000
  Mak135_V%radius_ak135(101) =   4868500.00000000
  Mak135_V%radius_ak135(102) =   4918000.00000000
  Mak135_V%radius_ak135(103) =   4967500.00000000
  Mak135_V%radius_ak135(104) =   5017000.00000000
  Mak135_V%radius_ak135(105) =   5066500.00000000
  Mak135_V%radius_ak135(106) =   5116000.00000000
  Mak135_V%radius_ak135(107) =   5165500.00000000
  Mak135_V%radius_ak135(108) =   5215000.00000000
  Mak135_V%radius_ak135(109) =   5264500.00000000
  Mak135_V%radius_ak135(110) =   5314000.00000000
  Mak135_V%radius_ak135(111) =   5363500.00000000
  Mak135_V%radius_ak135(112) =   5413000.00000000
  Mak135_V%radius_ak135(113) =   5462500.00000000
  Mak135_V%radius_ak135(114) =   5512000.00000000
  Mak135_V%radius_ak135(115) =   5561500.00000000
  Mak135_V%radius_ak135(116) =   5611000.00000000
  Mak135_V%radius_ak135(117) =   5661000.00000000
  Mak135_V%radius_ak135(118) =   5711000.00000000
  Mak135_V%radius_ak135(119) =   5711000.00000000
  Mak135_V%radius_ak135(120) =   5761000.00000000
  Mak135_V%radius_ak135(121) =   5811000.00000000
  Mak135_V%radius_ak135(122) =   5861000.00000000
  Mak135_V%radius_ak135(123) =   5911000.00000000
  Mak135_V%radius_ak135(124) =   5961000.00000000
  Mak135_V%radius_ak135(125) =   5961000.00000000
  Mak135_V%radius_ak135(126) =   6011000.00000000
  Mak135_V%radius_ak135(127) =   6061000.00000000
  Mak135_V%radius_ak135(128) =   6111000.00000000
  Mak135_V%radius_ak135(129) =   6161000.00000000
  Mak135_V%radius_ak135(130) =   6161000.00000000
  Mak135_V%radius_ak135(131) =   6206000.00000000
  Mak135_V%radius_ak135(132) =   6251000.00000000
  Mak135_V%radius_ak135(133) =   6291000.00000000
  Mak135_V%radius_ak135(134) =   6291000.00000000
  Mak135_V%radius_ak135(135) =   6328000.00000000
  Mak135_V%radius_ak135(136) =   6353000.00000000
  Mak135_V%radius_ak135(137) =   6353000.00000000
  Mak135_V%radius_ak135(138) =   6361000.00000000
  Mak135_V%radius_ak135(139) =   6361000.00000000
  Mak135_V%radius_ak135(140) =   6367700.00000000
  Mak135_V%radius_ak135(141) =   6367700.00000000
  Mak135_V%radius_ak135(142) =   6368000.00000000
  Mak135_V%radius_ak135(143) =   6368000.00000000
  Mak135_V%radius_ak135(144) =   6371000.00000000

  Mak135_V%density_ak135(  1) =   13.0122000000000
  Mak135_V%density_ak135(  2) =   13.0117000000000
  Mak135_V%density_ak135(  3) =   13.0100000000000
  Mak135_V%density_ak135(  4) =   13.0074000000000
  Mak135_V%density_ak135(  5) =   13.0036000000000
  Mak135_V%density_ak135(  6) =   12.9988000000000
  Mak135_V%density_ak135(  7) =   12.9929000000000
  Mak135_V%density_ak135(  8) =   12.9859000000000
  Mak135_V%density_ak135(  9) =   12.9779000000000
  Mak135_V%density_ak135( 10) =   12.9688000000000
  Mak135_V%density_ak135( 11) =   12.9586000000000
  Mak135_V%density_ak135( 12) =   12.9474000000000
  Mak135_V%density_ak135( 13) =   12.9351000000000
  Mak135_V%density_ak135( 14) =   12.9217000000000
  Mak135_V%density_ak135( 15) =   12.9072000000000
  Mak135_V%density_ak135( 16) =   12.8917000000000
  Mak135_V%density_ak135( 17) =   12.8751000000000
  Mak135_V%density_ak135( 18) =   12.8574000000000
  Mak135_V%density_ak135( 19) =   12.8387000000000
  Mak135_V%density_ak135( 20) =   12.8188000000000
  Mak135_V%density_ak135( 21) =   12.7980000000000
  Mak135_V%density_ak135( 22) =   12.7760000000000
  Mak135_V%density_ak135( 23) =   12.7530000000000
  Mak135_V%density_ak135( 24) =   12.7289000000000
  Mak135_V%density_ak135( 25) =   12.7037000000000
  Mak135_V%density_ak135( 26) =   12.1391000000000
  Mak135_V%density_ak135( 27) =   12.1133000000000
  Mak135_V%density_ak135( 28) =   12.0867000000000
  Mak135_V%density_ak135( 29) =   12.0593000000000
  Mak135_V%density_ak135( 30) =   12.0311000000000
  Mak135_V%density_ak135( 31) =   12.0001000000000
  Mak135_V%density_ak135( 32) =   11.9722000000000
  Mak135_V%density_ak135( 33) =   11.9414000000000
  Mak135_V%density_ak135( 34) =   11.9098000000000
  Mak135_V%density_ak135( 35) =   11.8772000000000
  Mak135_V%density_ak135( 36) =   11.8437000000000
  Mak135_V%density_ak135( 37) =   11.8092000000000
  Mak135_V%density_ak135( 38) =   11.7737000000000
  Mak135_V%density_ak135( 39) =   11.7373000000000
  Mak135_V%density_ak135( 40) =   11.6998000000000
  Mak135_V%density_ak135( 41) =   11.6612000000000
  Mak135_V%density_ak135( 42) =   11.6216000000000
  Mak135_V%density_ak135( 43) =   11.5809000000000
  Mak135_V%density_ak135( 44) =   11.5391000000000
  Mak135_V%density_ak135( 45) =   11.4962000000000
  Mak135_V%density_ak135( 46) =   11.4521000000000
  Mak135_V%density_ak135( 47) =   11.4069000000000
  Mak135_V%density_ak135( 48) =   11.3604000000000
  Mak135_V%density_ak135( 49) =   11.3127000000000
  Mak135_V%density_ak135( 50) =   11.2639000000000
  Mak135_V%density_ak135( 51) =   11.2137000000000
  Mak135_V%density_ak135( 52) =   11.1623000000000
  Mak135_V%density_ak135( 53) =   11.1095000000000
  Mak135_V%density_ak135( 54) =   11.0555000000000
  Mak135_V%density_ak135( 55) =   11.0001000000000
  Mak135_V%density_ak135( 56) =   10.9434000000000
  Mak135_V%density_ak135( 57) =   10.8852000000000
  Mak135_V%density_ak135( 58) =   10.8257000000000
  Mak135_V%density_ak135( 59) =   10.7647000000000
  Mak135_V%density_ak135( 60) =   10.7023000000000
  Mak135_V%density_ak135( 61) =   10.6385000000000
  Mak135_V%density_ak135( 62) =   10.5731000000000
  Mak135_V%density_ak135( 63) =   10.5062000000000
  Mak135_V%density_ak135( 64) =   10.4378000000000
  Mak135_V%density_ak135( 65) =   10.3679000000000
  Mak135_V%density_ak135( 66) =   10.2964000000000
  Mak135_V%density_ak135( 67) =   10.2233000000000
  Mak135_V%density_ak135( 68) =   10.1485000000000
  Mak135_V%density_ak135( 69) =   10.0722000000000
  Mak135_V%density_ak135( 70) =   9.99420000000000
  Mak135_V%density_ak135( 71) =   9.91450000000000
  Mak135_V%density_ak135( 72) =   5.77210000000000
  Mak135_V%density_ak135( 73) =   5.74580000000000
  Mak135_V%density_ak135( 74) =   5.71960000000000
  Mak135_V%density_ak135( 75) =   5.69340000000000
  Mak135_V%density_ak135( 76) =   5.43870000000000
  Mak135_V%density_ak135( 77) =   5.41760000000000
  Mak135_V%density_ak135( 78) =   5.39620000000000
  Mak135_V%density_ak135( 79) =   5.37480000000000
  Mak135_V%density_ak135( 80) =   5.35310000000000
  Mak135_V%density_ak135( 81) =   5.33130000000000
  Mak135_V%density_ak135( 82) =   5.30920000000000
  Mak135_V%density_ak135( 83) =   5.28700000000000
  Mak135_V%density_ak135( 84) =   5.26460000000000
  Mak135_V%density_ak135( 85) =   5.24200000000000
  Mak135_V%density_ak135( 86) =   5.21920000000000
  Mak135_V%density_ak135( 87) =   5.19630000000000
  Mak135_V%density_ak135( 88) =   5.17320000000000
  Mak135_V%density_ak135( 89) =   5.14990000000000
  Mak135_V%density_ak135( 90) =   5.12640000000000
  Mak135_V%density_ak135( 91) =   5.10270000000000
  Mak135_V%density_ak135( 92) =   5.07890000000000
  Mak135_V%density_ak135( 93) =   5.05480000000000
  Mak135_V%density_ak135( 94) =   5.03060000000000
  Mak135_V%density_ak135( 95) =   5.00620000000000
  Mak135_V%density_ak135( 96) =   4.98170000000000
  Mak135_V%density_ak135( 97) =   4.95700000000000
  Mak135_V%density_ak135( 98) =   4.93210000000000
  Mak135_V%density_ak135( 99) =   4.90690000000000
  Mak135_V%density_ak135(100) =   4.88170000000000
  Mak135_V%density_ak135(101) =   4.85620000000000
  Mak135_V%density_ak135(102) =   4.83070000000000
  Mak135_V%density_ak135(103) =   4.80500000000000
  Mak135_V%density_ak135(104) =   4.77900000000000
  Mak135_V%density_ak135(105) =   4.75280000000000
  Mak135_V%density_ak135(106) =   4.72660000000000
  Mak135_V%density_ak135(107) =   4.70010000000000
  Mak135_V%density_ak135(108) =   4.67350000000000
  Mak135_V%density_ak135(109) =   4.64670000000000
  Mak135_V%density_ak135(110) =   4.61980000000000
  Mak135_V%density_ak135(111) =   4.59260000000000
  Mak135_V%density_ak135(112) =   4.56540000000000
  Mak135_V%density_ak135(113) =   4.51620000000000
  Mak135_V%density_ak135(114) =   4.46500000000000
  Mak135_V%density_ak135(115) =   4.41180000000000
  Mak135_V%density_ak135(116) =   4.35650000000000
  Mak135_V%density_ak135(117) =   4.29860000000000
  Mak135_V%density_ak135(118) =   4.23870000000000
  Mak135_V%density_ak135(119) =   3.92010000000000
  Mak135_V%density_ak135(120) =   3.92060000000000
  Mak135_V%density_ak135(121) =   3.92180000000000
  Mak135_V%density_ak135(122) =   3.92330000000000
  Mak135_V%density_ak135(123) =   3.92730000000000
  Mak135_V%density_ak135(124) =   3.93170000000000
  Mak135_V%density_ak135(125) =   3.50680000000000
  Mak135_V%density_ak135(126) =   3.45770000000000
  Mak135_V%density_ak135(127) =   3.41100000000000
  Mak135_V%density_ak135(128) =   3.36630000000000
  Mak135_V%density_ak135(129) =   3.32430000000000
  Mak135_V%density_ak135(130) =   3.32430000000000
  Mak135_V%density_ak135(131) =   3.37110000000000
  Mak135_V%density_ak135(132) =   3.42680000000000
  Mak135_V%density_ak135(133) =   3.50200000000000
  Mak135_V%density_ak135(134) =   3.50200000000000
  Mak135_V%density_ak135(135) =   3.58010000000000
  Mak135_V%density_ak135(136) =   3.64100000000000
  Mak135_V%density_ak135(137) =   2.92000000000000
  Mak135_V%density_ak135(138) =   2.92000000000000
  Mak135_V%density_ak135(139) =   2.60000000000000
  Mak135_V%density_ak135(140) =   2.60000000000000
  Mak135_V%density_ak135(141) =   2.60000000000000
  Mak135_V%density_ak135(142) =   2.60000000000000
  Mak135_V%density_ak135(143) =   2.60000000000000
  Mak135_V%density_ak135(144) =   2.60000000000000

  Mak135_V%vp_ak135(  1) =   11.2622000000000
  Mak135_V%vp_ak135(  2) =   11.2618000000000
  Mak135_V%vp_ak135(  3) =   11.2606000000000
  Mak135_V%vp_ak135(  4) =   11.2586000000000
  Mak135_V%vp_ak135(  5) =   11.2557000000000
  Mak135_V%vp_ak135(  6) =   11.2521000000000
  Mak135_V%vp_ak135(  7) =   11.2477000000000
  Mak135_V%vp_ak135(  8) =   11.2424000000000
  Mak135_V%vp_ak135(  9) =   11.2364000000000
  Mak135_V%vp_ak135( 10) =   11.2295000000000
  Mak135_V%vp_ak135( 11) =   11.2219000000000
  Mak135_V%vp_ak135( 12) =   11.2134000000000
  Mak135_V%vp_ak135( 13) =   11.2041000000000
  Mak135_V%vp_ak135( 14) =   11.1941000000000
  Mak135_V%vp_ak135( 15) =   11.1832000000000
  Mak135_V%vp_ak135( 16) =   11.1715000000000
  Mak135_V%vp_ak135( 17) =   11.1590000000000
  Mak135_V%vp_ak135( 18) =   11.1457000000000
  Mak135_V%vp_ak135( 19) =   11.1316000000000
  Mak135_V%vp_ak135( 20) =   11.1166000000000
  Mak135_V%vp_ak135( 21) =   11.0983000000000
  Mak135_V%vp_ak135( 22) =   11.0850000000000
  Mak135_V%vp_ak135( 23) =   11.0718000000000
  Mak135_V%vp_ak135( 24) =   11.0585000000000
  Mak135_V%vp_ak135( 25) =   11.0427000000000
  Mak135_V%vp_ak135( 26) =   10.2890000000000
  Mak135_V%vp_ak135( 27) =   10.2854000000000
  Mak135_V%vp_ak135( 28) =   10.2745000000000
  Mak135_V%vp_ak135( 29) =   10.2565000000000
  Mak135_V%vp_ak135( 30) =   10.2329000000000
  Mak135_V%vp_ak135( 31) =   10.2049000000000
  Mak135_V%vp_ak135( 32) =   10.1739000000000
  Mak135_V%vp_ak135( 33) =   10.1415000000000
  Mak135_V%vp_ak135( 34) =   10.1095000000000
  Mak135_V%vp_ak135( 35) =   10.0768000000000
  Mak135_V%vp_ak135( 36) =   10.0439000000000
  Mak135_V%vp_ak135( 37) =   10.0103000000000
  Mak135_V%vp_ak135( 38) =   9.97610000000000
  Mak135_V%vp_ak135( 39) =   9.94100000000000
  Mak135_V%vp_ak135( 40) =   9.90510000000000
  Mak135_V%vp_ak135( 41) =   9.86820000000000
  Mak135_V%vp_ak135( 42) =   9.83040000000000
  Mak135_V%vp_ak135( 43) =   9.79140000000000
  Mak135_V%vp_ak135( 44) =   9.75130000000000
  Mak135_V%vp_ak135( 45) =   9.71000000000000
  Mak135_V%vp_ak135( 46) =   9.66730000000000
  Mak135_V%vp_ak135( 47) =   9.62320000000000
  Mak135_V%vp_ak135( 48) =   9.57770000000000
  Mak135_V%vp_ak135( 49) =   9.53060000000000
  Mak135_V%vp_ak135( 50) =   9.48140000000000
  Mak135_V%vp_ak135( 51) =   9.42970000000000
  Mak135_V%vp_ak135( 52) =   9.37600000000000
  Mak135_V%vp_ak135( 53) =   9.32050000000000
  Mak135_V%vp_ak135( 54) =   9.26340000000000
  Mak135_V%vp_ak135( 55) =   9.20420000000000
  Mak135_V%vp_ak135( 56) =   9.14260000000000
  Mak135_V%vp_ak135( 57) =   9.07920000000000
  Mak135_V%vp_ak135( 58) =   9.01380000000000
  Mak135_V%vp_ak135( 59) =   8.94610000000000
  Mak135_V%vp_ak135( 60) =   8.87610000000000
  Mak135_V%vp_ak135( 61) =   8.80360000000000
  Mak135_V%vp_ak135( 62) =   8.72830000000000
  Mak135_V%vp_ak135( 63) =   8.64960000000000
  Mak135_V%vp_ak135( 64) =   8.56920000000000
  Mak135_V%vp_ak135( 65) =   8.48610000000000
  Mak135_V%vp_ak135( 66) =   8.40010000000000
  Mak135_V%vp_ak135( 67) =   8.31220000000000
  Mak135_V%vp_ak135( 68) =   8.22130000000000
  Mak135_V%vp_ak135( 69) =   8.12830000000000
  Mak135_V%vp_ak135( 70) =   8.03820000000000
  Mak135_V%vp_ak135( 71) =   8.00000000000000
  Mak135_V%vp_ak135( 72) =   13.6601000000000
  Mak135_V%vp_ak135( 73) =   13.6570000000000
  Mak135_V%vp_ak135( 74) =   13.6533000000000
  Mak135_V%vp_ak135( 75) =   13.6498000000000
  Mak135_V%vp_ak135( 76) =   13.6498000000000
  Mak135_V%vp_ak135( 77) =   13.5899000000000
  Mak135_V%vp_ak135( 78) =   13.5311000000000
  Mak135_V%vp_ak135( 79) =   13.4741000000000
  Mak135_V%vp_ak135( 80) =   13.4156000000000
  Mak135_V%vp_ak135( 81) =   13.3584000000000
  Mak135_V%vp_ak135( 82) =   13.3017000000000
  Mak135_V%vp_ak135( 83) =   13.2465000000000
  Mak135_V%vp_ak135( 84) =   13.1895000000000
  Mak135_V%vp_ak135( 85) =   13.1337000000000
  Mak135_V%vp_ak135( 86) =   13.0786000000000
  Mak135_V%vp_ak135( 87) =   13.0226000000000
  Mak135_V%vp_ak135( 88) =   12.9663000000000
  Mak135_V%vp_ak135( 89) =   12.9093000000000
  Mak135_V%vp_ak135( 90) =   12.8524000000000
  Mak135_V%vp_ak135( 91) =   12.7956000000000
  Mak135_V%vp_ak135( 92) =   12.7384000000000
  Mak135_V%vp_ak135( 93) =   12.6807000000000
  Mak135_V%vp_ak135( 94) =   12.6226000000000
  Mak135_V%vp_ak135( 95) =   12.5638000000000
  Mak135_V%vp_ak135( 96) =   12.5030000000000
  Mak135_V%vp_ak135( 97) =   12.4427000000000
  Mak135_V%vp_ak135( 98) =   12.3813000000000
  Mak135_V%vp_ak135( 99) =   12.3181000000000
  Mak135_V%vp_ak135(100) =   12.2558000000000
  Mak135_V%vp_ak135(101) =   12.1912000000000
  Mak135_V%vp_ak135(102) =   12.1247000000000
  Mak135_V%vp_ak135(103) =   12.0571000000000
  Mak135_V%vp_ak135(104) =   11.9891000000000
  Mak135_V%vp_ak135(105) =   11.9208000000000
  Mak135_V%vp_ak135(106) =   11.8491000000000
  Mak135_V%vp_ak135(107) =   11.7768000000000
  Mak135_V%vp_ak135(108) =   11.7020000000000
  Mak135_V%vp_ak135(109) =   11.6265000000000
  Mak135_V%vp_ak135(110) =   11.5493000000000
  Mak135_V%vp_ak135(111) =   11.4704000000000
  Mak135_V%vp_ak135(112) =   11.3897000000000
  Mak135_V%vp_ak135(113) =   11.3068000000000
  Mak135_V%vp_ak135(114) =   11.2228000000000
  Mak135_V%vp_ak135(115) =   11.1355000000000
  Mak135_V%vp_ak135(116) =   11.0553000000000
  Mak135_V%vp_ak135(117) =   10.9222000000000
  Mak135_V%vp_ak135(118) =   10.7909000000000
  Mak135_V%vp_ak135(119) =   10.2000000000000
  Mak135_V%vp_ak135(120) =   10.0320000000000
  Mak135_V%vp_ak135(121) =   9.86400000000000
  Mak135_V%vp_ak135(122) =   9.69620000000000
  Mak135_V%vp_ak135(123) =   9.52800000000000
  Mak135_V%vp_ak135(124) =   9.36010000000000
  Mak135_V%vp_ak135(125) =   9.03020000000000
  Mak135_V%vp_ak135(126) =   8.84760000000000
  Mak135_V%vp_ak135(127) =   8.66500000000000
  Mak135_V%vp_ak135(128) =   8.48220000000000
  Mak135_V%vp_ak135(129) =   8.30070000000000
  Mak135_V%vp_ak135(130) =   8.30070000000000
  Mak135_V%vp_ak135(131) =   8.17500000000000
  Mak135_V%vp_ak135(132) =   8.05050000000000
  Mak135_V%vp_ak135(133) =   8.04500000000000
  Mak135_V%vp_ak135(134) =   8.04000000000000
  Mak135_V%vp_ak135(135) =   8.03790000000000
  Mak135_V%vp_ak135(136) =   8.03550000000000
  Mak135_V%vp_ak135(137) =   6.80000000000000
  Mak135_V%vp_ak135(138) =   6.80000000000000
  Mak135_V%vp_ak135(139) =   5.80000000000000
  Mak135_V%vp_ak135(140) =   5.80000000000000
  Mak135_V%vp_ak135(141) =   5.80000000000000
  Mak135_V%vp_ak135(142) =   5.80000000000000
  Mak135_V%vp_ak135(143) =   5.80000000000000
  Mak135_V%vp_ak135(144) =   5.80000000000000

  Mak135_V%vs_ak135(  1) =   3.66780000000000
  Mak135_V%vs_ak135(  2) =   3.66750000000000
  Mak135_V%vs_ak135(  3) =   3.66670000000000
  Mak135_V%vs_ak135(  4) =   3.66530000000000
  Mak135_V%vs_ak135(  5) =   3.66330000000000
  Mak135_V%vs_ak135(  6) =   3.66080000000000
  Mak135_V%vs_ak135(  7) =   3.65770000000000
  Mak135_V%vs_ak135(  8) =   3.65400000000000
  Mak135_V%vs_ak135(  9) =   3.64980000000000
  Mak135_V%vs_ak135( 10) =   3.64500000000000
  Mak135_V%vs_ak135( 11) =   3.63960000000000
  Mak135_V%vs_ak135( 12) =   3.63370000000000
  Mak135_V%vs_ak135( 13) =   3.62720000000000
  Mak135_V%vs_ak135( 14) =   3.62020000000000
  Mak135_V%vs_ak135( 15) =   3.61260000000000
  Mak135_V%vs_ak135( 16) =   3.60440000000000
  Mak135_V%vs_ak135( 17) =   3.59570000000000
  Mak135_V%vs_ak135( 18) =   3.58640000000000
  Mak135_V%vs_ak135( 19) =   3.57650000000000
  Mak135_V%vs_ak135( 20) =   3.56610000000000
  Mak135_V%vs_ak135( 21) =   3.55510000000000
  Mak135_V%vs_ak135( 22) =   3.54350000000000
  Mak135_V%vs_ak135( 23) =   3.53140000000000
  Mak135_V%vs_ak135( 24) =   3.51870000000000
  Mak135_V%vs_ak135( 25) =   3.50430000000000
  Mak135_V%vs_ak135( 26) =  0.000000000000000
  Mak135_V%vs_ak135( 27) =  0.000000000000000
  Mak135_V%vs_ak135( 28) =  0.000000000000000
  Mak135_V%vs_ak135( 29) =  0.000000000000000
  Mak135_V%vs_ak135( 30) =  0.000000000000000
  Mak135_V%vs_ak135( 31) =  0.000000000000000
  Mak135_V%vs_ak135( 32) =  0.000000000000000
  Mak135_V%vs_ak135( 33) =  0.000000000000000
  Mak135_V%vs_ak135( 34) =  0.000000000000000
  Mak135_V%vs_ak135( 35) =  0.000000000000000
  Mak135_V%vs_ak135( 36) =  0.000000000000000
  Mak135_V%vs_ak135( 37) =  0.000000000000000
  Mak135_V%vs_ak135( 38) =  0.000000000000000
  Mak135_V%vs_ak135( 39) =  0.000000000000000
  Mak135_V%vs_ak135( 40) =  0.000000000000000
  Mak135_V%vs_ak135( 41) =  0.000000000000000
  Mak135_V%vs_ak135( 42) =  0.000000000000000
  Mak135_V%vs_ak135( 43) =  0.000000000000000
  Mak135_V%vs_ak135( 44) =  0.000000000000000
  Mak135_V%vs_ak135( 45) =  0.000000000000000
  Mak135_V%vs_ak135( 46) =  0.000000000000000
  Mak135_V%vs_ak135( 47) =  0.000000000000000
  Mak135_V%vs_ak135( 48) =  0.000000000000000
  Mak135_V%vs_ak135( 49) =  0.000000000000000
  Mak135_V%vs_ak135( 50) =  0.000000000000000
  Mak135_V%vs_ak135( 51) =  0.000000000000000
  Mak135_V%vs_ak135( 52) =  0.000000000000000
  Mak135_V%vs_ak135( 53) =  0.000000000000000
  Mak135_V%vs_ak135( 54) =  0.000000000000000
  Mak135_V%vs_ak135( 55) =  0.000000000000000
  Mak135_V%vs_ak135( 56) =  0.000000000000000
  Mak135_V%vs_ak135( 57) =  0.000000000000000
  Mak135_V%vs_ak135( 58) =  0.000000000000000
  Mak135_V%vs_ak135( 59) =  0.000000000000000
  Mak135_V%vs_ak135( 60) =  0.000000000000000
  Mak135_V%vs_ak135( 61) =  0.000000000000000
  Mak135_V%vs_ak135( 62) =  0.000000000000000
  Mak135_V%vs_ak135( 63) =  0.000000000000000
  Mak135_V%vs_ak135( 64) =  0.000000000000000
  Mak135_V%vs_ak135( 65) =  0.000000000000000
  Mak135_V%vs_ak135( 66) =  0.000000000000000
  Mak135_V%vs_ak135( 67) =  0.000000000000000
  Mak135_V%vs_ak135( 68) =  0.000000000000000
  Mak135_V%vs_ak135( 69) =  0.000000000000000
  Mak135_V%vs_ak135( 70) =  0.000000000000000
  Mak135_V%vs_ak135( 71) =  0.000000000000000
  Mak135_V%vs_ak135( 72) =   7.28170000000000
  Mak135_V%vs_ak135( 73) =   7.27000000000000
  Mak135_V%vs_ak135( 74) =   7.25930000000000
  Mak135_V%vs_ak135( 75) =   7.24850000000000
  Mak135_V%vs_ak135( 76) =   7.24850000000000
  Mak135_V%vs_ak135( 77) =   7.22530000000000
  Mak135_V%vs_ak135( 78) =   7.20310000000000
  Mak135_V%vs_ak135( 79) =   7.18040000000000
  Mak135_V%vs_ak135( 80) =   7.15840000000000
  Mak135_V%vs_ak135( 81) =   7.13680000000000
  Mak135_V%vs_ak135( 82) =   7.11440000000000
  Mak135_V%vs_ak135( 83) =   7.09320000000000
  Mak135_V%vs_ak135( 84) =   7.07220000000000
  Mak135_V%vs_ak135( 85) =   7.05040000000000
  Mak135_V%vs_ak135( 86) =   7.02860000000000
  Mak135_V%vs_ak135( 87) =   7.00690000000000
  Mak135_V%vs_ak135( 88) =   6.98520000000000
  Mak135_V%vs_ak135( 89) =   6.96250000000000
  Mak135_V%vs_ak135( 90) =   6.94160000000000
  Mak135_V%vs_ak135( 91) =   6.91940000000000
  Mak135_V%vs_ak135( 92) =   6.89720000000000
  Mak135_V%vs_ak135( 93) =   6.87430000000000
  Mak135_V%vs_ak135( 94) =   6.85170000000000
  Mak135_V%vs_ak135( 95) =   6.82890000000000
  Mak135_V%vs_ak135( 96) =   6.80560000000000
  Mak135_V%vs_ak135( 97) =   6.78200000000000
  Mak135_V%vs_ak135( 98) =   6.75790000000000
  Mak135_V%vs_ak135( 99) =   6.73230000000000
  Mak135_V%vs_ak135(100) =   6.70700000000000
  Mak135_V%vs_ak135(101) =   6.68130000000000
  Mak135_V%vs_ak135(102) =   6.65540000000000
  Mak135_V%vs_ak135(103) =   6.62850000000000
  Mak135_V%vs_ak135(104) =   6.60090000000000
  Mak135_V%vs_ak135(105) =   6.57280000000000
  Mak135_V%vs_ak135(106) =   6.54310000000000
  Mak135_V%vs_ak135(107) =   6.51310000000000
  Mak135_V%vs_ak135(108) =   6.48220000000000
  Mak135_V%vs_ak135(109) =   6.45140000000000
  Mak135_V%vs_ak135(110) =   6.41820000000000
  Mak135_V%vs_ak135(111) =   6.38600000000000
  Mak135_V%vs_ak135(112) =   6.35190000000000
  Mak135_V%vs_ak135(113) =   6.31640000000000
  Mak135_V%vs_ak135(114) =   6.27990000000000
  Mak135_V%vs_ak135(115) =   6.24240000000000
  Mak135_V%vs_ak135(116) =   6.21000000000000
  Mak135_V%vs_ak135(117) =   6.08980000000000
  Mak135_V%vs_ak135(118) =   5.96070000000000
  Mak135_V%vs_ak135(119) =   5.61040000000000
  Mak135_V%vs_ak135(120) =   5.50470000000000
  Mak135_V%vs_ak135(121) =   5.39890000000000
  Mak135_V%vs_ak135(122) =   5.29220000000000
  Mak135_V%vs_ak135(123) =   5.18640000000000
  Mak135_V%vs_ak135(124) =   5.08060000000000
  Mak135_V%vs_ak135(125) =   4.87020000000000
  Mak135_V%vs_ak135(126) =   4.78320000000000
  Mak135_V%vs_ak135(127) =   4.69640000000000
  Mak135_V%vs_ak135(128) =   4.60940000000000
  Mak135_V%vs_ak135(129) =   4.51840000000000
  Mak135_V%vs_ak135(130) =   4.51840000000000
  Mak135_V%vs_ak135(131) =   4.50900000000000
  Mak135_V%vs_ak135(132) =   4.50000000000000
  Mak135_V%vs_ak135(133) =   4.49000000000000
  Mak135_V%vs_ak135(134) =   4.48000000000000
  Mak135_V%vs_ak135(135) =   4.48560000000000
  Mak135_V%vs_ak135(136) =   4.48390000000000
  Mak135_V%vs_ak135(137) =   3.90000000000000
  Mak135_V%vs_ak135(138) =   3.90000000000000
  Mak135_V%vs_ak135(139) =   3.20000000000000
  Mak135_V%vs_ak135(140) =   3.20000000000000
  Mak135_V%vs_ak135(141) =   3.20000000000000
  Mak135_V%vs_ak135(142) =   3.20000000000000
  Mak135_V%vs_ak135(143) =   3.20000000000000
  Mak135_V%vs_ak135(144) =   3.20000000000000

  if (SUPPRESS_CRUSTAL_MESH) then
    Mak135_V%vp_ak135(137:144) = Mak135_V%vp_ak135(136)
    Mak135_V%vs_ak135(137:144) = Mak135_V%vs_ak135(136)
    Mak135_V%density_ak135(137:144) = Mak135_V%density_ak135(136)
  endif

  Mak135_V%Qkappa_ak135(  1) =   601.270000000000
  Mak135_V%Qkappa_ak135(  2) =   601.320000000000
  Mak135_V%Qkappa_ak135(  3) =   601.460000000000
  Mak135_V%Qkappa_ak135(  4) =   601.700000000000
  Mak135_V%Qkappa_ak135(  5) =   602.050000000000
  Mak135_V%Qkappa_ak135(  6) =   602.490000000000
  Mak135_V%Qkappa_ak135(  7) =   603.040000000000
  Mak135_V%Qkappa_ak135(  8) =   603.690000000000
  Mak135_V%Qkappa_ak135(  9) =   604.440000000000
  Mak135_V%Qkappa_ak135( 10) =   605.280000000000
  Mak135_V%Qkappa_ak135( 11) =   606.260000000000
  Mak135_V%Qkappa_ak135( 12) =   607.310000000000
  Mak135_V%Qkappa_ak135( 13) =   608.480000000000
  Mak135_V%Qkappa_ak135( 14) =   609.740000000000
  Mak135_V%Qkappa_ak135( 15) =   611.120000000000
  Mak135_V%Qkappa_ak135( 16) =   612.620000000000
  Mak135_V%Qkappa_ak135( 17) =   614.210000000000
  Mak135_V%Qkappa_ak135( 18) =   615.930000000000
  Mak135_V%Qkappa_ak135( 19) =   617.780000000000
  Mak135_V%Qkappa_ak135( 20) =   619.710000000000
  Mak135_V%Qkappa_ak135( 21) =   621.500000000000
  Mak135_V%Qkappa_ak135( 22) =   624.080000000000
  Mak135_V%Qkappa_ak135( 23) =   626.870000000000
  Mak135_V%Qkappa_ak135( 24) =   629.890000000000
  Mak135_V%Qkappa_ak135( 25) =   633.260000000000
  Mak135_V%Qkappa_ak135( 26) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 27) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 28) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 29) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 30) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 31) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 32) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 33) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 34) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 35) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 36) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 37) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 38) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 39) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 40) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 41) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 42) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 43) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 44) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 45) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 46) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 47) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 48) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 49) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 50) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 51) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 52) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 53) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 54) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 55) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 56) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 57) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 58) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 59) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 60) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 61) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 62) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 63) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 64) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 65) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 66) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 67) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 68) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 69) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 70) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 71) =   57822.0000000000
  Mak135_V%Qkappa_ak135( 72) =   723.120000000000
  Mak135_V%Qkappa_ak135( 73) =   725.110000000000
  Mak135_V%Qkappa_ak135( 74) =   726.870000000000
  Mak135_V%Qkappa_ak135( 75) =   722.730000000000
  Mak135_V%Qkappa_ak135( 76) =   933.210000000000
  Mak135_V%Qkappa_ak135( 77) =   940.880000000000
  Mak135_V%Qkappa_ak135( 78) =   952.000000000000
  Mak135_V%Qkappa_ak135( 79) =   960.360000000000
  Mak135_V%Qkappa_ak135( 80) =   968.460000000000
  Mak135_V%Qkappa_ak135( 81) =   976.810000000000
  Mak135_V%Qkappa_ak135( 82) =   985.630000000000
  Mak135_V%Qkappa_ak135( 83) =   990.770000000000
  Mak135_V%Qkappa_ak135( 84) =   999.440000000000
  Mak135_V%Qkappa_ak135( 85) =   1008.79000000000
  Mak135_V%Qkappa_ak135( 86) =   1018.38000000000
  Mak135_V%Qkappa_ak135( 87) =   1032.14000000000
  Mak135_V%Qkappa_ak135( 88) =   1042.07000000000
  Mak135_V%Qkappa_ak135( 89) =   1048.09000000000
  Mak135_V%Qkappa_ak135( 90) =   1058.03000000000
  Mak135_V%Qkappa_ak135( 91) =   1064.23000000000
  Mak135_V%Qkappa_ak135( 92) =   1070.38000000000
  Mak135_V%Qkappa_ak135( 93) =   1085.97000000000
  Mak135_V%Qkappa_ak135( 94) =   1097.16000000000
  Mak135_V%Qkappa_ak135( 95) =   1108.58000000000
  Mak135_V%Qkappa_ak135( 96) =   1120.09000000000
  Mak135_V%Qkappa_ak135( 97) =   1127.02000000000
  Mak135_V%Qkappa_ak135( 98) =   1134.01000000000
  Mak135_V%Qkappa_ak135( 99) =   1141.32000000000
  Mak135_V%Qkappa_ak135(100) =   1148.76000000000
  Mak135_V%Qkappa_ak135(101) =   1156.04000000000
  Mak135_V%Qkappa_ak135(102) =   1163.16000000000
  Mak135_V%Qkappa_ak135(103) =   1170.53000000000
  Mak135_V%Qkappa_ak135(104) =   1178.19000000000
  Mak135_V%Qkappa_ak135(105) =   1186.06000000000
  Mak135_V%Qkappa_ak135(106) =   1193.99000000000
  Mak135_V%Qkappa_ak135(107) =   1202.04000000000
  Mak135_V%Qkappa_ak135(108) =   1210.02000000000
  Mak135_V%Qkappa_ak135(109) =   1217.91000000000
  Mak135_V%Qkappa_ak135(110) =   1226.52000000000
  Mak135_V%Qkappa_ak135(111) =   1234.54000000000
  Mak135_V%Qkappa_ak135(112) =   1243.02000000000
  Mak135_V%Qkappa_ak135(113) =   1251.69000000000
  Mak135_V%Qkappa_ak135(114) =   1260.68000000000
  Mak135_V%Qkappa_ak135(115) =   1269.44000000000
  Mak135_V%Qkappa_ak135(116) =   1277.93000000000
  Mak135_V%Qkappa_ak135(117) =   1311.17000000000
  Mak135_V%Qkappa_ak135(118) =   1350.54000000000
  Mak135_V%Qkappa_ak135(119) =   428.690000000000
  Mak135_V%Qkappa_ak135(120) =   425.510000000000
  Mak135_V%Qkappa_ak135(121) =   422.550000000000
  Mak135_V%Qkappa_ak135(122) =   419.940000000000
  Mak135_V%Qkappa_ak135(123) =   417.320000000000
  Mak135_V%Qkappa_ak135(124) =   413.660000000000
  Mak135_V%Qkappa_ak135(125) =   377.930000000000
  Mak135_V%Qkappa_ak135(126) =   366.340000000000
  Mak135_V%Qkappa_ak135(127) =   355.850000000000
  Mak135_V%Qkappa_ak135(128) =   346.370000000000
  Mak135_V%Qkappa_ak135(129) =   338.470000000000
  Mak135_V%Qkappa_ak135(130) =   200.970000000000
  Mak135_V%Qkappa_ak135(131) =   188.720000000000
  Mak135_V%Qkappa_ak135(132) =   182.570000000000
  Mak135_V%Qkappa_ak135(133) =   182.030000000000
  Mak135_V%Qkappa_ak135(134) =   1008.71000000000
  Mak135_V%Qkappa_ak135(135) =   972.770000000000
  Mak135_V%Qkappa_ak135(136) =   950.500000000000
  Mak135_V%Qkappa_ak135(137) =   1368.02000000000
  Mak135_V%Qkappa_ak135(138) =   1368.02000000000
  Mak135_V%Qkappa_ak135(139) =   1478.30000000000
  Mak135_V%Qkappa_ak135(140) =   1478.30000000000
  Mak135_V%Qkappa_ak135(141) =   1478.30000000000
  Mak135_V%Qkappa_ak135(142) =   1478.30000000000
  Mak135_V%Qkappa_ak135(143) =   1478.30000000000
  Mak135_V%Qkappa_ak135(144) =   1478.30000000000

  Mak135_V%Qmu_ak135(  1) =   85.0300000000000
  Mak135_V%Qmu_ak135(  2) =   85.0300000000000
  Mak135_V%Qmu_ak135(  3) =   85.0300000000000
  Mak135_V%Qmu_ak135(  4) =   85.0300000000000
  Mak135_V%Qmu_ak135(  5) =   85.0300000000000
  Mak135_V%Qmu_ak135(  6) =   85.0300000000000
  Mak135_V%Qmu_ak135(  7) =   85.0300000000000
  Mak135_V%Qmu_ak135(  8) =   85.0300000000000
  Mak135_V%Qmu_ak135(  9) =   85.0300000000000
  Mak135_V%Qmu_ak135( 10) =   85.0300000000000
  Mak135_V%Qmu_ak135( 11) =   85.0300000000000
  Mak135_V%Qmu_ak135( 12) =   85.0300000000000
  Mak135_V%Qmu_ak135( 13) =   85.0300000000000
  Mak135_V%Qmu_ak135( 14) =   85.0300000000000
  Mak135_V%Qmu_ak135( 15) =   85.0300000000000
  Mak135_V%Qmu_ak135( 16) =   85.0300000000000
  Mak135_V%Qmu_ak135( 17) =   85.0300000000000
  Mak135_V%Qmu_ak135( 18) =   85.0300000000000
  Mak135_V%Qmu_ak135( 19) =   85.0300000000000
  Mak135_V%Qmu_ak135( 20) =   85.0300000000000
  Mak135_V%Qmu_ak135( 21) =   85.0300000000000
  Mak135_V%Qmu_ak135( 22) =   85.0300000000000
  Mak135_V%Qmu_ak135( 23) =   85.0300000000000
  Mak135_V%Qmu_ak135( 24) =   85.0300000000000
  Mak135_V%Qmu_ak135( 25) =   85.0300000000000
  Mak135_V%Qmu_ak135( 26) =  0.000000000000000
  Mak135_V%Qmu_ak135( 27) =  0.000000000000000
  Mak135_V%Qmu_ak135( 28) =  0.000000000000000
  Mak135_V%Qmu_ak135( 29) =  0.000000000000000
  Mak135_V%Qmu_ak135( 30) =  0.000000000000000
  Mak135_V%Qmu_ak135( 31) =  0.000000000000000
  Mak135_V%Qmu_ak135( 32) =  0.000000000000000
  Mak135_V%Qmu_ak135( 33) =  0.000000000000000
  Mak135_V%Qmu_ak135( 34) =  0.000000000000000
  Mak135_V%Qmu_ak135( 35) =  0.000000000000000
  Mak135_V%Qmu_ak135( 36) =  0.000000000000000
  Mak135_V%Qmu_ak135( 37) =  0.000000000000000
  Mak135_V%Qmu_ak135( 38) =  0.000000000000000
  Mak135_V%Qmu_ak135( 39) =  0.000000000000000
  Mak135_V%Qmu_ak135( 40) =  0.000000000000000
  Mak135_V%Qmu_ak135( 41) =  0.000000000000000
  Mak135_V%Qmu_ak135( 42) =  0.000000000000000
  Mak135_V%Qmu_ak135( 43) =  0.000000000000000
  Mak135_V%Qmu_ak135( 44) =  0.000000000000000
  Mak135_V%Qmu_ak135( 45) =  0.000000000000000
  Mak135_V%Qmu_ak135( 46) =  0.000000000000000
  Mak135_V%Qmu_ak135( 47) =  0.000000000000000
  Mak135_V%Qmu_ak135( 48) =  0.000000000000000
  Mak135_V%Qmu_ak135( 49) =  0.000000000000000
  Mak135_V%Qmu_ak135( 50) =  0.000000000000000
  Mak135_V%Qmu_ak135( 51) =  0.000000000000000
  Mak135_V%Qmu_ak135( 52) =  0.000000000000000
  Mak135_V%Qmu_ak135( 53) =  0.000000000000000
  Mak135_V%Qmu_ak135( 54) =  0.000000000000000
  Mak135_V%Qmu_ak135( 55) =  0.000000000000000
  Mak135_V%Qmu_ak135( 56) =  0.000000000000000
  Mak135_V%Qmu_ak135( 57) =  0.000000000000000
  Mak135_V%Qmu_ak135( 58) =  0.000000000000000
  Mak135_V%Qmu_ak135( 59) =  0.000000000000000
  Mak135_V%Qmu_ak135( 60) =  0.000000000000000
  Mak135_V%Qmu_ak135( 61) =  0.000000000000000
  Mak135_V%Qmu_ak135( 62) =  0.000000000000000
  Mak135_V%Qmu_ak135( 63) =  0.000000000000000
  Mak135_V%Qmu_ak135( 64) =  0.000000000000000
  Mak135_V%Qmu_ak135( 65) =  0.000000000000000
  Mak135_V%Qmu_ak135( 66) =  0.000000000000000
  Mak135_V%Qmu_ak135( 67) =  0.000000000000000
  Mak135_V%Qmu_ak135( 68) =  0.000000000000000
  Mak135_V%Qmu_ak135( 69) =  0.000000000000000
  Mak135_V%Qmu_ak135( 70) =  0.000000000000000
  Mak135_V%Qmu_ak135( 71) =  0.000000000000000
  Mak135_V%Qmu_ak135( 72) =   273.970000000000
  Mak135_V%Qmu_ak135( 73) =   273.970000000000
  Mak135_V%Qmu_ak135( 74) =   273.970000000000
  Mak135_V%Qmu_ak135( 75) =   271.740000000000
  Mak135_V%Qmu_ak135( 76) =   350.880000000000
  Mak135_V%Qmu_ak135( 77) =   354.610000000000
  Mak135_V%Qmu_ak135( 78) =   359.710000000000
  Mak135_V%Qmu_ak135( 79) =   363.640000000000
  Mak135_V%Qmu_ak135( 80) =   367.650000000000
  Mak135_V%Qmu_ak135( 81) =   371.750000000000
  Mak135_V%Qmu_ak135( 82) =   375.940000000000
  Mak135_V%Qmu_ak135( 83) =   378.790000000000
  Mak135_V%Qmu_ak135( 84) =   383.140000000000
  Mak135_V%Qmu_ak135( 85) =   387.600000000000
  Mak135_V%Qmu_ak135( 86) =   392.160000000000
  Mak135_V%Qmu_ak135( 87) =   398.410000000000
  Mak135_V%Qmu_ak135( 88) =   403.230000000000
  Mak135_V%Qmu_ak135( 89) =   406.500000000000
  Mak135_V%Qmu_ak135( 90) =   411.520000000000
  Mak135_V%Qmu_ak135( 91) =   414.940000000000
  Mak135_V%Qmu_ak135( 92) =   418.410000000000
  Mak135_V%Qmu_ak135( 93) =   425.530000000000
  Mak135_V%Qmu_ak135( 94) =   431.030000000000
  Mak135_V%Qmu_ak135( 95) =   436.680000000000
  Mak135_V%Qmu_ak135( 96) =   442.480000000000
  Mak135_V%Qmu_ak135( 97) =   446.430000000000
  Mak135_V%Qmu_ak135( 98) =   450.450000000000
  Mak135_V%Qmu_ak135( 99) =   454.550000000000
  Mak135_V%Qmu_ak135(100) =   458.720000000000
  Mak135_V%Qmu_ak135(101) =   462.960000000000
  Mak135_V%Qmu_ak135(102) =   467.290000000000
  Mak135_V%Qmu_ak135(103) =   471.700000000000
  Mak135_V%Qmu_ak135(104) =   476.190000000000
  Mak135_V%Qmu_ak135(105) =   480.770000000000
  Mak135_V%Qmu_ak135(106) =   485.440000000000
  Mak135_V%Qmu_ak135(107) =   490.200000000000
  Mak135_V%Qmu_ak135(108) =   495.050000000000
  Mak135_V%Qmu_ak135(109) =   500.000000000000
  Mak135_V%Qmu_ak135(110) =   505.050000000000
  Mak135_V%Qmu_ak135(111) =   510.200000000000
  Mak135_V%Qmu_ak135(112) =   515.460000000000
  Mak135_V%Qmu_ak135(113) =   520.830000000000
  Mak135_V%Qmu_ak135(114) =   526.320000000000
  Mak135_V%Qmu_ak135(115) =   531.910000000000
  Mak135_V%Qmu_ak135(116) =   537.630000000000
  Mak135_V%Qmu_ak135(117) =   543.480000000000
  Mak135_V%Qmu_ak135(118) =   549.450000000000
  Mak135_V%Qmu_ak135(119) =   172.930000000000
  Mak135_V%Qmu_ak135(120) =   170.820000000000
  Mak135_V%Qmu_ak135(121) =   168.780000000000
  Mak135_V%Qmu_ak135(122) =   166.800000000000
  Mak135_V%Qmu_ak135(123) =   164.870000000000
  Mak135_V%Qmu_ak135(124) =   162.500000000000
  Mak135_V%Qmu_ak135(125) =   146.570000000000
  Mak135_V%Qmu_ak135(126) =   142.760000000000
  Mak135_V%Qmu_ak135(127) =   139.380000000000
  Mak135_V%Qmu_ak135(128) =   136.380000000000
  Mak135_V%Qmu_ak135(129) =   133.720000000000
  Mak135_V%Qmu_ak135(130) =   79.4000000000000
  Mak135_V%Qmu_ak135(131) =   76.5500000000000
  Mak135_V%Qmu_ak135(132) =   76.0600000000000
  Mak135_V%Qmu_ak135(133) =   75.6000000000000
  Mak135_V%Qmu_ak135(134) =   417.590000000000
  Mak135_V%Qmu_ak135(135) =   403.930000000000
  Mak135_V%Qmu_ak135(136) =   394.620000000000
  Mak135_V%Qmu_ak135(137) =   599.990000000000
  Mak135_V%Qmu_ak135(138) =   599.990000000000
  Mak135_V%Qmu_ak135(139) =   599.990000000000
  Mak135_V%Qmu_ak135(140) =   599.990000000000
  Mak135_V%Qmu_ak135(141) =   599.990000000000
  Mak135_V%Qmu_ak135(142) =   599.990000000000
  Mak135_V%Qmu_ak135(143) =   599.990000000000
  Mak135_V%Qmu_ak135(144) =   599.990000000000

! strip the crust and replace it by mantle
  if(USE_EXTERNAL_CRUSTAL_MODEL) then
    do i=NR_AK135-8,NR_AK135
      Mak135_V%density_ak135(i) = Mak135_V%density_ak135(NR_AK135-9)
      Mak135_V%vp_ak135(i) = Mak135_V%vp_ak135(NR_AK135-9)
      Mak135_V%vs_ak135(i) = Mak135_V%vs_ak135(NR_AK135-9)
      Mak135_V%Qkappa_ak135(i) = Mak135_V%Qkappa_ak135(NR_AK135-9)
      Mak135_V%Qmu_ak135(i) = Mak135_V%Qmu_ak135(NR_AK135-9)
    enddo
  endif

  end subroutine define_model_ak135

