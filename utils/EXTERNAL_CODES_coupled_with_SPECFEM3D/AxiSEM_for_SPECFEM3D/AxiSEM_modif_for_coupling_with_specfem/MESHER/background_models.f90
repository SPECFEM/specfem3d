!
!    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
!                    Simon Stahler, Kasra Hosseini, Stefanie Hempel
!
!    This file is part of AxiSEM.
!    It is distributed from the webpage <http://www.axisem.info>
!
!    AxiSEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    AxiSEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with AxiSEM.  If not, see <http://www.gnu.org/licenses/>.
!

!-----------------------------------------------------------------------------
!> This module is identical in the mesher and solver.
!! Function "velocity" retrieves the velocity (or density) of the given
!! background model (see different cases and respective routines below)
!! for radius r0 given its subdomain identification number idom
!! (determined respectively in mesher and solver).
!!
!! When adding new background models, one needs to define them both in terms
!! of subdomains for this module (e.g. as a polynomial), and in terms of the
!! discontinuities and their above/below elastic values for the mesher only
!! (see module model_discontinuities).
module background_models
  use global_parameters
  use interpolation
  implicit none

  public :: velocity, model_is_ani, model_is_anelastic !, arbitr_sub_solar_arr
  public :: read_ext_model, get_ext_disc, override_ext_q
  private
  character(len=6),           save  :: override_ext_q
  integer                   , save  :: nlayer = -1
  real(kind=dp), allocatable, save  :: vpv_layer(:)
  real(kind=dp), allocatable, save  :: vsv_layer(:)
  real(kind=dp), allocatable, save  :: vph_layer(:)
  real(kind=dp), allocatable, save  :: vsh_layer(:)
  real(kind=dp), allocatable, save  :: qka_layer(:)
  real(kind=dp), allocatable, save  :: qmu_layer(:)
  real(kind=dp), allocatable, save  :: rho_layer(:)
  real(kind=dp), allocatable, save  :: radius_layer(:)
  logical, save                     :: ext_model_is_ani, ext_model_is_anelastic
  type(interpolation_data), allocatable, save :: interp_vpv(:), interp_vsv(:), &
                                                 interp_vph(:), interp_vsh(:), &
                                                 interp_qka(:), interp_qmu(:), &
                                                 interp_rho(:)

contains

!-----------------------------------------------------------------------------
!> Wrapper function to call velocities upon different background models
!! for a given radius r0 [m], parameter type param (rho,vs,vp) and idom
real(kind=dp)  function velocity(r0, param, idom, bkgrdmodel2, lfbkgrdmodel2)

  real(kind=dp)   , intent(in)   :: r0
  integer, intent(in)            :: idom
  character(len=100), intent(in) :: bkgrdmodel2
  integer, intent(in)            :: lfbkgrdmodel2
  character(len=3)               :: param !rho, vs,vp

  select case(bkgrdmodel2(1:lfbkgrdmodel2))
     case('ak135')
        velocity = ak135(r0, param, idom)
     case('ak135f')
        velocity = ak135f(r0, param, idom)
     case('prem_iso')
        velocity = prem_sub(r0, param, idom)
     case('prem_ani')
        velocity = prem_ani_sub(r0, param, idom)
     case('prem_iso_solid')
        velocity = prem_solid_sub(r0, param, idom)
     case('prem_iso_light')
        velocity = prem_light_sub(r0, param, idom)
     case('prem_ani_light')
        velocity = prem_light_ani_sub(r0, param, idom)
     case('prem_iso_onecrust')
        velocity = prem_onecrust_sub(r0, param, idom)
     case('prem_ani_onecrust')
        velocity = prem_onecrust_ani_sub(r0, param, idom)
     case('prem_iso_solid_light')
        velocity = prem_solid_light_sub(r0, param, idom)
     case('iasp91')
        velocity = iasp91_sub(r0, param, idom)
!     case('solar')
!        velocity = arbitr_sub_solar(r0, param, idom, bkgrdmodel2)
     case('prem_iso_custom')
        velocity=prem_sub_vm(r0, param, idom)
     case('external')
        velocity = arbitr_sub_solar(r0, param, idom)
     case default
        write(6,*) 'Unknown background model: ', bkgrdmodel2
        stop
  end select

end function velocity
!=============================================================================

!-----------------------------------------------------------------------------
!> returns true if the model is radially anisotrpoic
logical function model_is_ani(bkgrdmodel2)

  character(len=100), intent(in) :: bkgrdmodel2

  select case(trim(bkgrdmodel2))
  case('prem_ani')
    model_is_ani = .true.
  case('prem_ani_onecrust')
    model_is_ani = .true.
  case('prem_ani_light')
    model_is_ani = .true.
  case('external')
    model_is_ani = ext_model_is_ani
  case default
    model_is_ani = .false.
  end select

end function model_is_ani
!=============================================================================

!-----------------------------------------------------------------------------
!> returns true if the model is radially anisotrpoic
logical function model_is_anelastic(bkgrdmodel2)

  character(len=100), intent(in) :: bkgrdmodel2

  select case(trim(bkgrdmodel2))
  case('ak135')
    model_is_anelastic = .true.
  case('ak135f')
    model_is_anelastic = .true.
  case('prem_iso')
    model_is_anelastic = .true.
  case('prem_ani')
    model_is_anelastic = .true.
  case('prem_ani_onecrust')
    model_is_anelastic = .true.
  case('prem_ani_light')
    model_is_anelastic = .true.
  case('prem_iso_light')
    model_is_anelastic = .true.
  case('prem_iso_custom')
    model_is_anelastic = .true.
  case('iasp91')
    model_is_anelastic = .true.
  case('external')
    model_is_anelastic = ext_model_is_anelastic
  case default
    model_is_anelastic = .false.
  end select

end function model_is_anelastic
!=============================================================================

!-----------------------------------------------------------------------------
!> from Montagner and Kennett (1995)
!! interpolated between discontinuities using matlab's polyfit, use radii!!!
!! use routine axisem_ak135_fitting.m and ak135/iasp whatever as nd-files
real(kind=dp) function ak135f(r0, param, idom)

    real(kind=dp)   , intent(in) :: r0
    integer, intent(in)          :: idom
    real(kind=dp)                :: r,x_ak
    real(kind=dp)                :: ro_ak, vp_ak, vs_ak
    real(kind=dp)                :: Qmu_ak, Qka_ak
    character(len=3), intent(in) :: param !rho, vs,vp

    r =r0 / 1000.
    x_ak = r / 6371. ! normalized


    select case(idom)
    case(12) ! depth: 6371.000000
        ro_ak = 13.012216 - 0.001140 * x_ak - 8.445249 * x_ak**2
        vp_ak = 11.264846 - 0.103927 * x_ak - 5.687562 * x_ak**2
        vs_ak = 3.667675 + 0.005479 * x_ak - 4.482579 * x_ak**2
        Qmu_ak = 85.03
        Qka_ak = 595.258179 + 166.678237 * x_ak
    case(11) ! depth: 5153.500000
        ro_ak = 12.277066 + 1.075439 * x_ak - 9.829445 * x_ak**2
        vp_ak = 10.134921 + 3.305589 * x_ak - 13.242147 * x_ak**2
        vs_ak = 0.0
        Qmu_ak = 0.0
        Qka_ak = 57822.0
    case(10) ! depth: 2891.500000
        ro_ak = 6.090533  + 2.034521 * x_ak - 4.792620 * x_ak**2
        vp_ak = 13.185166 + 2.121947 * x_ak - 2.292893 * x_ak**2
        vs_ak = 8.483334  - 2.972926 * x_ak + 1.414776 * x_ak**2
        Qmu_ak = 320.279248 -83.970345 * x_ak
        Qka_ak = 719.424861 + 9.016892 * x_ak
    case(9) ! depth: 2740.000000
        ro_ak = 5.880628 + 0.806448 * x_ak - 2.809695 * x_ak**2
        vp_ak = 15.048453 + 1.120394 * x_ak - 6.384131 * x_ak**2
        vs_ak = 7.454967 + 1.448974 * x_ak - 3.232585 * x_ak**2
        Qmu_ak = 6.940579 + 597.788236 * x_ak
        Qka_ak = 301.363871 + 1112.990919 * x_ak
    case(8) ! depth: 760.000000
        ro_ak = -1.851510 + 21.347947 * x_ak - 16.235856 * x_ak**2
        vp_ak = 37.426776 - 42.812610 * x_ak + 14.612271 * x_ak**2
        vs_ak = -36.840846 + 112.512879 * x_ak - 72.249561 * x_ak**2
        Qmu_ak = -125.610200 + 753.052200 * x_ak
        Qka_ak = -2797.238767 + 4625.983100 * x_ak
    case(7) ! depth: 660.000000
        ro_ak = 11.125709 - 16.016803 * x_ak + 8.900728 * x_ak**2
        vp_ak = 29.313708 - 21.244664 * x_ak - 0.086978 * x_ak**2
        vs_ak = 17.598241 - 13.243125 * x_ak - 0.144963 * x_ak**2
        Qmu_ak = 407.826350 - 262.048331 * x_ak
        Qka_ak = 762.534265 - 372.539674 * x_ak
    case(6) ! depth: 410.000000
        ro_ak = 25.946777 - 41.561564 * x_ak + 18.787205 * x_ak**2
        vp_ak = 32.879415 - 27.659605 * x_ak + 2.319408 * x_ak**2
        vs_ak =  6.725880 +  6.913882 * x_ak - 9.509573 * x_ak**2
        Qmu_ak = 528.635760 - 408.763360 * x_ak
        Qka_ak = 1555.736580 - 1260.056380 * x_ak
    case(5) ! depth: 210.000000
        ro_ak = 80.939817 - 166.517962 * x_ak + 89.196989 * x_ak**2
        vp_ak = 36.839365 - 41.141558  * x_ak + 12.026560 * x_ak**2
        vs_ak = 9.581677  - 9.112575   * x_ak +  4.008853 * x_ak**2
        Qmu_ak = 307.648222 - 236.434889 * x_ak
        Qka_ak = 1459.535556 - 1302.515556 * x_ak
    case(4) ! depth: 120.000000
        ro_ak = -8.325080 + 11.977480 * x_ak
        vp_ak =  8.910012 - 0.876012 * x_ak
        vs_ak =  6.062750 - 1.592750 * x_ak
        Qmu_ak = 147.946500 - 73.266500 * x_ak
        Qka_ak = 266.958500 - 86.008500 * x_ak
    case(3) ! depth: 80.000000
        ro_ak = 199.022923 - 408.226152 * x_ak + 212.892136 * x_ak**2
        vp_ak = -16.800554 +  50.525274 * x_ak - 25.691438 * x_ak**2
        vs_ak = -137.314994 + 285.398038 * x_ak - 143.603107 * x_ak**2
        Qmu_ak = 2747.697307 -2359.725421 * x_ak
        Qka_ak = 6930.368496 -5997.284866 * x_ak
    case(2) ! depth: 18.000000
        ro_ak = 2.920000
        vp_ak = 6.800000
        vs_ak = 3.900000
        Qmu_ak = 599.990000
        Qka_ak = 1368.020000
    case(1) ! depth: 10.000000
        ro_ak = 2.600000
        vp_ak = 5.800000
        vs_ak = 3.200000
        Qmu_ak = 599.990000
        Qka_ak = 1478.300000
    end select

    if (param=='rho') then
       ak135f = ro_ak * 1000.
    else if (param=='v_p') then
       ak135f = vp_ak * 1000.
    else if (param=='v_s') then
       ak135f = vs_ak * 1000.
    else if (param=='vpv') then
       ak135f = vp_ak * 1000.
    else if (param=='vsv') then
       ak135f = vs_ak * 1000.
    else if (param=='vph') then
       ak135f = vp_ak * 1000.
    else if (param=='vsh') then
       ak135f = vs_ak * 1000.
    else if (param=='eta') then
       ak135f = 1.
    else if (param=='Qmu') then
       ak135f = Qmu_ak
    else if (param=='Qka') then
       ak135f = Qka_ak
    else
       write(6,*)'ERROR IN AK135 FUNCTION:', param, 'NOT AN OPTION'
       stop
    endif

end function ak135f
!=============================================================================

!-----------------------------------------------------------------------------
!> from Kennett, Engdahl and Buland, 1995
!! interpolated between discontinuities using matlab's polyfit, use radii!!!
!! use routine axisem_ak135_fitting.m and ak135/iasp whatever as nd-files
real(kind=dp) function ak135(r0, param, idom)

  real(kind=dp)   , intent(in) :: r0
  integer, intent(in)          :: idom
  real(kind=dp)                :: r,x_ak
  real(kind=dp)                :: ro_ak, vp_ak, vs_ak
  real(kind=dp)                :: Qmu, Qkappa
  character(len=3), intent(in) :: param !rho, vs,vp

  r = r0 / 1000.
  x_ak = r / 6371. ! normalized


  if(idom==1) then
     ro_ak = 2.72
     vp_ak = 5.8
     vs_ak = 3.46
     Qmu = 600.0
     Qkappa = 57827.0
  else if(idom==2) then
     ro_ak = 2.92
     vp_ak = 6.5
     vs_ak = 3.85
     Qmu = 600.0
     Qkappa = 57827.0
  else if(idom==3) then
  ! moho -> 210
     ro_ak =  7.1576 - 3.859  * x_ak
     vp_ak = 17.4734 - 9.5332 * x_ak
     vs_ak =  5.8556 - 1.3825 * x_ak
     Qmu = 600.0
     Qkappa = 57827.0
  else if(idom==4) then
  ! 210 -> 410
     ro_ak =  7.1594 -  3.8608 * x_ak
     vp_ak = 30.7877 - 23.2542 * x_ak
     vs_ak = 15.2181 - 11.0601 * x_ak
     Qmu = 143.0
     Qkappa = 57827.0
  else if(idom==5) then
  ! 410 -> 660
     ro_ak = 11.1204 -  7.8713 * x_ak
     vp_ak = 29.389  - 21.4066 * x_ak
     vs_ak = 17.7173 - 13.5065 * x_ak
     Qmu = 143.0
     Qkappa = 57827.0
  else if(idom==6) then
  ! 660 -> D''
     ro_ak =  6.8294 - 1.7227  * x_ak -  1.1064 * x_ak**2 -  0.034409 * x_ak**3
     vp_ak = 26.8598 - 48.9644 * x_ak + 63.7326 * x_ak**2 - 32.4155   * x_ak**3
     vs_ak = 18.0019 - 43.6346 * x_ak + 60.4205 * x_ak**2 - 29.689    * x_ak**3
     Qmu = 312.0
     Qkappa = 57827.0
  else if(idom==7) then
  ! D'' -> CMB
     ro_ak = -65.8145 + 386.221  * x_ak - 691.6551 *x_ak**2 + 409.6742 * x_ak**3
     vp_ak =   3.4872 + 55.1872  * x_ak -  99.0089 *x_ak**2 +  58.7141 * x_ak**3
     vs_ak = -22.9553 + 164.0287 * x_ak - 294.2766 *x_ak**2 + 174.5113 * x_ak**3
     Qmu = 312.0
     Qkappa = 57827.0
  else if(idom==8) then
  ! CMB -> ICB
     ro_ak = 12.592  - 1.778  * x_ak - 1.6964 * x_ak**2 -  7.3524 * x_ak**3
     vp_ak = 10.7738 - 2.4831 * x_ak + 3.2584 * x_ak**2 - 14.9171 * x_ak**3
     vs_ak = 0
     Qmu = 0.0
     Qkappa = 57827.0
  else if(idom==9) then
  ! inner core
     ro_ak = 13.0122 - 0.0011863 *x_ak - 8.4449 * x_ak**2
     vp_ak = 11.2641 - 0.090247  *x_ak - 5.7431 * x_ak**2
     vs_ak =  3.6677 + 0.0049932 *x_ak - 4.4808 * x_ak**2
     Qmu = 84.6
     Qkappa = 1327.7
  endif

  if (param=='rho') then
     ak135 = ro_ak * 1000.
  else if (param=='v_p') then
     ak135 = vp_ak * 1000.
  else if (param=='v_s') then
     ak135 = vs_ak * 1000.
  else if (param=='vpv') then
     ak135 = vp_ak * 1000.
  else if (param=='vsv') then
     ak135 = vs_ak * 1000.
  else if (param=='vph') then
     ak135 = vp_ak * 1000.
  else if (param=='vsh') then
     ak135 = vs_ak * 1000.
  else if (param=='eta') then
     ak135 = 1.
  else if (param=='Qmu') then
     ak135 = Qmu
  else if (param=='Qka') then
     ak135 = Qkappa
  else
     write(6,*)'ERROR IN AK135 FUNCTION:', param, 'NOT AN OPTION'
     stop
  endif

end function ak135
!=============================================================================

!-----------------------------------------------------------------------------
!> isotropic prem model in terms of domains separated by discontinuities
real(kind=dp) function prem_sub(r0, param, idom)

  real(kind=dp)   , intent(in) :: r0
  integer, intent(in)          :: idom
  real(kind=dp)                :: r,x_prem
  real(kind=dp)                :: ro_prem,vp_prem,vs_prem
  real(kind=dp)                :: Qmu, Qkappa
  character(len=3), intent(in) :: param !rho, vs,vp

  r = r0 / 1000.

  x_prem = r / 6371.     ! Radius (normalized to x(surface)=1 )

  if(idom==1)then        ! upper crustal layer
     ro_prem = 2.6
     vp_prem = 5.8
     vs_prem = 3.2
     Qmu = 600.0
     Qkappa = 57827.0
  else if(idom==2)then
     ro_prem = 2.9                       ! lower crustal layer
     vp_prem = 6.8
     vs_prem = 3.9
     Qmu = 600.0
     Qkappa = 57827.0
  else if(idom==3)then
     ro_prem = 2.691  + 0.6924 * x_prem             ! upper mantle
     vp_prem = 4.1875 + 3.9382 * x_prem
     vs_prem = 2.1519 + 2.3481 * x_prem
     Qmu = 600.0
     Qkappa = 57827.0
  else if(idom==4)then
     ro_prem=  7.1089 -  3.8045 * x_prem
     vp_prem= 20.3926 - 12.2569 * x_prem
     vs_prem=  8.9496 -  4.4597 * x_prem
     Qmu = 143.0
     Qkappa = 57827.0
  else if(idom==5)then
     ro_prem = 11.2494 -  8.0298 * x_prem
     vp_prem = 39.7027 - 32.6166 * x_prem
     vs_prem = 22.3512 - 18.5856 * x_prem
     Qmu = 143.0
     Qkappa = 57827.0
  else if(idom==6)then
     ro_prem =  5.3197 - 1.4836 * x_prem
     vp_prem = 19.0957 - 9.8672 * x_prem
     vs_prem =  9.9839 - 4.9324 * x_prem
     Qmu = 143.0
     Qkappa = 57827.0
  else if(idom==7)then   !lower mantle
     ro_prem =  7.9565 -  6.4761 * x_prem + 5.5283 * x_prem**2 - 3.0807 * x_prem**3
     vp_prem = 29.2766 - 23.6027 * x_prem + 5.5242 * x_prem**2 - 2.5514 * x_prem**3
     vs_prem = 22.3459 - 17.2473 * x_prem - 2.0834 * x_prem**2 + 0.9783 * x_prem**3
     Qmu = 312.0
     Qkappa = 57827.0
  else if(idom==8)then
     ro_prem =  7.9565 -  6.4761 * x_prem +  5.5283 * x_prem**2 -  3.0807 * x_prem**3
     vp_prem = 24.9520 - 40.4673 * x_prem + 51.4832 * x_prem**2 - 26.6419 * x_prem**3
     vs_prem = 11.1671 - 13.7818 * x_prem + 17.4575 * x_prem**2 -  9.2777 * x_prem**3
     Qmu = 312.0
     Qkappa = 57827.0
  else if(idom==9)then
     ro_prem =  7.9565 - 6.4761 * x_prem + 5.5283 * x_prem**2 - 3.0807 * x_prem**3
     vp_prem = 15.3891 - 5.3181 * x_prem + 5.5242 * x_prem**2 - 2.5514 * x_prem**3
     vs_prem =  6.9254 + 1.4672 * x_prem - 2.0834 * x_prem**2 + 0.9783 * x_prem**3
     Qmu = 312.0
     Qkappa = 57827.0
  else if(idom==10)then  ! outer core
     ro_prem = 12.5815 - 1.2638 * x_prem - 3.6426 * x_prem**2 -  5.5281  * x_prem**3
     vp_prem = 11.0487 - 4.0362 * x_prem + 4.8023 * x_prem**2 - 13.5732  * x_prem**3
     vs_prem = 0.0
     Qmu = 0.0
     Qkappa = 57827.0
  else if(idom==11)then                        ! inner core
     ro_prem = 13.0885 - 8.8381 * x_prem**2
     vp_prem = 11.2622 - 6.3640 * x_prem**2
     vs_prem =  3.6678 - 4.4475 * x_prem**2
     Qmu = 84.6
     Qkappa = 1327.7
  endif

  if (param=='rho') then
     prem_sub = ro_prem * 1000.
  else if (param=='v_p') then
     prem_sub = vp_prem * 1000.
  else if (param=='v_s') then
     prem_sub = vs_prem * 1000.
  else if (param=='vpv') then
     prem_sub = vp_prem * 1000.
  else if (param=='vsv') then
     prem_sub = vs_prem * 1000.
  else if (param=='vph') then
     prem_sub = vp_prem * 1000.
  else if (param=='vsh') then
     prem_sub = vs_prem * 1000.
  else if (param=='eta') then
     prem_sub = 1.
  else if (param=='Qmu') then
     prem_sub = Qmu
  else if (param=='Qka') then
     prem_sub = Qkappa
  else
     write(6,*)'ERROR IN PREM_SUB FUNCTION:',param,'NOT AN OPTION'
     stop
  endif

end function prem_sub
!=============================================================================

!-----------------------------------------------------------------------------
!> prem model in terms of domains separated by discontinuities
real(kind=dp) function prem_ani_sub(r0, param, idom)

  real(kind=dp)   , intent(in) :: r0
  integer, intent(in)          :: idom
  real(kind=dp)                :: r,x_prem
  real(kind=dp)                :: ro_prem, vpv_prem, vsv_prem, vph_prem
  real(kind=dp)                :: vsh_prem, eta_aniso, Qmu, Qkappa
  character(len=3), intent(in) :: param !rho, vs,vp

  r = r0 / 1000.

  x_prem = r / 6371.     ! Radius (normalized to x(surface)=1 )
  eta_aniso = 1.

  IF(idom==1)THEN       ! upper crustal layer
     ro_prem  = 2.6
     vpv_prem = 5.8
     vsv_prem = 3.2
     vph_prem = vpv_prem
     vsh_prem = vsv_prem
     Qmu = 600.0
     Qkappa = 57827.0
  else if(idom==2)THEN   ! lower crustal layer
     ro_prem  = 2.9
     vpv_prem = 6.8
     vsv_prem = 3.9
     vph_prem = vpv_prem
     vsh_prem = vsv_prem
     Qmu = 600.0
     Qkappa = 57827.0
  else if(idom==3)THEN   ! upper mantle
     ro_prem   =  2.6910 + 0.6924 * x_prem
     vpv_prem  =  0.8317 + 7.2180 * x_prem
     vph_prem  =  3.5908 + 4.6172 * x_prem
     vsv_prem  =  5.8582 - 1.4678 * x_prem
     vsh_prem  = -1.0839 + 5.7176 * x_prem
     eta_aniso =  3.3687 - 2.4778 * x_prem
     Qmu = 600.0
     Qkappa = 57827.0
  else if(idom==4)THEN
     ro_prem  =  7.1089 -  3.8045 * x_prem
     vpv_prem = 20.3926 - 12.2569 * x_prem
     vsv_prem =  8.9496 -  4.4597 * x_prem
     vph_prem = vpv_prem
     vsh_prem = vsv_prem
     Qmu = 143.0
     Qkappa = 57827.0
  else if(idom==5)THEN
     ro_prem  = 11.2494 -  8.0298 * x_prem
     vpv_prem = 39.7027 - 32.6166 * x_prem
     vsv_prem = 22.3512 - 18.5856 * x_prem
     vph_prem = vpv_prem
     vsh_prem = vsv_prem
     Qmu = 143.0
     Qkappa = 57827.0
  else if(idom==6)THEN
     ro_prem  =  5.3197 - 1.4836 * x_prem
     vpv_prem = 19.0957 - 9.8672 * x_prem
     vsv_prem =  9.9839 - 4.9324 * x_prem
     vph_prem = vpv_prem
     vsh_prem = vsv_prem
     Qmu = 143.0
     Qkappa = 57827.0
  else if(idom==7)THEN   !lower mantle
     ro_prem  =  7.9565 - 6.4761 * x_prem + 5.5283 * x_prem**2 - 3.0807 * x_prem**3
     vpv_prem = 29.2766 -23.6027 * x_prem + 5.5242 * x_prem**2 - 2.5514 * x_prem**3
     vsv_prem = 22.3459 -17.2473 * x_prem - 2.0834 * x_prem**2 + 0.9783 * x_prem**3
     vph_prem = vpv_prem
     vsh_prem = vsv_prem
     Qmu = 312.0
     Qkappa = 57827.0
  else if(idom==8)THEN
     ro_prem  =  7.9565 -  6.4761 * x_prem +  5.5283 * x_prem**2 -  3.0807 * x_prem**3
     vpv_prem = 24.9520 - 40.4673 * x_prem + 51.4832 * x_prem**2 - 26.6419 * x_prem**3
     vsv_prem = 11.1671 - 13.7818 * x_prem + 17.4575 * x_prem**2 -  9.2777 * x_prem**3
     vph_prem = vpv_prem
     vsh_prem = vsv_prem
     Qmu = 312.0
     Qkappa = 57827.0
  else if(idom==9)THEN
     ro_prem  =  7.9565 - 6.4761 * x_prem + 5.5283 * x_prem**2 - 3.0807 * x_prem**3
     vpv_prem = 15.3891 - 5.3181 * x_prem + 5.5242 * x_prem**2 - 2.5514 * x_prem**3
     vsv_prem =  6.9254 + 1.4672 * x_prem - 2.0834 * x_prem**2 + 0.9783 * x_prem**3
     vph_prem = vpv_prem
     vsh_prem = vsv_prem
     Qmu = 312.0
     Qkappa = 57827.0
  else if(idom==10)THEN  ! outer core
     ro_prem  = 12.5815 - 1.2638 * x_prem - 3.6426 * x_prem**2 -  5.5281 * x_prem**3
     vpv_prem = 11.0487 - 4.0362 * x_prem + 4.8023 * x_prem**2 - 13.5732 * x_prem**3
     vsv_prem =  0.0
     vph_prem = vpv_prem
     vsh_prem = vsv_prem
     Qmu = 0.0
     Qkappa = 57827.0
  else if(idom==11)THEN                        ! inner core
     ro_prem  = 13.0885 - 8.8381 * x_prem**2
     vpv_prem = 11.2622 - 6.3640 * x_prem**2
     vsv_prem =  3.6678 - 4.4475 * x_prem**2
     vph_prem = vpv_prem
     vsh_prem = vsv_prem
     Qmu = 84.6
     Qkappa = 1327.7
  endif

  if (param=='rho') then
     prem_ani_sub = ro_prem * 1000.
  else if (param=='vpv') then
     prem_ani_sub = vpv_prem * 1000.
  else if (param=='vsv') then
     prem_ani_sub = vsv_prem * 1000.
  else if (param=='vph') then
     prem_ani_sub = vph_prem * 1000.
  else if (param=='vsh') then
     prem_ani_sub = vsh_prem * 1000.
  else if (param=='eta') then
     prem_ani_sub = eta_aniso
  else if (param=='Qmu') then
     prem_ani_sub = Qmu
  else if (param=='Qka') then
     prem_ani_sub = Qkappa
  !min/max velocities needed for the mesher:
  else if (param=='v_p') then
     prem_ani_sub = max(vpv_prem, vph_prem) * 1000.
  else if (param=='v_s') then
     prem_ani_sub = min(vsv_prem, vsh_prem) * 1000.
  else
     write(6,*)'ERROR IN PREM_ANI_SUB FUNCTION:',param,' NOT AN OPTION'
     stop
  endif

end function prem_ani_sub
!=============================================================================

!-----------------------------------------------------------------------------
!> isotropic prem model in terms of domains separated by discontinuities
!! No fluid outer core, but instead vs=vp/sqrt(3)
real(kind=dp) function prem_solid_sub(r0,param,idom)

  real(kind=dp)   , intent(in) :: r0
  integer, intent(in)          :: idom
  real(kind=dp)                :: r, x_prem
  real(kind=dp)                :: ro_prem, vp_prem, vs_prem
  character(len=3), intent(in) :: param !rho, vs,vp

  r = r0 / 1000.

  x_prem = r / 6371.     ! Radius (normalized to x(surface)=1 )

  if(idom==1)then        ! upper crustal layer
     ro_prem = 2.6
     vp_prem = 5.8
     vs_prem = 3.2
  else if(idom==2)then
     ro_prem = 2.9                       ! lower crustal layer
     vp_prem = 6.8
     vs_prem = 3.9
  else if(idom==3)then
     ro_prem = 2.691  + .6924  * x_prem             ! upper mantle
     vp_prem = 4.1875 + 3.9382 * x_prem
     vs_prem = 2.1519 + 2.3481 * x_prem
  else if(idom==4)then
     ro_prem =  7.1089 - 3.8045 * x_prem
     vp_prem = 20.3926 -12.2569 * x_prem
     vs_prem =  8.9496 - 4.4597 * x_prem
  else if(idom==5)then
     ro_prem = 11.2494 -  8.0298 * x_prem
     vp_prem = 39.7027 - 32.6166 * x_prem
     vs_prem = 22.3512 - 18.5856 * x_prem
  else if(idom==6)then
     ro_prem =  5.3197 - 1.4836 * x_prem
     vp_prem = 19.0957 - 9.8672 * x_prem
     vs_prem =  9.9839 - 4.9324 * x_prem
  else if(idom==7)then   !lower mantle
     ro_prem =  7.9565 -  6.4761 * x_prem + 5.5283 * x_prem**2 - 3.0807 * x_prem**3
     vp_prem = 29.2766 - 23.6027 * x_prem + 5.5242 * x_prem**2 - 2.5514 * x_prem**3
     vs_prem = 22.3459 - 17.2473 * x_prem - 2.0834 * x_prem**2 + 0.9783 * x_prem**3
  else if(idom==8)then
     ro_prem =  7.9565 -  6.4761 * x_prem +  5.5283 * x_prem**2 -  3.0807 * x_prem**3
     vp_prem = 24.9520 - 40.4673 * x_prem + 51.4832 * x_prem**2 - 26.6419 * x_prem**3
     vs_prem = 11.1671 - 13.7818 * x_prem + 17.4575 * x_prem**2 -  9.2777 * x_prem**3
  else if(idom==9)then
     ro_prem =  7.9565 - 6.4761 * x_prem + 5.5283 * x_prem**2 - 3.0807 * x_prem**3
     vp_prem = 15.3891 - 5.3181 * x_prem + 5.5242 * x_prem**2 - 2.5514 * x_prem**3
     vs_prem =  6.9254 + 1.4672 * x_prem - 2.0834 * x_prem**2 + 0.9783 * x_prem**3
  else if(idom==10)then  ! outer core
     ro_prem = 12.5815 - 1.2638 * x_prem - 3.6426 * x_prem**2 -  5.5281 * x_prem**3
     vp_prem = 11.0487 - 4.0362 * x_prem + 4.8023 * x_prem**2 - 13.5732 * x_prem**3
     vs_prem = vp_prem / sqrt(3.)
  else if(idom==11)then                        ! inner core
     ro_prem = 13.0885 - 8.8381 * x_prem**2
     vp_prem = 11.2622 - 6.3640 * x_prem**2
     vs_prem =  3.6678 - 4.4475 * x_prem**2
  endif

  if (param=='rho') then
     prem_solid_sub = ro_prem * 1000.
  else if (param=='v_p') then
     prem_solid_sub = vp_prem * 1000.
  else if (param=='v_s') then
     prem_solid_sub = vs_prem * 1000.
  else if (param=='vpv') then
     prem_solid_sub = vp_prem * 1000.
  else if (param=='vsv') then
     prem_solid_sub = vs_prem * 1000.
  else if (param=='vph') then
     prem_solid_sub = vp_prem * 1000.
  else if (param=='vsh') then
     prem_solid_sub = vs_prem * 1000.
  else if (param=='eta') then
     prem_solid_sub = 1.
  else
     write(6,*)'ERROR IN PREM_SUB FUNCTION:', param, 'NOT AN OPTION'
     stop
  endif

end function prem_solid_sub
!=============================================================================

!-----------------------------------------------------------------------------
!> isotropic prem model in terms of domains separated by discontinuities
!! but with lower crust extended to the surface
real(kind=dp) function prem_onecrust_sub(r0, param, idom)

  real(kind=dp)   , intent(in) :: r0
  integer, intent(in)          :: idom
  real(kind=dp)                :: r, x_prem
  real(kind=dp)                :: ro_prem, vp_prem, vs_prem
  character(len=3), intent(in) :: param !rho, vs,vp

  r = r0 / 1000.

  x_prem = r / 6371.     ! Radius (normalized to x(surface)=1 )

  if(idom==1)then        ! one crustal layer
     ro_prem = 2.6
     vp_prem = 5.8
     vs_prem = 3.2
  else if(idom==2)then      ! upper mantle
     ro_prem = 2.691  + 0.6924 * x_prem
     vp_prem = 4.1875 + 3.9382 * x_prem
     vs_prem = 2.1519 + 2.3481 * x_prem
  else if(idom==3)then
     ro_prem =  7.1089 -  3.8045 * x_prem
     vp_prem = 20.3926 - 12.2569 * x_prem
     vs_prem =  8.9496 -  4.4597 * x_prem
  else if(idom==4)then
     ro_prem = 11.2494 -  8.0298 * x_prem
     vp_prem = 39.7027 - 32.6166 * x_prem
     vs_prem = 22.3512 - 18.5856 * x_prem
  else if(idom==5)then
     ro_prem =  5.3197 - 1.4836 * x_prem
     vp_prem = 19.0957 - 9.8672 * x_prem
     vs_prem =  9.9839 - 4.9324 * x_prem
  else if(idom==6)then   !lower mantle
     ro_prem =  7.9565 -  6.4761 * x_prem + 5.5283 * x_prem**2 - 3.0807 * x_prem**3
     vp_prem = 29.2766 - 23.6027 * x_prem + 5.5242 * x_prem**2 - 2.5514 * x_prem**3
     vs_prem = 22.3459 - 17.2473 * x_prem - 2.0834 * x_prem**2 + 0.9783 * x_prem**3
  else if(idom==7)then
     ro_prem =  7.9565 - 6.4761  * x_prem +  5.5283 * x_prem**2 -  3.0807 * x_prem**3
     vp_prem = 24.9520 - 40.4673 * x_prem + 51.4832 * x_prem**2 - 26.6419 * x_prem**3
     vs_prem = 11.1671 - 13.7818 * x_prem + 17.4575 * x_prem**2 -  9.2777 * x_prem**3
  else if(idom==8)then
     ro_prem =  7.9565 - 6.4761 * x_prem + 5.5283 * x_prem**2 - 3.0807 * x_prem**3
     vp_prem = 15.3891 - 5.3181 * x_prem + 5.5242 * x_prem**2 - 2.5514 * x_prem**3
     vs_prem =  6.9254 + 1.4672 * x_prem - 2.0834 * x_prem**2 + 0.9783 * x_prem**3
  else if(idom==9)then  ! outer core
     ro_prem = 12.5815 - 1.2638 * x_prem - 3.6426 * x_prem**2 -  5.5281 * x_prem**3
     vp_prem = 11.0487 - 4.0362 * x_prem + 4.8023 * x_prem**2 - 13.5732 * x_prem**3
     vs_prem = 0.00
  else if(idom==10)then                        ! inner core
     ro_prem = 13.0885 - 8.8381 * x_prem**2
     vp_prem = 11.2622 - 6.3640 * x_prem**2
     vs_prem =  3.6678 - 4.4475 * x_prem**2
  endif

  if (param=='rho') then
     prem_onecrust_sub = ro_prem * 1000.
  else if (param=='v_p') then
     prem_onecrust_sub = vp_prem * 1000.
  else if (param=='v_s') then
     prem_onecrust_sub = vs_prem * 1000.
  else if (param=='vpv') then
     prem_onecrust_sub = vp_prem * 1000.
  else if (param=='vsv') then
     prem_onecrust_sub = vs_prem * 1000.
  else if (param=='vph') then
     prem_onecrust_sub = vp_prem * 1000.
  else if (param=='vsh') then
     prem_onecrust_sub = vs_prem*1000.
  else if (param=='eta') then
     prem_onecrust_sub = 1.
  else
     write(6,*)'ERROR IN PREM_SUB FUNCTION:', param, 'NOT AN OPTION'
     stop
  endif

end function prem_onecrust_sub
!=============================================================================

!-----------------------------------------------------------------------------
! prem model in terms of domains separated by discontinuities
real(kind=dp) function prem_onecrust_ani_sub(r0, param, idom)

  real(kind=dp)   , intent(in) :: r0
  integer, intent(in)          :: idom
  real(kind=dp)                :: r, x_prem
  real(kind=dp)                :: ro_prem, vpv_prem, vsv_prem, vph_prem
  real(kind=dp)                :: vsh_prem, eta_aniso, Qmu, Qkappa
  character(len=3), intent(in) :: param !rho, vs,vp

  r = r0 / 1000.

  x_prem = r / 6371.     ! Radius (normalized to x(surface)=1 )
  eta_aniso = 1.

  if(idom==1)then       ! upper crustal layer
     ro_prem  = 2.6
     vpv_prem = 5.8
     vsv_prem = 3.2
     vph_prem = vpv_prem
     vsh_prem = vsv_prem
     Qmu = 600.0
     Qkappa = 57827.0
  else if(idom==2)then   ! upper mantle
     ro_prem   =  2.6910 + 0.6924 * x_prem
     vpv_prem  =  0.8317 + 7.2180 * x_prem
     vph_prem  =  3.5908 + 4.6172 * x_prem
     vsv_prem  =  5.8582 - 1.4678 * x_prem
     vsh_prem  = -1.0839 + 5.7176 * x_prem
     eta_aniso =  3.3687 - 2.4778 * x_prem
     Qmu = 600.0
     Qkappa = 57827.0
  else if(idom==3)then
     ro_prem  =  7.1089 -  3.8045 * x_prem
     vpv_prem = 20.3926 - 12.2569 * x_prem
     vsv_prem =  8.9496 -  4.4597 * x_prem
     vph_prem = vpv_prem
     vsh_prem = vsv_prem
     Qmu = 143.0
     Qkappa = 57827.0
  else if(idom==4)then
     ro_prem  = 11.2494 -  8.0298 * x_prem
     vpv_prem = 39.7027 - 32.6166 * x_prem
     vsv_prem = 22.3512 - 18.5856 * x_prem
     vph_prem = vpv_prem
     vsh_prem = vsv_prem
     Qmu = 143.0
     Qkappa = 57827.0
  else if(idom==5)then
     ro_prem  =  5.3197 - 1.4836 * x_prem
     vpv_prem = 19.0957 - 9.8672 * x_prem
     vsv_prem =  9.9839 - 4.9324 * x_prem
     vph_prem = vpv_prem
     vsh_prem = vsv_prem
     Qmu = 143.0
     Qkappa = 57827.0
  else if(idom==6)then   !lower mantle
     ro_prem  =  7.9565 - 6.4761 * x_prem + 5.5283 * x_prem**2 - 3.0807 * x_prem**3
     vpv_prem = 29.2766 -23.6027 * x_prem + 5.5242 * x_prem**2 - 2.5514 * x_prem**3
     vsv_prem = 22.3459 -17.2473 * x_prem - 2.0834 * x_prem**2 + 0.9783 * x_prem**3
     vph_prem = vpv_prem
     vsh_prem = vsv_prem
     Qmu = 312.0
     Qkappa = 57827.0
  else if(idom==7)then
     ro_prem  =  7.9565 -  6.4761 * x_prem +  5.5283 * x_prem**2 -  3.0807 * x_prem**3
     vpv_prem = 24.9520 - 40.4673 * x_prem + 51.4832 * x_prem**2 - 26.6419 * x_prem**3
     vsv_prem = 11.1671 - 13.7818 * x_prem + 17.4575 * x_prem**2 -  9.2777 * x_prem**3
     vph_prem = vpv_prem
     vsh_prem = vsv_prem
     Qmu = 312.0
     Qkappa = 57827.0
  else if(idom==8)then
     ro_prem  =  7.9565 - 6.4761 * x_prem + 5.5283 * x_prem**2 - 3.0807 * x_prem**3
     vpv_prem = 15.3891 - 5.3181 * x_prem + 5.5242 * x_prem**2 - 2.5514 * x_prem**3
     vsv_prem =  6.9254 + 1.4672 * x_prem - 2.0834 * x_prem**2 + 0.9783 * x_prem**3
     vph_prem = vpv_prem
     vsh_prem = vsv_prem
     Qmu = 312.0
     Qkappa = 57827.0
  else if(idom==9)then  ! outer core
     ro_prem  = 12.5815 - 1.2638 * x_prem - 3.6426 * x_prem**2 -  5.5281 * x_prem**3
     vpv_prem = 11.0487 - 4.0362 * x_prem + 4.8023 * x_prem**2 - 13.5732 * x_prem**3
     vsv_prem = 0.0
     vph_prem = vpv_prem
     vsh_prem = vsv_prem
     Qmu = 0.0
     Qkappa = 57827.0
  else if(idom==10)then                        ! inner core
     ro_prem  = 13.0885 - 8.8381 * x_prem**2
     vpv_prem = 11.2622 - 6.3640 * x_prem**2
     vsv_prem =  3.6678 - 4.4475 * x_prem**2
     vph_prem = vpv_prem
     vsh_prem = vsv_prem
     Qmu = 84.6
     Qkappa = 1327.7
  endif

  if (param=='rho') then
     prem_onecrust_ani_sub = ro_prem * 1000.
  else if (param=='vpv') then
     prem_onecrust_ani_sub = vpv_prem * 1000.
  else if (param=='vsv') then
     prem_onecrust_ani_sub = vsv_prem * 1000.
  else if (param=='vph') then
     prem_onecrust_ani_sub = vph_prem * 1000.
  else if (param=='vsh') then
     prem_onecrust_ani_sub = vsh_prem * 1000.
  else if (param=='eta') then
     prem_onecrust_ani_sub = eta_aniso
  else if (param=='Qmu') then
     prem_onecrust_ani_sub = Qmu
  else if (param=='Qka') then
     prem_onecrust_ani_sub = Qkappa
  !min/max velocities needed for the mesher:
  else if (param=='v_p') then
     prem_onecrust_ani_sub = max(vpv_prem, vph_prem) * 1000.
  else if (param=='v_s') then
     prem_onecrust_ani_sub = min(vsv_prem, vsh_prem) * 1000.
  else
     write(6,*)'ERROR IN PREM_ANI_SUB FUNCTION:', param, ' NOT AN OPTION'
     stop
  endif

end function prem_onecrust_ani_sub
!=============================================================================

!-----------------------------------------------------------------------------
! isotropic prem_light model (crust removed) in terms of domains separated by disconts.
real(kind=dp) function prem_light_sub(r0, param, idom)

  real(kind=dp)   , intent(in) :: r0
  integer, intent(in)          :: idom
  real(kind=dp)                :: r, x_prem
  real(kind=dp)                :: ro_prem, vp_prem, vs_prem
  real(kind=dp)                :: Qmu, Qkappa
  character(len=3), intent(in) :: param !rho, vs,vp

  r = r0 / 1000.

  x_prem = r / 6371.     ! Radius (normalized to x(surface)=1 )

  if(idom==1)then
     ro_prem = 2.691  + 0.6924 * x_prem             ! upper mantle
     vp_prem = 4.1875 + 3.9382 * x_prem
     vs_prem = 2.1519 + 2.3481 * x_prem
     Qmu = 600.0
     Qkappa = 57827.0
  else if(idom==2)then
     ro_prem =  7.1089 -  3.8045 * x_prem
     vp_prem = 20.3926 - 12.2569 * x_prem
     vs_prem =  8.9496 -  4.4597 * x_prem
     Qmu = 143.0
     Qkappa = 57827.0
  else if(idom==3)then
     ro_prem = 11.2494 -  8.0298 * x_prem
     vp_prem = 39.7027 - 32.6166 * x_prem
     vs_prem = 22.3512 - 18.5856 * x_prem
     Qmu = 143.0
     Qkappa = 57827.0
  else if(idom==4)then
     ro_prem =  5.3197 - 1.4836 * x_prem
     vp_prem = 19.0957 - 9.8672 * x_prem
     vs_prem =  9.9839 - 4.9324 * x_prem
     Qmu = 143.0
     Qkappa = 57827.0
  else if(idom==5)then   !lower mantle
     ro_prem =  7.9565-  6.4761 * x_prem + 5.5283 * x_prem**2 - 3.0807 * x_prem**3
     vp_prem = 29.2766- 23.6027 * x_prem + 5.5242 * x_prem**2 - 2.5514 * x_prem**3
     vs_prem = 22.3459- 17.2473 * x_prem - 2.0834 * x_prem**2 + 0.9783 * x_prem**3
     Qmu = 312.0
     Qkappa = 57827.0
  else if(idom==6)then
     ro_prem =  7.9565 -  6.4761 * x_prem +  5.5283 * x_prem**2 -  3.0807 * x_prem**3
     vp_prem = 24.9520 - 40.4673 * x_prem + 51.4832 * x_prem**2 - 26.6419 * x_prem**3
     vs_prem = 11.1671 - 13.7818 * x_prem + 17.4575 * x_prem**2 -  9.2777 * x_prem**3
     Qmu = 312.0
     Qkappa = 57827.0
  else if(idom==7)then
     ro_prem =  7.9565 - 6.4761 * x_prem + 5.5283 * x_prem**2 - 3.0807 * x_prem**3
     vp_prem = 15.3891 - 5.3181 * x_prem + 5.5242 * x_prem**2 - 2.5514 * x_prem**3
     vs_prem =  6.9254 + 1.4672 * x_prem - 2.0834 * x_prem**2 + 0.9783 * x_prem**3
     Qmu = 312.0
     Qkappa = 57827.0
  else if(idom==8)then  ! outer core
     ro_prem = 12.5815 - 1.2638 * x_prem - 3.6426 * x_prem**2 -  5.5281 * x_prem**3
     vp_prem = 11.0487 - 4.0362 * x_prem + 4.8023 * x_prem**2 - 13.5732 * x_prem**3
     vs_prem = 0.0
     Qmu = 0.0
     Qkappa = 57827.0
  else if(idom==9)then                        ! inner core
     ro_prem = 13.0885 - 8.8381 * x_prem**2
     vp_prem = 11.2622 - 6.3640 * x_prem**2
     vs_prem =  3.6678 - 4.4475 * x_prem**2
     Qmu = 84.6
     Qkappa = 1327.7
  endif

  if (param=='rho') then
     prem_light_sub = ro_prem * 1000.
  else if (param=='v_p') then
     prem_light_sub = vp_prem * 1000.
  else if (param=='v_s') then
     prem_light_sub = vs_prem * 1000.
  else if (param=='vpv') then
     prem_light_sub = vp_prem * 1000.
  else if (param=='vsv') then
     prem_light_sub = vs_prem * 1000.
  else if (param=='vph') then
     prem_light_sub = vp_prem * 1000.
  else if (param=='vsh') then
     prem_light_sub = vs_prem * 1000.
  else if (param=='eta') then
     prem_light_sub = 1.
  else if (param=='Qmu') then
     prem_light_sub = Qmu
  else if (param=='Qka') then
     prem_light_sub = Qkappa
  else
     write(6,*)'ERROR IN PREM_LIGHT_SUB FUNCTION:', param, 'NOT AN OPTION'
     stop
  endif

end function prem_light_sub
!=============================================================================

!-----------------------------------------------------------------------------
! anisotropic prem_light model (crust removed) in terms of domains separated by disconts.
real(kind=dp) function prem_light_ani_sub(r0, param, idom)

  real(kind=dp)   , intent(in) :: r0
  integer, intent(in)          :: idom
  real(kind=dp)                :: r,x_prem
  real(kind=dp)                :: ro_prem, vpv_prem, vsv_prem, vph_prem
  real(kind=dp)                :: vsh_prem, eta_aniso, Qmu, Qkappa
  character(len=3), intent(in) :: param !rho, vs,vp

  r = r0 / 1000.

  x_prem = r / 6371.     ! Radius (normalized to x(surface)=1 )
  eta_aniso = 1.

  if(idom==1)then
     ro_prem  =  2.6910 + 0.6924 * x_prem
     vpv_prem =  0.8317 + 7.2180 * x_prem
     vph_prem =  3.5908 + 4.6172 * x_prem
     vsv_prem =  5.8582 - 1.4678 * x_prem
     vsh_prem = -1.0839 + 5.7176 * x_prem
     eta_aniso=  3.3687 - 2.4778 * x_prem
     Qmu = 600.0
     Qkappa = 57827.0
  else if(idom==2)then
     ro_prem  =  7.1089 - 3.8045 * x_prem
     vpv_prem = 20.3926 -12.2569 * x_prem
     vsv_prem =  8.9496 - 4.4597 * x_prem
     vph_prem = vpv_prem
     vsh_prem = vsv_prem
     Qmu = 143.0
     Qkappa = 57827.0
  else if(idom==3)then
     ro_prem  = 11.2494 -  8.0298 * x_prem
     vpv_prem = 39.7027 - 32.6166 * x_prem
     vsv_prem = 22.3512 - 18.5856 * x_prem
     vph_prem = vpv_prem
     vsh_prem = vsv_prem
     Qmu = 143.0
     Qkappa = 57827.0
  else if(idom==4)then
     ro_prem  =  5.3197 - 1.4836 * x_prem
     vpv_prem = 19.0957 - 9.8672 * x_prem
     vsv_prem =  9.9839 - 4.9324 * x_prem
     vph_prem = vpv_prem
     vsh_prem = vsv_prem
     Qmu = 143.0
     Qkappa = 57827.0
  else if(idom==5)then   !lower mantle
     ro_prem =  7.9565 -  6.4761 * x_prem + 5.5283 * x_prem**2 - 3.0807 * x_prem**3
     vpv_prem= 29.2766 - 23.6027 * x_prem + 5.5242 * x_prem**2 - 2.5514 * x_prem**3
     vsv_prem= 22.3459 - 17.2473 * x_prem - 2.0834 * x_prem**2 + 0.9783 * x_prem**3
     vph_prem = vpv_prem
     vsh_prem = vsv_prem
     Qmu = 312.0
     Qkappa = 57827.0
  else if(idom==6)then
     ro_prem  =  7.9565 -  6.4761 * x_prem +  5.5283 * x_prem**2 -  3.0807 * x_prem**3
     vpv_prem = 24.9520 - 40.4673 * x_prem + 51.4832 * x_prem**2 - 26.6419 * x_prem**3
     vsv_prem = 11.1671 - 13.7818 * x_prem + 17.4575 * x_prem**2 -  9.2777 * x_prem**3
     vph_prem = vpv_prem
     vsh_prem = vsv_prem
     Qmu = 312.0
     Qkappa = 57827.0
  else if(idom==7)then
     ro_prem  =  7.9565 - 6.4761 * x_prem + 5.5283 * x_prem**2 - 3.0807 * x_prem**3
     vpv_prem = 15.3891 - 5.3181 * x_prem + 5.5242 * x_prem**2 - 2.5514 * x_prem**3
     vsv_prem =  6.9254 + 1.4672 * x_prem - 2.0834 * x_prem**2 + 0.9783 * x_prem**3
     vph_prem = vpv_prem
     vsh_prem = vsv_prem
     Qmu = 312.0
     Qkappa = 57827.0
  else if(idom==8)then  ! outer core
     ro_prem  = 12.5815 - 1.2638 * x_prem - 3.6426 * x_prem**2 -  5.5281 * x_prem**3
     vpv_prem = 11.0487 - 4.0362 * x_prem + 4.8023 * x_prem**2 - 13.5732 * x_prem**3
     vsv_prem = 0.0
     vph_prem = vpv_prem
     vsh_prem = vsv_prem
     Qmu = 0.0
     Qkappa = 57827.0
  else if(idom==9)then                        ! inner core
     ro_prem  = 13.0885 - 8.8381 * x_prem**2
     vpv_prem = 11.2622 - 6.3640 * x_prem**2
     vsv_prem =  3.6678 - 4.4475 * x_prem**2
     vph_prem = vpv_prem
     vsh_prem = vsv_prem
     Qmu = 84.6
     Qkappa = 1327.7
  endif

  if (param=='rho') then
     prem_light_ani_sub = ro_prem * 1000.
  else if (param=='vpv') then
     prem_light_ani_sub = vpv_prem * 1000.
  else if (param=='vsv') then
     prem_light_ani_sub = vsv_prem * 1000.
  else if (param=='vph') then
     prem_light_ani_sub = vph_prem * 1000.
  else if (param=='vsh') then
     prem_light_ani_sub = vsh_prem * 1000.
  else if (param=='eta') then
     prem_light_ani_sub = eta_aniso
  else if (param=='Qmu') then
     prem_light_ani_sub = Qmu
  else if (param=='Qka') then
     prem_light_ani_sub = Qkappa
  !min/max velocities needed for the mesher:
  else if (param=='v_p') then
     prem_light_ani_sub = max(vpv_prem, vph_prem) * 1000.
  else if (param=='v_s') then
     prem_light_ani_sub = min(vsv_prem, vsh_prem) * 1000.
  else
     write(6,*)'ERROR IN PREM_ANI_SUB FUNCTION:', param, ' NOT AN OPTION'
     stop
  endif

end function prem_light_ani_sub
!=============================================================================

!-----------------------------------------------------------------------------
!> isotropic prem_light model (crust removed) in terms of domains separated by disconts.
!! No fluid outer core, but instead vs=vp/sqrt(3)
real(kind=dp) function prem_solid_light_sub(r0, param, idom)

  real(kind=dp)   , intent(in) :: r0
  integer, intent(in)          :: idom
  real(kind=dp)                :: r, x_prem
  real(kind=dp)                :: ro_prem, vp_prem, vs_prem
  character(len=3), intent(in) :: param !rho, vs,vp

  r = r0 / 1000.

  x_prem = r / 6371.     ! Radius (normalized to x(surface)=1 )

  if(idom==1)then
     ro_prem = 2.691  + 0.6924 * x_prem             ! upper mantle
     vp_prem = 4.1875 + 3.9382 * x_prem
     vs_prem = 2.1519 + 2.3481 * x_prem
  else if(idom==2)then
     ro_prem =  7.1089 -  3.8045 * x_prem
     vp_prem = 20.3926 - 12.2569 * x_prem
     vs_prem =  8.9496 -  4.4597 * x_prem
  else if(idom==3)then
     ro_prem = 11.2494 -  8.0298 * x_prem
     vp_prem = 39.7027 - 32.6166 * x_prem
     vs_prem = 22.3512 - 18.5856 * x_prem
  else if(idom==4)then
     ro_prem =  5.3197 - 1.4836 * x_prem
     vp_prem = 19.0957 - 9.8672 * x_prem
     vs_prem =  9.9839 - 4.9324 * x_prem
  else if(idom==5)then   !lower mantle
     ro_prem =  7.9565 -  6.4761 * x_prem + 5.5283 * x_prem**2 - 3.0807 * x_prem**3
     vp_prem = 29.2766 - 23.6027 * x_prem + 5.5242 * x_prem**2 - 2.5514 * x_prem**3
     vs_prem = 22.3459 - 17.2473 * x_prem - 2.0834 * x_prem**2 + 0.9783 * x_prem**3
  else if(idom==6)then
     ro_prem =  7.9565 -  6.4761 * x_prem +  5.5283 * x_prem**2 -  3.0807 * x_prem**3
     vp_prem = 24.9520 - 40.4673 * x_prem + 51.4832 * x_prem**2 - 26.6419 * x_prem**3
     vs_prem = 11.1671 - 13.7818 * x_prem + 17.4575 * x_prem**2 -  9.2777 * x_prem**3
  else if(idom==7)then
     ro_prem =  7.9565 - 6.4761 * x_prem + 5.5283 * x_prem**2 - 3.0807 * x_prem**3
     vp_prem = 15.3891 - 5.3181 * x_prem + 5.5242 * x_prem**2 - 2.5514 * x_prem**3
     vs_prem =  6.9254 + 1.4672 * x_prem - 2.0834 * x_prem**2 + 0.9783 * x_prem**3
  else if(idom==8)then  ! outer core
     ro_prem = 12.5815 - 1.2638 * x_prem - 3.6426 * x_prem**2 -  5.5281 * x_prem**3
     vp_prem = 11.0487 - 4.0362 * x_prem + 4.8023 * x_prem**2 - 13.5732 * x_prem**3
     vs_prem = vp_prem / sqrt(3.)
  else if(idom==9)then                        ! inner core
     ro_prem = 13.0885 - 8.8381 * x_prem**2
     vp_prem = 11.2622 - 6.3640 * x_prem**2
     vs_prem =  3.6678 - 4.4475 * x_prem**2
  endif

  if (param=='rho') then
     prem_solid_light_sub = ro_prem * 1000.
  else if (param=='v_p') then
     prem_solid_light_sub = vp_prem * 1000.
  else if (param=='v_s') then
     prem_solid_light_sub = vs_prem * 1000.
  else if (param=='vpv') then
     prem_solid_light_sub = vp_prem * 1000.
  else if (param=='vsv') then
     prem_solid_light_sub = vs_prem * 1000.
  else if (param=='vph') then
     prem_solid_light_sub = vp_prem * 1000.
  else if (param=='vsh') then
     prem_solid_light_sub = vs_prem * 1000.
  else if (param=='eta') then
     prem_solid_light_sub = 1.
  else
     write(6,*)'ERROR IN PREM_LIGHT_SUB FUNCTION:', param, 'NOT AN OPTION'
     stop
  endif

end function prem_solid_light_sub
!=============================================================================

!-----------------------------------------------------------------------------
!> iasp91 model in terms of domains separated by discontinuities
!! with PREM density and attenuation
real(kind=dp) function iasp91_sub(r0, param, idom)

  real(kind=dp)   , intent(in)    :: r0
  integer, intent(in)             :: idom
  character(len=3), intent(in)    :: param !rho, vs,vp
  real(kind=dp)                   :: r, x
  real(kind=dp)                   :: rho, vp, vs, Qmu, Qkappa

  real(kind=dp)                   :: REARTH
  real(kind=dp)                   :: R120, RMOHO
  real(kind=dp)                   :: x1, x2

  REARTH           = 6371000.
  ! compute real physical radius in meters
  r = r0
  x = r / REARTH     ! Radius (normalized to x(surface)=1 )


  ! IASP91
  RMOHO            = 6336000.
  R120             = 6251000.

  x1 = R120 / REARTH
  x2 = RMOHO / REARTH

  if(idom==1)then ! upper crust
     vp  = 5.8
     vs  = 3.36
     rho = 2.72
     Qmu = 600.0
     Qkappa = 57827.0
  else if(idom==2)then ! lower crust
     vp  = 6.5
     vs  = 3.75
     rho = 2.92
     Qmu = 600.0
     Qkappa = 57827.0
  else if(idom==3)then ! R120 < r <= RMOHO
     vp = 8.78541 - 0.74953  * x
     vs = 6.706231- 2.248585 * x
     !rho = 3.3713 + (3.3198 - 3.3713) * (x - x1) / (x2 - x1)
     !! VM VM modif
     rho =  2.69100 +  0.69240 * x
     Qmu = 600.0
     Qkappa = 57827.0

     ! MvD: keeping this test for old bug [18]
!!$     if(rho < 3.319 .or. rho > 3.372) then
!!$        write(6,*) R120 / 1000., RMOHO / 1000.
!!$        write(6,*) r0 / 1000.
!!$        write(6,*) idom
!!$        write(6,*) x, x1, x2
!!$        write(6,*) x - x1
!!$        write(6,*) x2 - x1
!!$        write(6,*) (x - x1) / (x2 - x1)
!!$        write(6,*) 'incorrect density computed for IASP91', rho
!!$        write(6,*) 'known bug, use other velocity model for now'
!!$        call abort()
!!$        stop 2
!!$     endif
  else if(idom==4)then ! R220 < r <= R120
     rho =  2.6910  +  0.6924  * x
     vp  = 25.41389 - 17.69722 * x
     vs  =  5.75020 -  1.2742  * x
     Qmu = 600.0
     Qkappa = 57827.0
  else if(idom==5)then ! R400 < r <= R220
     rho =  7.1089  -  3.8045  * x
     vp  = 30.78765 - 23.25415 * x
     vs  = 15.24213 - 11.08552 * x
     Qmu = 143.0
     Qkappa = 57827.0
  else if(idom==6)then ! R670 < r <= R400
     !rho =  5.3197  -  1.4836  * x
     rho =  11.24940  -8.02980 * x !! VM VM modif
     vp  = 29.38896 - 21.40656 * x
     vs  = 17.70732 - 13.50652 * x
     Qmu = 143.0
     Qkappa = 57827.0
  else if(idom==7)then ! R771 < r <= R670
     rho =  7.9565  -  6.4761  *x + 5.5283 * x**2 - 3.0807 * x**3
     vp  = 25.96984 - 16.93412 *x
     vs  = 20.76890 - 16.53147 *x
     Qmu = 312.0
     Qkappa = 57827.0
  else if(idom==8)then ! RTOPDDOUBLEPRIME < r <= R771
     rho =  7.9565 -  6.4761 * x +  5.5283 * x**2 -  3.0807 * x**3
     vp  = 25.1486 - 41.1538 * x + 51.9932 * x**2 - 26.6083 * x**3
     vs  = 12.9303 - 21.2590 * x + 27.8988 * x**2 - 14.1080 * x**3
     Qmu = 312.0
     Qkappa = 57827.0
  else if(idom==9)then ! RCMB < r <= RTOPDDOUBLEPRIME
     rho =  7.9565  - 6.4761  * x + 5.5283 * x**2 - 3.0807 * x**3
     vp  = 14.49470 - 1.47089 * x
     vs  =  8.16616 - 1.58206 * x
     Qmu = 312.0
     Qkappa = 57827.0
  else if(idom==10)then ! RICB < r <= RCMB
     rho = 12.5815  - 1.2638  * x -  3.6426  *x **2 - 5.5281 * x**3
     vp  = 10.03904 + 3.75665 * x - 13.67046 *x **2
     vs  =  0.0
     Qmu = 0.0
     Qkappa = 57827.0
  else if(idom==11)then ! 0. < r <= RICB
     rho = 13.0885  - 8.8381  * x**2
     vp  = 11.24094 - 4.09689 * x**2
     vs  =  3.56454 - 3.45241 * x**2
     Qmu = 84.6
     Qkappa = 1327.7
  else
     write(6,*) 'iasp91_sub: error with domain idom=', idom
     stop
  endif

  if (param=='rho') then
     iasp91_sub = rho * 1000.
  else if (param=='v_p') then
     iasp91_sub = vp * 1000.
  else if (param=='v_s') then
     iasp91_sub = vs * 1000.
  else if (param=='vpv') then
     iasp91_sub = vp * 1000.
  else if (param=='vsv') then
     iasp91_sub = vs * 1000.
  else if (param=='vph') then
     iasp91_sub = vp * 1000.
  else if (param=='vsh') then
     iasp91_sub = vs * 1000.
  else if (param=='eta') then
     iasp91_sub = 1.
  else if (param=='Qmu') then
     iasp91_sub = Qmu
  else if (param=='Qka') then
     iasp91_sub = Qkappa
  else
     write(6,*)'ERROR IN IASP91_SUB FUNCTION:',param,'NOT AN OPTION'
     stop
  endif

end function iasp91_sub
!=============================================================================

!-----------------------------------------------------------------------------
!> file-based, step-wise model in terms of domains separated by disconts.
!! format:
!! ndisc
!! r vp vs rho
!! ...
real(kind=dp) function arbitr_sub_solar(r0, param, idom)

  real(kind=dp), intent(in)      :: r0
  integer, intent(in)            :: idom
  character(len=3), intent(in)   :: param
  logical                        :: success
  type(interpolation_data)       :: interp

     !print *, 'R0: ', r0, ', idom:', idom
     select case(param)
     case('rho')
        call interpolate(interp_rho(idom), r0, arbitr_sub_solar, success)

     case('v_p', 'vph', 'vpv')
        call interpolate(interp_vpv(idom), r0, arbitr_sub_solar, success)

     case('v_s', 'vsh', 'vsv')
        call interpolate(interp_vsv(idom), r0, arbitr_sub_solar, success)

     case('eta')
        arbitr_sub_solar = 1
        success = .true.

     case('Qmu', 'Qka')
        select case(trim(override_ext_q))
        case('prem')
           arbitr_sub_solar = prem_sub(r0, param, idom)
           success = .true.

        case('ak135f')
           arbitr_sub_solar = ak135f(r0, param, idom)
           success = .true.

        case default
           if (ext_model_is_anelastic) then
              if (param=='Qmu') then
                 call interpolate(interp_qmu(idom), r0, arbitr_sub_solar, success)
              else
                 call interpolate(interp_qka(idom), r0, arbitr_sub_solar, success)
              endif
           else
              print *, 'ERROR: Parameter: ', trim(param), ' requested, but '
              print *, '       external model is purely elastic. Set OVERRIDE_EXT_Q in '
              print *, '       inparam_mesh, if you want to use PREM or AK135F Q.'
              stop
           endif

        end select

     case default
        print *, 'ERROR: Parameter ', trim(param), ' not implemented in external model'

     end select

     if (.not.success) stop 'Interpolation of parameter not successful'


end function arbitr_sub_solar
!=============================================================================

!=============================================================================
subroutine read_ext_model(fnam_ext_model, nlayer_out, rho_layer_out, &
                          vpv_layer_out, vsv_layer_out, radius_layer_out)

  character(len=*), intent(in)                       :: fnam_ext_model
  real(kind=dp), allocatable, intent(out), optional  :: vpv_layer_out(:)
  real(kind=dp), allocatable, intent(out), optional  :: vsv_layer_out(:)
  real(kind=dp), allocatable, intent(out), optional  :: rho_layer_out(:)
  real(kind=dp), allocatable, intent(out), optional  :: radius_layer_out(:)
  integer, intent(out), optional                     :: nlayer_out
  integer                          :: ilayer, ierr
  logical                          :: bkgrdmodelfile_exists, startatsurface
  character(len=128)               :: fmtstring

  ! Has the file already been read in?
  if (nlayer<1) then

     ! Does the file bkgrdmodel".bm" exist?
     inquire(file=trim(fnam_ext_model), exist=bkgrdmodelfile_exists)

     if (.not. bkgrdmodelfile_exists) then
        write(6,*)'ERROR: File: ', trim(fnam_ext_model),' does not exist!'
        stop
     endif

     open(unit=77, file=trim(fnam_ext_model), action='read')


     read(77,*) ext_model_is_ani, ext_model_is_anelastic
     read(77,*) nlayer
     fmtstring = '(A,A,A,I5,A)'
     write(6,fmtstring,advance='no') 'Model in file ', trim(fnam_ext_model), ' has ', nlayer, ' layers'
     fmtstring = '(A)'
     if (ext_model_is_ani) then
         write(6, fmtstring, advance='no') ' and is anisotropic'
     else
         write(6, fmtstring, advance='no') ' and is isotropic'
     endif
     if (ext_model_is_anelastic) then
         print *, 'and anelastic...'
     else
         print *, 'and elastic...'
     endif

     allocate(radius_layer(nlayer))
     allocate(vpv_layer(nlayer))
     allocate(vsv_layer(nlayer))
     allocate(rho_layer(nlayer))

     if (ext_model_is_anelastic) then
         allocate(qka_layer(nlayer))
         allocate(qmu_layer(nlayer))
     endif


     ! Read in first layer
     if (ext_model_is_anelastic) then
         read(77,*) radius_layer(1), rho_layer(1), vpv_layer(1), &
                    vsv_layer(1), qka_layer(1), qmu_layer(1)
     else
         read(77,*) radius_layer(1), rho_layer(1), vpv_layer(1), vsv_layer(1)
     endif

     ! Recognize order of layers
     if (radius_layer(1)==0) then
        print *, 'Layers in file ', trim(fnam_ext_model), ' start in the core'
        startatsurface = .false.
     else
        print *, 'Layers in file ', trim(fnam_ext_model), ' start at surface'
        startatsurface = .true.
     endif

     ! Read in all other layers
     do ilayer = 2, nlayer
        if (ext_model_is_anelastic) then
           read(77,*,iostat=ierr) radius_layer(ilayer), rho_layer(ilayer), vpv_layer(ilayer),&
                                   vsv_layer(ilayer), qka_layer(ilayer), qmu_layer(ilayer)
        else
           read(77,*,iostat=ierr) radius_layer(ilayer), rho_layer(ilayer), vpv_layer(ilayer),&
                                   vsv_layer(ilayer)
        endif

        if (ierr==IOSTAT_END) then
           print *, 'ERROR: File ', trim(fnam_ext_model), ' has only ', ilayer-1, ' layers'
           print *, '       Not ', nlayer, ' as specified in the header.'
           stop
        endif
        if (startatsurface) then
           if ((radius_layer(ilayer) - radius_layer(ilayer-1))>0.0d0) then
              print *, 'ERROR: Radius of layers in external model has to be monotonously (increasing)!'
              print *, 'Radius of layer:', ilayer-1, ' is:' , radius_layer(ilayer-1)
              print *, 'Radius of layer:', ilayer, ' is:' , radius_layer(ilayer)
              stop
           endif
        else
           if ((radius_layer(ilayer) - radius_layer(ilayer-1))<0.0d0) then
              print *, 'ERROR: Radius of layers in external model has to be monotonously (decreasing)!'
              print *, 'Radius of layer:', ilayer-1, ' is:' , radius_layer(ilayer-1)
              print *, 'Radius of layer:', ilayer, ' is:' , radius_layer(ilayer)
              stop
           endif
        endif
     enddo
     close(77)


     ! Reorder elements if the file is in the wrong order (Mesher always assumes from surface to core)
     if (.not.startatsurface) then
        print *, 'Reordering layers to start at the surface'
        radius_layer = radius_layer(nlayer:1:-1)
        vpv_layer    = vpv_layer(nlayer:1:-1)
        vsv_layer    = vsv_layer(nlayer:1:-1)
        rho_layer    = rho_layer(nlayer:1:-1)
        if (ext_model_is_anelastic) then
           qka_layer = qka_layer(nlayer:1:-1)
           qmu_layer = qmu_layer(nlayer:1:-1)
        endif
     endif

     ! Lowermost layer should be at the center of the earth.
     if (radius_layer(nlayer)>smallval_dble) radius_layer(nlayer) = 0

  endif

  if (present(nlayer_out)) then
     nlayer_out = nlayer
  endif
  if (present(radius_layer_out)) then
     allocate(radius_layer_out(nlayer))
     radius_layer_out = radius_layer
  endif
  if (present(rho_layer_out)) then
     allocate(rho_layer_out(nlayer))
     rho_layer_out = rho_layer
  endif
  if (present(vpv_layer_out)) then
     allocate(vpv_layer_out(nlayer))
     vpv_layer_out = vpv_layer
  endif
  if (present(vsv_layer_out)) then
     allocate(vsv_layer_out(nlayer))
     vsv_layer_out = vsv_layer
  endif

end subroutine read_ext_model
!-----------------------------------------------------------------------------------------

!=============================================================================
subroutine get_ext_disc(fnam_ext_model, ndisc_out, discont, vp, vs, rho)

  use global_parameters, only : smallval_dble

  character(len=*)                :: fnam_ext_model
  integer, intent(out), optional  :: ndisc_out
  real(kind=dp), allocatable, intent(out), optional :: discont(:), vp(:,:), vs(:,:), rho(:,:)
  integer :: idom, junk, ilayer
  real(kind=dp), allocatable :: grad_vp(:), grad_vs(:)

  integer, parameter         :: ndom_max = 100
  real(kind=dp), parameter   :: grad_threshold = 1.d-2
  real(kind=dp)              :: disc_tmp(ndom_max), vp_tmp(ndom_max,2), vs_tmp(ndom_max,2), rho_tmp(ndom_max,2)
  real(kind=dp)              :: vp_laststep, vs_laststep, rho_laststep, dx, dx_min
  integer                    :: upper_layer(ndom_max), lower_layer(ndom_max), ndisc
  integer, allocatable       :: isdisc(:)
  character(len=128)         :: fmtstring

  call read_ext_model(fnam_ext_model)

  allocate(grad_vp(nlayer-1))
  allocate(grad_vs(nlayer-1))

  ! Calculate gradient

  allocate(isdisc(nlayer))

  grad_vs = 0.0
  grad_vp = 0.0
  isdisc  = 0

  idom = 1

  upper_layer(1) = 1

  do ilayer = 2, nlayer-1
     if (abs(radius_layer(ilayer+1) - radius_layer(ilayer)) < smallval_dble) then
        ! First order discontinuity
        idom = idom + 1
        isdisc(ilayer) = 1
        lower_layer(idom-1) = ilayer
        upper_layer(idom)   = ilayer + 1
     else
        grad_vp(ilayer) = (vpv_layer(ilayer+1) - vpv_layer(ilayer)) / &
                          (radius_layer(ilayer+1) - radius_layer(ilayer))
        grad_vs(ilayer) = (vsv_layer(ilayer+1) - vsv_layer(ilayer)) / &
                          (radius_layer(ilayer+1) - radius_layer(ilayer))
        if (abs(grad_vp(ilayer))>grad_threshold.or. &
            abs(grad_vs(ilayer))>grad_threshold) then
           ! Second order discontinuity
           idom = idom + 1
           isdisc(ilayer) = 2
           lower_layer(idom-1) = ilayer
           upper_layer(idom)   = ilayer
        endif

     endif
  enddo

  lower_layer(idom) = nlayer

  ndisc = idom ! The first discontinuity is at the surface,
               ! the last at the ICB, above the last domain

  print *, '  External model has', idom, ' layers'
  ndisc = idom

  print *, '  Creating interpolation objects'

  ! Create interpolation objects for each domain
  allocate(interp_vpv(ndisc))
  allocate(interp_vsv(ndisc))
  allocate(interp_rho(ndisc))

  if (ext_model_is_anelastic) then
      allocate(interp_qka(ndisc))
      allocate(interp_qmu(ndisc))
  endif

  print *, '   idom, upper_layer, lower_layer,   r(ul),   r(ll)'
  fmtstring = '(I8, I13, I13, F9.1, F9.1)'
  do idom = 1, ndisc
     print fmtstring, idom, upper_layer(idom), lower_layer(idom), radius_layer(upper_layer(idom)), radius_layer(lower_layer(idom))
     interp_rho(idom) = interpolation_object(radius_layer(upper_layer(idom):lower_layer(idom)), &
                                             rho_layer(upper_layer(idom):lower_layer(idom)), &
                                             extrapolation_none)
     interp_vpv(idom) = interpolation_object(radius_layer(upper_layer(idom):lower_layer(idom)), &
                                             vpv_layer(upper_layer(idom):lower_layer(idom)), &
                                             extrapolation_none)
     interp_vsv(idom) = interpolation_object(radius_layer(upper_layer(idom):lower_layer(idom)), &
                                             vsv_layer(upper_layer(idom):lower_layer(idom)), &
                                             extrapolation_none)
     if (ext_model_is_anelastic) then
         interp_qka(idom) = interpolation_object(radius_layer(upper_layer(idom):lower_layer(idom)), &
                                                 qka_layer(upper_layer(idom):lower_layer(idom)), &
                                                 extrapolation_none)
         interp_qmu(idom) = interpolation_object(radius_layer(upper_layer(idom):lower_layer(idom)), &
                                                 qmu_layer(upper_layer(idom):lower_layer(idom)), &
                                                 extrapolation_none)
     endif
  enddo

  if (present(ndisc_out)) then
     ndisc_out = ndisc
     allocate(discont(ndisc))
     allocate(vp(ndisc,2))
     allocate(vs(ndisc,2))
     allocate(rho(ndisc,2))
     discont        = radius_layer(upper_layer(1:ndisc))
     vp(1:ndisc,1)  = vpv_layer(upper_layer(1:ndisc))
     vp(1:ndisc,2)  = vpv_layer(lower_layer(1:ndisc))
     vs(1:ndisc,1)  = vsv_layer(upper_layer(1:ndisc))
     vs(1:ndisc,2)  = vsv_layer(lower_layer(1:ndisc))
     rho(1:ndisc,1) = rho_layer(upper_layer(1:ndisc))
     rho(1:ndisc,2) = rho_layer(lower_layer(1:ndisc))
  endif


end subroutine

!-----------------------------------------------------------------------------
!> isotropic prem model in terms of domains separated by discontinuities
real(kind=dp) function prem_sub_vm(r0, param, idom)

  real(kind=dp)   , intent(in) :: r0
  integer, intent(in)          :: idom
  real(kind=dp)                :: r,x_prem
  real(kind=dp)                :: ro_prem,vp_prem,vs_prem
  real(kind=dp)                :: Qmu, Qkappa
  character(len=3), intent(in) :: param !rho, vs,vp

  r = r0 / 1000.

  x_prem = r / 6371.     ! Radius (normalized to x(surface)=1 )

  if(idom==1)then        ! upper crustal layer
     ro_prem = 2.6
     vp_prem = 5.8
     vs_prem = 3.2
     Qmu = 600.0
     Qkappa = 57827.0
  else if(idom==2)then
     ro_prem = 2.9                       ! lower crustal layer
     vp_prem = 6.8
     vs_prem = 3.9
     Qmu = 600.0
     Qkappa = 57827.0
  else if(idom==3)then
     ro_prem = 3.38   ! 2.691  + 0.6924 * x_prem             ! upper mantle
     vp_prem = 8.1106 ! 4.1875 + 3.9382 * x_prem
     vs_prem = 4.4910 ! 2.1519 + 2.3481 * x_prem
     Qmu = 600.0
     Qkappa = 57827.0
  else if(idom==4)then
     ro_prem=  3.435   !7.1089 -  3.8045 * x_prem
     vp_prem=  8.55895 !20.3926 - 12.2569 * x_prem
     vs_prem=  4.64390 !8.9496 -  4.4597 * x_prem
     Qmu = 143.0
     Qkappa = 57827.0
  else if(idom==5)then
     ro_prem = 3.85      !11.2494 -  8.0298 * x_prem
     vp_prem = 9.133917  !39.7027 - 32.6166 * x_prem
     vs_prem = 5.515931  !22.3512 - 18.5856 * x_prem
     Qmu = 143.0
     Qkappa = 57827.0
  else if(idom==6)then
     ro_prem = 3.975    ! 5.3197 - 1.4836 * x_prem
     vp_prem = 10.15776 ! 19.0957 - 9.8672 * x_prem
     vs_prem = 5.516017 ! 9.9839 - 4.9324 * x_prem
     Qmu = 143.0
     Qkappa = 57827.0
  else if(idom==7)then   !lower mantle
     ro_prem = 4.45    ! 7.9565 -  6.4761 * x_prem + 5.5283 * x_prem**2 - 3.0807 * x_prem**3
     vp_prem = 11.0656 ! 29.2766 - 23.6027 * x_prem + 5.5242 * x_prem**2 - 2.5514 * x_prem**3
     vs_prem =  6.2404 ! 22.3459 - 17.2473 * x_prem - 2.0834 * x_prem**2 + 0.9783 * x_prem**3
     Qmu = 312.0
     Qkappa = 57827.0
  else if(idom==8)then
     ro_prem =  4.45    !7.9565 -  6.4761 * x_prem +  5.5283 * x_prem**2 -  3.0807 * x_prem**3
     vp_prem = 11.0656 !24.9520 - 40.4673 * x_prem + 51.4832 * x_prem**2 - 26.6419 * x_prem**3
     vs_prem =  6.2404 !11.1671 - 13.7818 * x_prem + 17.4575 * x_prem**2 -  9.2777 * x_prem**3
     Qmu = 312.0
     Qkappa = 57827.0
  else if(idom==9)then
     ro_prem =   4.45  !7.9565 - 6.4761 * x_prem + 5.5283 * x_prem**2 - 3.0807 * x_prem**3
     vp_prem =  11.0656 !15.3891 - 5.3181 * x_prem + 5.5242 * x_prem**2 - 2.5514 * x_prem**3
     vs_prem =   6.2404 !6.9254 + 1.4672 * x_prem - 2.0834 * x_prem**2 + 0.9783 * x_prem**3
     Qmu = 312.0
     Qkappa = 57827.0
  else if(idom==10)then  ! outer core
     ro_prem = 12.5815 - 1.2638 * x_prem - 3.6426 * x_prem**2 -  5.5281  * x_prem**3
     vp_prem = 11.0487 - 4.0362 * x_prem + 4.8023 * x_prem**2 - 13.5732  * x_prem**3
     vs_prem = 0.0
     Qmu = 0.0
     Qkappa = 57827.0
  else if(idom==11)then                        ! inner core
     ro_prem = 13.0885 - 8.8381 * x_prem**2
     vp_prem = 11.2622 - 6.3640 * x_prem**2
     vs_prem =  3.6678 - 4.4475 * x_prem**2
     Qmu = 84.6
     Qkappa = 1327.7
  endif

 ! homog model
  if (idom < 10) then
     ro_prem = 4.5
     vp_prem = 11.
     vs_prem = 6.
  endif

  if (param=='rho') then
     prem_sub_vm = ro_prem * 1000.
  else if (param=='v_p') then
     prem_sub_vm = vp_prem * 1000.
  else if (param=='v_s') then
     prem_sub_vm = vs_prem * 1000.
  else if (param=='vpv') then
     prem_sub_vm = vp_prem * 1000.
  else if (param=='vsv') then
     prem_sub_vm = vs_prem * 1000.
  else if (param=='vph') then
     prem_sub_vm = vp_prem * 1000.
  else if (param=='vsh') then
     prem_sub_vm = vs_prem * 1000.
  else if (param=='eta') then
     prem_sub_vm = 1.
  else if (param=='Qmu') then
     prem_sub_vm = Qmu
  else if (param=='Qka') then
     prem_sub_vm = Qkappa
  else
     write(6,*)'ERROR IN PREM_SUB FUNCTION:',param,'NOT AN OPTION'
     stop
  endif

end function prem_sub_vm
!=============================================================================
end module background_models
