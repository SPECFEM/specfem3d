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

  subroutine prem_iso(myrank,x,rho,drhodr,vp,vs,Qkappa,Qmu,idoubling,CRUSTAL, &
      ONE_CRUST,check_doubling_flag,RICB,RCMB,RTOPDDOUBLEPRIME, &
      R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN)

  implicit none

  include "constants.h"

! given a normalized radius x, gives the non-dimesionalized density rho,
! speeds vp and vs, and the quality factors Qkappa and Qmu

  logical CRUSTAL,ONE_CRUST,check_doubling_flag

  integer idoubling,myrank

  double precision x,rho,drhodr,vp,vs,Qkappa,Qmu,RICB,RCMB,RTOPDDOUBLEPRIME, &
      R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN

  double precision r

! compute real physical radius in meters
  r = x * R_EARTH

! check flags to make sure we correctly honor the discontinuities
! we use strict inequalities since r has been slighly changed in mesher

 if(check_doubling_flag) then

!
!--- inner core
!

  if(r >= 0.d0 .and. r < RICB) then
    if(idoubling /= IFLAG_INNER_CORE_NORMAL .and. &
       idoubling /= IFLAG_MIDDLE_CENTRAL_CUBE .and. &
       idoubling /= IFLAG_BOTTOM_CENTRAL_CUBE .and. &
       idoubling /= IFLAG_TOP_CENTRAL_CUBE .and. &
       idoubling /= IFLAG_IN_FICTITIOUS_CUBE) &
         call exit_MPI(myrank,'wrong doubling flag for inner core point')
!
!--- outer core
!
  else if(r > RICB .and. r < RCMB) then
    if(idoubling /= IFLAG_OUTER_CORE_NORMAL) &
      call exit_MPI(myrank,'wrong doubling flag for outer core point')
!
!--- D" at the base of the mantle
!
  else if(r > RCMB .and. r < RTOPDDOUBLEPRIME) then
    if(idoubling /= IFLAG_MANTLE_NORMAL) &
      call exit_MPI(myrank,'wrong doubling flag for D" point')
!
!--- mantle: from top of D" to d670
!
  else if(r > RTOPDDOUBLEPRIME .and. r < R670) then
    if(idoubling /= IFLAG_MANTLE_NORMAL) &
      call exit_MPI(myrank,'wrong doubling flag for top D" -> d670 point')

!
!--- mantle: from d670 to d220
!
  else if(r > R670 .and. r < R220) then
    if(idoubling /= IFLAG_670_220) &
      call exit_MPI(myrank,'wrong doubling flag for d670 -> d220 point')

!
!--- mantle and crust: from d220 to MOHO and then to surface
!
  else if(r > R220) then
    if(idoubling /= IFLAG_220_80 .and. idoubling /= IFLAG_80_MOHO .and. idoubling /= IFLAG_CRUST) &
      call exit_MPI(myrank,'wrong doubling flag for d220 -> Moho -> surface point')

  endif

  endif

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
  if(CRUSTAL .and. .not. SUPPRESS_CRUSTAL_MESH) then
! fill with PREM mantle and later add CRUST2.0
    if(r > R80) then
      drhodr=0.6924d0
      rho=2.6910d0+0.6924d0*x
      vp=4.1875d0+3.9382d0*x
      vs=2.1519d0+2.3481d0*x
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


    else if (SUPPRESS_CRUSTAL_MESH) then
!! DK DK extend the Moho up to the surface instead of the crust
      drhodr=0.6924d0
      rho = 2.6910d0+0.6924d0*(RMOHO / R_EARTH)
      vp = 4.1875d0+3.9382d0*(RMOHO / R_EARTH)
      vs = 2.1519d0+2.3481d0*(RMOHO / R_EARTH)
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
  rho=rho*1000.0d0
  vp=vp*1000.0d0
  vs=vs*1000.0d0

  end subroutine prem_iso

!
!=====================================================================
!

  subroutine prem_aniso(myrank,x,rho,vpv,vph,vsv,vsh,eta_aniso,Qkappa,Qmu, &
      idoubling,CRUSTAL,ONE_CRUST,RICB,RCMB,RTOPDDOUBLEPRIME, &
      R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN)

  implicit none

  include "constants.h"

! given a normalized radius x, gives the non-dimesionalized density rho,
! speeds vp and vs, and the quality factors Qkappa and Qmu

  logical CRUSTAL,ONE_CRUST

  integer idoubling,myrank

  double precision x,rho,Qkappa,Qmu,vpv,vph,vsv,vsh,eta_aniso,RICB,RCMB, &
      RTOPDDOUBLEPRIME,R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN

  double precision r
  double precision scaleval

! compute real physical radius in meters
  r = x * R_EARTH

! check flags to make sure we correctly honor the discontinuities
! we use strict inequalities since r has been slighly changed in mesher

!
!--- inner core
!
  if(r >= 0.d0 .and. r < RICB) then
    if(idoubling /= IFLAG_INNER_CORE_NORMAL .and. &
       idoubling /= IFLAG_MIDDLE_CENTRAL_CUBE .and. &
       idoubling /= IFLAG_BOTTOM_CENTRAL_CUBE .and. &
       idoubling /= IFLAG_TOP_CENTRAL_CUBE .and. &
       idoubling /= IFLAG_IN_FICTITIOUS_CUBE) &
         call exit_MPI(myrank,'wrong doubling flag for inner core point')
!
!--- outer core
!
  else if(r > RICB .and. r < RCMB) then
    if(idoubling /= IFLAG_OUTER_CORE_NORMAL) &
      call exit_MPI(myrank,'wrong doubling flag for outer core point')
!
!--- D" at the base of the mantle
!
  else if(r > RCMB .and. r < RTOPDDOUBLEPRIME) then
    if(idoubling /= IFLAG_MANTLE_NORMAL) &
      call exit_MPI(myrank,'wrong doubling flag for D" point')
!
!--- mantle: from top of D" to d670
!
  else if(r > RTOPDDOUBLEPRIME .and. r < R670) then
    if(idoubling /= IFLAG_MANTLE_NORMAL) &
      call exit_MPI(myrank,'wrong doubling flag for top D" -> d670 point')

!
!--- mantle: from d670 to d220
!
  else if(r > R670 .and. r < R220) then
    if(idoubling /= IFLAG_670_220) &
      call exit_MPI(myrank,'wrong doubling flag for d670 -> d220 point')

!
!--- mantle and crust: from d220 to MOHO and then to surface
!
  else if(r > R220) then
    if(idoubling /= IFLAG_220_80 .and. idoubling /= IFLAG_80_MOHO .and. idoubling /= IFLAG_CRUST) &
      call exit_MPI(myrank,'wrong doubling flag for d220 -> Moho -> surface point')

  endif

! no anisotropy by default
  eta_aniso = 1.d0

!
!--- inner core
!
  if(r >= 0.d0 .and. r <= RICB) then
    rho=13.0885d0-8.8381d0*x*x
    vpv=11.2622d0-6.3640d0*x*x
    vsv=3.6678d0-4.4475d0*x*x
    vph=vpv
    vsh=vsv
    Qmu=84.6d0
    Qkappa=1327.7d0
!
!--- outer core
!
  else if(r > RICB .and. r <= RCMB) then
    rho=12.5815d0-1.2638d0*x-3.6426d0*x*x-5.5281d0*x*x*x
    vpv=11.0487d0-4.0362d0*x+4.8023d0*x*x-13.5732d0*x*x*x
    vsv=0.0d0
    vph=vpv
    vsh=vsv
    Qmu=0.0d0
    Qkappa=57827.0d0
!
!--- D" at the base of the mantle
!
  else if(r > RCMB .and. r <= RTOPDDOUBLEPRIME) then
    rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
    vpv=15.3891d0-5.3181d0*x+5.5242d0*x*x-2.5514d0*x*x*x
    vsv=6.9254d0+1.4672d0*x-2.0834d0*x*x+0.9783d0*x*x*x
    vph=vpv
    vsh=vsv
    Qmu=312.0d0
    Qkappa=57827.0d0
!
!--- mantle: from top of D" to d670
!
  else if(r > RTOPDDOUBLEPRIME .and. r <= R771) then
    rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
    vpv=24.9520d0-40.4673d0*x+51.4832d0*x*x-26.6419d0*x*x*x
    vsv=11.1671d0-13.7818d0*x+17.4575d0*x*x-9.2777d0*x*x*x
    vph=vpv
    vsh=vsv
    Qmu=312.0d0
    Qkappa=57827.0d0
  else if(r > R771 .and. r <= R670) then
    rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
    vpv=29.2766d0-23.6027d0*x+5.5242d0*x*x-2.5514d0*x*x*x
    vsv=22.3459d0-17.2473d0*x-2.0834d0*x*x+0.9783d0*x*x*x
    vph=vpv
    vsh=vsv
    Qmu=312.0d0
    Qkappa=57827.0d0
!
!--- mantle: above d670
!
  else if(r > R670 .and. r <= R600) then
    rho=5.3197d0-1.4836d0*x
    vpv=19.0957d0-9.8672d0*x
    vsv=9.9839d0-4.9324d0*x
    vph=vpv
    vsh=vsv
    Qmu=143.0d0
    Qkappa=57827.0d0
  else if(r > R600 .and. r <= R400) then
    rho=11.2494d0-8.0298d0*x
    vpv=39.7027d0-32.6166d0*x
    vsv=22.3512d0-18.5856d0*x
    vph=vpv
    vsh=vsv
    Qmu=143.0d0
    Qkappa=57827.0d0
  else if(r > R400 .and. r <= R220) then
    rho=7.1089d0-3.8045d0*x
    vpv=20.3926d0-12.2569d0*x
    vsv=8.9496d0-4.4597d0*x
    vph=vpv
    vsh=vsv
    Qmu=143.0d0
    Qkappa=57827.0d0
  else if(r > R220 .and. r <= R80) then

! anisotropy in PREM only above 220 km

    rho=2.6910d0+0.6924d0*x
    vpv=0.8317d0+7.2180d0*x
    vph=3.5908d0+4.6172d0*x
    vsv=5.8582d0-1.4678d0*x
    vsh=-1.0839d0+5.7176d0*x
    eta_aniso=3.3687d0-2.4778d0*x
    Qmu=80.0d0
    Qkappa=57827.0d0

  else
  if(CRUSTAL) then
! fill with PREM mantle and later add CRUST2.0
    if(r > R80) then
      rho=2.6910d0+0.6924d0*x
      vpv=0.8317d0+7.2180d0*x
      vph=3.5908d0+4.6172d0*x
      vsv=5.8582d0-1.4678d0*x
      vsh=-1.0839d0+5.7176d0*x
      eta_aniso=3.3687d0-2.4778d0*x
      Qmu=600.0d0
      Qkappa=57827.0d0
    endif
  else
! use PREM crust
    if(r > R80 .and. r <= RMOHO) then

! anisotropy in PREM only above 220 km

      rho=2.6910d0+0.6924d0*x
      vpv=0.8317d0+7.2180d0*x
      vph=3.5908d0+4.6172d0*x
      vsv=5.8582d0-1.4678d0*x
      vsh=-1.0839d0+5.7176d0*x
      eta_aniso=3.3687d0-2.4778d0*x
      Qmu=600.0d0
      Qkappa=57827.0d0

! no anisotropy in the crust in PREM

    else if(r > RMOHO .and. r <= RMIDDLE_CRUST) then
      rho=2.9d0
      vpv=6.8d0
      vsv=3.9d0
      vph=vpv
      vsh=vsv
      Qmu=600.0d0
      Qkappa=57827.0d0

! same properties everywhere in PREM crust (only one layer in the crust)
      if(ONE_CRUST) then
        rho=2.6d0
        vpv=5.8d0
        vsv=3.2d0
        vph=vpv
        vsh=vsv
        Qmu=600.0d0
        Qkappa=57827.0d0
      endif

    else if(r > RMIDDLE_CRUST .and. r <= ROCEAN) then
      rho=2.6d0
      vpv=5.8d0
      vsv=3.2d0
      vph=vpv
      vsh=vsv
      Qmu=600.0d0
      Qkappa=57827.0d0
    else if(r > ROCEAN) then
      rho=2.6d0
      vpv=5.8d0
      vsv=3.2d0
      vph=vpv
      vsh=vsv
      Qmu=600.0d0
      Qkappa=57827.0d0
    endif
  endif
  endif

! non-dimensionalize
! time scaling (s^{-1}) is done with scaleval
! do not scale anisotropy parameter eta_aniso, which is dimensionless
  scaleval=dsqrt(PI*GRAV*RHOAV)
  rho=rho*1000.0d0/RHOAV
  vpv=vpv*1000.0d0/(R_EARTH*scaleval)
  vsv=vsv*1000.0d0/(R_EARTH*scaleval)
  vph=vph*1000.0d0/(R_EARTH*scaleval)
  vsh=vsh*1000.0d0/(R_EARTH*scaleval)

  end subroutine prem_aniso

!
!=====================================================================
!

  subroutine prem_display_outer_core(myrank,x,rho,vp,vs,Qkappa,Qmu,idoubling)

! routine used for AVS or DX display of stability condition
! and number of points per wavelength only in the fluid outer core

  implicit none

  include "constants.h"

! given a normalized radius x, gives the non-dimesionalized density rho,
! speeds vp and vs, and the quality factors Qkappa and Qmu

  integer idoubling,myrank
  double precision x,rho,vp,vs,Qkappa,Qmu

  double precision scaleval

  if(idoubling /= IFLAG_OUTER_CORE_NORMAL) call exit_MPI(myrank,'wrong doubling flag for outer core point')

!
!--- outer core
!
  rho=12.5815d0-1.2638d0*x-3.6426d0*x*x-5.5281d0*x*x*x
  vp=11.0487d0-4.0362d0*x+4.8023d0*x*x-13.5732d0*x*x*x
  vs=0.0d0
  Qmu=0.0d0
  Qkappa=57827.0d0

! non-dimensionalize
! time scaling (s^{-1}) is done with scaleval
  scaleval = dsqrt(PI*GRAV*RHOAV)
  rho = rho*1000.0d0/RHOAV
  vp = vp*1000.0d0/(R_EARTH*scaleval)
  vs = vs*1000.0d0/(R_EARTH*scaleval)

  end subroutine prem_display_outer_core

!
!=====================================================================
!

  subroutine prem_density(x,rho,ONE_CRUST,RICB,RCMB,RTOPDDOUBLEPRIME, &
      R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN)

  implicit none

  include "constants.h"

  double precision x,rho,RICB,RCMB,RTOPDDOUBLEPRIME, &
      R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN

  logical ONE_CRUST

  double precision r

  r = x * R_EARTH

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

