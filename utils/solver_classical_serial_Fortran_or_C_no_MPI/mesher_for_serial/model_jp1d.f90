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

subroutine model_jp1d(myrank,x,rho,vp,vs,Qkappa,Qmu,idoubling, &
     check_doubling_flag,RICB,RCMB,RTOPDDOUBLEPRIME, &
     R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST)

  implicit none

  include "constants.h"

  ! given a normalized radius x, gives the non-dimesionalized density rho,
! speeds vp and vs, and the quality factors Qkappa and Qmu

  logical check_doubling_flag
  integer idoubling,myrank

  double precision x,rho,vp,vs,Qkappa,Qmu,RICB,RCMB,RTOPDDOUBLEPRIME, &
      R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST

  double precision r
  double precision scaleval

! compute real physical radius in meters
  r = x * R_EARTH

! check flags to make sure we correctly honor the discontinuities
! we use strict inequalities since r has been slighly changed in mesher

  if(check_doubling_flag) then

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
  if (r >= 0.d0 .and. r <= RICB) then
     rho=13.0885d0-8.8381d0*x*x
     vp=11.24094-4.09689*x**2
     vs=3.56454-3.45241*x**2
     Qmu=84.6d0
     Qkappa=1327.7d0
!
!--- outer core
!
  else if (r > RICB .and. r <= RCMB) then
     rho=12.5815d0-1.2638d0*x-3.6426d0*x*x-5.5281d0*x*x*x
     vp=10.03904+3.75665*x-13.67046*x**2
     vs=0.0d0
     Qmu=0.0d0
     Qkappa=57827.0d0
!
!--- D" at the base of the mantle
!
  else if (r > RCMB .and. r <= RTOPDDOUBLEPRIME) then
     rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
     vp=14.49470-1.47089*x
     vs=8.16616-1.58206*x
     Qmu=312.0d0
     Qkappa=57827.0d0
!
!--- mantle: from top of D" to d670
!
  else if(r > RTOPDDOUBLEPRIME .and. r <= R771) then
     rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
     vp=-355.58324*x**4 + 1002.03178*x**3 - 1057.3873425*x**2 + 487.0891011*x - 68.520645
     vs=-243.33862*x**4 + 668.06411*x**3 - 685.20113*x**2 + 308.04893*x - 43.737642
     Qmu=312.0d0
     Qkappa=57827.0d0
  else if(r > R771 .and. r <= R670) then
     rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
     vp=-174.468866*x**2 + 286.37769*x - 106.034798
     vs=-81.0865*x*x + 129.67095*x - 45.268933
     Qmu=312.0d0
     Qkappa=57827.0d0
!
!--- mantle: above d670
!
  else if(r > R670 .and. r <= 5871000.d0) then
     vp=-300.510146*x*x  + 511.17372648*x - 206.265832
     vs=-139.78275*x*x + 233.3097462*x - 91.0129372
     rho=3.3d0 + (vs-4.4d0)*0.7d0
     Qmu=143.0d0
     Qkappa=57827.0d0

  else if(r > 5871000.d0 .and. r <= R400) then
     vp=-601.0202917*x*x + 1063.3823*x - 459.9388738
     vs=-145.2465705*x*x + 243.2807524*x - 95.561877
     rho=3.3d0 + (vs - 4.4d0)*0.7d0
     Qmu=143.0d0
     Qkappa=57827.0d0

  else if(r > R400 .and. r <= R220) then
     vp=25.042512155*x*x - 68.8367583*x + 51.4120272
     vs=15.540158021*x*x - 40.2087657*x + 28.9578929
     rho=3.3d0 + (vs - 4.4d0)*0.7d0
     Qmu=143.0d0
     Qkappa=57827.0d0

  else if(r > R220 .and. r <= R80) then
     vp=27.0989608 - 19.473338*x
     vs=13.920596 - 9.6309917*x
     rho=3.3d0 + (vs - 4.4d0)*0.7d0
     Qmu=80.0d0
     Qkappa=57827.0d0

  else if(r > R80 .and. r <= RMOHO) then
     vp=26.7663028 - 19.13645*x
     vs=13.4601434 - 9.164683*x
     rho=3.3d0 + (vs - 4.4d0)*0.7d0
     Qmu=600.0d0
     Qkappa=57827.0d0

  else if(r > RMOHO .and. r <= RMIDDLE_CRUST) then
     rho=2.9d0
     vp = 6.7d0
     vs = 3.8d0
     Qmu=600.0d0
     Qkappa=57827.0d0
  else
     rho=2.6d0
     vp = 6.0d0
     vs = 3.5d0
     Qmu=600.0d0
     Qkappa=57827.0d0
  end if


! non-dimensionalize
! time scaling (s^{-1}) is done with scaleval
  scaleval=dsqrt(PI*GRAV*RHOAV)
  rho=rho*1000.0d0/RHOAV
  vp=vp*1000.0d0/(R_EARTH*scaleval)
  vs=vs*1000.0d0/(R_EARTH*scaleval)
end subroutine model_jp1d
