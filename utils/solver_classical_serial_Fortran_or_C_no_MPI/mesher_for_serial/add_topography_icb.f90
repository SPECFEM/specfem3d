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

  subroutine add_topography_icb(myrank,xelm,yelm,zelm,RICB,RCMB)

  implicit none

  include "constants.h"

  integer myrank

  double precision xelm(NGNOD)
  double precision yelm(NGNOD)
  double precision zelm(NGNOD)

  double precision RICB,RCMB

  integer ia

  double precision topoicb

  double precision r,theta,phi
  double precision gamma

! we loop on all the points of the element
  do ia = 1,NGNOD

! convert to r theta phi
    call xyz_2_rthetaphi_dble(xelm(ia),yelm(ia),zelm(ia),r,theta,phi)
    call reduce(theta,phi)

! compute topography on ICB; the routine subtopo_icb needs to be supplied by the user
!    call subtopo_icb(theta,phi,topoicb)
    topoicb = 0.0d0

! non-dimensionalize the topography, which is in km
! positive for a depression, so change the sign for a perturbation in radius
    topoicb = -topoicb / R_EARTH_KM

    gamma = 0.0d0
    if(r > 0.0d0 .and. r <= RICB/R_EARTH) then
! stretching between center and RICB
      gamma = r/(RICB/R_EARTH)
    elseif(r>= RICB/R_EARTH .and. r <= RCMB/R_EARTH) then
! stretching between RICB and RCMB
      gamma = (r - RCMB/R_EARTH) / (RICB/R_EARTH - RCMB/R_EARTH)
    endif
    if(gamma < -0.0001 .or. gamma > 1.0001) call exit_MPI(myrank,'incorrect value of gamma for CMB topography')

    xelm(ia) = xelm(ia)*(ONE + gamma * topoicb / r)
    yelm(ia) = yelm(ia)*(ONE + gamma * topoicb / r)
    zelm(ia) = zelm(ia)*(ONE + gamma * topoicb / r)

  enddo

  end subroutine add_topography_icb

