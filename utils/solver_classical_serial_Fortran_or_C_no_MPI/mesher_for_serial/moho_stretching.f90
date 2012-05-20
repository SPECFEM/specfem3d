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

  subroutine moho_stretching(myrank,xelm,yelm,zelm,RMOHO,R220)

  implicit none

  include "constants.h"

! ocean-continent function maximum spherical harmonic degree
  integer, parameter :: NL_OCEAN_CONTINENT = 12

! spherical harmonic coefficients of the ocean-continent function (km)
  double precision A_lm(0:NL_OCEAN_CONTINENT,0:NL_OCEAN_CONTINENT),B_lm(0:NL_OCEAN_CONTINENT,0:NL_OCEAN_CONTINENT)

  common /smooth_moho/ A_lm,B_lm

  integer myrank

  double precision xelm(NGNOD)
  double precision yelm(NGNOD)
  double precision zelm(NGNOD)

  double precision RMOHO,R220

  integer ia

  integer l,m
  double precision r,theta,phi
  double precision sint,cost,x(2*NL_OCEAN_CONTINENT+1),dx(2*NL_OCEAN_CONTINENT+1)
  double precision elevation
  double precision gamma

! we loop on all the points of the element
  do ia = 1,NGNOD

! convert to r theta phi
    call xyz_2_rthetaphi_dble(xelm(ia),yelm(ia),zelm(ia),r,theta,phi)
    call reduce(theta,phi)

    elevation = 0.0d0
    do l = 0,NL_OCEAN_CONTINENT
      sint = dsin(theta)
      cost = dcos(theta)
      call lgndr(l,cost,sint,x,dx)
      m = 0
      elevation = elevation + A_lm(l,m)*x(m+1)
      do m = 1,l
        elevation = elevation + (A_lm(l,m)*dcos(dble(m)*phi)+B_lm(l,m)*dsin(dble(m)*phi))*x(m+1)
      enddo
    enddo
    elevation = -0.25d0*elevation/R_EARTH_KM

    gamma = 0.0d0
    if(r >= RMOHO/R_EARTH) then
! stretching above the Moho
      gamma = (1.0d0 - r) / (1.0d0 - RMOHO/R_EARTH)
    elseif(r>= R220/R_EARTH .and. r< RMOHO/R_EARTH) then
! stretching between R220 and RMOHO
      gamma = (r - R220/R_EARTH) / (RMOHO/R_EARTH - R220/R_EARTH)
    endif
    if(gamma < -0.0001 .or. gamma > 1.0001) call exit_MPI(myrank,'incorrect value of gamma for Moho topography')

    xelm(ia) = xelm(ia)*(ONE + gamma * elevation / r)
    yelm(ia) = yelm(ia)*(ONE + gamma * elevation / r)
    zelm(ia) = zelm(ia)*(ONE + gamma * elevation / r)

  enddo

  end subroutine moho_stretching

  subroutine read_smooth_moho

  implicit none

! ocean-continent function maximum spherical harmonic degree
  integer, parameter :: NL_OCEAN_CONTINENT = 12

! spherical harmonic coefficients of the ocean-continent function (km)
  double precision A_lm(0:NL_OCEAN_CONTINENT,0:NL_OCEAN_CONTINENT),B_lm(0:NL_OCEAN_CONTINENT,0:NL_OCEAN_CONTINENT)

  common /smooth_moho/ A_lm,B_lm

!  integer l,m
!
! ocean-continent function (km)
!  open(unit=10,file='DATA/ocean_continent_function/ocean_continent_function.txt',status='old',action='read')
!  do l=0,NL_OCEAN_CONTINENT
!    read(10,*) A_lm(l,0),(A_lm(l,m),B_lm(l,m),m=1,l)
!  enddo
!  close(10)

  A_lm(0,0) = -3.8201999E-04
  B_lm(0,0) = 0.
  A_lm(1,0) = 13.88800
  B_lm(1,0) = 0.
  A_lm(1,1) = -15.24000
  B_lm(1,1) = -9.187200
  A_lm(2,0) = 11.21500
  B_lm(2,0) = 0.
  A_lm(2,1) = -6.754500
  B_lm(2,1) = -8.516700
  A_lm(2,2) = -8.327800
  B_lm(2,2) = -5.029200
  A_lm(3,0) = -3.614500
  B_lm(3,0) = 0.
  A_lm(3,1) = 5.394800
  B_lm(3,1) = -0.9220800
  A_lm(3,2) = -10.05100
  B_lm(3,2) = 13.98100
  A_lm(3,3) = -2.711200
  B_lm(3,3) = -13.57100
  A_lm(4,0) = 7.523300
  B_lm(4,0) = 0.
  A_lm(4,1) = 5.156100
  B_lm(4,1) = 2.184400
  A_lm(4,2) = -10.67300
  B_lm(4,2) = 2.640600
  A_lm(4,3) = -7.786300
  B_lm(4,3) = 0.3674500
  A_lm(4,4) = -3.076400
  B_lm(4,4) = 16.83000
  A_lm(5,0) = -9.681000
  B_lm(5,0) = 0.
  A_lm(5,1) = 0.5026800
  B_lm(5,1) = 2.111300
  A_lm(5,2) = -2.931000
  B_lm(5,2) = -4.329000
  A_lm(5,3) = -1.766800
  B_lm(5,3) = -3.621200
  A_lm(5,4) = 16.08200
  B_lm(5,4) = -4.493900
  A_lm(5,5) = -0.3705800
  B_lm(5,5) = -5.574500
  A_lm(6,0) = 4.407900
  B_lm(6,0) = 0.
  A_lm(6,1) = 0.3799000
  B_lm(6,1) = 1.589400
  A_lm(6,2) = -1.886400
  B_lm(6,2) = -0.5686300
  A_lm(6,3) = -0.9816800
  B_lm(6,3) = -5.827800
  A_lm(6,4) = 3.620600
  B_lm(6,4) = -2.713100
  A_lm(6,5) = 1.445600
  B_lm(6,5) = 3.964100
  A_lm(6,6) = 1.167400
  B_lm(6,6) = 2.134100
  A_lm(7,0) = -4.086100
  B_lm(7,0) = 0.
  A_lm(7,1) = 0.5462000
  B_lm(7,1) = -4.488100
  A_lm(7,2) = 3.116400
  B_lm(7,2) = 1.793600
  A_lm(7,3) = 2.594600
  B_lm(7,3) = -2.129100
  A_lm(7,4) = -5.445000
  B_lm(7,4) = 0.5381500
  A_lm(7,5) = -2.178100
  B_lm(7,5) = 1.766700
  A_lm(7,6) = -1.040000
  B_lm(7,6) = -5.541000
  A_lm(7,7) = 1.536500
  B_lm(7,7) = 3.700600
  A_lm(8,0) = -2.562200
  B_lm(8,0) = 0.
  A_lm(8,1) = 0.3736200
  B_lm(8,1) = 1.488000
  A_lm(8,2) = 1.347500
  B_lm(8,2) = 0.5288200
  A_lm(8,3) = -0.8493700
  B_lm(8,3) = -1.626500
  A_lm(8,4) = 0.2423400
  B_lm(8,4) = 4.202800
  A_lm(8,5) = 2.052200
  B_lm(8,5) = 0.6880400
  A_lm(8,6) = 2.838500
  B_lm(8,6) = 2.835700
  A_lm(8,7) = -4.981400
  B_lm(8,7) = -1.883100
  A_lm(8,8) = -1.102800
  B_lm(8,8) = -1.951700
  A_lm(9,0) = -1.202100
  B_lm(9,0) = 0.
  A_lm(9,1) = 1.020300
  B_lm(9,1) = 1.371000
  A_lm(9,2) = -0.3430100
  B_lm(9,2) = 0.8782800
  A_lm(9,3) = -0.4462500
  B_lm(9,3) = -0.3046100
  A_lm(9,4) = 0.7750700
  B_lm(9,4) = 2.351600
  A_lm(9,5) = -2.092600
  B_lm(9,5) = -2.377100
  A_lm(9,6) = 0.3126900
  B_lm(9,6) = 4.996000
  A_lm(9,7) = -2.284000
  B_lm(9,7) = 1.183700
  A_lm(9,8) = 1.445900
  B_lm(9,8) = 1.080000
  A_lm(9,9) = 1.146700
  B_lm(9,9) = 1.457800
  A_lm(10,0) = -2.516900
  B_lm(10,0) = 0.
  A_lm(10,1) = -0.9739500
  B_lm(10,1) = -0.7195500
  A_lm(10,2) = -2.846000
  B_lm(10,2) = -1.464700
  A_lm(10,3) = 2.720100
  B_lm(10,3) = 0.8241400
  A_lm(10,4) = -1.247800
  B_lm(10,4) = 1.220300
  A_lm(10,5) = -1.638500
  B_lm(10,5) = -1.099500
  A_lm(10,6) = 3.043000
  B_lm(10,6) = -1.976400
  A_lm(10,7) = -1.007300
  B_lm(10,7) = -1.604900
  A_lm(10,8) = 0.6620500
  B_lm(10,8) = -1.135000
  A_lm(10,9) = -3.576800
  B_lm(10,9) = 0.5554900
  A_lm(10,10) = 2.418700
  B_lm(10,10) = -1.482200
  A_lm(11,0) = 0.7158800
  B_lm(11,0) = 0.
  A_lm(11,1) = -3.694800
  B_lm(11,1) = 0.8491400
  A_lm(11,2) = 9.3208998E-02
  B_lm(11,2) = -1.276000
  A_lm(11,3) = 1.575600
  B_lm(11,3) = 0.1972100
  A_lm(11,4) = 0.8989600
  B_lm(11,4) = -1.063000
  A_lm(11,5) = -0.6301000
  B_lm(11,5) = -1.329400
  A_lm(11,6) = 1.389000
  B_lm(11,6) = 1.184100
  A_lm(11,7) = 0.5640700
  B_lm(11,7) = 2.286200
  A_lm(11,8) = 1.530300
  B_lm(11,8) = 0.7677500
  A_lm(11,9) = 0.8495500
  B_lm(11,9) = 0.7247500
  A_lm(11,10) = 2.106800
  B_lm(11,10) = 0.6588000
  A_lm(11,11) = 0.6067800
  B_lm(11,11) = 0.1366800
  A_lm(12,0) = -2.598700
  B_lm(12,0) = 0.
  A_lm(12,1) = -1.150500
  B_lm(12,1) = -0.8425700
  A_lm(12,2) = -0.1593300
  B_lm(12,2) = -1.241400
  A_lm(12,3) = 1.508600
  B_lm(12,3) = 0.3385500
  A_lm(12,4) = -1.941200
  B_lm(12,4) = 1.120000
  A_lm(12,5) = -0.4630500
  B_lm(12,5) = -6.4753003E-02
  A_lm(12,6) = 0.8967000
  B_lm(12,6) = 4.7417998E-02
  A_lm(12,7) = 4.5407999E-02
  B_lm(12,7) = 0.8876400
  A_lm(12,8) = -2.444400
  B_lm(12,8) = 1.172500
  A_lm(12,9) = -2.593400
  B_lm(12,9) = 0.1703700
  A_lm(12,10) = 0.5662700
  B_lm(12,10) = 0.7050800
  A_lm(12,11) = -0.1930000
  B_lm(12,11) = -2.008100
  A_lm(12,12) = -3.187900
  B_lm(12,12) = -1.672000

  end subroutine read_smooth_moho

