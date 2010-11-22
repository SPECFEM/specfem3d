!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 2
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and University of Pau / CNRS / INRIA
!         (c) California Institute of Technology July 2004
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

  double precision, parameter :: PI = 3.141592653589793d0

! small perturbation in distance (degrees)
! (0.0025 degrees corresponds to roughly 250 meters)
  double precision, parameter :: DDELTA = 0.0025

! small perturbation in depth (km)
  double precision, parameter :: DDEPTH = 1.0

! basic scalar moment
  double precision, parameter :: MOMENT = 1.0e+22

! maximum record length
  integer, parameter :: NDATAMAX = 25000

! maximum number of records
  integer, parameter :: NRECORDSMAX = 250

! number of inversion parameters:
! 6 for moment tensor only,
! 7 for moment tensor plus centroid time,
! 9 for moment tensor plus location,
! 10 for moment tensor plus location and centroid time,
  integer, parameter :: NPAR = 7

! weigh the data in the inversion using a NOISE_WINDOW before P
  logical, parameter :: WEIGHTING = .false.

! length of noise window before P in seconds if WEIGHTING is set to true
  double precision, parameter :: NOISE_WINDOW = 10.

! select data
  logical, parameter :: SELECT_DATA = .false.

! show the fit to the data before inversion (after offset removal)
  logical, parameter :: GRAPHICS = .true.

! rotate to get transverse and radial components
  logical, parameter :: ROTATE = .true.

! inversion for a zero-trace moment tensor
  logical, parameter :: ZERO_TRACE_INVERSION = .true.

