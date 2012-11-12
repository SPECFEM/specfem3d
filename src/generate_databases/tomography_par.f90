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

!--------------------------------------------------------------------------------------------------
!
! generic tomography file
!
! note: the idea is to use external tomography velocity models
!
! most of the routines here are place-holders, please add/implement your own routines
!
!--------------------------------------------------------------------------------------------------

  module tomography_par

  include "constants.h"

  ! for external tomography:
  ! (regular spaced, xyz-block file in ascii)

  ! number of external tomographic models
  integer :: NFILES_TOMO

  ! models dimensions
  double precision  :: END_X,END_Y,END_Z

  double precision, dimension(:), allocatable :: ORIG_X,ORIG_Y,ORIG_Z
  double precision, dimension(:), allocatable :: SPACING_X,SPACING_Y,SPACING_Z

  ! models parameter records
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: vp_tomography,vs_tomography,rho_tomography,z_tomography

  ! models entries
  integer, dimension(:), allocatable :: NX,NY,NZ
  integer, dimension(:), allocatable :: nrecord

  ! min/max statistics
  double precision, dimension(:), allocatable :: VP_MIN,VS_MIN,RHO_MIN,VP_MAX,VS_MAX,RHO_MAX

  ! process rank
  integer :: myrank_tomo

  end module tomography_par
