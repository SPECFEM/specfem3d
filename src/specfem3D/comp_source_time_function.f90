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

  double precision function comp_source_time_function(t,hdur)

  implicit none

  include "constants.h"

  double precision t,hdur

  double precision, external :: netlib_specfun_erf

  ! quasi Heaviside, small Gaussian moment-rate tensor with hdur
  comp_source_time_function = 0.5d0*(1.0d0 + netlib_specfun_erf(t/hdur))

  end function comp_source_time_function


!
!-------------------------------------------------------------------------------------------------
!

  double precision function comp_source_time_function_gauss(t,hdur)

  implicit none

  include "constants.h"

  double precision :: t,hdur
  double precision :: hdur_decay
  double precision,parameter :: SOURCE_DECAY_STRONG = 2.0d0/SOURCE_DECAY_MIMIC_TRIANGLE

  ! note: hdur given is hdur_gaussian = hdur/SOURCE_DECAY_MIMIC_TRIANGLE
  !           and SOURCE_DECAY_MIMIC_TRIANGLE ~ 1.68
  hdur_decay = hdur

  ! this here uses a stronger gaussian decay rate (empirical value) to avoid non-zero onset times;
  ! however, it should mimik a triangle source time function...
  !hdur_decay = hdur  / SOURCE_DECAY_STRONG

  ! note: a nonzero time to start the simulation with would lead to more high-frequency noise
  !          due to the (spatial) discretization of the point source on the mesh

  ! gaussian
  comp_source_time_function_gauss = exp(-(t/hdur_decay)**2)/(sqrt(PI)*hdur_decay)

  end function comp_source_time_function_gauss

!
!-------------------------------------------------------------------------------------------------
!

  double precision function comp_source_time_function_rickr(t,f0)

  implicit none

  include "constants.h"

  double precision t,f0

  ! ricker
  comp_source_time_function_rickr = (1.d0 - 2.d0*PI*PI*f0*f0*t*t ) &
                                    * exp( -PI*PI*f0*f0*t*t )

  !!! another source time function they have called 'ricker' in some old papers,
  !!! e.g., 'Finite-Frequency Kernels Based on Adjoint Methods' by Liu & Tromp, BSSA (2006)
  !!! in order to benchmark those simulations, the following formula is needed.
  ! comp_source_time_function_rickr = -2.d0*PI*PI*f0*f0*f0*t * exp(-PI*PI*f0*f0*t*t)

  end function comp_source_time_function_rickr

