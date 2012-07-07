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

! computes and displays the min or max angle in a mesh that give rise
! to a certain value of skewness

! Dimitri Komatitsch, University of Pau, France, February 2009.

  program convert_skewness_to_angle

  implicit none

  include "constants.h"

  double precision :: equiangle_skewness

  print *
  print *,'enter equiangle skewness to convert (0.=perfect, 1.=bad):'
  read(5,*) equiangle_skewness

  print *
  print *,'equiangle skewness = ',equiangle_skewness
  print *
  print *,'deviation angle from a right angle (90 degrees) is therefore = ',90.*equiangle_skewness
  print *
  print *,'worst angle in the mesh is therefore ',90.*(1. - equiangle_skewness)
  print *,'or ',180. - 90.*(1. - equiangle_skewness),' degrees'
  print *

  end program convert_skewness_to_angle

