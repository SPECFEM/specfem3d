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

subroutine stretching_function(r_top,r_bottom,ner,stretch_tab)

! define stretch_tab which contains r_top and r_bottom for each element layer in the crust for 3D models.

  implicit none

  include "constants.h"

  double precision :: r_top, r_bottom,value
  integer :: ner,i
  double precision, dimension (2,ner) :: stretch_tab
! for increasing execution speed but have less precision in stretching, increase step
! not very effective algorithm, but sufficient : used once per proc for meshing.
  double precision, parameter :: step = 0.001

! initialize array
  do i=1,ner
    stretch_tab(2,i)=(1.d0/ner)
  enddo

! fill with ratio of the layer one thickness for each element
  do while((stretch_tab(2,1) / stretch_tab(2,ner)) > MAX_RATIO_CRUST_STRETCHING)
    if (modulo(ner,2) /= 0) then
      value = -floor(ner/2.d0)*step
    else
      value = (0.5d0-floor(ner/2.d0))*step
    endif
    do i=1,ner
      stretch_tab(2,i) = stretch_tab(2,i) + value
      value = value + step
    enddo
  enddo

! deduce r_top and r_bottom
  ! r_top
  stretch_tab(1,1) = r_top
  do i=2,ner
    stretch_tab(1,i) = sum(stretch_tab(2,i:ner))*(r_top-r_bottom) + r_bottom
  enddo

  ! r_bottom
  stretch_tab(2,ner) = r_bottom
  do i=1,ner-1
    stretch_tab(2,i) = stretch_tab(1,i+1)
  enddo

end subroutine stretching_function

