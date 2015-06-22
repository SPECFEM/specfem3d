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

  subroutine count_number_of_sources(NSOURCES)

! count the total number of sources in the CMTSOLUTION file
! there are NLINES_PER_CMTSOLUTION_SOURCE lines per source in that file

  implicit none

  include "constants.h"

  integer, intent(out) :: NSOURCES

  integer ios,icounter

  character(len=150) CMTSOLUTION,dummystring

  call get_value_string(CMTSOLUTION, 'solver.CMTSOLUTION', 'DATA/CMTSOLUTION')

  open(unit=1,file=CMTSOLUTION,iostat=ios,status='old',action='read')
  if(ios /= 0) stop 'error opening CMTSOLUTION file'
  icounter = 0
  do while(ios == 0)
    read(1,"(a)",iostat=ios) dummystring
    if(ios == 0) icounter = icounter + 1
  enddo
  close(1)

  if(mod(icounter,NLINES_PER_CMTSOLUTION_SOURCE) /= 0) &
    stop 'total number of lines in CMTSOLUTION file should be a multiple of NLINES_PER_CMTSOLUTION_SOURCE'

  NSOURCES = icounter / NLINES_PER_CMTSOLUTION_SOURCE

  if(NSOURCES < 1) stop 'need at least one source in CMTSOLUTION file'

  end subroutine count_number_of_sources

