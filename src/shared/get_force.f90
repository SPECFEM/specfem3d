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

  subroutine get_force(tshift_force,hdur,lat,long,depth,NSOURCES,min_tshift_force_original,factor_force_source, &
                      comp_dir_vect_source_E,comp_dir_vect_source_N,comp_dir_vect_source_Z_UP)

  implicit none

  include "constants.h"

!--- input or output arguments of the subroutine below

  integer, intent(in) :: NSOURCES

  double precision, intent(out) :: min_tshift_force_original
  double precision, dimension(NSOURCES), intent(out) :: tshift_force,hdur,lat,long,depth,factor_force_source
  double precision, dimension(NSOURCES), intent(out) :: comp_dir_vect_source_E
  double precision, dimension(NSOURCES), intent(out) :: comp_dir_vect_source_N
  double precision, dimension(NSOURCES), intent(out) :: comp_dir_vect_source_Z_UP

!--- local variables below

  integer isource,dummyval
  double precision t_shift(NSOURCES)
  double precision length
  character(len=7) dummy
  character(len=256) string, FORCESOLUTION

  ! initializes
  lat(:) = 0.d0
  long(:) = 0.d0
  depth(:) = 0.d0
  t_shift(:) = 0.d0
  tshift_force(:) = 0.d0
  hdur(:) = 0.d0
  factor_force_source(:) = 0.d0
  comp_dir_vect_source_E(:) = 0.d0
  comp_dir_vect_source_N(:) = 0.d0
  comp_dir_vect_source_Z_UP(:) = 0.d0
!
!---- read info
!
  call get_value_string(FORCESOLUTION, 'solver.FORCESOLUTION', &
       IN_DATA_FILES_PATH(1:len_trim(IN_DATA_FILES_PATH))//'FORCESOLUTION')

  open(unit=1,file=trim(FORCESOLUTION),status='old',action='read')

! read source number isource
  do isource=1,NSOURCES

    read(1,"(a256)") string
    ! skips empty lines
    do while( len_trim(string) == 0 )
      read(1,"(a256)") string
    enddo

    ! read header with event information
    read(string,"(a6,i4)") dummy,dummyval

    ! read time shift
    read(1,"(a)") string
    read(string(12:len_trim(string)),*) t_shift(isource)

    ! read half duration
    read(1,"(a)") string
    read(string(6:len_trim(string)),*) hdur(isource)

    ! read latitude
    read(1,"(a)") string
    read(string(10:len_trim(string)),*) lat(isource)

    ! read longitude
    read(1,"(a)") string
    read(string(11:len_trim(string)),*) long(isource)

    ! read depth
    read(1,"(a)") string
    read(string(7:len_trim(string)),*) depth(isource)

    ! read magnitude
    read(1,"(a)") string
    read(string(21:len_trim(string)),*) factor_force_source(isource)

    ! read direction vector's East component
    read(1,"(a)") string
    read(string(29:len_trim(string)),*) comp_dir_vect_source_E(isource)

    ! read direction vector's North component
    read(1,"(a)") string
    read(string(29:len_trim(string)),*) comp_dir_vect_source_N(isource)

    ! read direction vector's vertical component
    read(1,"(a)") string
    read(string(32:len_trim(string)),*) comp_dir_vect_source_Z_UP(isource)

  enddo

  close(1)

  ! Sets tshift_force to zero to initiate the simulation!
  if(NSOURCES == 1)then
      tshift_force = 0.d0
      min_tshift_force_original = t_shift(1)
  else
      tshift_force(1:NSOURCES) = t_shift(1:NSOURCES)-minval(t_shift)
      min_tshift_force_original = minval(t_shift)
  endif

  do isource=1,NSOURCES
     ! checks half-duration
     ! half-duration is the dominant frequency of the source
     ! point forces use a Ricker source time function
     ! null half-duration indicates a very low-frequency source
     ! (see constants.h: TINYVAL = 1.d-9 )
     if( hdur(isource) < TINYVAL ) hdur(isource) = TINYVAL

     ! check (inclined) force source's direction vector
     length = sqrt( comp_dir_vect_source_E(isource)**2 + comp_dir_vect_source_N(isource)**2 + &
          comp_dir_vect_source_Z_UP(isource)**2 )
     if( length < TINYVAL) then
        print *, 'normal length: ', length
        print *, 'isource: ',isource
        stop 'error set force point normal length, make sure all forces have a non null direction vector'
     endif
  enddo

  end subroutine get_force
