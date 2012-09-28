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

  subroutine get_cmt(yr,jda,ho,mi,sec,tshift_cmt,hdur,lat,long,depth,moment_tensor,&
                    DT,NSOURCES,min_tshift_cmt_original)

  implicit none

  include "constants.h"

!--- input or output arguments of the subroutine below

  integer, intent(in) :: NSOURCES
  double precision, intent(in) :: DT

  integer, intent(out) :: yr,jda,ho,mi
  double precision, intent(out) :: sec,min_tshift_cmt_original
  double precision, dimension(NSOURCES), intent(out) :: tshift_cmt,hdur,lat,long,depth
  double precision, dimension(6,NSOURCES), intent(out) :: moment_tensor

!--- local variables below

  integer mo,da,julian_day,isource
  double precision t_shift(NSOURCES)
  character(len=5) datasource
  character(len=256) string, CMTSOLUTION

  ! initializes
  lat(:) = 0.d0
  long(:) = 0.d0
  depth(:) = 0.d0
  t_shift(:) = 0.d0
  tshift_cmt(:) = 0.d0
  hdur(:) = 0.d0
  moment_tensor(:,:) = 0.d0
  yr = 0
  jda = 0
  ho = 0
  mi = 0
  sec = 0.d0

!
!---- read hypocenter info
!
  call get_value_string(CMTSOLUTION, 'solver.CMTSOLUTION', &
       IN_DATA_FILES_PATH(1:len_trim(IN_DATA_FILES_PATH))//'CMTSOLUTION')

  open(unit=1,file=CMTSOLUTION,status='old',action='read')

! read source number isource
  do isource=1,NSOURCES

    read(1,"(a256)") string
    ! skips empty lines
    do while( len_trim(string) == 0 )
      read(1,"(a256)") string
    enddo

    ! read header with event information
    read(string,"(a4,i5,i3,i3,i3,i3,f6.2)") datasource,yr,mo,da,ho,mi,sec
    jda=julian_day(yr,mo,da)

    ! ignore line with event name
    read(1,"(a)") string

    ! read time shift
    read(1,"(a)") string
    !read(string(12:len_trim(string)),*) tshift_cmt(isource)
    read(string(12:len_trim(string)),*) t_shift(isource)

    ! read half duration
    read(1,"(a)") string
    read(string(15:len_trim(string)),*) hdur(isource)

    ! read latitude
    read(1,"(a)") string
    read(string(10:len_trim(string)),*) lat(isource)

    ! read longitude
    read(1,"(a)") string
    read(string(11:len_trim(string)),*) long(isource)

    ! read depth
    read(1,"(a)") string
    read(string(7:len_trim(string)),*) depth(isource)

    ! read Mrr
    read(1,"(a)") string
    read(string(5:len_trim(string)),*) moment_tensor(1,isource)

    ! read Mtt
    read(1,"(a)") string
    read(string(5:len_trim(string)),*) moment_tensor(2,isource)

    ! read Mpp
    read(1,"(a)") string
    read(string(5:len_trim(string)),*) moment_tensor(3,isource)

    ! read Mrt
    read(1,"(a)") string
    read(string(5:len_trim(string)),*) moment_tensor(4,isource)

    ! read Mrp
    read(1,"(a)") string
    read(string(5:len_trim(string)),*) moment_tensor(5,isource)

    ! read Mtp
    read(1,"(a)") string
    read(string(5:len_trim(string)),*) moment_tensor(6,isource)

    ! checks half-duration
    ! null half-duration indicates a Heaviside
    ! replace with very short error function
    if( hdur(isource) < 5. * DT ) hdur(isource) = 5. * DT

  enddo

  close(1)

  ! Sets tshift_cmt to zero to initiate the simulation!
  if(NSOURCES == 1)then
      tshift_cmt = 0.d0
      min_tshift_cmt_original = t_shift(1)
  else
      tshift_cmt(1:NSOURCES) = t_shift(1:NSOURCES)-minval(t_shift)
      min_tshift_cmt_original = minval(t_shift)
  endif


  !
  ! scale the moment tensor
  ! CMTSOLUTION file values are in dyne.cm
  ! 1 dyne is 1 gram * 1 cm / (1 second)^2
  ! 1 Newton is 1 kg * 1 m / (1 second)^2
  ! thus 1 Newton = 100,000 dynes
  ! therefore 1 dyne.cm = 1e-7 Newton.m
  !
  moment_tensor(:,:) = moment_tensor(:,:) * 1.d-7

  end subroutine get_cmt

! ------------------------------------------------------------------

  integer function julian_day(yr,mo,da)

  implicit none

  integer yr,mo,da

  integer mon(12)
  integer lpyr
  data mon /0,31,59,90,120,151,181,212,243,273,304,334/

  julian_day = da + mon(mo)
  if(mo>2) julian_day = julian_day + lpyr(yr)

  end function julian_day

! ------------------------------------------------------------------

  integer function lpyr(yr)

  implicit none

  integer yr
!
!---- returns 1 if leap year
!
  lpyr=0
  if(mod(yr,400) == 0) then
    lpyr=1
  else if(mod(yr,4) == 0) then
    lpyr=1
    if(mod(yr,100) == 0) lpyr=0
  endif

  end function lpyr

