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

  subroutine get_cmt(yr,jda,ho,mi,sec,t_cmt,hdur,lat,long,depth,moment_tensor,DT,NSOURCES)

  implicit none

  include "constants.h"

!--- input or output arguments of the subroutine below

  integer, intent(in) :: NSOURCES
  double precision, intent(in) :: DT

  integer, intent(out) :: yr,jda,ho,mi
  double precision, intent(out) :: sec
  double precision, dimension(NSOURCES), intent(out) :: t_cmt,hdur,lat,long,depth
  double precision, dimension(6,NSOURCES), intent(out) :: moment_tensor

!--- local variables below

  integer mo,da,julian_day,isource
  double precision scaleM
  character(len=5) datasource
  character(len=150) string, CMTSOLUTION

!
!---- read hypocenter info
!
  call get_value_string(CMTSOLUTION, 'solver.CMTSOLUTION', 'DATA/CMTSOLUTION')

  open(unit=1,file=CMTSOLUTION,status='old',action='read')

! read source number isource
  do isource=1,NSOURCES

! read header with event information
  read(1,"(a4,i5,i3,i3,i3,i3,f6.2)") datasource,yr,mo,da,ho,mi,sec
  jda=julian_day(yr,mo,da)

! ignore line with event name
  read(1,"(a)") string

! read time shift
  read(1,"(a)") string
  read(string(12:len_trim(string)),*) t_cmt(isource)

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

! null half-duration indicates a Heaviside
! replace with very short error function
  if(hdur(isource) < 5. * DT) hdur(isource) = 5. * DT

  enddo

  close(1)

!
! scale and non-dimensionalize the moment tensor
! CMTSOLUTION file values are in dyne.cm
! 1 dyne is 1 gram * 1 cm / (1 second)^2
! 1 Newton is 1 kg * 1 m / (1 second)^2
! thus 1 Newton = 100,000 dynes
! therefore 1 dyne.cm = 1e-7 Newton.m
!
  scaleM = 1.d7 * RHOAV * (R_EARTH**5) * PI*GRAV*RHOAV
  moment_tensor(:,:) = moment_tensor(:,:) / scaleM

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

! ------------------------------------------------------------------

! function to determine if year is a leap year
  logical function is_leap_year(yr)

  implicit none

  integer yr

  integer, external :: lpyr

!---- function lpyr above returns 1 if leap year
  if(lpyr(yr) == 1) then
    is_leap_year = .true.
  else
    is_leap_year = .false.
  endif

  end function is_leap_year

