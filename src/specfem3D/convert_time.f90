!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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

! open-source subroutines taken from the World Ocean Circulation Experiment (WOCE)
! web site at http://www.coaps.fsu.edu/woce/html/wcdtools.htm

! converted to Fortran90 by Dimitri Komatitsch, University of Pau, France, January 2008.
! Also converted "convtime" from a function to a subroutine.
! Also used a more complete test to detect leap years (the original version was incomplete).

! extended by Dimitri Komatitsch, University of Toulouse, France, April 2011,
! to go beyond the year 2020; I extended that to the year 3000 and thus had to write a loop to fill array "year()".

  subroutine convtime(timestamp,yr,mon,day,hr,minvalue)

! Originally written by Shawn Smith (ssmith AT coaps.fsu.edu)
! Updated Spring 1999 for Y2K compliance by Anthony Arguez (anthony AT coaps.fsu.edu).

! This subroutine will convert a given year, month, day, hour, and
! minutes to a minutes from 01 Jan 1980 00:00 time stamp.

  implicit none

  integer, parameter :: MAX_YEAR = 3000

  integer, intent(out) :: timestamp

  integer, intent(in) :: yr,mon,day,hr,minvalue

  integer :: year(1980:MAX_YEAR),month(12),leap_mon(12)

  integer ::  min_day,min_hr,iyr

! function to determine if year is a leap year
  logical, external :: is_leap_year

  data month /0, 44640, 84960, 129600, 172800, 217440, 260640, &
              305280, 349920, 393120, 437760, 480960/

  data leap_mon /0, 44640, 86400, 131040, 174240, 218880, 262080, &
                 306720, 351360, 394560, 439200, 482400/

  data min_day, min_hr /1440, 60/

! loop added by Dimitri Komatitsch to fill array "year()" automatically
  year(:) = 0
  do iyr = 1981,MAX_YEAR
    if (is_leap_year(iyr-1)) then
      year(iyr) = year(iyr-1) + 366*24*60 ! number of minutes in a year if leap year
    else
      year(iyr) = year(iyr-1) + 365*24*60 ! number of minutes in a year if not leap year
    endif
  enddo

! Test values to see if they fit valid ranges
  if (yr < 1980) stop 'Error in convtime: year out of range, must be >= 1980'

  if (mon < 1 .or. mon > 12) stop 'Error in convtime: month out of range (1-12)'

  if (mon == 2) then
   if (is_leap_year(yr) .and. (day < 1 .or. day > 29)) then
      stop 'Error in convtime: February day out of range (1-29)'
   else if (.not. is_leap_year(yr) .and. (day < 1 .or. day > 28)) then
      stop 'Error in convtime: February day out of range (1-28)'
   endif
  else if (mon == 4 .or. mon == 6 .or. mon == 9 .or. mon == 11) then
   if (day < 1 .or. day > 30) stop 'Error in convtime: day out of range (1-30)'
  else
   if (day < 1 .or. day > 31) stop 'Error in convtime: day out of range (1-31)'
  endif

  if (hr < 0 .or. hr > 23) stop 'Error in convtime: hour out of range (0-23)'

  if (minvalue < 0 .or. minvalue > 60) stop 'Error in convtime: minute out of range (0-60)'

! convert time (test if leap year)
  if (is_leap_year(yr)) then
   timestamp = year(yr)+leap_mon(mon)+((day-1)*min_day)+(hr*min_hr)+minvalue
  else
   timestamp = year(yr)+month(mon)+((day-1)*min_day)+(hr*min_hr)+minvalue
  endif

  end subroutine convtime

!
!----
!

  subroutine invtime(timestamp,yr,mon,day,hr,minvalue)

! This subroutine will convert a minutes timestamp to a year/month
! date. Based on the function convtime by Shawn Smith (COAPS).
!
! Written the spring of 1995, several iterations.
! James N. Stricherz (stricherz AT coaps.fsu.edu)
!
! Updated for Y2K compliance in July 1999.
! Shyam Lakshmin (lakshmin AT coaps.fsu.edu)
!
! This code returns correct results for the range of 01 Jan 1980 00:00
! through 31 Dec 2020 23:59. I know it does, because I tried each minute of that range.

  implicit none

  integer, parameter :: MAX_YEAR = 3000

  integer, intent(in) :: timestamp

  integer, intent(out) :: yr,mon,day,hr,minvalue

  integer :: year(1980:MAX_YEAR),month(13),leap_mon(13)

  integer :: min_day,min_hr,itime,tmon,ttime,thour,iyr,imon,iday,ihour

! function to determine if year is a leap year
  logical, external :: is_leap_year

  data month /0,  44640, 84960, 129600, 172800, 217440, 260640, &
            305280, 349920, 393120, 437760, 480960, 525600/

  data leap_mon /0,  44640,  86400, 131040, 174240, 218880, 262080, &
            306720, 351360, 394560, 439200, 482400, 527040/

  data min_day, min_hr /1440, 60/

! loop added by Dimitri Komatitsch to fill array "year()" automatically
  year(:) = 0
  do iyr = 1981,MAX_YEAR
    if (is_leap_year(iyr-1)) then
      year(iyr) = year(iyr-1) + 366*24*60 ! number of minutes in a year if leap year
    else
      year(iyr) = year(iyr-1) + 365*24*60 ! number of minutes in a year if not leap year
    endif
  enddo

! OK, let us invert the effects of the years: subtract off the
! number of minutes per year until it goes negative
! iyr then gives the year that the time (in minutes) occurs
  if (timestamp >= year(MAX_YEAR)) stop 'year too high in invtime'

  iyr = 1979
  itime=timestamp

 10 iyr=iyr+1
  ttime=itime-year(iyr)
  if (ttime <= 0) then
   if (iyr == 1980) iyr=iyr+1
   iyr=iyr-1
   itime=itime-year(iyr)
  else
   goto 10
  endif

! assign the return variable
  yr=iyr

! OK, the remaining time is less than one full year, so convert
! by the same method as above into months
  imon = 0

! if not leap year
  if (.not. is_leap_year(iyr)) then

! increment the month, and subtract off the minutes from the
! remaining time for a non-leap year
 20 imon=imon+1
   tmon=itime-month(imon)
   if (tmon > 0) then
      goto 20
   else if (tmon < 0) then
      imon=imon-1
      itime=itime-month(imon)
   else
      if (imon > 12) then
         imon=imon-12
         yr=yr+1
      endif
      mon=imon
      day = 1
      hr = 0
      minvalue = 0
      return
   endif

! if leap year
  else

! same thing, same code, but for a leap year
 30 imon=imon+1
   tmon=itime-leap_mon(imon)
   if (tmon > 0) then
      goto 30
   else if (tmon < 0) then
      imon=imon-1
      itime=itime-month(imon)
   else
      if (imon > 12) then
         imon=imon-12
         yr=yr+1
      endif
      mon=imon
      day = 1
      hr = 0
      minvalue = 0
      return
   endif
  endif

! assign the return variable
  mon=imon

! any remaining minutes will belong to day/hour/minutes
! OK, let us get the days
  iday = 0
 40 iday=iday+1
  ttime=itime-min_day
  if (ttime >= 0) then
   itime=ttime
   goto 40
  endif

! assign the return variable
  if (is_leap_year(iyr) .and. mon > 2) then
   day=iday-1
  else
   day=iday
  endif

! pick off the hours of the days...remember, hours can be 0, so we start at -1
  ihour=-1
 50 ihour=ihour+1
  thour=itime-min_hr
  if (thour >= 0) then
   itime=thour
   goto 50
  endif

! assign the return variables
  hr=ihour

! the remainder at this point is the minutes, so return them directly
  minvalue=itime

  end subroutine invtime

