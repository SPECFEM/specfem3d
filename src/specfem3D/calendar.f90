!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
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

  integer function julian_day(yr,mo,da)

  implicit none

  integer, intent(in) :: yr,mo,da

  integer :: mon(12)
  integer, external :: lpyr
  data mon /0,31,59,90,120,151,181,212,243,273,304,334/

  julian_day = da + mon(mo)
  if (mo > 2) julian_day = julian_day + lpyr(yr)

  end function julian_day

! ------------------------------------------------------------------

  integer function lpyr(yr)

  implicit none

  integer, intent(in) :: yr
!
!---- returns 1 if leap year
!
  lpyr=0
  if (mod(yr,400) == 0) then
    lpyr=1
  else if (mod(yr,4) == 0) then
    lpyr=1
    if (mod(yr,100) == 0) lpyr=0
  endif

  end function lpyr

! ------------------------------------------------------------------

! function to determine if year is a leap year
  logical function is_leap_year(yr)

  implicit none

  integer yr

  integer, external :: lpyr

!---- function lpyr above returns 1 if leap year
  if (lpyr(yr) == 1) then
    is_leap_year = .true.
  else
    is_leap_year = .false.
  endif

  end function is_leap_year


!----------------------------------------------------------------------------------------------
! open-source subroutines below taken from ftp://ftp.met.fsu.edu/pub/ahlquist/calendar_software
!----------------------------------------------------------------------------------------------

  integer function idaywk(jdayno)

! IDAYWK = compute the DAY of the WeeK given the Julian Day number,
!          version 1.0.

  implicit none

! Input variable
  integer, intent(in) :: jdayno
! jdayno = Julian Day number starting at noon of the day in question.

! Output of the function:
! idaywk = day of the week, where 0=Sunday, 1=Monday, ..., 6=Saturday.

!----------
! Compute the day of the week given the Julian Day number.
! You can find the Julian Day number given (day,month,year)
! using subroutine calndr below.
! Example: For the first day of the Gregorian calendar,
! Friday 15 October 1582, compute the Julian day number (option 3 of
! subroutine calndr) and compute the day of the week.
!     call calndr (3, 15, 10, 1582, jdayno)
!     write(*,*) jdayno, idaywk(jdayno)
! The numbers printed should be 2299161 and 5, where 5 refers to Friday.
!
! Copyright (C) 1999 Jon Ahlquist.
! Issued under the second GNU General Public License.
! See www.gnu.org for details.
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! If you find any errors, please notify:
! Jon Ahlquist
! Dept of Meteorology
! Florida State University
! Tallahassee, FL 32306-4520
! 15 March 1999.
!
!-----

! converted to Fortran90 by Dimitri Komatitsch,
! University of Pau, France, January 2008.

! jdSun is the Julian Day number starting at noon on any Sunday.
! I arbitrarily chose the first Sunday after Julian Day 1,
! which is Julian Day 6.
  integer, parameter :: jdSun = 6

  idaywk = mod(jdayno-jdSun,7)

! If jdayno-jdSun < 0, then we are taking the modulus of a negative
! number. Fortran's built-in mod function returns a negative value
! when the argument is negative.  In that case, we adjust the result
! to a positive value.
  if (idaywk < 0) idaywk = idaywk + 7

  end function idaywk

!
!----
!

  subroutine calndr(iday,month,iyear,idayct)

! CALNDR = CALeNDaR conversions, version 1.0

  implicit none

! specify the desired calendar conversion option.
! in order to return the julian day number, compatible with function idaywk from above,
! we choose option 3
! (tested with dates: Feb, 23 2010 -> idaywk = Tue
!                               Dec, 24 2009 -> idaywk = Thu
!                               Oct, 15 1582  -> idaywk = Fri ...which all look o.k. )
  integer, parameter :: ioptn = 3

! Input/Output variables
  integer, intent(inout) :: iday,month,iyear,idayct

!----------
!
! The subroutine calndr() performs calendar calculations using either
! the standard Gregorian calendar or the old Julian calendar.
! This subroutine extends the definitions of these calendar systems
! to any arbitrary year.  The algorithms in this subroutine
! will work with any date in the past or future,
! but overflows will occur if the numbers are sufficiently large.
! For a computer using a 32-bit integer, this routine can handle
! any date between roughly 5.8 million BC and 5.8 million AD
! without experiencing overflow during calculations.
!
! No external functions or subroutines are called.
!
!----------
!
! Input/output arguments for subroutine CALNDR()
!
! "ioptn" is the desired calendar conversion option explained below.
! Positive option values use the standard modern Gregorian calendar.
! Negative option values use the old Julian calendar which was the
! standard in Europe from its institution by Julius Caesar in 45 BC
! until at least 4 October 1582.  The Gregorian and Julian calendars
! are explained further below.
!
! (iday,month,iyear) is a calendar date where "iday" is the day of
! the month, "month" is 1 for January, 2 for February, etc.,
! and "iyear" is the year.  If the year is 1968 AD, enter iyear=1968,
! since iyear=68 would refer to 68 AD.
! For BC years, iyear should be negative, so 45 BC would be iyear=-45.
! By convention, there is no year 0 under the BC/AD year numbering
! scheme.  That is, years proceed as 2 BC, 1 BC, 1 AD, 2 AD, etc.,
! without including 0. The subroutine calndr() will print an error message
! and stop if you specify iyear = 0.
!
! "idayct" is a day count.  It is either the day number during the
! specified year or the Julian Day number, depending on the value
! of ioptn.  By day number during the specified year, we mean
! idayct=1 on 1 January, idayct=32 on 1 February, etc., to idayct=365
! or 366 on 31 December, depending on whether the specified year
! is a leap year.
!
! The values of input variables are not changed by this subroutine.
!
!
! ALLOWABLE VALUES FOR "IOPTN" and the conversions they invoke.
! Positive option values ( 1 to  5) use the standard Gregorian calendar.
! Negative option values (-1 to -5) use the old      Julian    calendar.
!
! Absolute
!  value
! of ioptn   Input variable(s)     Output variable(s)
!
!    1       iday,month,iyear      idayct
! Given a calendar date (iday,month,iyear), compute the day number
! (idayct) during the year, where 1 January is day number 1 and
! 31 December is day number 365 or 366, depending on whether it is
! a leap year.
!
!    2       idayct,iyear          iday,month
! Given the day number of the year (idayct) and the year (iyear),
! compute the day of the month (iday) and the month (month).
!
!    3       iday,month,iyear      idayct
! Given a calendar date (iday,month,iyear), compute the Julian Day
! number (idayct) that starts at noon of the calendar date specified.
!
!    4       idayct                iday,month,iyear
! Given the Julian Day number (idayct) that starts at noon,
! compute the corresponding calendar date (iday,month,iyear).
!
!    5       idayct                iday,month,iyear
! Given the Julian Day number (idayct) that starts at noon,
! compute the corresponding day number for the year (iday)
! and year (iyear).  On return from calndr(), "month" will always
! be set equal to 1 when ioptn=5.
!
! No inverse function is needed for ioptn=5 because it is
! available through option 3.  One simply calls calndr() with:
! ioptn = 3,
! iday  = day number of the year instead of day of the month,
! month = 1, and
! iyear = whatever the desired year is.
!
!----------
!
! EXAMPLES
! The first 6 examples are for the standard Gregorian calendar.
! All the examples deal with 15 October 1582, which was the first day
! of the Gregorian calendar.  15 October is the 288-th day of the year.
! Julian Day number 2299161 began at noon on 15 October 1582.
!
! Find the day number during the year on 15 October 1582
!     ioptn = 1
!     call calndr (ioptn, 15, 10, 1582,  idayct)
! calndr() should return idayct=288
!
! Find the day of the month and month for day 288 in year 1582.
!     ioptn = 2
!     call calndr (ioptn, iday, month, 1582, 288)
! calndr() should return iday=15 and month=10.
!
! Find the Julian Day number for 15 October 1582.
!     ioptn = 3
!     call calndr (ioptn, 15, 10, 1582, julian)
! calndr() should return julian=2299161
!
! Find the Julian Day number for day 288 during 1582 AD.
! When the input is day number of the year, one should specify month=1
!     ioptn = 3
!     call calndr (ioptn, 288, 1, 1582, julian)
! calndr() should return dayct=2299161
!
! Find the date for Julian Day number 2299161.
!     ioptn = 4
!     call calndr (ioptn, iday, month, iyear, 2299161)
! calndr() should return iday=15, month=10, and iyear=1582
!
! Find the day number during the year (iday) and year
! for Julian Day number 2299161.
!     ioptn = 5
!     call calndr (ioptn, iday, month, iyear, 2299161)
! calndr() should return iday=288, month = 1, iyear=1582
!
! Given 15 October 1582 under the Gregorian calendar,
! find the date (idayJ,imonthJ,iyearJ) under the Julian calendar.
! To do this, we call calndr() twice, using the Julian Day number
! as the intermediate value.
!     call calndr ( 3, 15,        10, 1582,    julian)
!     call calndr (-4, idayJ, monthJ, iyearJ,  julian)
! The first call to calndr() should return julian=2299161, and
! the second should return idayJ=5, monthJ=10, iyearJ=1582
!
!----------
!
! BASIC CALENDAR INFORMATION
!
! The Julian calendar was instituted by Julius Caesar in 45 BC.
! Every fourth year is a leap year in which February has 29 days.
! That is, the Julian calendar assumes that the year is exactly
! 365.25 days long.  Actually, the year is not quite this long.
! The modern Gregorian calendar remedies this by omitting leap years
! in years divisible by 100 except when the year is divisible by 400.
! Thus, 1700, 1800, and 1900 are leap years under the Julian calendar
! but not under the Gregorian calendar.  The years 1600 and 2000 are
! leap years under both the Julian and the Gregorian calendars.
! Other years divisible by 4 are leap years under both calendars,
! such as 1992, 1996, 2004, 2008, 2012, etc.  For BC years, we recall
! that year 0 was omitted, so 1 BC, 5 BC, 9 BC, 13 BC, etc., and 401 BC,
! 801 BC, 1201 BC, etc., are leap years under both calendars, while
! 101 BC, 201 BC, 301 BC, 501 BC, 601 BC, 701 BC, 901 BC, 1001 BC,
! 1101 BC, etc., are leap years under the Julian calendar but not
! the Gregorian calendar.
!
! The Gregorian calendar is named after Pope Gregory XIII.  He declared
! that the last day of the old Julian calendar would be Thursday,
! 4 October 1582 and that the following day, Friday, would be reckoned
! under the new calendar as 15 October 1582.  The jump of 10 days was
! included to make 21 March closer to the spring equinox.
!
! Only a few Catholic countries (Italy, Poland, Portugal, and Spain)
! switched to the Gregorian calendar on the day after 4 October 1582.
! It took other countries months to centuries to change to the
! Gregorian calendar.  For example, England's first day under the
! Gregorian calendar was 14 September 1752.  The same date applied to
! the entire British empire, including America.  Japan, Russia, and many
! eastern European countries did not change to the Gregorian calendar
! until the 20th century.  The last country to change was Turkey,
! which began using the Gregorian calendar on 1 January 1927.
!
! Therefore, between the years 1582 and 1926 AD, you must know
! the country in which an event was dated to interpret the date
! correctly.  In Sweden, there was even a year (1712) when February
! had 30 days.  Consult a book on calendars for more details
! about when various countries changed their calendars.
!
! DAY NUMBER DURING THE YEAR
! The day number during the year is simply a counter equal to 1 on
! 1 January, 32 on 1 February, etc., through 365 or 366 on 31 December,
! depending on whether the year is a leap year.  Sometimes this is
! called the Julian Day, but that term is better reserved for the
! day counter explained below.
!
! JULIAN DAY NUMBER
! The Julian Day numbering system was designed by Joseph Scaliger
! in 1582 to remove ambiguity caused by varying calendar systems.
! The name "Julian Day" was chosen to honor Scaliger's father,
! Julius Caesar Scaliger (1484-1558), an Italian scholar and physician
! who lived in France.  Because Julian Day numbering was especially
! designed for astronomers, Julian Days begin at noon so that the day
! counter does not change in the middle of an astronomer's observing
! period.  Julian Day 0 began at noon on 1 January 4713 BC under the
! Julian calendar.  A modern reference point is that 23 May 1968
! (Gregorian calendar) was Julian Day 2,440,000.
!
! JULIAN DAY NUMBER EXAMPLES
!
! The table below shows a few Julian Day numbers and their corresponding
! dates, depending on which calendar is used.  A negative 'iyear' refers
! to BC (Before Christ).
!
!                     Julian Day under calendar:
! iday  month   iyear     Gregorian   Julian
!  24     11   -4714            0        -38
!   1      1   -4713           38          0
!   1      1       1      1721426    1721424
!   4     10    1582      2299150    2299160
!  15     10    1582      2299161    2299171
!   1      3    1600      2305508    2305518
!  23      5    1968      2440000    2440013
!   5      7    1998      2451000    2451013
!   1      3    2000      2451605    2451618
!   1      1    2001      2451911    2451924
!
! From this table, we can see that the 10 day difference between the
! two calendars in 1582 grew to 13 days by 1 March 1900, since 1900 was
! a leap year under the Julian calendar but not under the Gregorian
! calendar.  The gap will widen to 14 days after 1 March 2100 for the
! same reason.
!
!----------
!
! PORTABILITY
!
! This subroutine is written in standard Fortran 90.
! It calls no external functions or subroutines and should run
! without problem on any computer having a 32-bit word or longer.
!
!----------
!
! ALGORITHM
!
! The goal in coding calndr() was clear, clean code, not efficiency.
! Calendar calculations usually take a trivial fraction of the time
! in any program in which dates conversions are involved.
! Data analysis usually takes the most time.
!
! Standard algorithms are followed in this subroutine.  Internal to
! this subroutine, we use a year counter "jyear" such that
!  jyear=iyear   when iyear is positive
!       =iyear+1 when iyear is negative.
! Thus, jyear does not experience a 1 year jump like iyear does
! when going from BC to AD.  Specifically, jyear = 0 when iyear=-1,
! i.e., when the year is 1 BC.
!
! For simplicity in dealing with February, inside this subroutine,
! we let the year begin on 1 March so that the adjustable month,
! February is the last month of the year.
! It is clear that the calendar used to work this way because the
! months September, October, November, and December refer to
! 7, 8, 9, and 10.  For consistency, jyear is incremented on 1 March
! rather than on 1 January.  Of course, everything is adjusted back to
! standard practice of years beginning on 1 January before answers
! are returned to the routine that calls calndr().
!
! Lastly, we use a trick to calculate the number of days from 1 March
! until the end of the month that precedes the specified month.
! That number of days is int(30.6001*(month+1))-122,
! where 30.6001 is used to avoid the possibility of round-off and
! truncation error.  For example, if 30.6 were used instead,
! 30.6*5 should be 153, but round-off error could make it 152.99999,
! which would then truncated to 152, causing an error of 1 day.
!
! Algorithm reference:
! Dershowitz, Nachum and Edward M. Reingold, 1990: Calendrical
! Calculations.  Software-Practice and Experience, vol. 20, number 9
! (September 1990), pp. 899-928.
!
! Copyright (C) 1999 Jon Ahlquist.
! Issued under the second GNU General Public License.
! See www.gnu.org for details.
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! If you find any errors, please notify:
! Jon Ahlquist
! Dept of Meteorology
! Florida State University
! Tallahassee, FL 32306-4520
! 15 March 1999.
!
!-----

! converted to Fortran90 by Dimitri Komatitsch,
! University of Pau, France, January 2008.

! Declare internal variables.
  integer jdref, jmonth, jyear, leap, n1yr, n4yr, n100yr, n400yr, ndays, ndy400, ndy100, nyrs, yr400, yrref
!
! Explanation of all internal variables.
! jdref   Julian Day on which 1 March begins in the reference year.
! jmonth  Month counter which equals month+1 if month > 2
!          or month+13 if month <= 2.
! jyear   Year index,  jyear=iyear if iyear > 0, jyear=iyear+1
!            if iyear < 0.  Thus, jyear does not skip year 0
!            like iyear does between BC and AD years.
! leap    =1 if the year is a leap year,  = 0 if not.
! n1yr    Number of complete individual years between iyear and
!            the reference year after all 4, 100,
!            and 400 year periods have been removed.
! n4yr    Number of complete 4 year cycles between iyear and
!            the reference year after all 100 and 400 year periods
!            have been removed.
! n100yr  Number of complete 100 year periods between iyear and
!            the reference year after all 400 year periods
!            have been removed.
! n400yr  Number of complete 400 year periods between iyear and
!            the reference year.
! ndays   Number of days since 1 March during iyear.  (In intermediate
!            steps, it holds other day counts as well.)
! ndy400  Number of days in 400 years.  Under the Gregorian calendar,
!            this is 400*365 + 100 - 3 = 146097.  Under the Julian
!            calendar, this is 400*365 + 100 = 146100.
! ndy100  Number of days in 100 years,  Under the Gregorian calendar,
!            this is 100*365 + 24 = 36524.   Under the Julian calendar,
!            this is 100*365 + 25 = 36525.
! nyrs    Number of years from the beginning of yr400
!              to the beginning of jyear.  (Used for option +/-3).
! yr400   The largest multiple of 400 years that is <= jyear.
!
!
!----------------------------------------------------------------
! Do preparation work.
!
! Look for out-of-range option values.
  if ((ioptn == 0) .or. (abs(ioptn) >= 6)) then
   write(*,*)'For calndr(), you specified ioptn = ', ioptn
   write(*,*) 'Allowable values are 1 to 5 for the Gregorian calendar'
   write(*,*) 'and -1 to -5 for the Julian calendar.'
   stop
  endif
!
! Options 1-3 have "iyear" as an input value.
! Internally, we use variable "jyear" that does not have a jump
! from -1 (for 1 BC) to +1 (for 1 AD).
  if (abs(ioptn) <= 3) then
   if (iyear > 0) then
      jyear = iyear
   else if (iyear == 0) then
      write(*,*) 'For calndr(), you specified the nonexistent year 0'
      stop
   else
      jyear = iyear + 1
   endif
!
!        Set "leap" equal to 0 if "jyear" is not a leap year
!        and equal to 1 if it is a leap year.
   leap = 0
   if ((jyear/4)*4 == jyear) then
      leap = 1
   endif
   if ((ioptn > 0) .and. ((jyear/100)*100 == jyear) .and. ((jyear/400)*400 /= jyear)) then
      leap = 0
   endif
  endif
!
! Options 3-5 involve Julian Day numbers, which need a reference year
! and the Julian Days that began at noon on 1 March of the reference
! year under the Gregorian and Julian calendars.  Any year for which
! "jyear" is divisible by 400 can be used as a reference year.
! We chose 1600 AD as the reference year because it is the closest
! multiple of 400 to the institution of the Gregorian calendar, making
! it relatively easy to compute the Julian Day for 1 March 1600
! given that, on 15 October 1582 under the Gregorian calendar,
! the Julian Day was 2299161.  Similarly, we need to do the same
! calculation for the Julian calendar.  We can compute this Julian
! Day knowing that on 4 October 1582 under the Julian calendar,
! the Julian Day number was 2299160.  The details of these calculations
! is next.
!    From 15 October until 1 March, the number of days is the remainder
! of October plus the days in November, December, January, and February:
! 17+30+31+31+28 = 137, so 1 March 1583 under the Gregorian calendar
! was Julian Day 2,299,298.  Because of the 10 day jump ahead at the
! switch from the Julian calendar to the Gregorian calendar, 1 March
! 1583 under the Julian calendar was Julian Day 2,299,308.  Making use
! of the rules for the two calendar systems, 1 March 1600 was Julian
! Day 2,299,298 + (1600-1583)*365 + 5 (due to leap years) =
! 2,305,508 under the Gregorian calendar and day 2,305,518 under the
! Julian calendar.
!    We also set the number of days in 400 years and 100 years.
! For reference, 400 years is 146097 days under the Gregorian calendar
! and 146100 days under the Julian calendar.  100 years is 36524 days
! under the Gregorian calendar and 36525 days under the Julian calendar.
  if (abs(ioptn) >= 3) then
!
!        Julian calendar values.
   yrref  =    1600
   jdref  = 2305518
!               = Julian Day reference value for the day that begins
!                 at noon on 1 March of the reference year "yrref".
   ndy400 = 400*365 + 100
   ndy100 = 100*365 +  25
!
!        Adjust for Gregorian calendar values.
   if (ioptn > 0) then
      jdref  = jdref  - 10
      ndy400 = ndy400 -  3
      ndy100 = ndy100 -  1
   endif
  endif
!
!----------------------------------------------------------------
! OPTIONS -1 and +1:
! Given a calendar date (iday,month,iyear), compute the day number
! of the year (idayct), where 1 January is day number 1 and 31 December
! is day number 365 or 366, depending on whether it is a leap year.
  if (abs(ioptn) == 1) then
!
!     Compute the day number during the year.
  if (month <= 2) then
   idayct = iday + (month-1)*31
  else
   idayct = iday + int(30.6001 * (month+1)) - 63 + leap
  endif
!
!----------------------------------------------------------------
! OPTIONS -2 and +2:
! Given the day number of the year (idayct) and the year (iyear),
! compute the day of the month (iday) and the month (month).
  else if (abs(ioptn) == 2) then
!
  if (idayct < 60+leap) then
   month  = (idayct-1)/31
   iday   = idayct - month*31
   month  = month + 1
  else
   ndays  = idayct - (60+leap)
!               = number of days past 1 March of the current year.
   jmonth = (10*(ndays+31))/306 + 3
!               = month counter, =4 for March, =5 for April, etc.
   iday   = (ndays+123) - int(30.6001*jmonth)
   month  = jmonth - 1
  endif
!
!----------------------------------------------------------------
! OPTIONS -3 and +3:
! Given a calendar date (iday,month,iyear), compute the Julian Day
! number (idayct) that starts at noon.
  else if (abs(ioptn) == 3) then
!
!     Shift to a system where the year starts on 1 March, so January
!     and February belong to the preceding year.
!     Define jmonth=4 for March, =5 for April, ..., =15 for February.
  if (month <= 2) then
    jyear  = jyear -  1
    jmonth = month + 13
  else
    jmonth = month +  1
  endif
!
!     Find the closest multiple of 400 years that is <= jyear.
  yr400 = (jyear/400)*400
!           = multiple of 400 years at or less than jyear.
  if (jyear < yr400) then
   yr400 = yr400 - 400
  endif
!
  n400yr = (yr400 - yrref)/400
!            = number of 400-year periods from yrref to yr400.
  nyrs   = jyear - yr400
!            = number of years from the beginning of yr400
!              to the beginning of jyear.
!
!     Compute the Julian Day number.
  idayct = iday + int(30.6001*jmonth) - 123 + 365*nyrs + nyrs/4 &
         + jdref + n400yr*ndy400
!
!     If we are using the Gregorian calendar, we must not count
!     every 100-th year as a leap year.  nyrs is less than 400 years,
!     so we do not need to consider the leap year that would occur if
!     nyrs were divisible by 400, i.e., we do not add nyrs/400.
  if (ioptn > 0) then
   idayct = idayct - nyrs/100
  endif
!
!----------------------------------------------------------------
! OPTIONS -5, -4, +4, and +5:
! Given the Julian Day number (idayct) that starts at noon,
! compute the corresponding calendar date (iday,month,iyear)
! (abs(ioptn)=4) or day number during the year (abs(ioptn)=5).
  else
!
!     Create a new reference date which begins on the nearest
!     400-year cycle less than or equal to the Julian Day for 1 March
!     in the year in which the given Julian Day number (idayct) occurs.
  ndays  = idayct - jdref
  n400yr = ndays / ndy400
!            = integral number of 400-year periods separating
!              idayct and the reference date, jdref.
  jdref  = jdref + n400yr*ndy400
  if (jdref > idayct) then
   n400yr = n400yr - 1
   jdref  = jdref  - ndy400
  endif
!
  ndays  = idayct - jdref
!            = number from the reference date to idayct.
!
  n100yr = min(ndays/ndy100, 3)
!            = number of complete 100-year periods
!              from the reference year to the current year.
!              The min() function is necessary to avoid n100yr=4
!              on 29 February of the last year in the 400-year cycle.
!
  ndays  = ndays - n100yr*ndy100
!            = remainder after removing an integral number of
!              100-year periods.
!
  n4yr   = ndays / 1461
!            = number of complete 4-year periods in the current century.
!              4 years consists of 4*365 + 1 = 1461 days.
!
  ndays  = ndays - n4yr*1461
!            = remainder after removing an integral number
!              of 4-year periods.
!
  n1yr   = min(ndays/365, 3)
!            = number of complete years since the last leap year.
!              The min() function is necessary to avoid n1yr=4
!              when the date is 29 February on a leap year,
!              in which case ndays=1460, and 1460/365 = 4.
!
  ndays  = ndays - 365*n1yr
!            = number of days so far in the current year,
!              where ndays = 0 on 1 March.
!
  iyear  = n1yr + 4*n4yr + 100*n100yr + 400*n400yr + yrref
!            = year, as counted in the standard way,
!              but relative to 1 March.
!
! At this point, we need to separate ioptn=abs(4), which seeks a
! calendar date, and ioptn=abs(5), which seeks the day number during
! the year.  First compute the calendar date if desired (abs(ioptn)=4).
  if (abs(ioptn) == 4) then
   jmonth = (10*(ndays+31))/306 + 3
!               = offset month counter.  jmonth=4 for March, =13 for
!                 December, =14 for January, =15 for February.
   iday   = (ndays+123) - int(30.6001*jmonth)
!               = day of the month, starting with 1 on the first day
!                 of the month.
!
!        Now adjust for the fact that the year actually begins
!        on 1 January.
   if (jmonth <= 13) then
      month = jmonth - 1
   else
      month = jmonth - 13
      iyear = iyear + 1
   endif
!
! This code handles abs(ioptn)=5, finding the day number during the year.
  else
!        ioptn=5 always returns month = 1, which we set now.
   month = 1
!
!        We need to determine whether this is a leap year.
   leap = 0
   if ((jyear/4)*4 == jyear) then
      leap = 1
   endif
   if ((ioptn > 0) .and. ((jyear/100)*100 == jyear) .and. ((jyear/400)*400 /= jyear)) then
      leap = 0
   endif
!
!        Now find the day number "iday".
!        ndays is the number of days since the most recent 1 March,
!        so ndays = 0 on 1 March.
   if (ndays <= 305) then
      iday  = ndays + 60 + leap
   else
      iday  = ndays - 305
      iyear = iyear + 1
   endif
  endif
!
!     Adjust the year if it is <= 0, and hence BC (Before Christ).
  if (iyear <= 0) then
   iyear = iyear - 1
  endif
!
! End the code for the last option, ioptn.
  endif

  end subroutine calndr

