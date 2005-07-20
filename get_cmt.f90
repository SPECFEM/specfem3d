!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 3
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology July 2005
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine get_cmt(yr,jda,ho,mi,sec,t_cmt,hdur,lat,long,depth,moment_tensor,DT,NSOURCES)

  implicit none

  include "constants.h"

  integer yr,jda,ho,mi,NSOURCES
  double precision sec
  double precision, dimension(NSOURCES) :: t_cmt,hdur,lat,long,depth
  double precision moment_tensor(6,NSOURCES)
  double precision DT

  integer mo,da,julian_day,isource
  character(len=5) datasource
  character(len=150) string

!
!---- read hypocenter info
!
  open(unit=1,file='DATA/CMTSOLUTION',status='old')

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

