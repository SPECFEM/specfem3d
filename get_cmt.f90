!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 1
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology October 2002
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine get_cmt(cmt_file,yr,jda,ho,mi,sec, &
                      t_cmt,hdur,elat,elon,depth,moment_tensor,DT)

  implicit none

  include "constants.h"

  integer yr,jda,ho,mi
  double precision sec,t_cmt,hdur,elat,elon,depth
  double precision moment_tensor(6)
  double precision DT
  character(len=150) cmt_file

  integer ios,lstr,mo,da,julian_day
  double precision scaleM
  character(len=5) datasource
  character(len=150) string

!
!---- read hypocenter info
!
  open(unit=1,file=cmt_file,iostat=ios,status='old')
  if(ios /= 0) stop 'error opening CMT file '

  read(1,"(a4,i5,i3,i3,i3,i3,f6.2)") datasource,yr,mo,da,ho,mi,sec

  jda=julian_day(yr,mo,da)
  t_cmt = 0.

  ios=0
  do while(ios == 0)
    read(1,"(a)",iostat=ios) string

    if(ios == 0) then
      lstr=len_trim(string)
      if(string(1:10) == 'time shift') then
        read(string(12:lstr),*) t_cmt
      else if(string(1:13) == 'half duration') then
        read(string(15:lstr),*) hdur
      else if(string(1:8) == 'latitude') then
        read(string(10:lstr),*) elat
      else if(string(1:9) == 'longitude') then
        read(string(11:lstr),*) elon
      else if(string(1:5) == 'depth') then
        read(string(7:lstr),*) depth
      else if(string(1:3) == 'Mrr') then
        read(string(5:lstr),*) moment_tensor(1)
      else if(string(1:3) == 'Mtt') then
        read(string(5:lstr),*) moment_tensor(2)
      else if(string(1:3) == 'Mpp') then
        read(string(5:lstr),*) moment_tensor(3)
      else if(string(1:3) == 'Mrt') then
        read(string(5:lstr),*) moment_tensor(4)
      else if(string(1:3) == 'Mrp') then
        read(string(5:lstr),*) moment_tensor(5)
      else if(string(1:3) == 'Mtp') then
        read(string(5:lstr),*) moment_tensor(6)
      endif
    endif
  enddo

  close(1)

! time shift not used in our code, but kept for compatibility of file format
  if(dabs(t_cmt) > DT) stop 't_cmt not implemented in current code'
  t_cmt = 0.

! null half-duration indicates a Heaviside
! replace with very short error function
  if(hdur < DT) hdur = 5. * DT

!
! scale the moment-tensor (dimensions dyn-cm)
!
  scaleM = 1.d7
  moment_tensor(:) = moment_tensor(:) / scaleM

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

