!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 2
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology July 2004
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  program make_cmts

! Jeroen Tromp, July 2001

  implicit none

  include "cmt.h"

  integer yr,jda,ho,mi
  double precision sec,t_cmt,hdur,elat,elon,depth
  double precision moment_tensor(6)
  character(len=150) cmt_file

  integer iu,i,ios,lstr,mo,da,julian_day
  double precision mb,ms
  double precision elatp,elonp
  character(len=24) reg
  character(len=5) datasource
  character(len=150) string

  open(unit=1,file='CMTSOLUTION',iostat=ios,status='old')
  if(ios /= 0) stop 'error opening CMT file '

  open(unit=2,file='CMTSOLUTION_latitude',iostat=ios,status='unknown')
  open(unit=3,file='CMTSOLUTION_longitude',iostat=ios,status='unknown')
  open(unit=4,file='CMTSOLUTION_depth',iostat=ios,status='unknown')
  open(unit=5,file='CMTSOLUTION_Mrr',iostat=ios,status='unknown')
  open(unit=6,file='CMTSOLUTION_Mtt',iostat=ios,status='unknown')
  open(unit=7,file='CMTSOLUTION_Mpp',iostat=ios,status='unknown')
  open(unit=8,file='CMTSOLUTION_Mrt',iostat=ios,status='unknown')
  open(unit=9,file='CMTSOLUTION_Mrp',iostat=ios,status='unknown')
  open(unit=10,file='CMTSOLUTION_Mtp',iostat=ios,status='unknown')

  read(1,"(a4,i5,i3,i3,i3,i3,f6.2,f9.4,f10.4,f6.1,f4.1,f4.1,1x,a)") &
           datasource,yr,mo,da,ho,mi,sec,elat,elon,depth,mb,ms,reg

  do iu=2,10
    write(iu,"(a3,i5,i3,i3,i3,i3,f6.2,f9.4,f10.4,f6.1,f4.1,f4.1,1x,a)") &
             datasource,yr,mo,da,ho,mi,sec,elat,elon,depth,mb,ms,reg
  enddo

  jda=julian_day(yr,mo,da)

  ios=0
  do while(ios == 0)

    read(1,"(a)",iostat=ios) string

    if(ios == 0) then

      lstr=len_trim(string)

      if(string(1:10) == 'event name') then
        do iu=2,10
          write(iu,"(a)") string(1:lstr)
        enddo
      else if(string(1:10) == 'time shift') then
        read(string(12:lstr),*) t_cmt
        do iu=2,10
          write(iu,"(a)") string(1:lstr)
        enddo
      else if(string(1:13) == 'half duration') then
        read(string(15:lstr),*) hdur
        do iu=2,10
          write(iu,"(a)") string(1:lstr)
        enddo
      else if(string(1:8) == 'latitude') then
        read(string(10:lstr),*) elat
        elatp = elat + DDELTA
        if(elatp > 90.0) elatp = 180.0 - elatp
        write(2,"(a9,5x,f9.4)") string(1:9),elatp
        do iu=3,10
          write(iu,"(a)") string(1:lstr)
        enddo
      else if(string(1:9) == 'longitude') then
        read(string(11:lstr),*) elon
        write(2,"(a)") string(1:lstr)
        elonp = elon + DDELTA
        if(elonp > 180.0) elonp = elonp - 360.0
        write(3,"(a10,4x,f9.4)") string(1:10),elonp
        do iu=4,10
          write(iu,"(a)") string(1:lstr)
        enddo
      else if(string(1:5) == 'depth') then
        read(string(7:lstr),*) depth
        write(2,"(a)") string(1:lstr)
        write(3,"(a)") string(1:lstr)
        write(4,"(a6,8x,f9.4)") string(1:6),depth+DDEPTH
        do iu=5,10
          write(iu,"(a)") string(1:lstr)
        enddo
      else if(string(1:3) == 'Mrr') then
        read(string(5:lstr),*) moment_tensor(1)
        do iu=2,4
          write(iu,"(a)") string(1:lstr)
        enddo
        write(5,"(a4,4x,e15.6)") string(1:4),MOMENT
        write(6,"(a4,4x,e15.6)") string(1:4),0.0
        write(7,"(a4,4x,e15.6)") string(1:4),0.0
        write(8,"(a4,4x,e15.6)") string(1:4),0.0
        write(9,"(a4,4x,e15.6)") string(1:4),0.0
        write(10,"(a4,4x,e15.6)") string(1:4),0.0
      else if(string(1:3) == 'Mtt') then
        read(string(5:lstr),*) moment_tensor(2)
        do iu=2,4
          write(iu,"(a)") string(1:lstr)
        enddo
        write(5,"(a4,4x,e15.6)") string(1:4),0.0
        write(6,"(a4,4x,e15.6)") string(1:4),MOMENT
        write(7,"(a4,4x,e15.6)") string(1:4),0.0
        write(8,"(a4,4x,e15.6)") string(1:4),0.0
        write(9,"(a4,4x,e15.6)") string(1:4),0.0
        write(10,"(a4,4x,e15.6)") string(1:4),0.0
      else if(string(1:3) == 'Mpp') then
        read(string(5:lstr),*) moment_tensor(3)
        do iu=2,4
          write(iu,"(a)") string(1:lstr)
        enddo
        write(5,"(a4,4x,e15.6)") string(1:4),0.0
        write(6,"(a4,4x,e15.6)") string(1:4),0.0
        write(7,"(a4,4x,e15.6)") string(1:4),MOMENT
        write(8,"(a4,4x,e15.6)") string(1:4),0.0
        write(9,"(a4,4x,e15.6)") string(1:4),0.0
        write(10,"(a4,4x,e15.6)") string(1:4),0.0
      else if(string(1:3) == 'Mrt') then
        read(string(5:lstr),*) moment_tensor(4)
        do iu=2,4
          write(iu,"(a)") string(1:lstr)
        enddo
        write(5,"(a4,4x,e15.6)") string(1:4),0.0
        write(6,"(a4,4x,e15.6)") string(1:4),0.0
        write(7,"(a4,4x,e15.6)") string(1:4),0.0
        write(8,"(a4,4x,e15.6)") string(1:4),MOMENT
        write(9,"(a4,4x,e15.6)") string(1:4),0.0
        write(10,"(a4,4x,e15.6)") string(1:4),0.0
      else if(string(1:3) == 'Mrp') then
        read(string(5:lstr),*) moment_tensor(5)
        do iu=2,4
          write(iu,"(a)") string(1:lstr)
        enddo
        write(5,"(a4,4x,e15.6)") string(1:4),0.0
        write(6,"(a4,4x,e15.6)") string(1:4),0.0
        write(7,"(a4,4x,e15.6)") string(1:4),0.0
        write(8,"(a4,4x,e15.6)") string(1:4),0.0
        write(9,"(a4,4x,e15.6)") string(1:4),MOMENT
        write(10,"(a4,4x,e15.6)") string(1:4),0.0
      else if(string(1:3) == 'Mtp') then
        read(string(5:lstr),*) moment_tensor(6)
        do iu=2,4
          write(iu,"(a)") string(1:lstr)
        enddo
        write(5,"(a4,4x,e15.6)") string(1:4),0.0
        write(6,"(a4,4x,e15.6)") string(1:4),0.0
        write(7,"(a4,4x,e15.6)") string(1:4),0.0
        write(8,"(a4,4x,e15.6)") string(1:4),0.0
        write(9,"(a4,4x,e15.6)") string(1:4),0.0
        write(10,"(a4,4x,e15.6)") string(1:4),MOMENT
      endif

    endif

  enddo

  close(1)

  close(2)
  close(3)
  close(4)
  close(5)
  close(6)
  close(7)
  close(8)
  close(9)
  close(10)

  end program make_cmts

! -------------------------------------------------------

  integer function julian_day(yr,mo,da)
  integer yr,mo,da

  integer mon(12)
  integer lpyr
  data mon /0,31,59,90,120,151,181,212,243,273,304,334/

  julian_day=da+mon(mo)
  if(mo>2) julian_day=julian_day+lpyr(yr)

  end function julian_day

! -------------------------------------------------------

  integer function lpyr(yr)
  integer yr
!
!---- returns 1 if yr is a leap year
!
  lpyr=0
  if(mod(yr,400) == 0) then
    lpyr=1
  else if(mod(yr,4) == 0) then
    lpyr=1
    if(mod(yr,100) == 0) then
      lpyr=0
    endif
  endif

  end function lpyr

