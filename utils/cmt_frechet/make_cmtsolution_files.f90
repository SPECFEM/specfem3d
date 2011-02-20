!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 2
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and University of Pau / CNRS / INRIA
!         (c) California Institute of Technology July 2004
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

  program make_cmts

! Jeroen Tromp, July 2001

  implicit none

  include "cmt.h"

  integer yr,jda,ho,mi
  double precision sec,tshift_cmt,hdur,lat,long,depth
  double precision moment_tensor(6)
  character(len=256) cmt_file

  integer iu,i,ios,lstr,mo,da,julian_day
  double precision mb,ms
  double precision latp,longp
  character(len=24) reg
  character(len=5) datasource
  character(len=256) string

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
           datasource,yr,mo,da,ho,mi,sec,lat,long,depth,mb,ms,reg

  do iu=2,10
    write(iu,"(a3,i5,i3,i3,i3,i3,f6.2,f9.4,f10.4,f6.1,f4.1,f4.1,1x,a)") &
             datasource,yr,mo,da,ho,mi,sec,lat,long,depth,mb,ms,reg
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
        read(string(12:lstr),*) tshift_cmt
        do iu=2,10
          write(iu,"(a)") string(1:lstr)
        enddo
      else if(string(1:13) == 'half duration') then
        read(string(15:lstr),*) hdur
        do iu=2,10
          write(iu,"(a)") string(1:lstr)
        enddo
      else if(string(1:8) == 'latitude') then
        read(string(10:lstr),*) lat
        latp = lat + DDELTA
        if(latp > 90.0) latp = 180.0 - latp
        write(2,"(a9,5x,f9.4)") string(1:9),latp
        do iu=3,10
          write(iu,"(a)") string(1:lstr)
        enddo
      else if(string(1:9) == 'longitude') then
        read(string(11:lstr),*) long
        write(2,"(a)") string(1:lstr)
        longp = long + DDELTA
        if(longp > 180.0) longp = longp - 360.0
        write(3,"(a10,4x,f9.4)") string(1:10),longp
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

