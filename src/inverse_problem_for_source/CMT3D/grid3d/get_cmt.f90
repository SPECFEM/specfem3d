! read cmt file

subroutine get_cmt(cmt_file,yr,mo,jda,ho,mi,sec, &
     t_cmt,hdur,elat,elon,depth,moment_tensor)

  implicit none

  character(len=*) :: cmt_file
  integer :: yr,jda,ho,mi
  real*8 :: sec,t_cmt,hdur,elat,elon,depth
  real*8 :: moment_tensor(6)

  integer :: i,ios,lstr,mo,da
  real*8 :: mb,ms
  character(len=24) :: reg
  character(len=5) :: datasource
  character(len=150) :: string

  integer :: julian_day

  !
  !---- first read hypocenter info
  !
  open(unit=1,file=cmt_file,iostat=ios,status='old')
  if (ios /= 0) stop 'Error opening CMT file '

  read(1,"(a4,i5,i3,i3,i3,i3,f6.2,f9.4,f10.4,f6.1,f4.1,f4.1,1x,a)",iostat=ios) &
       datasource,yr,mo,da,ho,mi,sec,elat,elon,depth,mb,ms,reg
  if (ios /= 0) stop 'Error reading the information line of the CMT file'
  jda=julian_day(yr,mo,da)
  ios=0
  do while(ios == 0)
     read(1,"(a)",iostat=ios) string
     lstr=len_trim(string)

     if (string(1:10) == 'event name') then
     else if (string(1:10) == 'time shift') then
        read(string(12:lstr),*) t_cmt
     else if (string(1:13) == 'half duration') then
        read(string(15:lstr),*) hdur
     else if (string(1:8) == 'latitude') then
        read(string(10:lstr),*) elat
     else if (string(1:9) == 'longitude') then
        read(string(11:lstr),*) elon
     else if (string(1:5) == 'depth') then
        read(string(7:lstr),*) depth
     else if (string(1:3) == 'Mrr') then
        read(string(5:lstr),*) moment_tensor(1)
     else if (string(1:3) == 'Mtt') then
        read(string(5:lstr),*) moment_tensor(2)
     else if (string(1:3) == 'Mpp') then
        read(string(5:lstr),*) moment_tensor(3)
     else if (string(1:3) == 'Mrt') then
        read(string(5:lstr),*) moment_tensor(4)
     else if (string(1:3) == 'Mrp') then
        read(string(5:lstr),*) moment_tensor(5)
     else if (string(1:3) == 'Mtp') then
        read(string(5:lstr),*) moment_tensor(6)
     endif
  enddo

  close(1)

end subroutine get_cmt

! --------

integer function julian_day(yr,mo,da)

  implicit none

  integer yr,mo,da

  integer mon(12)
  data mon/0,31,59,90,120,151,181,212,243,273,304,334/

  integer :: lpyr

  julian_day=da+mon(mo)
  if (mo > 2) julian_day=julian_day+lpyr(yr)

end function julian_day

! -------

integer function lpyr(yr)

  implicit none

  integer :: yr
  !
  !---- returns 1 if yr is a leap year
  !
  lpyr=0
  if (mod(yr,400) == 0) then
     lpyr=1
  else if (mod(yr,4) == 0) then
     lpyr=1
     if (mod(yr,100) == 0) then
        lpyr=0
     endif
  endif

end function lpyr
