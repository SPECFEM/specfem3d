
  program convert_raw_Wald

  implicit none

! number of values per line in Rob Graves' data files
  integer, parameter :: NVAL = 6

  integer nstep,hour,minutes,nlines,iline,it,ival

  double precision a(NVAL)
  double precision seconds,dt,initial_time

  open(unit=11,file='tutu',status='old')
  open(unit=12,file='tutu2',status='unknown')

! skip title
  read(11,*)

! read timing info
  read(11,*) nstep,dt,hour,minutes,seconds

! compute number of lines of NVAL values
  nlines = nstep / NVAL

! initialize time step number
  it = 0

! initial time
  initial_time = 3600 * hour + 60 * minutes + seconds

  print *,'initial time = ',initial_time
  print *

! loop on all the lines
  do iline = 1,nlines

    read(11,*) (a(ival),ival=1,NVAL)

! unit is cm/s/s, convert to m/s/s (i.e., SI)

    do ival = 1,NVAL
      it = it + 1
      write(12,*) sngl((it-1)*dt + initial_time),sngl(a(ival)/100.d0)
    enddo

  enddo

  close(11)
  close(12)

  end program convert_raw_Wald

