!------------------------------------------------

subroutine rasc(ascfile,dat,npts,b,dt,nmax,nerr)

  implicit none
  character(len=*) ascfile
  real dat(1), b, dt
  integer npts, nmax, nerr

  real tv
  integer nline, ios, IOUNIT

  IOUNIT=99

  open(IOUNIT,file=trim(ascfile),status='old',iostat=ios)
  if (ios /= 0) then
     print 'Error opening '//trim(ascfile); nerr=1
     return
  endif

  nline=0
  do while(ios == 0)
     read(IOUNIT,*,iostat=ios) tv,dat(nline+1)
     if (ios == 0) nline=nline+1
     if (nline==1) then
        b=tv
     else if (nline == 2) then
        dt=tv-b
     endif
  enddo

  if (nline > nmax) then
     print *, 'npts exceeding nmax '; nerr=2
     return
  endif

  npts=nline
  nerr=0
  close(IOUNIT)

end subroutine rasc

subroutine wasc(ascfile,dat,npts,b,dt,nerr)

  implicit none

  character(len=*) :: ascfile
  real :: dat(1), b, dt
  integer :: npts, nerr, i, ios, IOUNIT

  IOUNIT=99
  open(IOUNIT,file=trim(ascfile),status='unknown',iostat=ios)
  if (ios /= 0) then
     print 'Error opening '//trim(ascfile)//' for writing'
     nerr = 1; return
  endif
  do i = 1, npts
     write(IOUNIT,'(2g15.7)') b+(i-1)*dt, dat(i)
  enddo
  close(IOUNIT)

  nerr=0

end subroutine wasc


