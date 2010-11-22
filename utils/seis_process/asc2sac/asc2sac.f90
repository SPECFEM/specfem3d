program asc2sac

  implicit none

  integer, parameter :: NDATAMAX = 40000
  integer :: npts, nl, nerr
  character(len=250) :: cnl, ascfile, sacfile
  real b, dt, data(NDATAMAX)

  call getarg(1,ascfile)
  call getarg(2,cnl)
  call getarg(3,sacfile)
  if (trim(ascfile) == '' .or. trim(cnl) == '' .or. trim(sacfile) == '') &
    stop 'Usage: asc2sac ascfile npts sacfile'
  read(cnl,*) nl

  call rasc(ascfile,data,npts,b,dt,NDATAMAX,nerr)
  if (npts > NDATAMAX .or. npts /= nl) stop 'Check npts'

  call wsac1(sacfile,data,npts,b,dt,nerr)


end program asc2sac
