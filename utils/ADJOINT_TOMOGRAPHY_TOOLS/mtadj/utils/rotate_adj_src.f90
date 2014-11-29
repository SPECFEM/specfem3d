program rotate_adj_src

  ! this program rotates the three component seismograms in the sac format
  ! from (R, T) to (N, E)

  implicit none

  integer, parameter :: NDIM = 40000  ! from mt_constants.f90
  real, parameter :: EPS = 1.0e-2  ! check constistency

  character(len=150) :: bazch, zfile,tfile,rfile,efile,nfile
  integer :: nptst, nptsr, npts, nerr
  logical :: r_exist, t_exist, z_exist
  real :: baz,t0t,dtt,t0r,dtr, t0,dt, costh, sinth
  real, dimension(NDIM) :: rdata, tdata, edata, ndata


  ! read in arguments
  call getarg(1,bazch)
  call getarg(2,zfile)
  call getarg(3,tfile)
  call getarg(4,rfile)
  call getarg(5,efile)
  call getarg(6,nfile)

  if (trim(bazch)=='' .or. trim(zfile) =='' .or. trim(tfile)=='' .or. &
       trim(rfile) == '' .or. trim(efile) =='' .or. trim(nfile) =='') then
     stop 'rotate_adj_src baz(radians) zfile tfile rfile efile nfile'
  endif

  read(bazch,*) baz

  ! check existence of Z, T, R files
  inquire(file=trim(zfile),exist=z_exist)
  inquire(file=trim(tfile),exist=t_exist)
  inquire(file=trim(rfile),exist=r_exist)

  ! only Z component -> t0, dt, npts
  if (.not. t_exist .and. .not. r_exist) then
     if (.not. z_exist) stop 'At least one file should exist: zfile, tfile, rfile'
     call rsac1(zfile,tdata,npts,t0,dt,NDIM,nerr)
     if (nerr > 0) stop 'Error reading Z file'
     tdata = 0.
     rdata = 0.
  endif

  ! read R/T components
  if (t_exist) then
     call rsac1(tfile,tdata,nptst,t0t,dtt,NDIM,nerr)
     if (nerr > 0) stop 'Error reading T file'
     if (.not. r_exist) rdata = 0.
  endif
  if (r_exist) then
     call rsac1(rfile,rdata,nptsr,t0r,dtr,NDIM,nerr)
     if (nerr > 0) stop 'Error reading R file'
     if (.not. t_exist) tdata = 0.
  endif
  ! check consistency of T and R components -> t0,dt,npts
  if (t_exist .and. r_exist) then
     if (abs(t0t-t0r)>EPS .or. abs(dtt-dtr)>EPS .or. nptst /= nptsr) &
          stop 'check t0 and npts'
  endif
  if (t_exist) then
     t0 =t0t; dt = dtt; npts = nptst
  else if (r_exist) then
     t0 =t0r; dt = dtr; npts = nptsr
  endif
  if (NDIM < npts) stop 'Increase NDIM'

  ! [N E]' = [ -costh sinth; -sinth costh ] [R T]'
  costh = cos(baz)
  sinth = sin(baz)
  ndata(1:npts) = - costh * rdata(1:npts) + sinth * tdata(1:npts)
  edata(1:npts) = - sinth * rdata(1:npts) - costh * tdata(1:npts)

  ! write N/E adjoint source
  call wsac1(nfile,ndata,npts,t0,dt,nerr)
  if (nerr > 0) stop 'Error writing N file'
  call wsac1(efile,edata,npts,t0,dt,nerr)
  if (nerr > 0) stop 'Error writing E file'

  if (.not. z_exist) then
     tdata = 0.
     call wsac1(zfile,tdata,npts,t0,dt,nerr)
     if (nerr > 0) stop 'Error writing Z file'
  endif

end program rotate_adj_src
