program rotate_adj_src

  use ascii_rw   ! ascii read and write module

  implicit none

  character(len=150) :: bazch, zfile,tfile,rfile,efile,nfile

  double precision :: baz
  integer :: nptst, nptsr, npts
  logical :: r_exist, t_exist, z_exist
  double precision :: t0t,dtt,t0r,dtr, t0,dt, costh, sinth
  integer, parameter :: NDIM = 40000  ! check ma_constants.f90
  double precision, dimension(NDIM) :: rdata, tdata, edata, ndata

  call getarg(1,bazch)
  call getarg(2,zfile)
  call getarg(3,tfile)
  call getarg(4,rfile)
  call getarg(5,efile)
  call getarg(6,nfile)

  if (trim(bazch)=='' .or. trim(zfile) =='' .or. trim(tfile)=='' .or. &
             trim(rfile) == '' .or. trim(efile) =='' .or. trim(nfile) =='') then
    stop 'rotate_adj_src baz(radian!) zfile tfile rfile efile nfile'
  endif

  read(bazch,*) baz

  inquire(file=trim(tfile),exist=t_exist)
  inquire(file=trim(rfile),exist=r_exist)
  inquire(file=trim(zfile),exist=z_exist)

  ! initialization
  rdata = 0; tdata = 0

  ! at least one file (T,R,Z) should be present
  if (.not. t_exist .and. .not. r_exist) then
    if (.not. z_exist) stop 'At least one file should exist: zfile, tfile, rfile'
  ! need to read Z comp adjoint source for [to,dt,npts]
    call drascii(zfile,tdata,npts,t0,dt)
  endif

  ! read in T file
  if (t_exist) then
    call drascii(tfile,tdata,nptst,t0t,dtt)
  endif
  ! read in R file
  if (r_exist) then
    call drascii(rfile,rdata,nptsr,t0r,dtr)
  endif

  ! check consistency of t0,dt,npts
  if (t_exist .and. r_exist) then
    if (abs(t0t-t0r)>1.0e-2 .or. abs(dtt-dtr)>1.0e-2 .or. nptst /= nptsr) &
               stop 'check t0 and npts'
  endif
  if (t_exist) then
    t0 =t0t; dt = dtt; npts = nptst
  else if (r_exist) then
    t0 =t0r; dt = dtr; npts = nptsr
  endif
  if (NDIM < npts) stop 'Increase NDIM'

  ! rotate T,R to E,N based on baz (note in radian!)
  costh = cos(baz)
  sinth = sin(baz)
  edata(1:npts) = -costh * tdata(1:npts) - sinth * rdata(1:npts)
  ndata(1:npts) =  sinth * tdata(1:npts) - costh * rdata(1:npts)

  ! write E,N files
  call dwascii(efile,edata,npts,t0,dt)
  call dwascii(nfile,ndata,npts,t0,dt)

  ! write Z file if did not exist
  if (.not. z_exist) then
    tdata = 0.
    call dwascii(zfile,tdata,npts,t0,dt)
  endif

end program rotate_adj_src
