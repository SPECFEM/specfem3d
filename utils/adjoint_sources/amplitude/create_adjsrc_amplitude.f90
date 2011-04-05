program create_adjsrc_amplitude

! this program cuts a certain portion of displacement seismograms and
! converts them into adjoint sources for generating classical amplitude
! kernels following Tromp et al. (2005) eq.67.  
! Modified from cut_velocity.f90 (Qinya Liu, Caltech, May 2007)
! by Ebru, Princeton, March 2011.
!
! call by: ./xcreate_adjsrc_amplitude t1 t2 ifile[0-5] E/N/Z-ascii-files [baz]
!
  implicit none

  integer :: i, is, ie, nstep, j, itime ,ifile,ios, i1, i2, nstep_old
  character(len=256) :: arg(100), file(100)
  character(len=256) :: filename
  integer,parameter :: NMAX = 30000
  real*8, parameter :: EPS = 1.0d-17
  real*8, parameter :: PI = 3.1415926d0
  real*8 :: ts, te, data(5,NMAX), out(NMAX), adj(NMAX), tw(NMAX), norm
  real*8 :: dt, t0, t0_old, dt_old, costh, sinth, th, baz
  logical :: lrot

  i = 1
  lrot = .false.

  ! reads in file arguments
  do while (1 == 1) 
    call getarg(i,arg(i))
    if (i < 6 .and. trim(arg(i)) == '') then
      print*,'Usage: '
      print*,'  xcreate_adjsrc_amplitude t1 t2 ifile[0-5] E/N/Z-ascii-files [baz]'
      print*,'with'
      print*,'  t1: window start time'
      print*,'  t2: window end time'
      print*,'  ifile: 0 = adjoint source calculated for each seismogram component'     
      print*,'  ifile: 1 = adjoint source given by East component only'     
      print*,'  ifile: 2 = adjoint source given by North component'     
      print*,'  ifile: 3 = adjoint source given by Z component'     
      print*,'  ifile: 4 = adjoint source given by rotated transversal component (requires baz)'     
      print*,'  ifile: 5 = adjoint source given by rotated radial component (requires baz)'     
      print*,'  E/N/Z-ascii-files : displacement traces stored as ascii files'
      print*,'  [baz]: (optional) back-azimuth, requires ifile = 4 or ifile = 5'
      stop 'create_adjsrc_amplitude t1 t2 ifile[0-5] E/N/Z-ascii-files [baz]'
    endif
    if (trim(arg(i)) == '') exit
    if (i == 1) then 
      read(arg(i),*,iostat=ios) ts
      if (ios /= 0) stop 'Error reading ts'
    else if (i == 2) then
      read(arg(i),*,iostat=ios) te
      if (ios /= 0) stop 'Error reading te'
    else if (i == 3) then
      read(arg(i),*) ifile
      if (ios /= 0) stop 'Error reading ifile'
    else if (i == 4 .or. i == 5 .or. i == 6) then
      file(i-3) = trim(arg(i))
    else if (i == 7) then
      read(arg(i),*,iostat=ios) baz
      if (ios /= 0) stop 'Error reading baz'
      lrot = .true.
    else if (i > 7) then
      stop 'Error: create_adjsrc_amplitude t1 t2 ifile[0-5] E/N/Z-ascii-files [baz]'
    endif
    i = i + 1
  enddo

  ! checks rotation baz and ifile parameter
  i = i - 1
  if (lrot) then
    if (i /= 7) stop 'create_adjsrc_amplitude t1 t2 ifile[0-5] E/N/Z-ascii-files [baz]'
    if (ifile /= 4 .and. ifile /= 5) stop 'ifile = 4 or 5 when baz is present'
    th = (baz - 180.0) / 180.0 * PI
    costh = cos(th)
    sinth = sin(th)
  else
    if (ifile > 3 .or. ifile < 0) stop 'Error ifile should be between 0 - 3 when baz is not present'
    if (i /= 6) stop 'create_adjsrc_amplitude t1 t2 ifile[0-5] E/N/Z-ascii-files [baz]'
  endif
  
  ! user output
  print *, 'ifile = ', ifile, '  lrot = ', lrot
  print *, ' '

  ! reads seismograms (ascii format)
  do i = 1, 3
    filename = trim(file(i))
    print *, 'reading asc file '//trim(filename)//' ...'
    call dread_ascfile_c(trim(filename)//char(0),t0,dt,nstep,data(i,:))
    if (nstep > NMAX) stop 'Change the data array range limit'
    if (i == 1) then
      t0_old = t0; dt_old = dt; nstep_old = nstep
    else
      if (i > 1 .and. abs(t0_old - t0) > EPS &
               .and. abs(dt_old - dt) > EPS &
               .and. nstep_old /= nstep) &
                 stop 'Error different t0, dt, nstep'
    endif
  enddo
  print *, ' '  
  print *, 'start time:',t0
  print *, 'time step:',dt
  print *, 'number of steps:',nstep   
  print *, ' '
 
  ! component rotation
  if (lrot) then
    data(4,:) = costh * data(1,:) - sinth * data(2,:)
    data(5,:) = sinth * data(1,:) + costh * data(2,:)
    call dwrite_ascfile_c(trim('t.txt')//char(0),t0,dt,nstep,data(4,:))
    call dwrite_ascfile_c(trim('r.txt')//char(0),t0,dt,nstep,data(5,:))
    i1 = 3; i2 = 5
  else
    i1 = 1; i2 = 3
  endif
    
  ! loops over seismogram components 
  do i = i1, i2
    ! start and end index
    is = (ts - t0) / dt + 1
    ie = (te - t0) / dt + 1
    if (is < 1 .or. ie <= is .or. ie > nstep) then
      print *, 'Error in ts, te'; stop
    endif
    
    ! time window (parabola shaped)
    tw(1:nstep) = 0.
    if( i == i1 ) open(44,file='plot_time_window.txt',status='unknown')
    do j = is, ie
      tw(j) = 1 - (2 * (dble(j) - is)/(ie - is) - 1) ** 2
      if( i == i1 ) write(44,*) j,tw(j)
    enddo
    if( i == i1 ) close(44)
    
    ! displacement array
    do itime = 1, nstep
       out(itime) =  data(i,itime) 
    enddo
    
    ! normalization factor
    norm = dt * sum( tw(1:nstep) * out(1:nstep) * out(1:nstep))
    print *, 'i = ', i, 'norm = ', norm
    if (ifile /= 0 .and. ifile /= i) norm = 0.0

    ! adjoint source
    if (abs(norm) > EPS) then
      adj(1:nstep) =  out(1:nstep) * tw(1:nstep) / norm
    else
      print *, 'norm < EPS for file '//trim(file(i))
      adj(:) = 0.
    endif
    data(i,:) = adj(:)
    
  enddo
  print *, ' '
  
  ! component rotation back to cartesian x-y-z
  if (lrot) then
    call dwrite_ascfile_c(trim('t-cut.txt')//char(0),t0,dt,nstep,data(4,:))
    call dwrite_ascfile_c(trim('r-cut.txt')//char(0),t0,dt,nstep,data(5,:))
    data(1,:) = costh * data(4,:) + sinth * data(5,:)
    data(2,:) = -sinth * data(4,:) + costh * data(5,:)
  endif

  ! file output for component BHE/BHN/BHZ
  do i = 1, 3
    filename = trim(file(i))//'.adj'  
    print *, 'write to asc file '//trim(filename)
    call dwrite_ascfile_c(trim(filename)//char(0),t0,dt,nstep,data(i,:))
  enddo

end program create_adjsrc_amplitude


