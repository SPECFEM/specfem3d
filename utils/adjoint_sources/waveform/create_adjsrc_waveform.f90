program create_adjsrc_waveform

! this program cuts certain portion of the seismograms and converts them into
! the adjoint sources for generating banana-dougnut kernels.
!
!
! call by: ./xcreate_adjsrc_waveform t1 t2 ifile[0-5] E/N/Z-ascii-file DATA_DIR/ [baz]
!
  implicit none

  integer,parameter :: NMAX = 30000
  real*8, parameter :: EPS = 1.0d-17
  real*8, parameter :: PI = 3.1415926d0

  real*8 :: trace_syn(5,NMAX),trace_dat(5,NMAX), out(NMAX), adj(NMAX), tw(NMAX)
  real*8 :: ts, te, norm
  real*8 :: dt, t0, t0_old, dt_old, costh, sinth, th, baz

  integer :: i, is, ie, nstep, j, itime ,ifile, ios, i1, i2, nstep_old
  integer :: idx, idx2, nargs, ncomp,itaper_length

  character(len=256) :: arg(100), file(3)
  character(len=256) :: filename
  character(len=256) :: data_dir

  character(len=256) :: basename,basename2,syn_dir
  character(len=64) :: net_sta
  character(len=3) :: channel
  character(len=1) :: comp(3) = (/ 'X', 'Y', 'Z' /)

  logical :: lrot,single_file

  ! initializes
  lrot = .false.
  single_file = .false.

  ! reads in file arguments
  nargs = command_argument_count()
  ! arguments: single file input can have 5 arguments
  !            3 files (for each component) can have 7 or 8 (with baz)
  if (nargs /= 5 .and. nargs /= 7 .and. nargs /= 8) then
    print *,'Usage: '
    print *,'  xcreate_adjsrc_waveform t1 t2 ifile[0-5] E/N/Z-ascii-files DATA_DIR [baz]'
    print *,'with'
    print *,'  t1: window start time'
    print *,'  t2: window end time'
    print *,'  ifile: 0 = adjoint source calculated for each seismogram component'
    print *,'  ifile: 1 = adjoint source given by East component only'
    print *,'  ifile: 2 = adjoint source given by North component'
    print *,'  ifile: 3 = adjoint source given by Z component'
    print *,'  ifile: 4 = adjoint source given by rotated transversal component (requires baz)'
    print *,'  ifile: 5 = adjoint source given by rotated radial component (requires baz)'
    print *,'  E/N/Z-ascii-files : displacement traces (synthetics) stored as ascii files'
    print *,'  DATA_DIR: directory holding (true) data (for waveform misfit; must have same filenames as synthetics)'
    print *,'  [baz]: (optional) back-azimuth, requires ifile = 4 or ifile = 5'

    stop 'create_adjsrc_waveform t1 t2 ifile[0-5] E/N/Z-ascii-files DATA_DIR [baz]'
  endif

  i = 1
  do while (1 == 1)
    call get_command_argument(i,arg(i))
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
    endif

    ! input files
    if (nargs == 5) then
      ! single trace input
      if (i == 4) then
        file(1) = trim(arg(i))
        file(2) = ''
        file(3) = ''
        single_file = .true.
      else if (i == 5) then
        data_dir = trim(arg(i))
        exit
      endif
    else
      ! 3-components
      if (i == 4 .or. i == 5 .or. i == 6) then
        file(i-3) = trim(arg(i))
      else if (i == 7) then
        data_dir = trim(arg(i))
      else if (i == 8) then
        read(arg(i),*,iostat=ios) baz
        if (ios /= 0) stop 'Error reading baz'
        lrot = .true.
      else if (i > 8) then
        stop 'Error: create_adjsrc_waveform t1 t2 ifile[0-5] E/N/Z-ascii-files DATA_DIR [baz]'
      endif
    endif
    i = i + 1
  enddo

  ! checks rotation baz and ifile parameter
  if (lrot) then
    if (ifile /= 4 .and. ifile /= 5) stop 'ifile = 4 or 5 when baz is present'
    th = (baz - 180.0) / 180.0 * PI
    costh = cos(th)
    sinth = sin(th)
  else
    if (ifile > 3 .or. ifile < 0) stop 'Error ifile should be between 0 - 3 when baz is not present'
  endif

  ! user output
  print *, 'xcreate_adjsrc_waveform:'
  print *, '  measurement window start/end = ',ts,'/',te
  print *, '  component ifile = ', ifile, '  lrot = ', lrot
  print *, '  single file = ',single_file

  ! gets directory from filename
  filename = trim(file(1))
  ! gets directory name
  ! example: filename = ../OUTPUT_FILES/DB.X1.MXX.semd  -> basename = DB.X1.MXX.semd
  idx = index(filename,'/',.true.)  ! position of '/' reading from back
  if (idx == 0) then
    ! no / found
    syn_dir = ''
  else if (idx < len_trim(filename)) then
    syn_dir = filename(1:idx)
  else
    stop 'error getting basename, please check file names...'
  endif
  print *, '  data directory : ',trim(data_dir)
  print *, '  synthetics directory (derived from input filenames): ',trim(syn_dir)
  print *, ' '

  ! number of components to read
  if (single_file) then
    ncomp = 1
  else
    ncomp = 3
  endif

  ! reads seismograms (ascii format)
  trace_syn(:,:) = 0.0
  trace_dat(:,:) = 0.0
  do i = 1, ncomp
    filename = trim(file(i))
    ! gets basename
    ! example: filename = ../OUTPUT_FILES/DB.X1.MXX.semd  -> basename = DB.X1.MXX.semd
    idx = index(filename,'/',.true.)
    if (idx == 0) then
      ! no / found
      basename = trim(filename)
    else if (idx < len_trim(filename)) then
      basename = filename(idx+1:len_trim(filename))
    else
      stop 'error getting basename, please check file names...'
    endif
    ! gets channel/component from basename
    ! example: basename = DB.X1.MXX.semd -> format net.sta.comp.ending has channel MXX
    ! checks if basename has '.'
    idx = index(basename,'.')  ! DB.***
    if (idx == 0) stop 'error getting component from basename, file name must have format "net.sta.comp.ending"'
    ! second occurrence
    basename2 = basename(idx+1:len_trim(basename))
    idx2 = index(basename2,'.')  ! X1.***
    channel = basename2(idx2+1:idx2+3)
    !print *,'debug: ',trim(basename),' ',trim(basename2),' ',trim(channel),' component: ',channel(3:3)

    ! synthetics
    print *, 'reading synthetics asc file '//trim(filename)//' ...'
    call dread_ascfile_c(trim(filename)//char(0),t0,dt,nstep,trace_syn(i,:))
    if (nstep > NMAX) stop 'Change the synthetic trace array range limit'
    if (i == 1) then
      t0_old = t0; dt_old = dt; nstep_old = nstep
    else
      if (i > 1 .and. abs(t0_old - t0) > EPS &
                .and. abs(dt_old - dt) > EPS &
                .and. nstep_old /= nstep) &
        stop 'Error synthetics with different t0, dt, nstep'
    endif

    ! data
    filename = trim(data_dir)//'/'//trim(basename)
    print *, 'reading data asc file '//trim(filename)//' ...'
    call dread_ascfile_c(trim(filename)//char(0),t0,dt,nstep,trace_dat(i,:))
    if (nstep > NMAX) stop 'Change the data trace array range limit'
    if (abs(t0_old - t0) > EPS &
        .and. abs(dt_old - dt) > EPS &
        .and. nstep_old /= nstep) then
      print *,'Error: data trace has different t0,dt,nstep ',t0,dt,nstep
      print *,'       must match synthetic trace with t0,dt,nstep ',t0_old,dt_old,nstep_old
      stop 'Error data has different t0, dt, nstep; length and sampling rate must match between data and synthetics'
    endif
  enddo
  print *, '  start time     :',t0
  print *, '  time step      :',dt
  print *, '  number of steps:',nstep

  ! component rotation
  if (lrot) then
    ! synthetics
    trace_syn(4,:) = costh * trace_syn(1,:) - sinth * trace_syn(2,:)
    trace_syn(5,:) = sinth * trace_syn(1,:) + costh * trace_syn(2,:)
    call dwrite_ascfile_c(trim('syn_t.txt')//char(0),t0,dt,nstep,trace_syn(4,:))
    call dwrite_ascfile_c(trim('syn_r.txt')//char(0),t0,dt,nstep,trace_syn(5,:))
    ! data
    trace_dat(4,:) = costh * trace_dat(1,:) - sinth * trace_dat(2,:)
    trace_dat(5,:) = sinth * trace_dat(1,:) + costh * trace_dat(2,:)
    call dwrite_ascfile_c(trim('dat_t.txt')//char(0),t0,dt,nstep,trace_dat(4,:))
    call dwrite_ascfile_c(trim('dat_r.txt')//char(0),t0,dt,nstep,trace_dat(5,:))
    i1 = 3; i2 = 5
  else
    i1 = 1; i2 = ncomp
  endif

  ! loops over seismogram components
  do i = i1, i2
    ! start and end index
    is = int((ts - t0) / dt) + 1
    ie = int((te - t0) / dt) + 1
    if (is < 1 .or. ie <= is .or. ie > nstep) then
      print *, 'Error in ts, te'; stop
    endif

    ! taper
    ! time window (parabola shaped)
    tw(1:nstep) = 0.0
    if ( i == i1 ) open(44,file='plot_time_window.txt',status='unknown')
    ! select taper
    if (1 == 0) then
      ! parabola taper
      do j = is, ie
        tw(j) = 1 - (2 * (dble(j) - is)/(ie - is) - 1) ** 2
        if ( i == i1 ) write(44,*) j,tw(j)
      enddo
    else
      ! cosine taper
      itaper_length = int( 0.1*(ie - is)) ! 10 percent taper length
      do j = is, ie
        if (j < is + itaper_length) then
          tw(j) = (1.0-cos(pi*2.0*(j-is)/(itaper_length-1)))/2.0
        else if (j <= ie - itaper_length) then
          tw(j) = 1.0
        else
          tw(j) = (1.0-cos(pi*2.0*(ie-j)/(itaper_length-1)))/2.0
        endif
      enddo
    endif
    if ( i == i1 ) close(44)

    ! waveform misfit
    ! note: Tromp et al. 2005, eq. (9), defines [ s_i(T-t) - d_i(T-t) ], i.e., synthetics minus data.
    !       in exploration, this often switched to data minus synthetics.
    !
    !       also, the adjoint source in SEM/ is not time-reversed yet as in eq. (9).
    !       to account for this, SPECFEM will read in the trace in reverse order.
    out(:) = 0.0
    do itime = 1, nstep
       out(itime) =  trace_syn(i,itime) - trace_dat(i,itime)
    enddo

    ! normalization factor
    ! (will not be used to scale adjoint source, only for information here)
    norm = dt * sum( tw(1:nstep) * out(1:nstep) * out(1:nstep))
    print *, 'i = ', i, 'component = ',channel(3:3), ' norm = ', norm
    ! zero out unwanted component
    if (ifile /= 0) then
      if (single_file) then
        if (comp(ifile) /= channel(3:3)) then
          norm = 0.0
          out(:) = 0.0
        endif
      else
        if (ifile /= i) then
          norm = 0.0
          out(:) = 0.0
        endif
      endif
    endif

    ! adjoint source (windowed)
    adj(:) = 0.0
    if (abs(norm) > EPS) then
      adj(1:nstep) = out(1:nstep) * tw(1:nstep)
    else
      ! user output
      if (single_file) then
        if (ifile /= 0 .and. comp(ifile) /= channel(3:3)) then
          print *,'  component set to zero'
        else
          print *, '  norm < EPS for file '//trim(file(i))
        endif
      else
        if (ifile /= 0 .and. ifile /= i) then
          print *,'  component set to zero'
        else
          print *, '  norm < EPS for file '//trim(file(i))
        endif
      endif
      ! setting to zero
      adj(:) = 0.0
    endif

    ! store adjoint source in trace array
    trace_syn(i,:) = adj(:)

  enddo

  ! component rotation back to Cartesian x-y-z
  if (lrot) then
    call dwrite_ascfile_c(trim('syn_t-cut.txt')//char(0),t0,dt,nstep,trace_syn(4,:))
    call dwrite_ascfile_c(trim('syn_r-cut.txt')//char(0),t0,dt,nstep,trace_syn(5,:))
    trace_syn(1,:) = costh * trace_syn(4,:) + sinth * trace_syn(5,:)
    trace_syn(2,:) = -sinth * trace_syn(4,:) + costh * trace_syn(5,:)
  endif

  ! file output for component BHE/BHN/BHZ
  do i = 1, ncomp
    filename = trim(file(i))
    ! gets basename
    ! example: filename = ../OUTPUT_FILES/DB.X1.MXX.semd  -> basename = DB.X1.MXX.semd
    idx = index(filename,'/',.true.)  ! position of '/' reading from back
    if (idx == 0) then
      ! no / found
      syn_dir = ''
      basename = trim(filename)
    else if (idx < len_trim(filename)) then
      syn_dir = filename(1:idx)
      basename = filename(idx+1:len_trim(filename))
    endif
    ! gets channel/component from basename
    ! example: basename = DB.X1.MXX.semd -> format net.sta.comp.ending has channel MXX
    ! checks if basename has '.'
    idx = index(basename,'.')  ! DB.***
    if (idx == 0) stop 'error getting component from basename, file name must have format "net.sta.comp.ending"'
    ! second occurrence
    basename2 = basename(idx+1:len_trim(basename))
    idx2 = index(basename2,'.')  ! X1.***
    net_sta = basename(1:(idx+idx2)-1)     ! DB.X1
    channel = basename2(idx2+1:idx2+3) ! MXX

    !filename = trim(file(i))//'.adj'
    filename = trim(syn_dir)//trim(net_sta)//'.'//trim(channel)//'.adj'

    print *, 'write to adjoint source asc file '//trim(filename)
    call dwrite_ascfile_c(trim(filename)//char(0),t0,dt,nstep,trace_syn(i,:))
  enddo

end program create_adjsrc_waveform


