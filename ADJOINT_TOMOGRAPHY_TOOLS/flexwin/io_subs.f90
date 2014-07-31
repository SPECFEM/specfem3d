!----------------------------------------------------------------------

  subroutine write_seis_gmt(basename)
  use seismo_variables
  implicit none

  character*120 :: basename
  character*240 :: file_obs, file_syn
  double precision :: max_value_obs, max_value_syn, max_value

  integer :: i,ier

  ! create the filenames
  file_obs=trim(basename)//'.obs'
  file_syn=trim(basename)//'.syn'

  ! find the maximum value
  !max_value_obs = max (maxval(obs_lp), abs(minval(obs_lp)) )
  !max_value_syn = max (maxval(synt_lp), abs(minval(synt_lp)) )
  max_value_obs = maxval( abs(obs_lp) )
  max_value_syn = maxval( abs(synt_lp) )
  max_value = max(max_value_obs, max_value_syn)

  ! open the files for writing
  open(unit=11, file=file_obs,iostat=ier)
  if( ier /= 0 ) then
    print*,'Error opening: ',trim(file_obs)
    stop
  endif
  open(unit=12, file=file_syn,iostat=ier)
  if( ier /= 0 ) then
    print*,'Error opening: ',trim(file_syn)
    stop
  endif

  ! write the header - observed
  write(11,'("# NPTS = ",i10)') npts
  write(11,'("# PLOT_MAX = ",e12.6)') max_value
  write(11,'("# T_START = ",f10.2)') b
  write(11,'("# T_END = ",f10.2)') b+(npts-1)*dt

  ! write the header - synthetic
  write(12,'("# NPTS = ",i10)') npts
  write(12,'("# PLOT_MAX = ",e12.6)') max_value
  write(12,'("# T_START = ",f10.2)') b
  write(12,'("# T_END = ",f10.2)') b+(npts-1)*dt

  ! write the time series
  do i = 1, npts
    write(11,'(f10.2,2x,e12.6)') b+(i-1)*dt, obs_lp(i)
    write(12,'(f10.2,2x,e12.6)') b+(i-1)*dt, synt_lp(i)
  end do

  ! close the files
  close(11)
  close(12)

  end subroutine

!----------------------------------------------------------------------

  subroutine write_env_gmt(basename)
  use seismo_variables
  implicit none

  character*120 :: basename
  character*240 :: file_obs_env, file_syn_env
  double precision :: max_value_obs, max_value_syn, max_value

  integer :: i

  ! create the filenames
  file_obs_env=trim(basename)//'.env.obs'
  file_syn_env=trim(basename)//'.env.syn'

  ! find the maximum value
  max_value_obs = maxval(env_obs_lp)
  max_value_syn = maxval(env_synt_lp)
  max_value = max(max_value_obs, max_value_syn)

  ! open the files for writing
  open(unit=13, file=file_obs_env)
  open(unit=14, file=file_syn_env)

  ! write the header - observed envelope
  write(13,'("# NPTS = ",i10)') npts
  write(13,'("# PLOT_MAX = ",e12.6)') max_value
  write(13,'("# T_START = ",f10.2)') b
  write(13,'("# T_END = ",f10.2)') b+(npts-1)*dt

  ! write the header - synthetic envelope
  write(14,'("# NPTS = ",i10)') npts
  write(14,'("# PLOT_MAX = ",e12.6)') max_value
  write(14,'("# T_START = ",f10.2)') b
  write(14,'("# T_END = ",f10.2)') b+(npts-1)*dt

  ! write the time series
  do i = 1, npts
    write(13,'(f10.2,2x,e12.6)') b+(i-1)*dt, env_obs_lp(i)
    write(14,'(f10.2,2x,e12.6)') b+(i-1)*dt, env_synt_lp(i)
  end do

  ! close the files
  close(13)
  close(14)

  end subroutine

!----------------------------------------------------------------------

  subroutine write_win_gmt(basename)
  use seismo_variables
  implicit none

  character*120 :: basename
  character*240 :: file_win
  integer, dimension(MAX_PHASES) :: phase_indexes
  integer :: i, k, n_phases

  ! set filename
  file_win=trim(basename)//'.win'

  open(unit=11, file=file_win)
  write(11,'("# NUM_WIN = ",i10)') num_win
  do i = 1, num_win
    call phases_in_window(win_start(i), win_end(i), n_phases, phase_indexes)
    write(11,'(i4,1x,2f10.2,1x,i4,60(1x,a8))') i, win_start(i), win_end(i), n_phases, (ph_names(phase_indexes(k)), k=1,n_phases)
  end do
  close(11)

  end subroutine

!----------------------------------------------------------------------

  subroutine write_win_qual_gmt(basename)
  use seismo_variables
  implicit none

  character*120 :: basename
  character*240 :: file_win
  integer :: i

  ! set filename
  file_win=trim(basename)//'.win.qual'

  open(unit=11, file=file_win)
  write(11,'("# NUM_WIN = ",i10)') num_win
  write(11,*) "# i win_start win_end Tshift CC dlnA "
  do i = 1, num_win
    write(11,'(i4,5(1x,f10.5))') i, win_start(i), win_end(i), &
                                 Tshift(i),CC(i),dlnA(i)
  end do
  close(11)

  end subroutine

!----------------------------------------------------------------------

  subroutine read_win_qual_gmt(basename)
  use seismo_variables
  implicit none

  character*120 :: basename
  character*240 :: file_win
  integer :: i, i_dummy
  character*120 :: c_dummy

  ! set filename
  file_win=trim(basename)//'.win.qual'

  open(unit=11, file=file_win, status='old')
!  read(11,'("# NUM_WIN = ",i10)') num_win
  read(11,'(a12,i10)') c_dummy, num_win
  read(11,'(a)') c_dummy
  do i = 1, num_win
  ! read file, ignoring input "quality" values
    read(11,'(i4,5(1x,f10.5))') i_dummy, win_start(i), win_end(i), &
                                 Tshift(i),CC(i),dlnA(i)
    ! calculate indexes for start and and of windows
    i_start(i)=1+int((win_start(i)-b)/dt)
    i_end(i)  =1+int((win_end(i)  -b)/dt)
  end do
  close(11)

  end subroutine

!----------------------------------------------------------------------

  subroutine write_phases_gmt(basename)
  use seismo_variables
  implicit none

  character*120 :: basename
  character*240 :: file_phases
  integer :: i

  ! set filename
  file_phases=trim(basename)//'.phases'

  open(unit=11, file=file_phases)
  write(11,'("# NUM_PHASES = ",i10)') num_phases
  write(11,'(a)') "# PHASE    T_TIME(s)"
  do i = 1, num_phases
    write(11,'(2x,a8,1x,f9.2)') ph_names(i), ph_times(i)
  enddo

  close(11)

  end subroutine

!----------------------------------------------------------------------

  subroutine write_stalta_gmt(basename)
  use seismo_variables
  implicit none

  character*120 :: basename
  character*240 :: file_stalta
  double precision :: max_qual, max_tshift, max_dlna, max_sta_lta, max_s2n
  integer :: i

  ! set filename
  file_stalta=trim(basename)//'.stalta'

  ! get the max amplitude of envelopes
  max_qual = maxval(CC_LIMIT)
  max_tshift = maxval(TSHIFT_LIMIT)
  max_dlna = maxval(DLNA_LIMIT)
  max_sta_lta = maxval(STA_LTA)
  max_s2n = maxval(S2N_LIMIT)

  ! open the files
  open(unit=15, file=file_stalta)
  ! write the header - with f2 quality
  write(15,'("# NPTS = ",i10)') npts
  write(15,'("# CC_MAX = ",f10.4)') max_qual
  write(15,'("# TSHIFT_LIMIT = ",f10.4)') max_tshift
  write(15,'("# DLNA_MAX = ",f10.4)') max_dlna
  write(15,'("# STA_LTA_MAX = ",f10.4)') max_sta_lta
  write(15,'("# S2N_LIMIT = ",f10.4)') max_s2n
  write(15,'("# T_START = ",f10.2)') b
  write(15,'("# T_END = ",f10.2)') b+(npts-1)*dt

  do i = 1, npts
    write(15,'(f10.2,6(2x,e12.6))') b+(i-1)*dt, STA_LTA(i), STALTA_W_LEVEL(i), &
    CC_LIMIT(i), TSHIFT_LIMIT(i), DLNA_LIMIT(i), S2N_LIMIT(i)
  end do

  close(15)

  end subroutine

!----------------------------------------------------------------------

  subroutine write_info_gmt(basename)
  use seismo_variables
  implicit none

  character*120 :: basename
  character*240 :: file_info

  ! set filename
  file_info=trim(basename)//'.info'

  ! set filenames for output (use full path)
  file_info=trim(basename)//'.info'

  ! write the seismograms and envelopes and f1f2
  ! open the files
  open(unit=16, file=file_info)
  ! write the header - info
  write(16,'("# NPTS = ",i10)') npts
  write(16,'("# DT = ",f10.2)') dt
  write(16,'("# B = ",f10.2)') b
  write(16,'("# EVLA = ",f10.4)') evla
  write(16,'("# EVLO = ",f10.4)') evlo
  write(16,'("# STLA = ",f10.4)') stla
  write(16,'("# STLO = ",f10.4)') stlo
  write(16,'("# EVDP = ",f10.4)') evdp
  write(16,'("# AZ = ",f10.2)') azimuth
  write(16,'("# BAZ = ",f10.2)') backazimuth
  write(16,'("# DDG = ",f10.4)') dist_deg
  write(16,'("# DKM = ",f12.4)') dist_km
  write(16,'("# NUM_WIN = ",i10)') num_win

  close(16)

  end subroutine

!----------------------------------------------------------------------

  subroutine write_seismos_gmt(basename)
  use seismo_variables
  implicit none

  character*120 :: basename

  ! write the seismograms in sac
  if (RUN_BANDPASS)  call write_lp_sac(basename)

  ! write the seismograms for gmt
  call write_seis_gmt(basename)

  ! write the envelopes
  call write_env_gmt(basename)

  ! write sta lta time series
  call write_stalta_gmt(basename)

  ! write information
  call write_info_gmt(basename)

  ! write the windows (with phase information)
  call write_win_gmt(basename)

  ! write the windows with quality information
  call write_win_qual_gmt(basename)

  ! write the phases
  call write_phases_gmt(basename)

  end subroutine write_seismos_gmt

! -----------------------------------------------------------------------
  subroutine write_lp_sac(basename)
  use seismo_variables
  implicit none

  character*120 :: basename
  character*240 :: file_obs, file_syn
  integer :: io_error, i, nerr
  real, dimension(NDIM) :: obs_filt, syn_filt, time

  obs_filt(1:npts) = sngl(obs_lp(1:npts))
  syn_filt(1:npts) = sngl(synt_lp(1:npts))
  do i = 1, npts
    !time(i) = b + (i-1)*dt       ! compilation warning
    time(i) = sngl(b + (i-1)*dt)
  enddo

  file_obs=trim(basename)//'.obs_lp.sac'
  file_syn=trim(basename)//'.syn_lp.sac'

  call newhdr()
  call setkhv('kstnm',trim(kstnm),nerr)
  call setkhv('knetwk',trim(knetwk),nerr)
  call setkhv('kcmpnm',trim(kcmpnm),nerr)
  call setfhv('b',sngl(b),nerr)
  call setfhv('delta',sngl(dt),nerr)
  call setnhv('npts',npts,nerr)
  call setihv('iftype','ITIME',nerr)
  call setihv('iztype','IB',nerr)
  call setlhv('leven',1,nerr)

  call wsac0(file_obs,time(1:npts),obs_filt(1:npts),io_error)
  call wsac0(file_syn,time(1:npts),syn_filt(1:npts),io_error)

  end subroutine

!----------------------------------------------------------------------

  subroutine slash_index(filename,spos)
  implicit none
  character*120, intent(in) :: filename
  integer, intent(out) :: spos

  integer :: pos, prev_pos, k

  spos=0
  k=0
  pos=index(filename,'/')
  if (pos.eq.0) then
    spos = 0
  else
   do while (pos .ne. 0)
     k=k+1
     prev_pos = pos
     pos = index(filename(prev_pos+1:len_trim(filename)), '/')
   enddo
   spos = prev_pos
  endif

  end subroutine slash_index

! -----------------------------------------------------------------------

  subroutine write_mt_input(basename)
    use seismo_variables
    implicit none

    character*120 :: basename
    character*240 :: file_s, file_o
    character*240 :: file_mt

    integer :: i

    file_o=trim(basename)//'.obs_lp.sac'
    file_s=trim(basename)//'.syn_lp.sac'

    file_mt=trim(basename)//'.mt_input'
    open(unit=11, file=file_mt)

    !write(11,*) 1
    write(11,'(a)') trim(file_o)
    write(11,'(a)') trim(file_s)
    write(11,*) num_win
    do i = 1, num_win
       write(11,'(2(f10.4,1x))') win_start(i), win_end(i)
    end do
    close(11)

  end subroutine write_mt_input
  ! -----------------------------------------------------------------------

  subroutine write_mt_input_2(basename,file_o,file_s)
    use seismo_variables
    implicit none

    character*120 :: basename
    character*240 :: file_o, file_s
    character*240 :: file_mt

    integer :: i

    if(num_win > 0) then

       file_mt=trim(basename)//'.mt_input'
       open(unit=11, file=file_mt)

       write(11,'(a)') trim(file_o)
       write(11,'(a)') trim(file_s)
       write(11,*) num_win
       do i = 1, num_win
          write(11,'(2(f10.4,1x))') win_start(i), win_end(i)
       enddo
       close(11)

    endif

  end subroutine write_mt_input_2

! -----------------------------------------------------------------------
