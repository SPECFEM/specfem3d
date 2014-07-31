
!----------------------------------------------------------------------

  module seismo_variables
  use user_parameters
  implicit none

  !=========================================
  ! Declaration of parameters from PAR_FILE
  !=========================================

  ! -------------------------------------------------------------
  ! boolean parameters
  logical :: DEBUG
  logical :: MAKE_SEISMO_PLOTS
  logical :: MAKE_WINDOW_FILES
  logical :: BODY_WAVE_ONLY

  ! -------------------------------------------------------------
  ! period min/max for filtering
  logical :: RUN_BANDPASS
  double precision :: WIN_MIN_PERIOD
  double precision :: WIN_MAX_PERIOD
  double precision :: FSTART
  double precision :: FEND

  ! -------------------------------------------------------------
  ! E(t) water level
    double precision :: STALTA_BASE

  ! -------------------------------------------------------------
  ! TSHIFT
    double precision :: TSHIFT_BASE, TSHIFT_REFERENCE

  ! -------------------------------------------------------------
  ! limit on dlnA (dA/A) for window acceptance
    double precision :: DLNA_BASE, DLNA_REFERENCE

  ! -------------------------------------------------------------
  ! limit on CC for window acceptance
    double precision :: CC_BASE

  ! -------------------------------------------------------------
  ! limit on signal-to-noise on the observed data

  ! boolean switch for check_data_quality
    logical :: DATA_QUALITY

  ! if DATA_QUALITY = .true. and if two different measurements of
  ! signal-to-noise ratios exceeds these two base levels,
  ! then the data time series (and syn) is kept
    double precision :: SNR_INTEGRATE_BASE
    double precision :: SNR_MAX_BASE

  ! -------------------------------------------------------------
  ! limit on signal to noise ratio in a particular window.
    double precision :: WINDOW_S2N_BASE

  ! -------------------------------------------------------------
  ! Fine tuning constants
    double precision :: C_0
    double precision :: C_1
    double precision :: C_2
    double precision :: C_3a
    double precision :: C_3b
    double precision :: C_4a
    double precision :: C_4b

    double precision :: WEIGHT_SPACE_COVERAGE
    double precision :: WEIGHT_AVERAGE_CC
    double precision :: WEIGHT_N_WINDOWS

  !=========================================
  ! Global variables and arrays
  !=========================================

  ! un-windowed synthetic and observed seismograms
  double precision, dimension (NDIM) :: synt, obs

  ! un-windowed synthetic and observed seismograms filtered
  double precision, dimension (NDIM) :: synt_lp, obs_lp

  ! envelopes of un-windowed, filtered, synthetic and observed seismograms
  double precision, dimension (NDIM) :: env_synt_lp, env_obs_lp

  ! STA/LTA array
  double precision, dimension (NDIM) :: STA_LTA

  ! time-dependent criteria arrays for selection parameters
  double precision, dimension (NDIM) :: STALTA_W_LEVEL, TSHIFT_LIMIT, CC_LIMIT, DLNA_LIMIT, S2N_LIMIT

  ! seismogram parameters
  double precision :: dt, b
  character*8 :: kstnm,knetwk,kcmpnm
  integer :: npts
  real :: evla, evlo, stla, stlo, evdp
  real :: azimuth, backazimuth, dist_deg, dist_km
  real :: P_pick, S_pick
  integer :: n_p ! index of the first arrival (filtered)

  ! window parameters
  integer :: num_win
  integer, dimension(NWINDOWS) :: i_start, i_end
  double precision, dimension(NWINDOWS) :: win_start, win_end
  double precision, dimension(NWINDOWS) :: F1, F2, Tshift, CC, dlnA
  double precision, dimension(NWINDOWS) :: F1_after, F2_after, Tshift_after, CC_after, dlnA_after
  double precision, dimension(NWINDOWS) :: F1_aux, F2_aux, Tshift_aux, CC_aux, dlnA_aux

  ! phase parameters for ttimes
  integer, parameter :: MAX_PHASES = 60
  integer :: num_phases
  character*8, dimension(MAX_PHASES) :: ph_names
  double precision, dimension(MAX_PHASES) :: ph_times

  ! start/end times for signal and noise
  double precision :: noise_start, noise_end, signal_start, signal_end

  ! indices of start/end times for signal and noise
  integer :: in1, in2, is1, is2

  end module seismo_variables

!----------------------------------------------------------------------

  subroutine read_parameter_file()
  use seismo_variables
  implicit none

  integer, parameter :: IIN = 11
  integer, parameter :: NHEAD = 12

  integer :: idummy
  character(len=34) junk

  open(unit=IIN,file='PAR_FILE',status='old')

  ! ignore header
  do idummy=1,NHEAD
    read(IIN,*)
  enddo

  !--------------------------------------------------------
  ! read parameters, skipping empty lines and comment lines

  ! boolean parameters
  read(IIN,*)
  read(IIN,*)
  read(IIN,*)
  read(IIN,3) junk,DEBUG
  read(IIN,3) junk,MAKE_SEISMO_PLOTS
  read(IIN,3) junk,MAKE_WINDOW_FILES
  read(IIN,3) junk,BODY_WAVE_ONLY
  if (DEBUG) then
    write(*,*) 'DEBUG: ------------------------------------:'
    write(*,*) 'DEBUG: PARAMETERS read in from PAR_FILE are:'
    write(*,*) '       DEBUG',DEBUG
    write(*,*) '       MAKE_SEISMO_PLOTS',MAKE_SEISMO_PLOTS
    write(*,*) '       MAKE_WINDOW_FILES',MAKE_WINDOW_FILES
    write(*,*) '       BODY_WAVE_ONLY',BODY_WAVE_ONLY
  endif

  ! period min/max for filtering
  read(IIN,*)
  read(IIN,*)
  read(IIN,*)
  read(IIN,3) junk,RUN_BANDPASS
  read(IIN,2) junk,WIN_MIN_PERIOD
  read(IIN,2) junk,WIN_MAX_PERIOD
  FSTART = 1./WIN_MAX_PERIOD
  FEND   = 1./WIN_MIN_PERIOD
  if (DEBUG) then
    write(*,*) '       WIN_MIN_PERIOD',WIN_MIN_PERIOD
    write(*,*) '       WIN_MAX_PERIOD',WIN_MAX_PERIOD
  endif

  ! E(t) water level
  read(IIN,*)
  read(IIN,*)
  read(IIN,*)
  read(IIN,2) junk,STALTA_BASE
  if (DEBUG) write(*,*) '       STALTA_BASE',STALTA_BASE

  ! Tshift
  read(IIN,*)
  read(IIN,*)
  read(IIN,*)
  read(IIN,2) junk,TSHIFT_BASE
  read(IIN,2) junk,TSHIFT_REFERENCE
  if (DEBUG) write(*,*) '       TSHIFT_BASE',TSHIFT_BASE
  if (DEBUG) write(*,*) '       TSHIFT_REFERENCE',TSHIFT_REFERENCE

  ! limit on dlnA for window acceptance
  read(IIN,*)
  read(IIN,*)
  read(IIN,*)
  read(IIN,2) junk,DLNA_BASE
  read(IIN,2) junk,DLNA_REFERENCE
  if (DEBUG) write(*,*) '       DLNA_BASE',DLNA_BASE
  if (DEBUG) write(*,*) '       DLNA_REFERENCE',DLNA_REFERENCE

  ! limit on CC
  read(IIN,*)
  read(IIN,*)
  read(IIN,*)
  read(IIN,2) junk,CC_BASE
  if (DEBUG) write(*,*) '       CC_BASE',CC_BASE

  ! boolean switch for check_data_quality
  read(IIN,*)
  read(IIN,*)
  read(IIN,*)
  read(IIN,3) junk,DATA_QUALITY

  ! if DATA_QUALITY = .true. and if two different measurements of
  ! signal-to-noise ratios exceeds these two base levels,
  ! then the data time series (and syn) is kept
  read(IIN,*)
  read(IIN,*)
  read(IIN,*)
  read(IIN,*)
  read(IIN,2) junk,SNR_INTEGRATE_BASE
  read(IIN,2) junk,SNR_MAX_BASE

  ! limit on signal to noise ratio in a particular window.
  read(IIN,*)
  read(IIN,*)
  read(IIN,*)
  read(IIN,2) junk,WINDOW_S2N_BASE

  ! Fine tuning constants
  read(IIN,*)
  read(IIN,*)
  read(IIN,*)
  read(IIN,2) junk,C_0
  read(IIN,2) junk,C_1
  read(IIN,2) junk,C_2
  read(IIN,2) junk,C_3a
  read(IIN,2) junk,C_3b
  read(IIN,2) junk,C_4a
  read(IIN,2) junk,C_4b
  read(IIN,*)
  read(IIN,2) junk,WEIGHT_SPACE_COVERAGE
  read(IIN,2) junk,WEIGHT_AVERAGE_CC
  read(IIN,2) junk,WEIGHT_N_WINDOWS

  if (DEBUG) then
    write(*,*) 'DEBUG: ------------------------------------:'
  endif

  !--------------------------------------------------------
  ! close parameter file
  close(IIN)


  !--------------------------------------------------------
  ! line formats
2 format(a,f20.8)
3 format(a,l20)

  ! unused formats
  ! 1 format(a,i20)
  ! 4 format(a,a)


  end subroutine
!----------------------------------------------------------------------

  subroutine read_sac_files(file_s, file_o,ier)
  use seismo_variables
  implicit none

  character (len=240) :: file_s, file_o
  integer :: ier
  ! local parameters
  real :: b1, dt1, b2, dt2, synt_sngl(NDIM), obs_sngl(NDIM)
  integer :: npts1, npts2, nerr

  ! initialize
  ier = 0

! read synthetic
  call rsac1(file_s,synt_sngl,npts1,b1,dt1,NDIM,nerr)
  if (nerr.ne.0) then
    write(*,*)
    write(*,*)' !!! Error reading file ', file_s
    !write(*,*)'     Program stop !!!'
    write(*,*)
    write(*,*)'  npts1:',npts1,'b1:',b1,'dt1:',dt1,'NDIM:',NDIM
    write(*,*)'  error: ',nerr
    !stop
    ier = 1
    return
  endif
  synt(:) = 0.
  synt(1:npts1) = dble(synt_sngl(1:npts1))

! read observed
  call rsac1(file_o,obs_sngl, npts2,b2,dt2,NDIM,nerr)
  if (nerr.ne.0) then
    write(*,*)
    write(*,*)' !!! Error reading file ', file_o
    !write(*,*)'     Program stop !!!'
    write(*,*)
    write(*,*)'  npts2:',npts2,'b2:',b2,'dt2:',dt2,'NDIM:',NDIM
    write(*,*)'  error: ',nerr
    write(*,*)
    !stop
    ier = 1
    return
  endif
  obs(:) = 0.
  obs(1:npts1) = dble(obs_sngl(1:npts1))

! check sample rates are equal
  if( abs(dt1-dt2).gt.1e-05) then
    write(*,*)
    !write(*,*) 'dt1-syn, dt2-dat,  abs(dt1-dt2):'
    write(*,*) 'error: dt1-syn, dt2-dat,  abs(dt1-dt2):'
    write(*,*) dt1, dt2,  abs(dt1-dt2)
    write(*,*)' !!! sampling rates differ, program stop !!!'
    write(*,*)
    !stop
    ier = 1
    return
  endif

! set global variable dt
  dt = dble(dt1)
  !write(*,*)'sampling rate dt=',dt

! check start times are equal
  if(abs(b1-b2).gt.2.0*dt) then
    write(*,*)
    write(*,*) 'b1, b2, dt, abs(b1-b2), 2*dt:'
    write(*,*) b1, b2, dt, abs(b1-b2), 2.0*dt
    write(*,*)' !!! start times differ, program stop !!!'
    write(*,*)
    !stop
    ier = 1
    return
  endif

! set global variable b
  b=dble(b1)
  !write(*,*) 'seismogram start time b=',b

! set global npts to the npts of the shortest seismogram
  npts = npts1
  if (npts2 .lt. npts1) npts = npts2
!  write(*,*)'number of points =',npts

! DEBUG
  if (DEBUG) write(*,'(a,f10.1,f10.4,i10)') ' DEBUG : b, dt, npts ', b, dt, npts


! read event and station header parameters from observation file
  call getfhv('evla', evla, nerr)
  call getfhv('evlo', evlo, nerr)
  call getfhv('stla', stla, nerr)
  call getfhv('stlo', stlo, nerr)
  call getfhv('evdp', evdp, nerr)
  call getkhv('kstnm', kstnm, nerr)
  call getkhv('kcmpnm', kcmpnm, nerr)
  call getkhv('knetwk', knetwk, nerr)

  if (BODY_WAVE_ONLY) then
     call getfhv('t1', P_pick, nerr)
     call getfhv('t2', S_pick, nerr)
  endif

  ! LQY fixed -- CHT: why does this output as:  BZN    ^@BHR    ^@AZ     ^@
  if (DEBUG) write(*,*) 'DEBUG : sta,net,comp = ', kstnm(1:7), kcmpnm(1:7), knetwk(1:7)

  ! calculate distances and azimuths
  call distaz(evla,evlo,stla,stlo,azimuth,backazimuth,dist_deg,dist_km)

  ! Frequency limits may be conditional on station or event information
  ! so call user function to modify them if required
  call modify_T0_T1_on_condition

  if(RUN_BANDPASS) then
     ! clean up the seismograms
     call detrend(obs,npts)
     call t_taper(obs, npts,HANNING,0.05)

     call detrend(synt,npts)
     call t_taper(synt, npts,HANNING,0.05)

     !call prepare_lp_seis
     call prepare_bp_seis

  else
     obs_lp(:)  = obs(:)
     synt_lp(:) = synt(:)

  endif

  end subroutine read_sac_files

!----------------------------------------------------------------------

  subroutine calculate_windows_before_quality
  use seismo_variables
  implicit none
  integer i

  do i = 1, num_win
    call calc_criteria(obs_lp,synt_lp,npts,i_start(i),i_end(i),dt,Tshift(i),CC(i),dlnA(i))
  enddo
  end subroutine

!----------------------------------------------------------------------

  subroutine calculate_windows_aux_quality
  use seismo_variables
  implicit none
  integer i

  do i = 1, num_win
    call calc_criteria(obs_lp,synt_lp,npts,i_start(i),i_end(i),dt,Tshift_aux(i),CC_aux(i),dlnA_aux(i))

  enddo
  end subroutine

!----------------------------------------------------------------------

! ------------------------------------------------------------------
! function costaper(ipoint, ndata, tas)
! ------------------------------------------------------------------
  real function costaper(npoints, ipoint, iexp)
  integer, intent (in) :: npoints, ipoint, iexp
  real pi

  pi = asin(1.0d0)*2
  costaper = sngl(1.0d0 - cos(pi*(ipoint-1)/npoints)**iexp)
  return
  end

!----------------------------------------------------------------------

  subroutine t_taper(x, n, t_type, width)

  use user_parameters
  implicit none

  integer, intent(in) :: n, t_type
  double precision, dimension(*), intent(inout) :: x
  real, intent(in) :: width

  integer i, n_width
  double precision omega, f0, f1

  ! set the number of points corresponding to the taper width
  n_width=int(floor(n*width))

  ! set the taper properties according to type
  select case (t_type)
    case (HANNING)
      omega = PI/dble(n_width)
      f0 = 0.5
      f1 = 0.5
    case (HAMMING)
      omega = PI/dble(n_width)
      f0 = 0.54
      f1 = 0.46
    case (COSINE)
      omega = PI/(2*dble(n_width))
      f0 = 1.0
      f1 = 1.0
  end select

  ! apply the taper symmetrically to both ends of the data
  do i = 1, n_width
    x(i) = x(i) * (f0-f1*cos(omega*(i-1)))
    x(n+1-i) = x(n+1-i) * (f0-f1*cos(omega*(i-1)))
  end do

  end subroutine t_taper

!----------------------------------------------------------------------
! calls sac function envelope (single precision version)
  subroutine sac_envelope(n, seis, env)

  ! AM : added this to fix memory bug
  use user_parameters
  implicit none

  integer, intent(in) :: n
  double precision, intent(in), dimension(*) :: seis
  double precision, intent(out), dimension(*) :: env

  ! AM : removed allocation here to fix a memory bug
  !real, dimension(:), allocatable :: dummy_seis1, dummy_seis2
  real, dimension(NDIM) :: dummy_seis1, dummy_seis2
  external :: envelope

  ! AM : removed allocation here to fix memory bug
  !allocate(dummy_seis1(n))
  !allocate(dummy_seis2(n))

  if( NDIM < n ) stop 'error envelope size'

  dummy_seis1(:)=0 ; dummy_seis2(:)=0
  dummy_seis1(1:n)=sngl(seis(1:n))
  call envelope(n,dummy_seis1(1:n),dummy_seis2(1:n))
  env(1:n)=dble(dummy_seis2(1:n))

  ! AM : removed deallocation here to fix memory bug
  !deallocate(dummy_seis1)
  !deallocate(dummy_seis2)

  end subroutine sac_envelope

!----------------------------------------------------------------------

  subroutine bandpass(x, n, delta_t, x_filt)
  use seismo_variables
  implicit none

  integer, intent(in) :: n
  double precision, intent(in),  dimension(*) :: x
  double precision, intent(out), dimension(*) :: x_filt
  double precision, intent(in) :: delta_t

  real, dimension(:), allocatable :: x_sngl
!  real delta_t_sngl

  allocate(x_sngl(n))

  x_sngl(1:n) = sngl(x(1:n))
!  delta_t_sngl = sngl(delta_t)

  ! old version - uses old SacLib
  ! does band-pass filter
  !call xapiir(x_sngl,n,'BU',sngl(TRBDNDW),sngl(APARM),IORD,'BP',sngl(FSTART),sngl(FEND),delta_t_sngl,PASSES)


  ! new version, uses subroutines in libsac.a
  ! does band-pass filter
  call xapiir(x_sngl,n,'BU',TRBDNDW,APARM,IORD,'BP',FSTART,FEND,delta_t,PASSES)

  x_filt(1:n) = dble(x_sngl(1:n))

  deallocate(x_sngl)

  end subroutine bandpass

!----------------------------------------------------------------------
  subroutine prepare_bp_seis
  use seismo_variables
  implicit none

  ! make filtered seismograms
  obs_lp(:) = 0.
  synt_lp(:) = 0.
  call bandpass(obs,npts,dt,obs_lp)
  call bandpass(synt,npts,dt,synt_lp)

  end subroutine prepare_bp_seis

!----------------------------------------------------------------------

  subroutine prepare_lp_env
  use seismo_variables
  implicit none

  ! make filtered envelopes
  env_obs_lp(:) = 0.
  env_synt_lp(:) = 0.
  if (DEBUG) write(*,*) 'DEBUG first point of obs/synt envelope before ',env_obs_lp(1), env_synt_lp(1)
  call sac_envelope(npts,synt_lp,env_synt_lp)
  if (DEBUG) write(*,*) 'DEBUG first point of obs/synt envelope middle ',env_obs_lp(1), env_synt_lp(1)
  call sac_envelope(npts,obs_lp,env_obs_lp)
  if (DEBUG) write(*,*) 'DEBUG first point of obs/synt envelope after ',env_obs_lp(1), env_synt_lp(1)

  end subroutine prepare_lp_env

!----------------------------------------------------------------------
  subroutine prepare_STA_LTA
  use seismo_variables
  implicit none

  double precision :: TOL, Cs, Cl, sta, lta, noise
  integer :: i, n_extend
  double precision, dimension(:), allocatable :: extended_syn

  ! the sta_lta value peaks at the begining of the trace.
  ! there is no phase incoming at that time, but there is some small trend to positive values.
  ! tweak only sets the sta/lta value when a threshold value in the envelope is reached
  logical,parameter:: TWEAK = .false.
  integer :: thres_noise
  double precision :: lta_max, lta_org

! set the Cs Cl for STA/LTA ratio using the Bai & Kennett (2001) expressions

  Cs=10**(-dt/WIN_MIN_PERIOD)
  Cl=10**(-dt/(12*WIN_MIN_PERIOD))
  TOL=1e-9

! set pre-extension for synthetic data and allocate extended_syn
!  n_extend=5*12*WIN_MIN_PERIOD/dt
  n_extend=1000*int(WIN_MIN_PERIOD/dt)
  allocate(extended_syn(npts+n_extend))

! set noise level
!  n_sample=12*WIN_MIN_PERIOD/dt
!  noise=sum(env_synt_lp(1:n_sample))/n_sample
  noise=maxval(env_synt_lp)/10**5

  ! determines the onset of the envelope (above some threshold value)
  if( TWEAK ) then
    do i=1,npts
      if( env_synt_lp(i) >= noise*1000. ) then
        thres_noise = i
        exit
      endif
    enddo
  endif


! copy the original synthetic into the extended array, right justified
!  call random_number(extended_syn)
!  extended_syn=noise*extended_syn
  extended_syn=noise
  extended_syn(n_extend+1:n_extend+npts)=env_synt_lp(1:npts)+noise

  if (DEBUG) write(*,*) 'DEBUG : Cs, Cl = ', Cs, Cl
  if (DEBUG) write(*,*) 'Number of points used to pre-extend synthetics ', n_extend

! initialise sta/lta variables
  STA_LTA(:)=0.0
  sta=0.0 ; lta = 0.0

! warm up the sta and lta
  do i = 1, n_extend
    sta=(Cs*sta+extended_syn(i))
    lta=(Cl*lta+extended_syn(i))
  enddo

  ! determine long term average maximum value
  if( TWEAK ) then
    lta_org = lta
    lta_max = 0.0
    do i = 1, npts
      lta=(Cl*lta+extended_syn(n_extend+i))
      if( lta > lta_max ) lta_max = lta
    enddo
    if (DEBUG) write(*,*) 'DEBUG : lta_max = ', lta_max
    lta = lta_org
  endif

! calculate sta/lta for the envelope
  do i = 1, npts
    sta=(Cs*sta+extended_syn(n_extend+i))
    lta=(Cl*lta+extended_syn(n_extend+i))
    if (lta.gt.TOL) then
      if( TWEAK ) then
        ! additional envelope criteria
        if( lta .gt. TOL*lta_max ) then
          if( i >= thres_noise) then
            STA_LTA(i)=sta/lta
          endif
        endif
      else
        STA_LTA(i)=sta/lta
      endif
    endif

  enddo

  deallocate(extended_syn)


  end subroutine

!----------------------------------------------------------------------

  subroutine detrend(x,n)
  implicit none

  integer, intent(in) :: n
  double precision, dimension(*) :: x

  double precision :: ds1,ds2,dan,davei,davex,dslope,dai
  integer :: i, an

  an = n
  dan=n
  ds1=0
  ds2=0

  do i=1,n
    ds1 = ds1+ x(i)
    ds2 = ds2 + ds1
  enddo
  davei = 0.5 * (dan+1.0)
  davex = ds1/dan
  dslope = -12.0*(ds2-davei*ds1)/(dan*(dan*dan-1.0))
  do i=1,n
    dai = i-1
    x(i) = x(i)- davex - dslope*(dai-davei)
  enddo

  end subroutine detrend

!------------------------------------------------------------------------
  subroutine int2string(i,istring)
  implicit none

  integer, intent(in) :: i
  character*8, intent(out) :: istring

  if(i.lt.0) stop 'Complete subroutine for negative integers'

  if(abs(i).lt.10) then
    write(istring,'(i1)') i
  elseif(abs(i).lt.100) then
    write(istring,'(i2)') i
  elseif(abs(i).lt.1000) then
    write(istring,'(i3)') i
  elseif(abs(i).lt.10000) then
    write(istring,'(i4)') i
  else
    stop 'Complete subroutine for integers over 9999 '
  endif

  end subroutine

!------------------------------------------------------------------------
  subroutine phases_in_window(t_start, t_end, n, indexes)
  use seismo_variables
  implicit none

  double precision, parameter :: TOL = 1

  double precision, intent(in) :: t_start, t_end
  integer, intent(out) :: n
  integer, dimension(*), intent(out) :: indexes

  integer :: i

  n=0
  do i = 1, num_phases
   if (ph_times(i) > t_start-TOL .and. ph_times(i) < t_end+TOL) then
     n=n+1
     indexes(n)=i
   endif
  enddo

  end subroutine phases_in_window

!------------------------------------------------------------------------
  subroutine xcorr_calc(d,s,npts,i1,i2,ishift,cc_max)
  implicit none

  ! inputs:
  ! s(npts) = synthetic
  ! d(npts) = data (or observed)
  ! i1, i2 = start and stop indexes of window within s and d

  double precision, dimension(*), intent(in) :: s,d
  integer, intent(in) :: npts, i1, i2

  ! outputs:
  ! ishift = index lag (d-s) for max cross correlation
  ! cc_max = maximum of cross correlation (normalised by sqrt(synthetic*data))
  integer, intent(out) :: ishift
  double precision, intent(out) :: cc_max

  ! local variables
  integer :: nlen
  integer :: i_left, i_right, i, j, id_left, id_right
  double precision :: cc, norm, norm_s

  ! initialise shift and cross correlation to zero
  ishift = 0
  cc_max = 0.0d0

  if (i1.lt.1 .or. i1.gt.i2 .or. i2.gt.npts) then
    write(*,*) 'Error with window limits: i1, i2, npts ', i1, i2, npts
    return
  endif

  ! length of window (number of points, including ends)
  nlen = i2 - i1 + 1

  ! power of synthetic signal in window
  norm_s = sqrt(sum(s(i1:i2)*s(i1:i2)))

  ! left and right limits of index (time) shift search
  ! NOTE: This looks OUTSIDE the time window of interest to compute TSHIFT and CC.
  !       How far to look outside, in theory, should be another parameter.
  !       However, it does not matter as much if the data and synthetics are
  !          zeroed outside the windows, as currently done in calc_criteria.
  i_left = -1*int(nlen/2.0)
  i_right = int(nlen/2.0)

  ! i is the index to shift to be applied to DATA (d)
  do i = i_left, i_right

    ! normalization factor varies as you take different windows of d
    id_left = max(1,i1+i)      ! left index for data window
    id_right = min(npts,i2+i)  ! right index for data window
    norm = norm_s * sqrt(sum(d(id_left:id_right)*(d(id_left:id_right))))

    ! cc as a function of i
    cc = 0.
    do j = i1, i2   ! loop over full window length
      if((j+i).ge.1 .and. (j+i).le.npts) cc = cc + s(j)*d(j+i)  ! d is shifted by i
    enddo
    cc = cc/norm

    ! keeping cc-max only
    if (cc .gt. cc_max) then
      cc_max = cc
      ishift = i
    endif
  enddo

  ! EXAMPLE: consider the following indexing:
  ! Two records are from 1 to 100, window is i1=20 to i2=41.
  !    --> nlen = 22, i_left = -11, i_right = 11
  !    i   i1+i   i2+i  id_left  id_right
  !  -11     9     30      9        30
  !   -5    15     36     15        36
  !    0    20     41     20        41
  !    5    25     46     25        46
  !   10    31     52     31        52

end subroutine xcorr_calc

!------------------------------------------------------------------------

subroutine calc_criteria(d,s,npts,i1,i2,dt,tshift,cc_max,dlnA)
  use user_parameters  ! for NDIM in d_win and s_win
  implicit none

  ! CHT: modified version of calc_criteria_original
  ! The point here is that we only look EXACTLY within the time window of interest,
  ! rather than looking outside to line up the signal that is selected in the synthetics.

  double precision, dimension(*), intent(in) :: d, s
  integer, intent(in) :: npts,i1,i2
  double precision, intent(in) :: dt
  double precision, intent(out) :: tshift,cc_max,dlnA

  double precision, dimension(NDIM) :: d_win, s_win
  integer :: ishift

  ! do cross-correlation
  ! CHT: zero the data and synthetics outside the window (see comments in xcorr_calc)
  d_win(:) = 0. ; d_win(i1:i2) = d(i1:i2)
  s_win(:) = 0. ; s_win(i1:i2) = s(i1:i2)
  call xcorr_calc(d_win,s_win,npts,i1,i2,ishift,cc_max)

  ! cross-correlation time-shift
  tshift = ishift*dt

  ! calculate dlnA : definition of Dahlen and Baig (2002), Eq. 3,17,18 : dlnA = Aobs/Asyn - 1
  ! NOTE 1: these records are unshifted
  ! NOTE 2: this measurement will reflect any noise in the data, too
  !dlnA = sqrt( ( sum( d_win(i1:i2)*d_win(i1:i2) )) / (sum( s_win(i1:i2)*s_win(i1:i2) )) ) - 1.0

  ! CHT revision 15-oct-2008
  ! The previously used expression for dlnA is a first-order perturbation of ln(A1/A2) = (A1-A2)/A2
  ! The new expression is better suited to getting values between -1 and 1,
  ! with dlnA = 0 indicating perfect fit, as before.
  dlnA = 0.5 * log( sum(d_win(i1:i2)*d_win(i1:i2)) / sum(s_win(i1:i2)*s_win(i1:i2)) )

end subroutine calc_criteria

!------------------------------------------------------------------------

subroutine calc_criteria_original(d,s,npts,i1,i2,dt,tshift,cc_max,dlnA)
  implicit none

  double precision, parameter :: TOL = 1e-7

  double precision, dimension(*), intent(in) :: d, s
  integer, intent(in) :: npts,i1,i2
  double precision, intent(in) :: dt
  double precision, intent(out) ::  tshift,cc_max,dlnA

  double precision, dimension(:), allocatable :: s_cor,d_loc
  integer :: ishift, npts_win, i

  npts_win = i2-i1+1

  ! allocate memory for s_cor (the corrected synthetic)
  allocate(s_cor(npts_win))
  allocate(d_loc(npts_win))

  d_loc(1:npts_win) = d(i1:i2)

  ! do cross-correlation
  call xcorr_calc(d,s,npts,i1,i2,ishift,cc_max)
  tshift = ishift*dt

  ! apply time shift to synthetic seismogram
  do i = 1, npts_win
     s_cor(i) = 0.
     if( (i1-1+i-ishift) .gt. 1 .and. (i1-1+i-ishift) .lt. npts ) s_cor(i) = s(i1-1+i-ishift)
  enddo

  ! calculate dlnA
  dlnA = sqrt( ( sum( d(i1:i2) * d(i1:i2) )) / (sum( s_cor(1:npts_win) * s_cor(1:npts_win) )) ) - 1.0

end subroutine calc_criteria_original

!------------------------------------------------------------------------

  subroutine prepare_for_window_sel_stalta
  use seismo_variables

  implicit none

! get the travel times
  if (DEBUG) print *, "Calling ttimes"
  call ttimes(dist_deg,evdp,num_phases,ph_names,ph_times)

  ! find the index of the first arrival (global, filtered)
  n_p = 1+int((ph_times(1)-WIN_MIN_PERIOD/2.0-b)/dt)
  if (n_p .le. 0) n_p = 1

  call prepare_lp_env

  call prepare_sta_lta

! set up the selection criteria arrays (subroutine is in user_functions.f90
! file)
  call set_up_criteria_arrays

  end subroutine

!------------------------------------------------------------------------
