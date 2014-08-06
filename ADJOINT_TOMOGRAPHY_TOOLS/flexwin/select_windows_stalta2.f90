!
! $Id:$
!
!
!----------------------------------------------------------------------

  subroutine select_windows_stalta2()
  use seismo_variables
  implicit none

  integer, dimension(NWINDOWS) :: maxima_lp, minima_lp
  integer, dimension(NWINDOWS*NWINDOWS) :: iM, iL, iR
  double precision, dimension(NWINDOWS) :: CC_local, Tshift_local, dlnA_local
  integer :: nmin, nmax, nwin
  integer :: k
  double precision :: w_level_local

  logical :: data_ok
  double precision :: noise_int, signal_int, snr_int, signal_amp, noise_amp, snr_amp

! initialise global arrays
  num_win = 0
  i_start(:) = 0 ; i_end(:) = 0
  win_start(:) = 0 ; win_end(:) = 0
  Tshift(:) = 0 ; CC(:) = 0 ; dlnA(:) = 0
  Tshift_after(:) = 0 ; CC_after(:) = 0 ; dlnA_after(:) = 0

! initialise local arrays
  maxima_lp(:) = 0; minima_lp(:) = 0
  iM(:) = 0; iL(:) = 0; iR(:) = 0
  CC_local(:) = 0; dlnA_local(:) = 0; Tshift_local(:) = 0

! make all preparatory steps (create envelopes and sta_lta time series)
  if(DEBUG) print *, "Preparing windows for stalta"
  call prepare_for_window_sel_stalta

  ! check data quality (signal-to-noise ratios)
  if(DATA_QUALITY) then
    if(DEBUG) print *, "Checking data quality"
    ! returns data_ok = .true. if two signal-to-noise criteria are met
    call check_data_quality(data_ok,signal_int,noise_int,snr_int,signal_amp,noise_amp,snr_amp)
    ! if data quality is not ok, then return without selecting windows
    if (.not. data_ok) then
      num_win = 0
      return
    endif
  endif

!!$  if (1==1) then
!!$     ! keep the manually selected window that has met the signal-to-noise criteria;
!!$     ! this option is used for time-reversing coda surface waves
!!$     nwin = 1
!!$     iL(1) = is1
!!$     iM(1) = int((is1+is2)/2)
!!$     iR(1) = is2
!!$
!!$     ! reject window if it does not satisfy F2, tshift and dlnA criteria;
!!$     ! here we are only interested in using the dlnA criterion
!!$     call reject_on_fit_criteria(nwin,iM,iL,iR,CC_local,Tshift_local,dlnA_local)
!!$     if (nwin.eq.0) then
!!$        write(*,*) 'NO WINDOWS SELECTED FOR THIS TRACE'
!!$        num_win = 0
!!$        return
!!$     endif
!!$     if(DEBUG) call display_windows_full(nwin,iM,iL,iR,CC_local,Tshift_local,dlnA_local)
!!$
!!$  else

     ! find list of maxima and minima, ignoring water level
     w_level_local = 0.0
     if(DEBUG) print *, "Finding max min "
     call find_maxima(STA_LTA,npts,maxima_lp,nmax,w_level_local)
     call find_minima(STA_LTA,npts,minima_lp,nmin,w_level_local)

     ! find all possible windows within seismogram, with central maxima above the water level
     if(DEBUG) print *, "Making all windows "
     call setup_M_L_R(nmax,maxima_lp,nmin,minima_lp,nwin,iM,iL,iR)
     if (nwin.eq.0) then
        write(*,*) 'NO WINDOWS SELECTED FOR THIS TRACE'
        num_win = 0
        return
     endif

     ! remove windows with internal minima below the water level
     if(DEBUG) print *, "Rejecting on water level "
     call reject_on_water_level(nwin,iM,iL,iR,nmin,minima_lp,C_0)
     if (nwin.eq.0) then
        write(*,*) 'NO WINDOWS SELECTED FOR THIS TRACE'
        num_win = 0
        return
     endif
     ! sort windows
     call sort_on_start_time(nwin,iM,iL,iR)
     !if(DEBUG) call display_windows(nwin,iM,iL,iR)

     ! remove small windows
     if(DEBUG) print *, "Rejecting on size "
     call reject_on_window_width(nwin,iM,iL,iR,C_1)
     if (nwin.eq.0) then
        write(*,*) 'NO WINDOWS SELECTED FOR THIS TRACE'
        num_win = 0
        return
     endif
     !if(DEBUG) call display_windows(nwin,iM,iL,iR)

     ! reject on prominence of central maximum
     if(DEBUG) print *, "Rejecting on prominence "
     call reject_on_prominence(nwin,iM,iL,iR,nmin,minima_lp,C_2)
     if (nwin.eq.0) then
        write(*,*) 'NO WINDOWS SELECTED FOR THIS TRACE'
        num_win = 0
        return
     endif
     !if(DEBUG) call display_windows(nwin,iM,iL,iR)

     ! reject windows containing more than one distinct phase
     call reject_on_phase_separation(nwin,iM,iL,iR,nmin,minima_lp,nmax,maxima_lp,C_3a, C_3b)
     if (nwin.eq.0) then
        write(*,*) 'NO WINDOWS SELECTED FOR THIS TRACE'
        num_win = 0
        return
     endif
     !if(DEBUG) call display_windows(nwin,iM,iL,iR)

     ! curtail window ends by time decay
     call curtail_window_length(nwin,iL,iR,nmax,maxima_lp,C_4a,C_4b)
     !if(DEBUG) call display_windows(nwin,iM,iL,iR)

     ! REPEAT: remove small windows, since curtailing may have shrunk them
     ! NOTE: perhaps this only needs to be done once (here)
     if(DEBUG) print *, "Rejecting on size (REPEAT)"
     call reject_on_window_width(nwin,iM,iL,iR,C_1)
     if (nwin.eq.0) then
        write(*,*) 'NO WINDOWS SELECTED FOR THIS TRACE'
        num_win = 0
        return
     endif
     if(DEBUG) call display_windows(nwin,iM,iL,iR)

     ! check window quality
     if(DATA_QUALITY) call check_window_s2n(nwin,iM,iL,iR)
     if (nwin.eq.0) then
        write(*,*) 'NO WINDOWS SELECTED FOR THIS TRACE'
        num_win = 0
        return
     endif
     if(DEBUG) call display_windows(nwin,iM,iL,iR)

     ! reject windows that do not satisfy F2, tshift and dlnA criteria
     call reject_on_fit_criteria(nwin,iM,iL,iR,CC_local,Tshift_local, dlnA_local)
     if (nwin.eq.0) then
        write(*,*) 'NO WINDOWS SELECTED FOR THIS TRACE'
        num_win = 0
        return
     endif
     if(DEBUG) call display_windows_full(nwin,iM,iL,iR,CC_local,Tshift_local,dlnA_local)

     ! now that all criteria are satisfied, reject any duplicate windows
     call reject_on_duplicate(nwin,iM,iL,iR,CC_local,Tshift_local,dlnA_local)
     if(DEBUG) call display_windows_full(nwin,iM,iL,iR,CC_local,Tshift_local,dlnA_local)

     ! resolve the overlaps
     call resolve_overlaps(nwin,iM,iL,iR,CC_local,Tshift_local,dlnA_local)
     if(DEBUG) call display_windows_full(nwin,iM,iL,iR,CC_local,Tshift_local,dlnA_local)

!!$  endif

  ! set global number of windows
  num_win = nwin
  do k = 1,num_win
    win_start(k) = b+(iL(k)-1)*dt
    win_end(k)   = b+(iR(k)-1)*dt
  enddo

  write(*,*) "Selected windows, start and end time, CC, Tshift, dlnA  "
  do k = 1,num_win
    write(*,'(i4,1x,5f10.3)') k, win_start(k), win_end(k), CC_local(k), Tshift_local(k), dlnA_local(k)

    ! sanity check for the end of the record
    if(win_start(k) .lt. b) win_start(k) = b
    if(win_end(k) .gt. b+(npts-1)*dt) win_end(k) = b+(npts-1)*dt

    ! calculate indexes for start and and of windows
    i_start(k) = 1+int((win_start(k)-b)/dt)
    i_end(k)   = 1+int((win_end(k)  -b)/dt)

  enddo

  ! setup quality criteria
  CC(1:num_win) = CC_local(1:num_win)
  dlnA(1:num_win) = dlnA_local(1:num_win)
  Tshift(1:num_win) = Tshift_local(1:num_win)

  end subroutine select_windows_stalta2


  subroutine check_window_s2n(nwin,iM,iL,iR)
    use seismo_variables
    implicit none

    integer, intent(inout) :: nwin
    integer, dimension(*), intent(inout) :: iM, iL, iR

    integer :: iwin, i
    integer :: nwin_new
    integer, dimension(NWINDOWS) :: iM_new, iL_new, iR_new
    logical :: accept
    double precision :: time, noise_amp, signal_amp, amp_ratio
    double precision :: noise_pow, signal_pow, pow_ratio

    ! Determine the max amplitude of the noise, as defined by the
    ! global variable noise_end.
    do i = 1, NDIM
       time = b+(i-1)*dt
       if (time > noise_end) exit
    enddo
    ! i is the index corresponding to noise_end
    noise_amp = maxval( abs(obs_lp(1:i)) )
    noise_pow = sum( (obs_lp(1:i))**2 )/(i-1)

    ! Determine the max amplitude of the windows that have their
    ! left bound time greater than the noise_end time.
    ! Window will be rejected based on the amp_ratio.

    nwin_new = 0
    do iwin = 1, nwin
       accept = .true.
       signal_amp = maxval( abs(obs_lp(iL(iwin):iR(iwin))) )
       signal_pow = sum( (obs_lp(iL(iwin):iR(iwin)))**2 )/(iR(iwin)-iL(iwin))
       amp_ratio = signal_amp / noise_amp
       pow_ratio = signal_pow / noise_pow

       if(DEBUG) then
           write (*,*) 'DEBUG : amp_ratio,signal_amp,noise_amp : ',amp_ratio,signal_amp,noise_amp
           write (*,*) 'DEBUG : pow_ratio,signal_pow,noise_pow : ',pow_ratio,signal_pow,noise_pow
           write (*,*) 'DEBUG : iwin, amp_ratio : ', iwin, amp_ratio, S2N_LIMIT(iM(iwin))
       endif

       if (amp_ratio < S2N_LIMIT(iM(iwin))) then
       !if (pow_ratio < S2N_LIMIT(iM(iwin))) then
          accept = .false.
       endif

       if (accept) then
          nwin_new = nwin_new + 1
          iM_new(nwin_new) = iM(iwin)
          iR_new(nwin_new) = iR(iwin)
          iL_new(nwin_new) = iL(iwin)
       endif
    enddo

    ! update the iM, iR, iL arrays to contain only accepted proto-windows
    nwin=nwin_new
    iM(1:nwin) = iM_new(1:nwin_new)
    iL(1:nwin) = iL_new(1:nwin_new)
    iR(1:nwin) = iR_new(1:nwin_new)

    if (DEBUG) write(*,*) 'DEBUG : window noise quality rejection retained ', nwin, ' acceptable windows'

  end subroutine check_window_s2n


  subroutine check_data_quality(data_ok,signal_int,noise_int,snr_int,signal_amp,noise_amp,snr_amp)
    use seismo_variables
    implicit none

    ! data_ok is a logical value to indicate whether or not
    ! to continue on with the window selection algorithm.
    ! If noise to signal ratio is too high, data_ok is false,
    ! and that particular station will be rejected.

    logical, intent(out) :: data_ok
    double precision, intent(out) :: signal_int,noise_int,snr_int,signal_amp,noise_amp,snr_amp

    double precision :: time_obs_noise, time_obs_signal
    integer :: i, j, k, m

    ! initialize values
    data_ok = .true.
    noise_int = 0.0
    signal_int = 0.0
    snr_int = 0.0
    signal_amp = 0.0
    noise_amp = 0.0
    snr_amp = 0.0

    if(DEBUG) then
       write(*,*) 'DEBUG : noise_start = ', sngl(noise_start)
       write(*,*) 'DEBUG : noise_end = ', sngl(noise_end)
       write(*,*) 'DEBUG : signal_start = ', sngl(signal_start)
       write(*,*) 'DEBUG : signal_end = ', sngl(signal_end)
    endif

    ! compute noise
    ! see user_functions.f90 (default: noise_start = b)
    do m = 1, npts
       time_obs_noise = b+(m-1)*dt
       if (time_obs_noise > noise_start) exit
    enddo
    in1 = m-1       ! index corresponding to noise_start
    do i = 1, npts
       time_obs_noise = b+(i-1)*dt
       if (time_obs_noise > noise_end ) exit
    enddo
    in2 = i-1       ! index corresponding to noise_end
    noise_int = sum( (obs_lp(in1:in2))**2 )/(in2-in1)
    noise_amp = maxval( abs(obs_lp(in1:in2)) )

    ! compute signal
    ! see user_functions.f90 (default: signal_start = noise_end)
    do j = 1, npts
       time_obs_signal = b+(j-1)*dt
       if (time_obs_signal > signal_start) exit
    enddo
    is1 = j     ! index corresponding to signal_start
    do k = 1, npts
       time_obs_signal = b+(k-1)*dt
       if (time_obs_signal > signal_end) exit
    enddo
    is2 = k-1     ! index corresponding to signal_end
    print *, 'is1 = ', is1, ', is2 = ', is2
    signal_int = sum( (obs_lp(is1:is2))**2 )/(is2-is1)
    signal_amp = maxval( abs(obs_lp(is1:is2)) )

    ! Calculate signal to noise ratio and amplitude ratio.
    snr_int = signal_int / noise_int
    snr_amp = signal_amp / noise_amp

    ! If snr_int or snr_amp is less than a certain value, then data_ok is .false.
    if (snr_int < SNR_INTEGRATE_BASE .or. snr_amp < SNR_MAX_BASE) then
       write(*,*) 'DATA QUALITY NOT OK'
       if(DEBUG) write(*,*) 'snr_int < SNR_INTEGRATE_BASE: ', &
          sngl(snr_int), sngl(SNR_INTEGRATE_BASE), snr_int < SNR_INTEGRATE_BASE
       if(DEBUG) write(*,*) 'snr_amp < SNR_MAX_BASE: ', &
          sngl(snr_amp), sngl(SNR_MAX_BASE), snr_amp < SNR_MAX_BASE
       data_ok = .false.
    endif

    if (DEBUG) then
       print *, 'DEBUG : in1,in2,is1,is2,npts :', in1,in2,is1,is2,npts
       print *, 'DEBUG : data_ok :', data_ok
       print *, 'DEBUG : signal_int : ', sngl(signal_int)
       print *, 'DEBUG : noise_int :',  sngl(noise_int)
       print *, 'DEBUG : snr_int :', sngl(snr_int)
       print *, 'DEBUG : signal_amp : ', sngl(signal_amp)
       print *, 'DEBUG : noise_amp :',  sngl(noise_amp)
       print *, 'DEBUG : snr_amp :', sngl(snr_amp)
       print *, 'DEBUG : SNR_INTEGRATE_BASE :', sngl(SNR_INTEGRATE_BASE)
       print *, 'DEBUG : SNR_MAX_BASE :', sngl(SNR_MAX_BASE)
    endif

  end subroutine check_data_quality


  subroutine setup_M_L_R(nmax,maxima_lp,nmin,minima_lp,nwin,iM,iL,iR)
  use seismo_variables
  implicit none
  integer, intent(in) :: nmax, nmin
  integer, dimension(*), intent(in) :: maxima_lp, minima_lp
  integer, intent(out) :: nwin
  integer, dimension(NWINDOWS*NWINDOWS), intent(out) :: iM,iL,iR

  integer :: i, R, L
  integer :: ii, i_left, i_right

  nwin = 0
  do i = 1, nmax
    ! only continue if there are available minima on either side
    if (maxima_lp(i) > minima_lp(1) .and. maxima_lp(i) < minima_lp(nmin) ) then
      ! only continue if this maximum is above the water level
      if (STA_LTA(maxima_lp(i)) .gt. STALTA_W_LEVEL(maxima_lp(i))) then

        ! find the first minimum right of this maximum
        R = 0
        do ii = 1,nmin
          if (minima_lp(ii) > maxima_lp(i)) then
            R = ii
            exit
          endif
        enddo
        ! find the first minimum left of this maximum
        L = 0
        do ii = nmin,1,-1
          if (minima_lp(ii) < maxima_lp(i)) then
            L = ii
            exit
          endif
        enddo

        ! iterate over minima to get all possible windows around this maximum
        do i_left = L,1,-1
          do i_right = R,nmin,1
            ! checks number of windows so far
            if (nwin .ge. NWINDOWS*NWINDOWS) then
               print *, 'setup_M_L_R: ', nwin, NWINDOWS, NWINDOWS*NWINDOWS
               print *, 'setup_M_L_R: limit number (nwin = NWINDOWS*NWINDOWS)'
               !print *, 'setup_M_L_R: ignoring trace (nwin > NWINDOWS*NWINDOWS)'
               !nwin = 0
               !exit do-loop
               exit
               !return
               !stop 'setup_M_L_R: Increase NWINDOWS'
            endif

            ! set the index of the maximum, left and right boundaries
            nwin = nwin+1
            iM(nwin) = maxima_lp(i)
            iL(nwin) = minima_lp(i_left)
            iR(nwin) = minima_lp(i_right)
          enddo
        enddo

      endif
    endif
    if (nwin .ge. NWINDOWS*NWINDOWS) exit
  end do
  if (DEBUG) write(*,*) 'DEBUG : window generation found ', nwin, ' possible windows'

  end subroutine setup_M_L_R

  subroutine reject_on_water_level(nwin,iM,iL,iR,nmin,minima_lp,c_value)
  use seismo_variables
  implicit none

  integer, intent(inout) :: nwin
  integer, dimension(*), intent(inout) :: iM, iL, iR
  integer, intent(in) :: nmin
  integer, dimension(*), intent(in) :: minima_lp
  double precision, intent(in) :: c_value

  integer :: iwin, imin, min_index
  integer :: nwin_new
  integer, dimension(NWINDOWS) :: iM_new, iL_new, iR_new
  logical :: accept

  nwin_new = 0
  iM_new(:) = 0
  iL_new(:) = 0
  iR_new(:) = 0

  ! for each proto-window, check that no internal minima fall below the
  ! local water level of the window center
  do iwin = 1, nwin
    accept = .true.
    do imin = 1, nmin
       ! if the current minimum is within our window AND it falls beneath
       ! the water level, then reject this proto-window
       min_index = minima_lp(imin)
       if ( min_index .gt. iL(iwin) &
          .and. min_index .lt. iR(iwin) &
          .and. STA_LTA(min_index) .le. c_value*STALTA_W_LEVEL(iM(iwin)) ) then
         accept = .false.
         exit
       end if
    enddo
    ! if the proto-window is acceptable, then accept it
    if (accept) then
      ! checks windows
      if (nwin_new .ge. NWINDOWS) then
         print *, 'reject_on_water_level: limit windows number (nwin = NWINDOWS)',nwin_new
         !print *, 'reject_on_water_level: ignoring trace (nwin > NWINDOWS*NWINDOWS)'
         !nwin = 0
         ! exit do-loop
         exit
         !return
         !stop 'reject_on_water_level : Increase NWINDOWS'
      endif
      nwin_new = nwin_new + 1
      iM_new(nwin_new) = iM(iwin)
      iR_new(nwin_new) = iR(iwin)
      iL_new(nwin_new) = iL(iwin)
    endif
  enddo

  ! update the iM, iR, iL arrays to contain only accepted proto-windows
  nwin = nwin_new
  iM(1:nwin) = iM_new(1:nwin_new)
  iL(1:nwin) = iL_new(1:nwin_new)
  iR(1:nwin) = iR_new(1:nwin_new)

  if (DEBUG) write(*,*) 'DEBUG : water level rejection retained ', nwin, ' acceptable windows'

  end subroutine

  subroutine reject_on_phase_separation(nwin,iM,iL,iR,nmin,minima_lp,nmax,maxima_lp,c_value,c_value2)
  use seismo_variables
  implicit none

  integer, intent(inout) :: nwin
  integer, dimension(*), intent(inout) :: iM, iL, iR
  integer, intent(in) :: nmin, nmax
  integer, dimension(*), intent(in) :: minima_lp, maxima_lp
  double precision, intent(in) :: c_value, c_value2

  integer :: iwin, imin, min_index, imax, max_index
  integer :: nwin_new
  integer, dimension(NWINDOWS) :: iM_new, iL_new, iR_new
  double precision :: stalta_min, d_stalta, d_stalta_center, d_time, f_time
  logical :: accept

  nwin_new = 0
  ! for each proto-window, check that all included maxima are part of the
  ! same phase group as the central maximum iM
  do iwin = 1, nwin
    accept = .true.

    ! find the lowest minimum within the window
    stalta_min = STA_LTA(iL(iwin))
    do imin = 1,nmin
       min_index = minima_lp(imin)
       if ( min_index .gt. iL(iwin) &
           .and. min_index .le. iR(iwin) &
           .and. STA_LTA(min_index) .lt. stalta_min ) then
         stalta_min = STA_LTA(min_index)
       endif
    end do
    if(stalta_min .lt. STALTA_W_LEVEL(iM(iwin)) ) stalta_min = STALTA_W_LEVEL(iM(iwin))

    ! find the height of the central maximum above this minimum value
    d_stalta_center = STA_LTA(iM(iwin)) - stalta_min

    ! run check on maxima on either side of the central maximum
    do imax = 1, nmax
       max_index = maxima_lp(imax)
       if ( max_index .gt. iL(iwin) &
           .and. max_index .lt. iR(iwin) &
           .and. max_index .ne. iM(iwin) ) then
          ! find height of current maximum above lowest minimum
          d_stalta = (STA_LTA(max_index) - stalta_min)
          ! find scaled time between current maximum and central maximum
          d_time = abs(iM(iwin)-max_index)*dt / WIN_MIN_PERIOD
          ! find value of time decay function
          if (d_time .ge. c_value2) then
            f_time = exp(-((d_time-c_value2)/c_value2)**2)
          else
            f_time = 1.0
          endif
          ! check condition
          if (d_stalta .gt. c_value*d_stalta_center*f_time) then
            accept = .false.
            exit
          endif
       end if
    end do

    ! if the proto-window is acceptable, then accept it
    if (accept) then
      nwin_new = nwin_new + 1
      iM_new(nwin_new) = iM(iwin)
      iR_new(nwin_new) = iR(iwin)
      iL_new(nwin_new) = iL(iwin)
    endif
  enddo

  ! update the iM, iR, iL arrays to contain only accepted proto-windows
  nwin = nwin_new
  iM(1:nwin) = iM_new(1:nwin_new)
  iL(1:nwin) = iL_new(1:nwin_new)
  iR(1:nwin) = iR_new(1:nwin_new)

  if (DEBUG) write(*,*) 'DEBUG : single phase group rejection retained ', nwin, ' acceptable windows'

  end subroutine

  subroutine reject_on_window_width(nwin,iM,iL,iR,c_value)
  use seismo_variables
  implicit none

  integer, intent(inout) :: nwin
  integer, dimension(*), intent(inout) :: iM, iL, iR
  double precision, intent(in) :: c_value

  integer :: iwin
  integer :: nwin_new
  integer, dimension(NWINDOWS) :: iM_new, iL_new, iR_new
  logical :: accept
  double precision :: window_length

  window_length = c_value * WIN_MIN_PERIOD / dt

  nwin_new = 0
  do iwin = 1,nwin
    accept = .true.
    ! check window length against minimum acceptable length
    if (iR(iwin) - iL(iwin) .lt. window_length  ) then
      accept = .false.
    end if

    if (accept) then
      ! checks windows
      if (nwin_new .ge. NWINDOWS) then
         print *, 'reject_on_window_width: limit (nwin = NWINDOWS)'
         exit
         !print *, 'reject_on_window_width: ignoring trace (nwin > NWINDOWS*NWINDOWS)'
         !nwin = 0
         !return
         !stop 'Increase NWINDOWS'
      endif
      nwin_new = nwin_new + 1
      iM_new(nwin_new) = iM(iwin)
      iR_new(nwin_new) = iR(iwin)
      iL_new(nwin_new) = iL(iwin)
    endif
  enddo

  ! update the iM, iR, iL arrays to contain only accepted proto-windows
  nwin = nwin_new
  iM(1:nwin) = iM_new(1:nwin_new)
  iL(1:nwin) = iL_new(1:nwin_new)
  iR(1:nwin) = iR_new(1:nwin_new)

  if (DEBUG) write(*,*) 'DEBUG : window length rejection retained ', nwin, ' acceptable windows'

  end subroutine

!=============================================================

!!$  subroutine reject_on_duplicate(nwin,iM,iL,iR)
!!$  use seismo_variables
!!$  implicit none
!!$  integer, intent(inout) :: nwin
!!$  integer, dimension(*), intent(inout) :: iM, iL, iR
!!$
!!$  integer :: iwin, nwin_new
!!$  integer, dimension(NWINDOWS) :: iM_new, iL_new, iR_new
!!$  logical :: accept
!!$  logical, dimension(NWINDOWS) :: duplicate
!!$
!!$  duplicate(1:nwin) = .false.
!!$
!!$  nwin_new = 0
!!$  do iwin = 1,nwin
!!$    accept = .true.
!!$    if (duplicate(iwin)) accept = .false.
!!$    do iwin2 = iwin,nwin
!!$      ! check if other windows are duplicates of this one, and if so set them as duplicates
!!$      if (iwin2.ne.iwin &
!!$        .and. iL(iwin).eq.iL(iwin2) &
!!$        .and. iR(iwin).eq.iR(iwin2) ) then
!!$        duplicate(iwin2) = .true.
!!$        continue
!!$      endif
!!$    enddo
!!$
!!$    if (accept) then
!!$      nwin_new = nwin_new + 1
!!$      iM_new(nwin_new) = iM(iwin)
!!$      iR_new(nwin_new) = iR(iwin)
!!$      iL_new(nwin_new) = iL(iwin)
!!$    endif
!!$  enddo
!!$
!!$  ! update the iM, iR, iL arrays to contain only accepted proto-windows
!!$  nwin = nwin_new
!!$  iM(1:nwin) = iM_new(1:nwin_new)
!!$  iL(1:nwin) = iL_new(1:nwin_new)
!!$  iR(1:nwin) = iR_new(1:nwin_new)
!!$
!!$  if (DEBUG) write(*,*) 'DEBUG : duplicate rejection retained ', nwin, ' acceptable windows'
!!$
!!$  end subroutine

!=============================================================

  subroutine reject_on_duplicate(nwin,iM,iL,iR,CC_local,Tshift_local,dlnA_local)
  use seismo_variables
  implicit none

  integer, intent(inout) :: nwin
  integer, dimension(*), intent(inout) :: iM, iL, iR
  double precision, dimension(*), intent(inout) :: CC_local,Tshift_local,dlnA_local

  integer :: iwin, nwin_new, ikeep, iwin2
  integer, dimension(NWINDOWS) :: iwin_keep
  logical :: accept
  logical, dimension(NWINDOWS) :: duplicate

  duplicate(1:nwin) = .false.
  iwin_keep(1:nwin) = 0

  nwin_new = 0
  do iwin = 1,nwin
    accept = .true.
    if (duplicate(iwin)) accept = .false.
    do iwin2 = iwin,nwin
      ! check if other windows are duplicates of this one, and if so set them as duplicates
      if (iwin2.ne.iwin &
        .and. iL(iwin).eq.iL(iwin2) &
        .and. iR(iwin).eq.iR(iwin2) ) then
        duplicate(iwin2) = .true.
        continue
      endif
    enddo

    if (accept) then
      nwin_new = nwin_new + 1
      iwin_keep(nwin_new) = iwin
      !write(*,*) 'Accept ', nwin_new, ' : ', iwin_keep(nwin_new)
    endif
  enddo

  ! update the iM, iR, iL arrays to contain only accepted proto-windows
  ! KEY: the ikeep index is ALWAYS greater than or equal to the iwin index
  nwin = nwin_new
  do iwin = 1,nwin
     ikeep = iwin_keep(iwin)
     iM(iwin) = iM(ikeep)
     iR(iwin) = iR(ikeep)
     iL(iwin) = iL(ikeep)

     CC_local(iwin) = CC_local(ikeep)
     Tshift_local(iwin) = Tshift_local(ikeep)
     dlnA_local(iwin) = dlnA_local(ikeep)
  enddo

  if (DEBUG) write(*,*) 'DEBUG : duplicate rejection retained ', nwin, ' acceptable windows'

  end subroutine

!=============================================================

  subroutine reject_on_prominence(nwin,iM,iL,iR,nmin,minima_lp,c_value)
  use seismo_variables
  implicit none

  integer, intent(inout) :: nwin
  integer, dimension(*), intent(inout) :: iM, iL, iR
  integer, intent(in) :: nmin
  integer, dimension(*), intent(in) :: minima_lp
  double precision, intent(in) :: c_value


  integer :: nwin_new
  integer, dimension(NWINDOWS) :: iM_new, iL_new, iR_new
  logical :: accept
  integer :: iwin, imin, i_left, i_right
  double precision :: delta_left, delta_right

  nwin_new=0
  do iwin = 1,nwin
    accept = .true.

    i_right = 0
    i_left = 0
    ! find index of minimum on the right of the central maximum
    do imin = 1, nmin
      if (minima_lp(imin) .gt. iM(iwin)) then
        i_right = minima_lp(imin)
        exit
      endif
    enddo
    ! find index of minimum on the left of the central maximum
    do imin = nmin, 1, -1
      if (minima_lp(imin) .lt. iM(iwin)) then
        i_left = minima_lp(imin)
        exit
      endif
    enddo

    delta_left = STA_LTA(iM(iwin)) - STA_LTA(i_left)
    delta_right = STA_LTA(iM(iwin)) - STA_LTA(i_right)

    ! check prominence condition
    ! if the maximum is not high enough above EITHER of the two adjecent
    ! minima, then reject the window
    if ( delta_left .lt. c_value * STALTA_W_LEVEL(iM(iwin)) &
        .or. delta_right .lt. c_value * STALTA_W_LEVEL(iM(iwin)) ) then
      accept = .false.
    end if

    if (accept) then
      nwin_new = nwin_new + 1
      iM_new(nwin_new) = iM(iwin)
      iR_new(nwin_new) = iR(iwin)
      iL_new(nwin_new) = iL(iwin)
    endif
  enddo

  ! update the iM, iR, iL arrays to contain only accepted proto-windows
  nwin = nwin_new
  iM(1:nwin) = iM_new(1:nwin_new)
  iL(1:nwin) = iL_new(1:nwin_new)
  iR(1:nwin) = iR_new(1:nwin_new)

  if (DEBUG) write(*,*) 'DEBUG : prominence of central maximum rejection retained ', &
                         nwin, ' acceptable windows'

  end subroutine

  ! =================================================

  subroutine curtail_window_length(nwin,iL,iR,nmax,maxima_lp,cl_value,cr_value)
  use seismo_variables
  implicit none

  integer, intent(inout) :: nwin
  integer, dimension(*), intent(inout) :: iL, iR
  integer, intent(in) :: nmax
  integer, dimension(*), intent(in) :: maxima_lp
  double precision, intent(in) :: cl_value, cr_value


  integer :: iwin, i_left, i_right, imax, n_left, n_right
  double precision :: time_decay_left, time_decay_right, delta_left, delta_right

  time_decay_left = WIN_MIN_PERIOD * cl_value / dt
  time_decay_right = WIN_MIN_PERIOD * cr_value / dt

  n_left = 0
  n_right = 0
  do iwin = 1,nwin

    i_right = 0
    i_left = 0
    ! find index of maximum on the right of left boundary
    do imax = 1, nmax
      if (maxima_lp(imax) .gt. iL(iwin)) then
        i_left = maxima_lp(imax)
        exit
      endif
    enddo
    ! find index of maximum on the left of right boundary
    do imax = nmax, 1, -1
      if (maxima_lp(imax) .lt. iR(iwin)) then
        i_right = maxima_lp(imax)
        exit
      endif
    enddo

    delta_left = i_left - iL(iwin)
    delta_right = iR(iwin) - i_right

    ! check condition
    if (delta_left .gt. time_decay_left) then
      n_left = n_left+1
      iL(iwin) = int(i_left - time_decay_left)
    end if
    if (delta_right .gt. time_decay_right) then
      n_right = n_right+1
      iR(iwin) = int(i_right + time_decay_right)
    end if

  enddo

  if (DEBUG) write(*,*) 'DEBUG : curtailed ', n_left, ' windows on the left'
  if (DEBUG) write(*,*) 'DEBUG : curtailed ', n_right, ' windows on the right'

  end subroutine

  !------------------------------------------------

  subroutine reject_on_fit_criteria(nwin,iM,iL,iR,CC_local,Tshift_local, dlnA_local)
  use seismo_variables
  implicit none

  integer, intent(inout) :: nwin
  integer, dimension(*), intent(inout) :: iM, iL, iR
  double precision, dimension(*), intent(out) :: CC_local, Tshift_local, dlnA_local

  integer :: iwin, nwin_new
  integer, dimension(NWINDOWS) :: iM_new, iL_new, iR_new
  double precision :: CC_temp, Tshift_temp, dlnA_temp
  double precision :: tshift_min, tshift_max, dlnA_min, dlnA_max
  logical :: accept

  nwin_new = 0
  ! for each proto-window, check that all included maxima are part of the
  ! same phase group as the central maximum iM
  do iwin = 1, nwin
    accept=.true.

    ! check the conditions, and set accept to false if rejection
    call calc_criteria(obs_lp,synt_lp,npts,iL(iwin),iR(iwin),dt,Tshift_temp,CC_temp,dlnA_temp)

    ! here we allow for a systematic shift in the time shift and amplitude measurements
    tshift_min = TSHIFT_REFERENCE - TSHIFT_LIMIT(iM(iwin))
    tshift_max = TSHIFT_REFERENCE + TSHIFT_LIMIT(iM(iwin))
    dlnA_min   = DLNA_REFERENCE - DLNA_LIMIT(iM(iwin))
    dlnA_max   = DLNA_REFERENCE + DLNA_LIMIT(iM(iwin))

    if ( (Tshift_temp .lt. tshift_min) .or. (Tshift_temp .gt.tshift_max ) ) then
       write(*,*) iwin, ' : rejection based on not satisfying TSHIFT_MIN < TSHIFT < TSHIFT_MAX'
       write(*,*) 'Tshift : ',tshift_min , Tshift_temp, tshift_max
       accept = .false.
    endif

    if ( accept .and. ( (dlnA_temp .lt. dlnA_min) .or. (dlnA_temp .gt. dlnA_max ) )) then
       write(*,*) iwin, ' : rejection based on not satisfying DLNA_MIN < DLNA < DLNA_MAX'
       write(*,*) 'dlnA : ', dlnA_min , dlnA_temp, dlnA_max
       accept = .false.
    endif

    if( accept .and. (CC_temp .lt. CC_LIMIT(iM(iwin))) ) then
       write(*,*) iwin, ' : rejection based on CC < CC_MIN'
       write(*,*) 'CC : ', CC_temp, CC_LIMIT(iM(iwin))
       accept = .false.
    endif

    ! if the proto-window is acceptable, then accept it
    if (accept) then
      nwin_new = nwin_new + 1
      iM_new(nwin_new) = iM(iwin)
      iR_new(nwin_new) = iR(iwin)
      iL_new(nwin_new) = iL(iwin)
      CC_local(nwin_new) = CC_temp
      Tshift_local(nwin_new) = Tshift_temp
      dlnA_local(nwin_new) = dlnA_temp
    endif
  enddo

  ! update the iM, iR, iL arrays to contain only accepted proto-windows
  nwin = nwin_new
  iM(1:nwin) = iM_new(1:nwin_new)
  iL(1:nwin) = iL_new(1:nwin_new)
  iR(1:nwin) = iR_new(1:nwin_new)

  if (DEBUG) write(*,*) 'DEBUG : fit criteria rejection retained ', nwin, ' acceptable windows'

!!$  ! CHT: this is done by display_windows_full
!!$  if(DEBUG) then
!!$    write(*,'(" DEBUG : iwin, iL, iM, iR,  CC,  Tshift,  dlnA")')
!!$    do iwin = 1, nwin
!!$      write(*,'(" DEBUG : ",i2,1x,3(i4,1x),3(f6.2,1x))') iwin, iL(iwin), iM(iwin), iR(iwin), &
!!$            CC_local(iwin), Tshift_local(iwin), dlnA_local(iwin)
!!$    enddo
!!$  endif

  end subroutine



  subroutine sort_on_start_time(nwin,iM,iL,iR)
  use seismo_variables
  implicit none

  integer, intent(inout) :: nwin
  integer, dimension(*), intent(inout) :: iM, iL, iR

  integer :: iwin, jwin, iL_early, iL_early_index
  integer, dimension(NWINDOWS) :: iM_new, iL_new, iR_new, order_new
  logical, dimension(NWINDOWS) :: dealt_with

  ! initialise dealt_with to false
  dealt_with(1:nwin) = .false.

  ! iterate of windows
  do iwin = 1, nwin
    ! initialise iL_early to a large number
    iL_early = NDIM
    ! iterate over all non dealt with windows to find earliest time
    do jwin = 1, nwin
      if ((.not. dealt_with(jwin)) .and. iL(jwin) .lt. iL_early) then
        iL_early_index = jwin
        iL_early = iL(jwin)
      endif
    enddo
    ! set this earliest window to dealt with
    dealt_with(iL_early_index) = .true.
    ! put this earliest window in the ordered list
    order_new(iwin) = iL_early_index
  enddo

  ! reconstruct the sorted arrays using order_new
  do iwin = 1,nwin
    iL_new(iwin) = iL(order_new(iwin))
    iM_new(iwin) = iM(order_new(iwin))
    iR_new(iwin) = iR(order_new(iwin))
  enddo

  ! update the iM, iR, iL arrays to contain only accepted proto-windows
  iM(1:nwin) = iM_new(1:nwin)
  iL(1:nwin) = iL_new(1:nwin)
  iR(1:nwin) = iR_new(1:nwin)

  if (DEBUG) write(*,*) 'DEBUG : sorted ', nwin, ' windows'

  end subroutine


  subroutine resolve_overlaps(nwin, iM, iL, iR, CC_local, Tshift_local, dlnA_local)
  use seismo_variables
  implicit none
  integer, intent(inout) :: nwin
  integer, dimension(*), intent(inout) :: iM, iL, iR
  double precision, dimension(*), intent(inout) :: CC_local, Tshift_local, dlnA_local

  integer :: iwin, iwin2
  integer :: nwin_new, n_groups
  integer, dimension(NWINDOWS) :: iM_new, iL_new, iR_new
  double precision, dimension(NWINDOWS) :: CC_new, Tshift_new, dlnA_new
  integer, dimension(NWINDOWS,NWINDOWS) :: groups
  integer, dimension(NWINDOWS) :: iL_group, iR_group
  integer, dimension(NWINDOWS) :: group_size
  logical, dimension(NWINDOWS) :: assigned

  ! for combinations of overlapping windows
  logical :: iterate
  integer :: ii, igroup, igroup_window
  integer :: n_comb, icomb, win_length, group_length
  integer :: prev_comb_start, prev_comb_end
  integer, dimension(NWINDOWS) :: iL_comb, iR_comb
  integer, dimension(NWINDOWS) :: comb_size
  integer :: chosen_combination
  double precision, dimension(NWINDOWS) :: space_coverage
  double precision, dimension(NWINDOWS) :: average_CC
  double precision, dimension(NWINDOWS) :: n_window_fraction
  double precision, dimension(NWINDOWS) :: score
  integer, dimension(NWINDOWS,NWINDOWS) :: comb

  nwin_new = 0

  ! initialise assigned for each window to false
  assigned(1:nwin) = .false.

  ! create groups
  n_groups = 0
  do iwin = 1, nwin
    if (.not.assigned(iwin)) then
      ! start a group with this window
      assigned(iwin) = .true.
      n_groups = n_groups+1
      !checks bounds
      if( n_groups > NWINDOWS ) then
        write(*,*) '  too many groups: set NWINDOWS higher!'
        stop 'error too many groups (n_groups)'
      endif
      group_size(n_groups) = 1
      ! initialise group listing
      groups(n_groups,group_size(n_groups)) = iwin
      ! set left and right limits of group
      iL_group(n_groups) = iL(iwin)
      iR_group(n_groups) = iR(iwin)
      ! check other windows for any that overlap with this one,
      ! ignoring those that have already been assigned to groups
      do iwin2 = 1, nwin
        if (.not.assigned(iwin2) &
            .and. iR(iwin2).gt.iL_group(n_groups) &
            .and. iL(iwin2).lt.iR_group(n_groups)) then
          ! assign this window to the current overlap list
          assigned(iwin2) = .true.
          group_size(n_groups) = group_size(n_groups) + 1
          groups(n_groups,group_size(n_groups)) = iwin2
          ! update the group boundaries
          if ( iL(iwin2).lt.iL_group(n_groups) ) iL_group(n_groups) = iL(iwin2)
          if ( iR(iwin2).gt.iR_group(n_groups) ) iR_group(n_groups) = iR(iwin2)
        endif
      enddo
    endif
  enddo

  if (DEBUG) write(*,*) 'DEBUG : windows fall into ', n_groups, ' groups'

  if (DEBUG) then
    do igroup = 1, n_groups
      write(*,*) 'DEBUG : ', igroup, ' : ', groups(igroup,1:group_size(igroup))
    enddo
  endif

  ! create combinations for each group
  do igroup = 1, n_groups
    if(DEBUG) write(*,*) 'DEBUG : Analysing group ', igroup
    n_comb = 0
    iterate = .true.
    ! add the single windows in this group to the combination
    do igroup_window = 1, group_size(igroup)
      ! extract the window from the group
      iwin = groups(igroup,igroup_window)
      ! increment the number of combinations
      ! checks bounds
      if (n_comb .ge. NWINDOWS) then
        print*,'limiting number of combinations: n_comb=',n_comb
        !stop 'Too many window combinations'
        exit
      endif
      n_comb = n_comb+1
      comb_size(n_comb) = 1
      ! set the combination to this window
      comb(n_comb,comb_size(n_comb)) = iwin
      ! set the combination limits
      iL_comb(n_comb) = iL(iwin)
      iR_comb(n_comb) = iR(iwin)
    enddo
    ! set the counters which keep track of the additions made to the combination
    ! list in the previous iteration
    prev_comb_start = 1
    prev_comb_end   = n_comb
    do while (iterate)
      ! for each of the current combinations
      do icomb = prev_comb_start, prev_comb_end
        ! check if we can add another window to this combination (only on the
        ! right, so only if this combination does not already reach the group
        ! right limit)
        if (iR_comb(icomb).lt.iR_group(igroup)) then
          do igroup_window = 1, group_size(igroup)
            iwin=groups(igroup,igroup_window)
            ! if the window is to the right of (or adjacent to) the current combination
            if(iL(iwin).ge.iR_comb(icomb)) then
              ! make a new combination with this window added
              ! checks bounds
              if( n_comb >= NWINDOWS ) then
                write(*,*) '  limits combinations: set NWINDOWS higher!'
                exit
                !stop 'error too many combinations (n_comb)'
              endif
              n_comb = n_comb+1
              ! copy the old combination to the new combination
              comb(n_comb,1:comb_size(icomb)) = comb(icomb,1:comb_size(icomb))
              ! extend it by one window
              comb_size(n_comb) = comb_size(icomb)+1
              comb(n_comb,comb_size(n_comb)) = iwin
              ! set its limits
              iL_comb(n_comb) = iL_comb(icomb)
              iR_comb(n_comb) = iR(iwin)
            endif
          enddo
        endif
        if(n_comb >= NWINDOWS) exit
      enddo
      ! check if we have added some new combinations
      if (n_comb.gt.prev_comb_end) then
        ! update the counters
        prev_comb_start = prev_comb_end+1
        prev_comb_end = n_comb
      else
        ! have not added any new combinations - stop iteration
        iterate = .false.
      endif
    enddo ! end of iterative combination construction loop

    ! start the scoring loop
    group_length = iR_group(igroup) - iL_group(igroup)
    do icomb = 1, n_comb

      ! score the combination length
      win_length = 0
      do ii = 1,comb_size(icomb)
        iwin = comb(icomb,ii)
        ! add the window length to the total
        win_length = win_length + (iR(iwin)-iL(iwin))
      enddo
      space_coverage(icomb) = float(win_length)/float(group_length)

      ! score the average CC
      ! If the longer window has a better cross-correlation than the
      ! average of the constituent windows, then you want the longer window.
      average_CC(icomb) = 0.
      do ii = 1,comb_size(icomb)
        iwin = comb(icomb,ii)
        ! add the window length to the total
        average_CC(icomb) = average_CC(icomb) + CC_local(iwin)
      enddo
      average_CC(icomb) = average_CC(icomb)/float(comb_size(icomb))

      ! score the number of windows in combination
      n_window_fraction(icomb) = 1.0 - dble(comb_size(icomb))/dble(group_size(igroup))

      ! make the summary score
      score(icomb) = ( WEIGHT_SPACE_COVERAGE*space_coverage(icomb) &
                    + WEIGHT_AVERAGE_CC*average_CC(icomb) &
                    + WEIGHT_N_WINDOWS*n_window_fraction(icomb) ) / &
                    (WEIGHT_SPACE_COVERAGE + WEIGHT_AVERAGE_CC + WEIGHT_N_WINDOWS)

    enddo ! end the scoring loop

    ! Sanity check debug output
    if (DEBUG) then
      write(*,*) 'checking the scoring for all combinations'
      do icomb = 1, n_comb
        write(*,'(a,i6,a,1f8.4,a,3f8.4,a)') ' DEBUG : ', icomb, &
        ' (',score(icomb),' = ',space_coverage(icomb),average_CC(icomb),n_window_fraction(icomb),') : '
        write(*,*) 'DEBUG : Combination is :', comb(icomb,1:comb_size(icomb))
      enddo
    endif

    ! find the best-scoring combination
    chosen_combination = 1
    if (n_comb.gt.1) then
      do icomb = 2, n_comb
        if (score(icomb).gt.score(chosen_combination)) chosen_combination = icomb
      enddo
    endif

    if (DEBUG) write(*,*) 'DEBUG : Chosen combination is ', chosen_combination

    ! loop over chosen combination, adding its windows to the final list of windows
    do ii = 1,comb_size(chosen_combination)
      iwin = comb(chosen_combination,ii)
      ! checks bounds
      if( nwin_new >= NWINDOWS ) exit
      nwin_new = nwin_new + 1
      iL_new(nwin_new) = iL(iwin)
      iR_new(nwin_new) = iR(iwin)
      iM_new(nwin_new) = iM(iwin)
      CC_new(nwin_new) = CC_local(iwin)
      Tshift_new(nwin_new) = Tshift_local(iwin)
      dlnA_new(nwin_new) = dlnA_local(iwin)
    enddo

  enddo ! end of loop over groups

  ! update the iM, iR, iL arrays to contain only accepted proto-windows
  ! update also the CC, dlnA and Tshift arrays
  nwin = nwin_new
  iM(1:nwin) = iM_new(1:nwin_new)
  iL(1:nwin) = iL_new(1:nwin_new)
  iR(1:nwin) = iR_new(1:nwin_new)
  CC_local(1:nwin) = CC_new(1:nwin_new)
  dlnA_local(1:nwin) = dlnA_new(1:nwin_new)
  Tshift_local(1:nwin) = Tshift_new(1:nwin_new)

  if (DEBUG) write(*,*) 'DEBUG : overlap resolution retained ', nwin, ' preferred windows'

  end subroutine

!----------------------------

  subroutine display_windows_full(nwin,iM,iL,iR,CC_local,Tshift_local,dlnA_local)
    use seismo_variables
    implicit none

    integer, intent(in) :: nwin
    integer, dimension(*), intent(in) :: iM, iL, iR
    double precision, dimension(NWINDOWS), intent(in) :: CC_local, Tshift_local, dlnA_local
    integer :: iwin

    if(DEBUG) then
       write(*,'(" DEBUG : iwin, iL, iM, iR, t_start, t_end, Tshift, dlnA, CC")')
       do iwin = 1, nwin
          write(*,'(" DEBUG : ",i4,1x,3(i4,1x),2(f8.2,1x),3(f8.2,1x))') iwin, iL(iwin), iM(iwin), &
               iR(iwin), b+(iL(iwin)-1)*dt, b+(iR(iwin)-1)*dt, &
               Tshift_local(iwin), dlnA_local(iwin), CC_local(iwin)
       enddo
    endif

  end subroutine display_windows_full

!----------------------------

 subroutine display_windows(nwin,iM,iL,iR)
   use seismo_variables
   implicit none

   integer, intent(in) :: nwin
   integer, dimension(*), intent(in) :: iM, iL, iR
   integer :: iwin

   if(DEBUG) then
      write(*,'(" DEBUG : iwin, iL, iM, iR, t_start, t_end")')
      do iwin = 1, nwin
         write(*,'(" DEBUG : ",i4,1x,3(i4,1x),2(f8.2,1x))') iwin, iL(iwin), iM(iwin), iR(iwin), b+(iL(iwin)-1)*dt, b+(iR(iwin)-1)*dt
      enddo
   endif

 end subroutine display_windows
