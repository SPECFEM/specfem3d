! -------------------------------------------------------------
! edit here to change T0 and T1 on some condition 
! Note, this function is called AFTER the seismogram has been 
! read but before it is filtered.
! -------------------------------------------------------------

subroutine modify_T0_T1_on_condition
  use seismo_variables

  ! do nothing

  ! adjust fstart and fend accordingly
  !FSTART=1./WIN_MAX_PERIOD
  !FEND=1./WIN_MIN_PERIOD

end subroutine modify_T0_T1_on_condition

! -------------------------------------------------------------
! Edit here to change the time dependent properties of the selection criteria
! Note, this function is called AFTER the seismogram has been read.
! -------------------------------------------------------------
subroutine set_up_criteria_arrays
  use seismo_variables 

  integer :: i
  double precision :: time

  double precision :: Pnl_start, S_end, Sw_start, Sw_end
  double precision :: Nlam, dtresh, vref
 
!===========================

! -----------------------------------------------------------------
! This is the basic version of the subroutine - no variation with time
! -----------------------------------------------------------------
   do i = 1, npts
     time = b+(i-1)*dt
     DLNA_LIMIT(i) = DLNA_BASE
     CC_LIMIT(i) = CC_BASE
     TSHIFT_LIMIT(i) = TSHIFT_BASE       ! WIN_MIN_PERIOD/2.0
     STALTA_W_LEVEL(i) = STALTA_BASE
     S2N_LIMIT(i) = WINDOW_S2N_BASE
   enddo

!!$  if (.not. BODY_WAVE_ONLY) then
!!$     Pnl_start =  -5.0 + dist_km/7.8
!!$     Sw_start  = -15.0 + dist_km/3.5
!!$     Sw_end    =  35.0 + dist_km/3.1
!!$  else
!!$     Pnl_start =  P_pick - 5.0
!!$     S_end     =  S_pick + 5.0
!!$     Sw_start  = -15.0 + dist_km/3.5
!!$     Sw_end    =  35.0 + dist_km/3.1
!!$  endif

  ! regional (Qinya's formulation):
  ! -------------------------------------------------------------
  ! see Liu et al. (2004), p. 1755, but note that the PARENTHESES
  ! that are listed in the publication should not be there
  ! THESE ARE PROBABLY NOT ACCURATE ENOUGH FOR LONGER PATHS.

  Sw_start  = -15.0 + dist_km/3.5
  Sw_end    =  35.0 + dist_km/3.1

  if (BODY_WAVE_ONLY) then
     !Pnl_start =  P_pick - 5.0
     !S_end     =  S_pick + 5.0
     Pnl_start =  P_pick - 2.5*WIN_MIN_PERIOD
     S_end     =  S_pick + 2.5*WIN_MIN_PERIOD

  else
     Pnl_start =  -5.0 + dist_km/7.8
     S_end     =  Sw_start
  endif

  ! variables for signal to noise ratio criteria.
  noise_start  = b
  noise_end    = Pnl_start
  signal_start = noise_end
  signal_end   = Sw_end

  if(DEBUG) then
     if (BODY_WAVE_ONLY) then
         write(*,*) 'DEBUG : P_pick = ', P_pick
         write(*,*) 'DEBUG : S_pick = ', S_pick
     endif
     write(*,*) 'DEBUG : signal_end = ', sngl(signal_end)
     write(*,*) 'DEBUG : noise_end = ', sngl(noise_end)
  endif

 ! --------------------------------
 ! modulate criteria in time
  do i = 1, npts
     time = b+(i-1)*dt     ! time

     ! raises STA/LTA water level before P wave arrival.
     if(time.lt.Pnl_start) then
        STALTA_W_LEVEL(i) = 10.*STALTA_BASE
     endif

     ! raises STA/LTA water level after surface wave arrives
     if (BODY_WAVE_ONLY) then
        if(time.gt.S_end) then
        !if(time.gt.Sw_end) then
           STALTA_W_LEVEL(i) = 10.*STALTA_BASE
        endif
        
     else
!!$        ! set time- and distance-specific Tshift and DlnA to mimic Qinya's criteria
!!$        ! (see Liu et al., 2004, p. 1755; note comment above)
!!$        if(time.ge.Pnl_start .and. time.lt.Sw_start) then
!!$           !DLNA_LIMIT(i) = 1.5  ! ratio is 2.5, and dlna is ratio-1
!!$           TSHIFT_LIMIT(i) = 3.0 + dist_km/80.0
!!$        endif
!!$        if(time.ge.Sw_start .and. time.le.Sw_end) then
!!$           !DLNA_LIMIT(i) = 1.5  ! ratio is 2.5, and dlna is ratio-1
!!$           TSHIFT_LIMIT(i) = 3.0 + dist_km/50.0
!!$        endif

        ! double the STA/LTA water level after the surface waves
        if(time.gt.Sw_end) then
           STALTA_W_LEVEL(i) = 10.0*STALTA_BASE
        endif

!!$        ! allow 100 seconds to possibly capture additional phases
!!$        if(time.gt. (Sw_end+100.0) ) then
!!$           STALTA_W_LEVEL(i) = 10.*STALTA_BASE
!!$        endif

     endif

  enddo

 ! --------------------------------
 ! if the distance to the station is less than N wavelengths, then reject records
 ! by reasing the entire water level

  Nlam = 1.7    ! number of wavelengths
  vref = 2.0    ! reference velocity, km/s
  dtresh = Nlam*WIN_MIN_PERIOD*vref
  if (dist_km .le. dtresh ) then
     if(DEBUG) then
         write(*,*) 'REJECT by raising water level: station is too close for this period range'
         write(*,*) 'dist_km, dtresh = Nlam*WIN_MIN_PERIOD, Nlam, WIN_MIN_PERIOD :'
         write(*,'(4f12.4)') dist_km, dtresh, Nlam, WIN_MIN_PERIOD
     endif
     do i = 1,npts
        STALTA_W_LEVEL(i) = 10.*STALTA_BASE
     enddo
  endif

! The following is for check_window quality_s2n

! -----------------------------------------------------------------
! Start of user-dependent portion

! This is where you modulate the time dependence of the selection
! criteria.  You have access to the following parameters from the 
! seismogram itself:
!
! dt, b, kstnm, knetwk, kcmpnm
! evla, evlo, stla, stlo, evdp, azimuth, backazimuth, dist_deg, dist_km
! num_phases, ph_names, ph_times
!
! Example of modulation:
!-----------------------
! To increase s2n limit after arrival of R1 try
!
! R_vel=3.2
! R_time=dist_km/R_vel
! do i = 1, npts
!   time=b+(i-1)*dt
!   if (time.gt.R_time) then
!     S2N_LIMIT(i)=2*WINDOW_S2N_BASE
!   endif
! enddo
!
! End of user-dependent portion
! -----------------------------------------------------------------

end subroutine set_up_criteria_arrays
! -------------------------------------------------------------
