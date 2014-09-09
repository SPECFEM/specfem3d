
! -------------------------------------------------------------
! edit here to change T0 and T1 on some condition 
! Note, this function is called AFTER the seismogram has been 
! read but before it is filtered.
! -------------------------------------------------------------

  subroutine modify_T0_T1_on_condition
  use seismo_variables
  implicit none

  ! do nothing

  ! adjust fstart and fend accordingly
  FSTART=1./WIN_MAX_PERIOD
  FEND=1./WIN_MIN_PERIOD

  end subroutine

! -------------------------------------------------------------
! edit here to change the time dependent properties of the 
! selection criteria
! Note, this function is called AFTER the seismogram has been 
! read and filtered.
! -------------------------------------------------------------

  subroutine set_up_criteria_arrays
  use seismo_variables 
  implicit none

  integer :: i
  double precision :: time
  double precision :: R_vel, R_time
  double precision :: Q_vel, Q_time


! -----------------------------------------------------------------
! This is the basic version of the subroutine - no variation with time
! -----------------------------------------------------------------
  do i = 1, npts
    time=b+(i-1)*dt
    DLNA_LIMIT(i)=DLNA_BASE
    CC_LIMIT(i)=CC_BASE
    TSHIFT_LIMIT(i)=TSHIFT_BASE
    STALTA_W_LEVEL(i)=STALTA_BASE
    S2N_LIMIT(i)=WINDOW_S2N_BASE
  enddo

  ! these values will be used for signal2noise calculations
  ! if DATA_QUALITY=.true.
  if (DATA_QUALITY) then
    noise_start=b
    noise_end=max(ph_times(1)-WIN_MIN_PERIOD,b+dt)
    signal_start=noise_end
    signal_end=b+(npts-1)*dt
  endif


! -----------------------------------------------------------------
! Start of user-dependent portion

! This is where you reset the signal_end and noise_end values
! if you want to modify them from the defaults.
! This is also where you modulate the time dependence of the selection
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

 ! --------------------------------
 ! Set approximate end of rayleigh wave arrival
 R_vel=3.2
 R_time=dist_km/R_vel
 ! --------------------------------
 ! Set approximate start of love wave arrival
 Q_vel=4.2
 Q_time=dist_km/Q_vel

 ! reset the signal_end time to be the end of the Rayleigh waves
 if (DATA_QUALITY) then
   signal_end=R_time
 endif

 ! --------------------------------
 ! modulate criteria in time
 do i = 1, npts
   time=b+(i-1)*dt
   ! --------------------------------
   ! if we are beyond the Rayleigh wave, then make the all criteria stronger
   ! ratio criterion stronger
   if (time.gt.R_time) then
     S2N_LIMIT(i)=10*WINDOW_S2N_BASE    ! only pick big signals
     CC_LIMIT(i)= 0.95                  ! only pick very similar signals
     TSHIFT_LIMIT(i)=TSHIFT_BASE/3.0    ! only pick small timeshifts
     DLNA_LIMIT(i)=DLNA_BASE/3.0        ! only pick small amplitude anomalies
     STALTA_W_LEVEL(i)=STALTA_BASE*2.0     ! pick only distinctive arrivals
   endif
   ! --------------------------------
   ! if we are in the surface wave times, then make the cross-correlation
   ! criterion less severe
   if (time.gt.Q_time .and. time.lt.R_time) then
     CC_LIMIT(i)=0.9*CC_LIMIT(i)
   endif
   ! --------------------------------
   ! modulate criteria according to event depth
   !
   ! if an intermediate depth event
   if (evdp.ge.70 .and. evdp.lt.300) then
     TSHIFT_LIMIT(i)=TSHIFT_BASE*1.4
   ! if a deep event
   elseif (evdp.ge.300) then
     TSHIFT_LIMIT(i)=TSHIFT_BASE*1.7
   endif
 enddo




!
! End of user-dependent portion
! -----------------------------------------------------------------

  end subroutine
  ! -------------------------------------------------------------
