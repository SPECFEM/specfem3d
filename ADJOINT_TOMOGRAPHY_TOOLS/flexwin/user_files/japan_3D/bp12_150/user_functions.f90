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
  ! read.
  ! -------------------------------------------------------------
  subroutine set_up_criteria_arrays
  use seismo_variables 

  integer :: i
  double precision :: time

  ! for qinya's scsn picking
  double precision :: Pnl_start, S_end, Sw_start, Sw_end

  ! for global
!  double precision :: Rw_end

  ! global:
  ! -------------------------------------------------------------
!  Rw_end=dist_km/3.6
!  do i = 1, npts
!    time=b+(i-1)*dt
!    DLNA_LIMIT(i)=DLNA_BASE
!    CC_LIMIT(i)=CC_BASE
!    TSHIFT_LIMIT(i)=WIN_MIN_PERIOD/2.0
!    if (time.le.Rw_end) then
!      STALTA_W_LEVEL(i)=STALTA_BASE
!    else
!      STALTA_W_LEVEL(i)=2*STALTA_BASE
!    endif
!    if (time.lt. ph_times(1)- 2*WIN_MIN_PERIOD) then
!      STALTA_W_LEVEL(i)=10*STALTA_BASE
!    endif
!  enddo

  ! regional (Qinya's formulation):
  ! -------------------------------------------------------------
  ! see Liu et al. (2004), p. 1755, but note that the PARENTHESES
  ! that are listed in the publication should not be there

     Pnl_start =  P_pick - 60.0
     S_end     =  S_pick + 60.0
     Sw_start  = -15.0 + dist_km/3.5
     Sw_end    =  15 + dist_km/3.2
 

  ! loop over all seismogram points
  do i = 1, npts
     ! calculate time
     time = b+(i-1)*dt
     ! set the values to base ones by default
     DLNA_LIMIT(i) = DLNA_BASE
     TSHIFT_LIMIT(i) = WIN_MIN_PERIOD/2.0
     CC_LIMIT(i) = CC_BASE
     STALTA_W_LEVEL(i) = STALTA_BASE

!!$     ! set time- and distance-specific Tshift and DlnA to mimic Qinya's criteria
!!$     ! (see Liu et al., 2004, p. 1755; note comment above)
!!$     if(time.ge.Pnl_start .and. time.lt.Sw_start) then
!!$        DLNA_LIMIT(i) = 1.5  ! ratio is 2.5, and dlna is ratio-1
!!$        TSHIFT_LIMIT(i) = 3.0 + dist_km/80.0
!!$     endif
!!$     if(time.ge.Sw_start .and. time.le.Sw_end) then
!!$        DLNA_LIMIT(i) = 1.5  ! ratio is 2.5, and dlna is ratio-1
!!$        TSHIFT_LIMIT(i) = 3.0 + dist_km/50.0
!!$     endif
!!$
!!$     double the STA/LTA water level after the surface waves
!!$     if(time.gt.Sw_end) then
!!$       STALTA_W_LEVEL(i) = 2.*STALTA_BASE
!!$     endif
!!$
!!$     ! raises STA/LTA water level after surface wave arrives. 
!!$     if (BODY_WAVE_ONLY) then
!!$        if(time.gt.S_end) then
!!$           STALTA_W_LEVEL(i) = 10.*STALTA_BASE
!!$        endif
!!$     endif

     ! double the STA/LTA water level after the surface waves
    ! if(time.lt.Pnl_start.or.time.gt.S_end) then
    !   STALTA_W_LEVEL(i)=5*STALTA_BASE
    ! endif

  enddo

  end subroutine
  ! -------------------------------------------------------------
