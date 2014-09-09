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
  double precision :: Pnl_start, Sw_start, Sw_end

  ! for Min's japan picking
  double precision :: P_start

  ! for global
!  double precision :: Rw_end

  ! global:
  ! -------------------------------------------------------------
!  Rw_end=dist_km/3.6
!  do i = 1, npts
!    time=b+(i-1)*dt
!    DLNA_LIMIT(i)=DLNA_BASE
!    CC_LIMIT(i)=CC_BASE
!    TSHIFT_LIMIT(i)=WIN_LP_PERIOD/2.0
!    if (time.le.Rw_end) then
!      STALTA_W_LEVEL(i)=STALTA_BASE
!    else
!      STALTA_W_LEVEL(i)=2*STALTA_BASE
!    endif
!    if (time.lt. ph_times(1)- 2*WIN_LP_PERIOD) then
!      STALTA_W_LEVEL(i)=10*STALTA_BASE
!    endif
!  enddo

  ! regional (Qinya's formulation):
  ! -------------------------------------------------------------
   Pnl_start=-5.0+dist_km/7.5
   Sw_start=-15.0+dist_km/3.5
   Sw_end=35.0+dist_km/3.2
   P_start=ph_times(1)-60
   ! loop over all seismogram points
   do i = 1, npts
     ! calculate time
     time=b+(i-1)*dt
     ! set the values to base ones by default
     DLNA_LIMIT(i)=DLNA_BASE
     TSHIFT_LIMIT(i)=WIN_LP_PERIOD/2.0
     CC_LIMIT(i)=CC_BASE
     STALTA_W_LEVEL(i)=STALTA_BASE
     ! set time and distance specific Tshift and DlnA to mimic Qinya's criteria
     !if(time.ge.Pnl_start .and. time.lt.Sw_start) then
     !  DLNA_LIMIT(i)=1.5  ! ratio is 2.5, and dlna is ratio-1
     !  TSHIFT_LIMIT(i)=3.0+dist_km/80.0
     !endif
     !if(time.ge.Sw_start .and. time.le.Sw_end) then
     !  DLNA_LIMIT(i)=1.0  ! ratio is 2.5, and dlna is ratio-1
     !  TSHIFT_LIMIT(i)=3.0+dist_km/50.0
     !endif
     ! double the STA/LTA water level after the surface waves
     if(time.lt.P_start.or.time.gt.Sw_end) then
       STALTA_W_LEVEL(i)=5*STALTA_BASE
     endif
   enddo

  end subroutine
  ! -------------------------------------------------------------
