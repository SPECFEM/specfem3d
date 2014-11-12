module ma_variables

  use ma_constants
!
! multi-taper measurements
!
! Ying Zhou: The fit between the recovered data and the data can be improved
! by either increasing the window width (HWIN above) or by decreasing NPI.
! In her experience, NPI = 2.5 is good for noisy data.
! For synthetic data, we can use a lower NPI.
! number of tapers should be fixed as twice NPI -- see Latex notes
!
! See write_par_file.pl and measure_adj.f90

  character(len=150) :: OUT_DIR

  double precision :: TLONG, TSHORT
  double precision :: WTR, NPI, DT_FAC, ERR_FAC, DT_MAX_SCALE, NCYCLE_IN_WINDOW
  !double precision :: BEFORE_QUALITY, AFTER_QUALITY, BEFORE_TSHIFT, AFTER_TSHIFT
  double precision :: TSHIFT_MIN, TSHIFT_MAX, DLNA_MIN, DLNA_MAX, CC_MIN
  double precision :: DT_SIGMA_MIN, DLNA_SIGMA_MIN

  integer :: ntaper, ipwr_t, ipwr_w, ERROR_TYPE
  integer :: imeas0, imeas, itaper, is_mtm0, is_mtm

  logical :: DISPLAY_DETAILS,OUTPUT_MEASUREMENT_FILES,RUN_BANDPASS,COMPUTE_ADJOINT_SOURCE,USE_PHYSICAL_DISPERSION

end module ma_variables


module ma_weighting

! module for weighting/normalizing measurements

  logical,parameter :: DO_WEIGHTING = .false.

  ! transverse, radial and vertical weights
  double precision :: weight_T, weight_R, weight_Z
  ! body waves: number of picks on vertical, radial and transverse component
  double precision :: num_P_SV_V,num_P_SV_R,num_SH_T
  ! surface waves: number of pick on vertical, radial and transverse
  double precision :: num_Rayleigh_V,num_Rayleigh_R,num_Love_T

  ! typical surface wave speed in km/s, to calculate surface wave arrival times
  ! Love waves faster than Rayleigh
  double precision, parameter :: surface_vel = 4.0

  ! wave type pick
  integer, parameter :: P_SV_V = 1
  integer, parameter :: P_SV_R = 2
  integer, parameter :: SH_T = 3
  integer, parameter :: Rayleigh_V = 4
  integer, parameter :: Rayleigh_R = 5
  integer, parameter :: Love_T = 6

end module ma_weighting
