  module user_parameters
 
  ! -------------------------------------------------------------
  ! DEBUG flag
  logical, parameter :: DEBUG = .true.
 
  ! -------------------------------------------------------------
  ! option to pick windows only around body wave arrivals.
  logical, parameter :: BODY_WAVE_ONLY = .false.

  ! -------------------------------------------------------------
  ! period min/max for filtering

  ! SCSN regional events
    double precision, parameter :: WIN_MIN_PERIOD = 6
    double precision, parameter :: WIN_MAX_PERIOD = 50

  ! global events
  ! double precision, parameter :: WIN_MIN_PERIOD = 50
  ! double precision, parameter :: WIN_MAX_PERIOD = 500

  ! -------------------------------------------------------------
  ! E(t) water level
  ! SCSN regional events
    double precision, parameter :: STALTA_BASE = 0.12
  ! global events
  ! double precision, parameter :: STALTA_BASE = 0.07

  ! -------------------------------------------------------------
  ! TSHIFT
  ! SCSN regional events
    double precision, parameter :: TSHIFT_BASE = 3.0
  ! global events
  ! double precision, parameter :: TSHIFT_BASE = 20.0

  ! -------------------------------------------------------------
  ! limit on CC for window acceptance
    double precision, parameter :: CC_BASE = 0.8

  ! -------------------------------------------------------------
  ! limit on dlnA (dA/A) for window acceptance
    double precision, parameter :: DLNA_BASE = 1.0

  ! -------------------------------------------------------------
  ! limit on signal-to-noise on the observed data

  ! boolean switch for check_data_quality
    logical, parameter :: DATA_QUALITY = .true.

  ! if DATA_QUALITY = .true. and if two different measurements of
  ! signal-to-noise ratios exceeds these two base levels,
  ! then the data time series (and syn) is kept
    double precision, parameter :: SNR_INTEGRATE_BASE = 3.0  
    double precision, parameter :: SNR_MAX_BASE = 3.5

  ! -------------------------------------------------------------
  ! Fine tuning constants 
    double precision, parameter :: C_0  = 1.0
    double precision, parameter :: C_1  = 3.0
    double precision, parameter :: C_2  = 0.8
    double precision, parameter :: C_3a = 1.0 ! height condition
    double precision, parameter :: C_3b = 1.0 ! time separation / T0
    double precision, parameter :: C_4a  = 5.0 ! curtail on left
    double precision, parameter :: C_4b  = 5.0 ! curtail on right

    double precision, parameter :: WEIGHT_SPACE_COVERAGE = 1.0
    double precision, parameter :: WEIGHT_AVERAGE_CC = 1.0

  !===================================================================
  ! filter parameters for xapiir subroutine (filter type is BP)
  double precision, parameter :: TRBDNDW = 0.3
  double precision, parameter :: APARM = 30.
  integer, parameter :: IORD = 4
  integer, parameter :: PASSES = 2

  double precision, parameter :: FSTART = 1./WIN_MAX_PERIOD
  double precision, parameter :: FEND = 1./WIN_MIN_PERIOD

  ! -------------------------------------------------------------
  ! array dimensions
  integer, parameter :: NDIM = 32000
  integer, parameter :: NWINDOWS = 2000

  ! -------------------------------------------------------------
  ! miscellaneous

  ! mathematical constants
  double precision, parameter :: PI = 3.1415926535897
  double precision, parameter :: E  = 2.7182818284590

  ! filter types
  integer, parameter :: HANNING = 1
  integer, parameter :: HAMMING = 2
  integer, parameter :: COSINE  = 3

  ! -------------------------------------------------------------

  end module user_parameters

