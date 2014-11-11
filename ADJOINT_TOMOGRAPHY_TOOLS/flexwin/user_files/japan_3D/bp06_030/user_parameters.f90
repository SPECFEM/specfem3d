  module user_parameters
 
  ! -------------------------------------------------------------
  ! boolean parameters

  logical, parameter :: DEBUG = .true.
  logical, parameter :: MAKE_SEISMO_PLOTS = .true.
  logical, parameter :: MAKE_WINDOW_FILES = .true.
  logical, parameter :: BODY_WAVE_ONLY = .true.

  ! -------------------------------------------------------------
  ! period min/max for filtering
    double precision, parameter :: WIN_MIN_PERIOD = 6
    double precision, parameter :: WIN_MAX_PERIOD = 30

  ! -------------------------------------------------------------
  ! E(t) water level
    double precision, parameter :: STALTA_BASE = 0.08

  ! -------------------------------------------------------------
  ! TSHIFT
    double precision, parameter :: TSHIFT_BASE = 2.0

  ! -------------------------------------------------------------
  ! limit on CC for window acceptance
    double precision, parameter :: CC_BASE = 0.7

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
    double precision, parameter :: C_1  = 4.0
    double precision, parameter :: C_2  = 1.0
    double precision, parameter :: C_3a = 3.0 ! height condition
    double precision, parameter :: C_3b = 1.4 ! time separation / T0
    double precision, parameter :: C_4a  = 2.0 ! curtail on left
    double precision, parameter :: C_4b  = 6.0 ! curtail on right

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
  integer, parameter :: NWINDOWS = 3000

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

