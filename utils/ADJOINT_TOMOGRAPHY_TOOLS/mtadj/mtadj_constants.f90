module mtadj_constants

  ! adjoint source types

  integer, parameter :: IKER_WF = 0
  integer, parameter :: IKER_CC = 1
  integer, parameter :: IKER_FD = 2
  character(len=2), parameter :: CKER(3) = (/ 'wf', 'cc', 'fd' /)

  ! taper types
  integer, parameter :: ITAP_BC = 0
  integer, parameter :: ITAP_CS = 1
  integer, parameter :: ITAP_MT = 2
  character(len=2), parameter :: CTAP(3) = (/ 'bc', 'cs', 'mt' /)

  ! maximum number of tapers
  integer, parameter :: NMAX_TAPER = 10

  ! debug flags
  logical, parameter :: DEBUG = .true.
  logical, parameter :: OUTPUT_MEASURE = .true.
  logical, parameter :: OUTPUT_ADJSRC = .true.

  ! constants
  real, parameter :: PI = 3.141592653
  real, parameter :: TWOPI = 2.0 * PI
  complex ,parameter :: CCI = cmplx(0.,1.)
  real, parameter :: LARGE_VAL = 1.0d8
  real, parameter :: EPS_dt=1.0e-4

  ! FFT parameters
  integer, parameter :: LNPT = 13, NPT = 2**LNPT, NDIM = 40000
  real*8, parameter :: FORWARD_FFT = 1.0
  real*8, parameter :: REVERSE_FFT = -1.0

  ! phase correction control parameters, set this between (PI, 2PI),
  ! use a higher value for conservative phase wrapping
  real, parameter :: PHASE_STEP = 1.5 * PI

  ! filter parameters for xapiir subroutine (filter type is BP)
  ! TRBDNDW and APARM are used only for chebyshev filters
  ! 2-passes guarantees a zero-phase filter (non-causal,which is ok)
  real*8, parameter :: TRBDNDW = 0.3
  real*8, parameter :: APARM = 30.
  integer, parameter :: IORD = 4
  integer, parameter :: PASSES = 2


end module mtadj_constants

