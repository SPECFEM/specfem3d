module ma_constants

  ! number of entries in window_chi output file
  integer, parameter :: N_MEASUREMENT = 5
  integer, parameter :: NCHI = 3*(N_MEASUREMENT-1) + 8

  ! constants
  double precision, parameter :: PI = 3.141592653589793d+00
  double precision, parameter :: TWOPI = 2.0 * PI
  complex*16, parameter :: CCI = cmplx(0.,1.)
  double precision, parameter :: LARGE_VAL = 1.0d8

  ! FFT parameters
  integer, parameter :: LNPT = 15, NPT = 2**LNPT, NDIM = 80000
  double precision, parameter :: FORWARD_FFT = 1.0  
  double precision, parameter :: REVERSE_FFT = -1.0   

  ! phase correction control parameters, set this between (PI, 2PI),
  ! use a higher value for conservative phase wrapping
  double precision, parameter :: PHASE_STEP = 1.5 * PI

  ! filter parameters for xapiir bandpass subroutine (filter type is BP)
  ! (These should match the filter used in pre-processing.)
  double precision, parameter :: TRBDNDW = 0.3
  double precision, parameter :: APARM = 30.
  integer, parameter :: IORD = 4
  integer, parameter :: PASSES = 2

  ! takes waveform of first trace dat_dtw, without taking the difference waveform to the second trace syn_dtw
  ! this is useful to cissor out later reflections which appear in data (no synthetics needed)
  logical, parameter :: NO_WAVEFORM_DIFFERENCE = .false. 

  ! constructs adjoint sources for a "ray density" kernel, where all misfits are equal to one
  logical, parameter :: DO_RAY_DENSITY_SOURCE = .false.

end module ma_constants
