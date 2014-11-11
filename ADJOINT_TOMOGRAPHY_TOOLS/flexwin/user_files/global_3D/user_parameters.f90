  module user_parameters
 
  !===================================================================
  ! filter parameters for xapiir subroutine (filter type is BP)
  double precision, parameter :: TRBDNDW = 0.3
  double precision, parameter :: APARM = 30.
  integer, parameter :: IORD = 4
  integer, parameter :: PASSES = 2


  ! -------------------------------------------------------------
  ! array dimensions
  integer, parameter :: NDIM = 32000
  integer, parameter :: NWINDOWS = 8000

  ! -------------------------------------------------------------
  ! miscellaneous - do not modify!
  ! -------------------------------------------------------------

  ! mathematical constants
  double precision, parameter :: PI = 3.1415926535897
  double precision, parameter :: E  = 2.7182818284590

  ! filter types
  integer, parameter :: HANNING = 1
  integer, parameter :: HAMMING = 2
  integer, parameter :: COSINE  = 3

  ! -------------------------------------------------------------

  end module user_parameters

