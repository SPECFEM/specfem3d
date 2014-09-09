  module user_parameters
 
  !===================================================================
  ! filter parameters for xapiir subroutine (filter type is BP)
  double precision, parameter :: TRBDNDW = 0.3
  double precision, parameter :: APARM = 30.0
  integer, parameter :: IORD = 5
  integer, parameter :: PASSES = 2

  ! -------------------------------------------------------------
  ! array dimensions
  ! note that some integer arrays (iM,iL,iR) are NWINDOWS * NWINDOWS
  ! THESE SHOULD PROBABLY BE USER PARAMETERS, SINCE THEY WILL AFFECT
  ! THE SPEED OF THE PROGRAM (ALTHOUGH NOT THE OUTPUT).
  integer, parameter :: NDIM = 10000
  integer, parameter :: NWINDOWS = 2500

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

