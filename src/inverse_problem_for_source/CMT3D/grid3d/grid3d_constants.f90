module grid3d_constants

  ! mathematical constant pi
  real, parameter :: PI = 3.14159265

  ! radian-to-degree conversion
  real, parameter :: R2D = 180./PI
  real, parameter :: D2R = PI/180.

  ! maximum array dimensions
  integer, parameter :: NMW_MAX = 5
  integer, parameter :: NST_MAX =  20
  integer, parameter :: NDIP_MAX = 20
  integer, parameter :: NRAKE_MAX = 40

  ! maximum total array dimension 80000
  integer, parameter :: NMEC_MAX=NMW_MAX*NST_MAX*NDIP_MAX*NRAKE_MAX

  ! maximum number of records (NRECMAX < NWINMAX)
  integer, parameter :: NRECMAX = 1200

  ! maximum number of files/windows
  integer, parameter :: NWINMAX = 1800

  ! maximum number of data points
  integer, parameter :: NDATAMAX=30000

  ! number of pars for moment only
  integer, parameter :: NM = 6

  ! moment element names
  character(len=3), parameter :: PAR_NAME(NM) =  &
       (/'Mrr','Mtt','Mpp','Mrt','Mrp','Mtp'/)

  ! small numbers
  real, parameter :: EPS2 = 1.0e-2
  real, parameter :: EPS5 = 1.0e-5
  real, parameter :: SMALL = -huge(1.0)

  ! io unit for parameter files
  integer, parameter :: IOPAR = 10

  ! io unit for flexwin output files
  integer, parameter :: IOWIN = 20

  ! io unit for input file
  integer, parameter :: IOGRD = 30

  ! io unit for cmt solution
  integer, parameter :: IOCMT = 40

  ! debug boolean
  logical, parameter :: DEBUG = .true.

  ! number of regions for azimuthal weighting
  integer, parameter :: NREGIONS = 10

  ! reference distance for Pnl, Rayleigh and Love wave weighting
  real, parameter :: REF_DIST = 100.0

end module grid3d_constants
