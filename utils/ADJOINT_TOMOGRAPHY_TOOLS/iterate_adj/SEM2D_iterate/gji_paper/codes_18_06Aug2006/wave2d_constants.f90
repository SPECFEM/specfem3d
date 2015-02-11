module wave2d_constants

!
! GRID, TIME-STEP, AND SOURCE PARAMETERS
!

! NFRAME : number of frames to save
! NSAVE  : timestep increment to save the wavefield
! NSTEP  : number of timesteps
  integer, parameter :: NFRAME = 10    ! 10,12,17
  integer, parameter :: NSAVE  = 400   ! 400
  integer, parameter :: NSTEP  = NFRAME*NSAVE

! time step in seconds
  double precision, parameter :: DT = 6.0d-02 ! (0.02)

! temporal properties of source (source time function)
  integer, parameter :: ISRC_TIME = 1                   ! type (1)
  double precision, parameter :: hdur = 10.0            ! HALF-duration (s)
  double precision, parameter :: tshift = 2.*DT*NSAVE   ! time shift (s)
  !double precision, parameter :: tshift = 8.*hdur
  logical, parameter :: SRC_TAPER  = .true.

! spatial properties of sources
! (1) point source, (2) finite segment, (3) CA shelf boundary, (4) CA coast, (5) finite circle
! (6) a point source EVENT
  integer, parameter :: ISRC_SPACE = 6  ! see wave2d.f90

! spatial properties of receivers
! IREC_SPACE
!   (1) individual station(s)
!   (2) SoCal (used for GJI paper)
!   (3) regular mesh on land
!   (4) regular mesh
! NMESH_REC : determines the number of receivers in a regular mesh (IREC_SPACE=3)
! STATION_GRID_BUFFER : exclude stations within this distance from edge of grid
! STATION_COAST_BUFFER : exclude stations within this distance from edge of coast
  integer, parameter          :: IREC_SPACE = 2 ! see wave2d.f90
  integer, parameter          :: NMESH_REC = 10
  double precision, parameter :: SOURCE_GRID_BUFFER   = 4.0d+03  ! m
  double precision, parameter :: STATION_GRID_BUFFER  = 15.0d+03 ! m
  double precision, parameter :: STATION_COAST_BUFFER = 0.0d+03  ! m

! model specification for c(th,ph)
! (0) het map, (1) homo map, (2) checkerboard, (3) read in
  !integer, parameter          :: IMODEL = 3

! bounds for bandpass filter (in seconds), see also below (fmin,etc)
  double precision, parameter :: hwid = 3.0  ! HALF-width of window
  double precision, parameter :: tmin = 2.*hdur-hwid
  double precision, parameter :: tmax = 2.*hdur+hwid

! mesh specifications
  double precision, parameter :: LENGTH = 480.0d+03 ! m (200)
  double precision, parameter :: HEIGHT = 480.0d+03 ! m (80)
  integer, parameter :: NEX = 40 !40
  integer, parameter :: NEZ = 40 !40
  double precision, parameter :: LAT_MIN = 32.0d0
  double precision, parameter :: LON_MIN = -120.d0

! boolean parameters
! IKER: (0) waveform
!       (1) traveltime, cross-correlation, misfit
!       (2) amplitude, cross-correlation, misfit
!       (3) traveltime, multitaper
!       (4) amplitude, multitaper
!       (5) traveltime, cross-correlation, sampling
!       (6) amplitude, cross-correlation, sampling
  integer, parameter :: IKER = 1
  integer, parameter :: ISURFACE = 1, NCOMP = 1, NABSORB = 4   ! surface waves
!  integer, parameter :: ISURFACE = 0, NCOMP = 3, NABSORB = 3   ! body waves

! iteration and smoothing parameters
!  logical, parameter :: IUPDATE = .false.
  integer, parameter :: NITERATION = 0
  integer, parameter :: POLY_ORDER = 2             ! 2 (equally good) or 3
  !double precision, parameter :: SIGMA = 10.0d+03 ! m

! parameters controlling what to write to file
! NOTE: for the tomography simulations, ALL of these can be .false.

  logical, parameter :: WRITE_STF_F           = .false.
  logical, parameter :: WRITE_SEISMO_F        = .false.    ! true
  logical, parameter :: WRITE_SPECTRA_F       = .false.
  logical, parameter :: WRITE_SPECTRAL_MAP_F  = .false.

  logical, parameter :: WRITE_STF_A           = .false.
  logical, parameter :: WRITE_SEISMO_A        = .false.
  logical, parameter :: WRITE_SPECTRA_A       = .false.
  logical, parameter :: WRITE_SPECTRAL_MAP_A  = .false.

  logical, parameter :: WRITE_KERNELS = .false.    ! kernel snapshots
  logical, parameter :: WRITE_SNAPSHOTS = .false.  ! wavefield snapshots


! MODEL (S.I. units)
  double precision, parameter :: DENSITY = 2.6d+03 ! kg/m^3
  double precision, parameter :: INCOMPRESSIBILITY = 5.2d+10 ! Pa
  double precision, parameter :: RIGIDITY = 2.66d+10 ! Pa

!---------------------------------------------------------------
! CHT: do not change these

! UTM zone for Southern California region
!  integer, parameter :: UTM_PROJECTION_ZONE = 11

! to suppress UTM projection for SCEC benchmarks
  logical, parameter :: SUPPRESS_UTM_PROJECTION = .false.

! flag for projection from latitude/longitude to UTM, and back
  integer, parameter :: ILONGLAT2UTM = 0, IUTM2LONGLAT = 1

! flag for projection from latitude/longitude to mesh-UTM, and back
  integer, parameter :: ILONLAT2MESH = 0, IMESH2LONLAT = 1

! max number of fake receivers
  integer, parameter :: MAX_SR_FAKE = 1000

! max number of events, receivers, and phass
  integer, parameter :: MAX_EVENT = 50
  integer, parameter :: MAX_SR    = 1400
  integer, parameter :: MAX_PHASE = 1
  integer, parameter :: MAX_COMP  = NCOMP

! parameter for FFTW
  integer, parameter :: NOUT = NSTEP/2 + 1

! filter parameters for bandpass
  double precision, parameter :: fmin = 1./tmax, fmax = 1./tmin
  double precision, parameter :: trbdndw = 0.3, a = 30.
  integer, parameter :: passes = 2, iord = 4

!---------------------------------------------------------------
!
! GRID AND GLL POINTS
!
  integer, parameter :: NELE = MAX(NEX,NEZ)
  integer, parameter :: NSPEC = NEX*NEZ

! number of GLL points (polynomial degree plus one)
  integer, parameter :: NGLLX = 5
  integer, parameter :: NGLLZ = 5
  integer, parameter :: NGLL = MAX(NGLLX,NGLLZ)

! number of points per surface element
  integer, parameter :: NGLLSQUARE = NGLLX * NGLLZ

! number of global points
  integer, parameter :: NGLOB = ((NGLLX-1)*NEX + 1)*((NGLLZ-1)*NEZ +1)

! number of local points
  integer, parameter :: NLOCAL = NGLLX * NGLLZ * NSPEC

! number of nodes for 2D and 3D shape functions for hexahedra
! we use 8-node mesh bricks, which are more stable than 27-node elements
  integer, parameter :: NGNOD = 8, NGNOD2D = 4

! number of iterations to solve the system for xi and eta
  integer, parameter :: NUM_ITER = 1

! very large and very small values
  double precision, parameter :: HUGEVAL = 1.d+30, TINYVAL = 1.d-9

! for the Gauss-Lobatto-Legendre points and weights
  double precision, parameter :: GAUSSALPHA = 0.d0,GAUSSBETA = 0.d0

!
! CONSTANTS
!
! pi
  double precision, parameter :: PI = 3.141592653589793d+00
  double precision, parameter :: FOUR_THIRDS = 4.d0/3.d0
  double precision, parameter :: ONE_THIRD = 1.d0/3.d0
  double precision, parameter :: ONEOVERTWO = 0.5d0
  double precision, parameter :: EPS = 1.0d-35
  double precision, parameter :: DEG = 180./PI

! normalization factor of point source force
  double precision, parameter :: FNORM = 1.0d10

! factors from the multitaper method
  integer, parameter :: MAXTAPER=5, NDIM=8000*4, lnpt=14, npt=2**lnpt
  double precision, parameter :: wtr=0.02, ZZIGN=-1.0

end module wave2d_constants
