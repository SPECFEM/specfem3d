!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 1
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology October 2002
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

!
! solver in single or double precision depending on the machine (4 or 8 bytes)
!
!  ALSO CHANGE FILE  precision.h ACCORDINGLY
!
  integer, parameter :: SIZE_REAL = 4
  integer, parameter :: SIZE_DOUBLE = 8

! uncomment this to run in single precision
! integer, parameter :: CUSTOM_REAL = SIZE_REAL
! uncomment this to run in double precision (increases memory size by 2)
  integer, parameter :: CUSTOM_REAL = SIZE_DOUBLE

!----------- parameters that can be changed by the user -----------

!! DK DK UGLY for tsurf source, multiply moment by mu = rho cs^2
  logical, parameter :: MULTIPLY_MU_TSURF = .false.

! set to .true.  if running on a Beowulf-type machine with local disks
! set to .false. if running on a shared-memory machine with common file system
! if running on a Beowulf, also modify name of nodes in filter_machine_file.f90
  logical, parameter :: USE_LOCAL_PATH_BEOWULF = .true.

! on some processors (e.g. Pentiums) it is necessary to suppress underflows
! by using a small initial field instead of zero
  logical, parameter :: FIX_UNDERFLOW_PROBLEM = .true.

! save AVS or OpenDX files in mesher or not
! do not use if you do not plan to use AVS or OpenDX to visualize the mesh
  logical, parameter :: SAVE_AVS_DX_MESH_FILES = .true.

! minimum thickness in meters to include the effect of the oceans
! to avoid taking into account spurious oscillations in topography model
  double precision, parameter :: MINIMUM_THICKNESS_3D_OCEANS = 10.d0

! min and max density in the model
  double precision, parameter :: DENSITY_MAX = 3000.d0
  double precision, parameter :: DENSITY_MIN = 2000.d0

! density of sea water
  real(kind=CUSTOM_REAL), parameter :: RHO_OCEANS = 1020.0

! extend basin model below threshold and above topography to make sure
! there is no small gap between interpolated maps and sediments
  logical, parameter :: EXTEND_VOXET_BELOW_BASEMENT = .true.
  logical, parameter :: EXTEND_VOXET_ABOVE_TOPO = .true.
  double precision, parameter :: DISTMAX_ASSUME_SEDIMENTS = 210.d0
  integer, parameter :: NCELLS_EXTEND = 8

! depth at which we start to honor the basement interface
  double precision, parameter :: Z_THRESHOLD_HONOR_BASEMENT = -4700.d0

! apply heuristic rule to modify doubling regions to balance angles
  logical, parameter :: APPLY_HEURISTIC_RULE = .true.

! to suppress UTM projection for SCEC benchmarks
  logical, parameter :: SUPPRESS_UTM_PROJECTION = .false.

! interval at which we output time step info and max of norm of displacement
  integer, parameter :: ITAFF_TIME_STEPS = 100

! number of GLL points in each direction of an element (degree plus one)
  integer, parameter :: NGLLX = 5
  integer, parameter :: NGLLY = NGLLX
  integer, parameter :: NGLLZ = NGLLX

! number of points per surface element
  integer, parameter :: NGLLSQUARE = NGLLX * NGLLY

! number of points per spectral element
  integer, parameter :: NGLLCUBE = NGLLX * NGLLY * NGLLZ

! input, output and main MPI I/O files
  integer, parameter :: ISTANDARD_OUTPUT = 6
  integer, parameter :: IIN = 40,IOUT = 41
! uncomment this to write messages to a text file
  integer, parameter :: IMAIN = 42
! uncomment this to write messages to the screen
! integer, parameter :: IMAIN = ISTANDARD_OUTPUT

! flag to print the source time function and spectrum
  logical, parameter :: PRINT_SOURCE_TIME_FUNCTION = .false.

!----------- do not modify anything below -------------

! some useful constants
  double precision, parameter :: PI = 3.141592653589793d0,TWO_PI = 2.d0 * PI

! 3-D simulation
  integer, parameter :: NDIM = 3

! dimension of the boundaries of the slices
  integer, parameter :: NDIM2D = 2

! number of nodes for 2D and 3D shape functions for hexahedra
! we use 8-node mesh bricks, which are more stable than 27-node elements
  integer, parameter :: NGNOD = 8, NGNOD2D = 4

  double precision, parameter :: ZERO = 0.d0,ONE = 1.d0,TWO = 2.d0,HALF = 0.5d0

  real(kind=CUSTOM_REAL), parameter :: &
    ONE_THIRD   = 1._CUSTOM_REAL/3._CUSTOM_REAL, &
    FOUR_THIRDS = 4._CUSTOM_REAL/3._CUSTOM_REAL

  double precision, parameter :: HUGEVAL = 1.d+30,TINYVAL = 1.d-9

! declare real value independently of the machine
  real(kind=CUSTOM_REAL), parameter :: HUGEVAL_SNGL = 1.e+30_CUSTOM_REAL

  integer, parameter :: HUGEINT = 100000000

! define flag for elements
  integer, parameter :: IFLAG_ONE_LAYER_TOPOGRAPHY = 1
  integer, parameter :: IFLAG_BASEMENT_TOPO = 2
  integer, parameter :: IFLAG_16km_BASEMENT = 3
  integer, parameter :: IFLAG_MOHO_16km = 4
  integer, parameter :: IFLAG_HALFSPACE_MOHO = 5

! define flag for regions for attenuation
  integer, parameter :: NUM_REGIONS_ATTENUATION = 2

  integer, parameter :: IATTENUATION_SEDIMENTS = 1
  integer, parameter :: IATTENUATION_BEDROCK = 2

! number of standard linear solids in parallel for attenuation
  integer, parameter :: N_SLS = 3

! flag for the four edges of each slice and for the bottom edge
  integer, parameter :: XI_MIN = 1
  integer, parameter :: XI_MAX = 2
  integer, parameter :: ETA_MIN = 3
  integer, parameter :: ETA_MAX = 4
  integer, parameter :: BOTTOM = 5

! flag for projection from latitude/longitude to UTM, and back
  integer, parameter :: ILONGLAT2UTM = 0, IUTM2LONGLAT = 1

! smallest real number on the Pentium and the SGI =  1.1754944E-38
! largest real number on the Pentium and the SGI  =  3.4028235E+38
! small negligible initial value to avoid very slow underflow trapping
! but not too small to avoid trapping on velocity and acceleration in Newmark
  real(kind=CUSTOM_REAL), parameter :: VERYSMALLVAL = 1.E-24_CUSTOM_REAL

! displacement threshold above which we consider the code became unstable
  real(kind=CUSTOM_REAL), parameter :: STABILITY_THRESHOLD = 1.E+25_CUSTOM_REAL

! geometrical tolerance for boundary detection
  double precision, parameter :: SMALLVAL = 0.00001d0

! do not use tags for MPI messages, use dummy tag instead
  integer, parameter :: itag = 0,itag2 = 0

! for the Gauss-Lobatto-Legendre points and weights
  double precision, parameter :: GAUSSALPHA = 0.d0,GAUSSBETA = 0.d0

! number of iterations to solve the system for xi and eta
  integer, parameter :: NUM_ITER = 4

! size of topography and bathymetry file for Southern California
  integer, parameter :: NX_TOPO = 1401,NY_TOPO = 1001
  double precision, parameter :: ORIG_LAT_TOPO = 32.d0
  double precision, parameter :: ORIG_LONG_TOPO = -121.d0
  double precision, parameter :: DEGREES_PER_CELL_TOPO = 5.d0 / 1000.d0

! size of Lupei Zhu's Moho map file for Southern California
  integer, parameter :: NX_MOHO = 71,NY_MOHO = 51
  double precision, parameter :: ORIG_LAT_MOHO = 32.d0
  double precision, parameter :: ORIG_LONG_MOHO = -121.d0
  double precision, parameter :: DEGREES_PER_CELL_MOHO = 0.1d0

! size of basement map file
  integer, parameter :: NX_BASEMENT = 161,NY_BASEMENT = 144
  double precision, parameter :: ORIG_X_BASEMENT = 316000.
  double precision, parameter :: ORIG_Y_BASEMENT = 3655000.
  double precision, parameter :: SPACING_X_BASEMENT = 1000.
  double precision, parameter :: SPACING_Y_BASEMENT = 1000.

!
! new Gocad Voxets Peter July 29, 2002 - high-res and medium-res blocks
!

! size of the medium-resolution Gocad voxet

  integer, parameter :: NX_GOCAD_MR = 194, NY_GOCAD_MR = 196, NZ_GOCAD_MR = 100

  double precision, parameter :: ORIG_X_GOCAD_MR = 283000.
  double precision, parameter :: ORIG_Y_GOCAD_MR = 3655000.
  double precision, parameter :: ORIG_Z_GOCAD_MR = -15000.

  double precision, parameter :: SPACING_X_GOCAD_MR = 1000.
  double precision, parameter :: SPACING_Y_GOCAD_MR = 1000.
  double precision, parameter :: SPACING_Z_GOCAD_MR = 200.

! maximum size of model for tapering of transition between Hauksson and MR
  double precision, parameter :: END_X_GOCAD_MR = ORIG_X_GOCAD_MR + SPACING_X_GOCAD_MR * (NX_GOCAD_MR - 1)
  double precision, parameter :: END_Y_GOCAD_MR = ORIG_Y_GOCAD_MR + SPACING_Y_GOCAD_MR * (NY_GOCAD_MR - 1)

! size of the high-resolution Gocad voxet

  integer, parameter :: NX_GOCAD_HR = 185, NY_GOCAD_HR = 196, NZ_GOCAD_HR = 100

  double precision, parameter :: ORIG_X_GOCAD_HR = 371052.25
  double precision, parameter :: ORIG_Y_GOCAD_HR = 3725250.
  double precision, parameter :: ORIG_Z_GOCAD_HR = -9500.

  double precision, parameter :: SPACING_X_GOCAD_HR = 250.
  double precision, parameter :: SPACING_Y_GOCAD_HR = 250.
  double precision, parameter :: SPACING_Z_GOCAD_HR = 100.

! maximum size of model for tapering of transition between HR and MR
  double precision, parameter :: END_X_GOCAD_HR = ORIG_X_GOCAD_HR + SPACING_X_GOCAD_HR * (NX_GOCAD_HR - 1)
  double precision, parameter :: END_Y_GOCAD_HR = ORIG_Y_GOCAD_HR + SPACING_Y_GOCAD_HR * (NY_GOCAD_HR - 1)

! implement smooth transition between Hauksson, HR and MR Gocad blocks
  logical, parameter :: TAPER_GOCAD_TRANSITIONS = .true.

!
!--- larger Hauksson model for entire So-Cal, 15 km resolution
!

! number of non-constant layers
  integer, parameter :: NLAYERS_HAUKSSON = 9
! depth of layers
  double precision, parameter :: Z_HAUKSSON_LAYER_1 =  -1000.d0
  double precision, parameter :: Z_HAUKSSON_LAYER_2 =  -4000.d0
  double precision, parameter :: Z_HAUKSSON_LAYER_3 =  -6000.d0
  double precision, parameter :: Z_HAUKSSON_LAYER_4 = -10000.d0
  double precision, parameter :: Z_HAUKSSON_LAYER_5 = -15000.d0
  double precision, parameter :: Z_HAUKSSON_LAYER_6 = -17000.d0
  double precision, parameter :: Z_HAUKSSON_LAYER_7 = -22000.d0
  double precision, parameter :: Z_HAUKSSON_LAYER_8 = -31000.d0
  double precision, parameter :: Z_HAUKSSON_LAYER_9 = -33000.d0

  integer, parameter :: NGRID_NEW_HAUKSSON = 201

! corners of new Hauksson's interpolated grid
  double precision, parameter :: UTM_X_ORIG_HAUKSSON = 122035.012d0
  double precision, parameter :: UTM_X_END_HAUKSSON  = 766968.628d0
  double precision, parameter :: UTM_Y_ORIG_HAUKSSON = 3547232.986d0
  double precision, parameter :: UTM_Y_END_HAUKSSON  = 4098868.501d0

  double precision, parameter :: SPACING_UTM_X_HAUKSSON = (UTM_X_END_HAUKSSON - UTM_X_ORIG_HAUKSSON) / (NGRID_NEW_HAUKSSON-1.d0)
  double precision, parameter :: SPACING_UTM_Y_HAUKSSON = (UTM_Y_END_HAUKSSON - UTM_Y_ORIG_HAUKSSON) / (NGRID_NEW_HAUKSSON-1.d0)

! layers in the So-Cal regional model
  double precision, parameter :: DEPTH_5p5km_SOCAL = -5500.d0
  double precision, parameter :: DEPTH_16km_SOCAL = -16000.d0
  double precision, parameter :: DEPTH_MOHO_SOCAL = -35000.d0

! reference surface of the model before adding topography
  double precision, parameter :: Z_SURFACE = 0.d0

! for vectorization of loops
  integer, parameter :: NGLLSQUARE_NDIM = NGLLSQUARE * NDIM
  integer, parameter :: NGLLCUBE_NDIM = NGLLCUBE * NDIM

! number of points in each AVS or OpenDX quadrangular cell for movies
  integer, parameter :: NGNOD2D_AVS_DX = 4

! magic ratio for heuristic rule
! this gives 120 degree angles in doubling
! standard value 0.5 gives 135-135-90, which is not optimal
  double precision, parameter :: MAGIC_RATIO = 0.6056d0

! type of elements for heuristic rule
  integer, parameter :: ITYPE_UNUSUAL_1  = 1
  integer, parameter :: ITYPE_UNUSUAL_1p = 2
  integer, parameter :: ITYPE_UNUSUAL_4  = 3
  integer, parameter :: ITYPE_UNUSUAL_4p = 4

