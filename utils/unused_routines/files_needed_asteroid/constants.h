!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and University of Pau / CNRS / INRIA
!         (c) California Institute of Technology July 2005
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

! constants.h.  Generated from constants.h.in by configure.

!
! solver in single or double precision depending on the machine (4 or 8 bytes)
!
! ALSO CHANGE FILE precision.h ACCORDINGLY
!
  integer, parameter :: SIZE_REAL = 4
  integer, parameter :: SIZE_DOUBLE = 8

! set to SIZE_REAL to run in single precision
! set to SIZE_DOUBLE to run in double precision (increases memory size by 2)
  integer, parameter :: CUSTOM_REAL = SIZE_REAL

!----------- parameters that can be changed by the user -----------

! set to .false.  if running on a Beowulf-type machine with local disks
! set to .true. if running on a shared-memory machine with common file system
! if running on a Beowulf, also modify name of nodes in filter_machine_file.f90
  logical, parameter :: LOCAL_PATH_IS_ALSO_GLOBAL = .false.

! apply heuristic rule to modify doubling regions to balance angles
  logical, parameter :: APPLY_HEURISTIC_RULE = .true.

! number of GLL points in each direction of an element (degree plus one)
  integer, parameter :: NGLLX = 5
  integer, parameter :: NGLLY = NGLLX
  integer, parameter :: NGLLZ = NGLLX

! input, output and main MPI I/O files
  integer, parameter :: ISTANDARD_OUTPUT = 6
  integer, parameter :: IIN = 40,IOUT = 41
! uncomment this to write messages to a text file
  integer, parameter :: IMAIN = 42
! uncomment this to write messages to the screen
! integer, parameter :: IMAIN = ISTANDARD_OUTPUT
! I/O unit for source and receiver vtk file
  integer, parameter :: IOVTK = 98

! minimum thickness in meters to include the effect of the oceans
! to avoid taking into account spurious oscillations in topography model
  double precision, parameter :: MINIMUM_THICKNESS_3D_OCEANS = 10.d0

! min and max density in the model
  double precision, parameter :: DENSITY_MAX = 3000.d0
  double precision, parameter :: DENSITY_MIN = 2000.d0

! density of sea water
  real(kind=CUSTOM_REAL), parameter :: RHO_OCEANS = 1020.0

! extend model below threshold and above topography to make sure
! there is no small gap between interpolated maps and sediments
  logical, parameter :: EXTEND_VOXET_BELOW_BASEMENT = .true.
  logical, parameter :: EXTEND_VOXET_ABOVE_TOPO = .true.
  double precision, parameter :: DISTMAX_ASSUME_SEDIMENTS = 210.d0
  integer, parameter :: NCELLS_EXTEND = 8

! depth at which we start to honor the basement interface
  double precision, parameter :: Z_THRESHOLD_HONOR_BASEMENT = -4700.d0

! flag to print the details of source location
  logical, parameter :: SHOW_DETAILS_LOCATE_SOURCE = .false.

! maximum length of station and network name for receivers
  integer, parameter :: MAX_LENGTH_STATION_NAME = 32
  integer, parameter :: MAX_LENGTH_NETWORK_NAME = 8

! number of sources to be gathered by MPI_Gather
  integer, parameter :: NGATHER_SOURCES = 10000

! we mimic a triangle of half duration equal to half_duration_triangle
! using a Gaussian having a very close shape, as explained in Figure 4.2
! of the manual. This source decay rate to mimic an equivalent triangle
! was found by trial and error
  double precision, parameter :: SOURCE_DECAY_MIMIC_TRIANGLE = 1.628d0

! ---------------------------------------------------------------------------------------
! LQY -- Following 3 variables stays here temporarily,
!        we need to move them to Par_file at a proper time
! ---------------------------------------------------------------------------------------
! save moho mesh and compute Moho boundary kernels
  logical, parameter :: SAVE_MOHO_MESH = .false.

! number of steps to save the state variables in the forward simulation,
! to be used in the backward reconstruction in the presence of attenuation
  integer, parameter :: NSTEP_Q_SAVE = 200

! the scratch disk to save the state variables saved in the forward
! simulation, this can be a global scratch disk in case you run out of
! space on the local scratch disk
  character(len=150), parameter :: LOCAL_PATH_Q = '/ibrixfs1/scratch/lqy/DATABASES_MPI_Q/'

!------------------------------------------------------
! nlegoff -- Variables that should be read/computed elsewhere.
!            Temporarily declared here.
!------------------------------------------------------
! whether or not an external mesh is used (provided by CUBIT for example)
  logical, parameter :: USE_EXTERNAL_MESH = .true.

! no lagrange interpolation on seismograms (we take the value on one NGLL point)
  logical, parameter :: FASTER_RECEIVERS_POINTS_ONLY = .true.
  logical, parameter :: FASTER_SOURCES_POINTS_ONLY = .true.

! the receivers can be located inside the model
  logical, parameter :: RECVS_CAN_BE_BURIED_EXT_MESH = .false.
  logical, parameter :: SOURCES_CAN_BE_BURIED_EXT_MESH = .false.

! the seismograms are normal to surface
  logical, parameter :: EXT_MESH_RECV_NORMAL = .false.

!
  logical, parameter :: EXTERNAL_MESH_MOVIE_SURFACE = .false.
  logical, parameter :: EXTERNAL_MESH_CREATE_SHAKEMAP = .true.

! deltat
  double precision, parameter :: DT_ext_mesh = 0.001d0

! NSTEP
  integer, parameter :: NSTEP_ext_mesh = 100

! number of nodes per element as provided by the external mesh
  integer, parameter :: ESIZE = 8

! geometry tolerance parameter to calculate number of independent grid points
! sensitive to actual size of model, assumes reference sphere of radius 1
! this is an absolute value for normalized coordinates in the Earth
  double precision, parameter :: SMALLVAL_TOL = 1.d-10

!------------------------------------------------------
!----------- do not modify anything below -------------
!------------------------------------------------------

! on some processors (e.g. Pentiums) it is necessary to suppress underflows
! by using a small initial field instead of zero
  logical, parameter :: FIX_UNDERFLOW_PROBLEM = .true.

! some useful constants
  double precision, parameter :: PI = 3.141592653589793d0
  double precision, parameter :: TWO_PI = 2.d0 * PI

! 3-D simulation
  integer, parameter :: NDIM = 3

! dimension of the boundaries of the slices
  integer, parameter :: NDIM2D = 2

! number of nodes for 2D and 3D shape functions for hexahedra
! we use 8-node mesh bricks, which are more stable than 27-node elements
  integer, parameter :: NGNOD = 8, NGNOD2D = 4

! a few useful constants
  double precision, parameter :: ZERO = 0.d0,ONE = 1.d0,TWO = 2.d0,HALF = 0.5d0

  real(kind=CUSTOM_REAL), parameter :: &
    ONE_THIRD   = 1._CUSTOM_REAL/3._CUSTOM_REAL, &
    FOUR_THIRDS = 4._CUSTOM_REAL/3._CUSTOM_REAL

! very large and very small values
  double precision, parameter :: HUGEVAL = 1.d+30,TINYVAL = 1.d-9

! very large real value declared independently of the machine
  real(kind=CUSTOM_REAL), parameter :: HUGEVAL_SNGL = 1.e+30_CUSTOM_REAL

! very large integer value
  integer, parameter :: HUGEINT = 100000000

! define flag for elements
  integer, parameter :: IFLAG_ONE_LAYER_TOPOGRAPHY = 1
  integer, parameter :: IFLAG_BASEMENT_TOPO = 2
  integer, parameter :: IFLAG_16km_BASEMENT = 3
  integer, parameter :: IFLAG_MOHO_16km = 4
  integer, parameter :: IFLAG_HALFSPACE_MOHO = 5

! define flag for regions for attenuation
  integer, parameter :: NUM_REGIONS_ATTENUATION = 13

  integer, parameter :: IATTENUATION_SEDIMENTS_40 = 1
  integer, parameter :: IATTENUATION_SEDIMENTS_50 = 2
  integer, parameter :: IATTENUATION_SEDIMENTS_60 = 3
  integer, parameter :: IATTENUATION_SEDIMENTS_70 = 4
  integer, parameter :: IATTENUATION_SEDIMENTS_80 = 5
  integer, parameter :: IATTENUATION_SEDIMENTS_90 = 6
  integer, parameter :: IATTENUATION_SEDIMENTS_100 = 7
  integer, parameter :: IATTENUATION_SEDIMENTS_110 = 8
  integer, parameter :: IATTENUATION_SEDIMENTS_120 = 9
  integer, parameter :: IATTENUATION_SEDIMENTS_130 = 10
  integer, parameter :: IATTENUATION_SEDIMENTS_140 = 11
  integer, parameter :: IATTENUATION_SEDIMENTS_150 = 12
  integer, parameter :: IATTENUATION_BEDROCK = 13

! Olsen's constant for Q_mu = constant * v_s attenuation rule
  real, parameter :: OLSEN_ATTENUATION_RATIO = 0.05

! number of standard linear solids in parallel for attenuation
  integer, parameter :: N_SLS = 3

! flag for the four edges of each slice and for the bottom edge
  integer, parameter :: XI_MIN = 1
  integer, parameter :: XI_MAX = 2
  integer, parameter :: ETA_MIN = 3
  integer, parameter :: ETA_MAX = 4
  integer, parameter :: BOTTOM = 5

! number of points per surface element
  integer, parameter :: NGLLSQUARE = NGLLX * NGLLY

! number of points per spectral element
  integer, parameter :: NGLLCUBE = NGLLX * NGLLY * NGLLZ

! for vectorization of loops
  integer, parameter :: NGLLSQUARE_NDIM = NGLLSQUARE * NDIM
  integer, parameter :: NGLLCUBE_NDIM = NGLLCUBE * NDIM

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

! number of lines per source in CMTSOLUTION file
  integer, parameter :: NLINES_PER_CMTSOLUTION_SOURCE = 13

! number of iterations to solve the system for xi and eta
  integer, parameter :: NUM_ITER = 4

! size of topography and bathymetry file for Southern California
  integer, parameter :: NX_TOPO_SOCAL = 1401,NY_TOPO_SOCAL = 1001
  double precision, parameter :: ORIG_LAT_TOPO_SOCAL = 32.d0
  double precision, parameter :: ORIG_LONG_TOPO_SOCAL = -121.d0
  double precision, parameter :: DEGREES_PER_CELL_TOPO_SOCAL = 5.d0 / 1000.d0
  character(len=100), parameter :: TOPO_FILE_SOCAL = 'DATA/la_topography/topo_bathy_final.dat'

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

!  Salton Sea Gocad voxet
  integer, parameter :: GOCAD_ST_NU = 638, GOCAD_ST_NV = 219, GOCAD_ST_NW = 76
  double precision, parameter :: GOCAD_ST_O_X = 720844.0, GOCAD_ST_O_Y = 3401799.250, &
    GOCAD_ST_O_Z =      -6354.334
  double precision, parameter :: GOCAD_ST_U_X = -209197.89, GOCAD_ST_U_Y =  320741.71
  double precision, parameter :: GOCAD_ST_V_X = 109670.74, GOCAD_ST_V_Y = 71530.72
  double precision, parameter :: GOCAD_ST_W_Z =  7666.334
  double precision, parameter :: GOCAD_ST_NO_DATA_VALUE = -99999

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
! DEPTH_MOHO_SOCAL = -35 km was based on Dreger and Helmberger (1990)
! and is (July 2007) the preferred Moho depth for Dreger.
! The depth of 32 km is used in the standard processing (Wald et al., 1995)
! of SoCal events and is the value in the original Kanamori-Hadley (1975) model.
  double precision, parameter :: DEPTH_5p5km_SOCAL = -5500.d0
  double precision, parameter :: DEPTH_16km_SOCAL = -16000.d0
  double precision, parameter :: DEPTH_MOHO_SOCAL = -32000.d0

! reference surface of the model before adding topography
  double precision, parameter :: Z_SURFACE = 0.d0

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

