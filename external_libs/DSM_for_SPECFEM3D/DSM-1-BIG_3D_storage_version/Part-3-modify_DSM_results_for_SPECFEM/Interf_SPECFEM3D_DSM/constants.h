!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and University of Pau / CNRS / INRIA
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            April 2011
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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


! src/shared/constants.h.  Generated from constants.h.in by configure.

!
! solver in single or double precision depending on the machine (4 or 8 bytes)
!
! ALSO CHANGE FILE precision.h ACCORDINGLY
!
  integer, parameter :: SIZE_REAL = 4
  integer, parameter :: SIZE_DOUBLE = 8

! usually the size of integer and logical variables is the same as regular single-precision real variable
  integer, parameter :: SIZE_INTEGER = SIZE_REAL
  integer, parameter :: SIZE_LOGICAL = SIZE_REAL

! set to SIZE_REAL to run in single precision
! set to SIZE_DOUBLE to run in double precision (increases memory size by 2)
  integer, parameter :: CUSTOM_REAL = SIZE_REAL

!----------- parameters that can be changed by the user -----------

! set to .false.  if running on a Beowulf-type machine with local disks
! set to .true. if running on a shared-memory machine with common file system
! if running on a Beowulf, also modify name of nodes in filter_machine_file.f90
  logical, parameter :: LOCAL_PATH_IS_ALSO_GLOBAL = .true.

! apply heuristic rule to modify doubling regions to balance angles
  logical, parameter :: APPLY_HEURISTIC_RULE = .true.

! adds/superimposes velocity model values from 'model_external_values.f90'
  logical, parameter :: USE_MODEL_EXTERNAL_VALUES = .false.

! use inlined products of Deville et al. (2002) to speedup the calculations to compute internal forces
  logical, parameter :: USE_DEVILLE_PRODUCTS = .true.

! number of GLL points in each direction of an element (degree plus one)
  integer, parameter :: NGLLX = 5
  integer, parameter :: NGLLY = NGLLX
  integer, parameter :: NGLLZ = NGLLX

! number of points per surface element
  integer, parameter :: NGLLSQUARE = NGLLX * NGLLY

! number of points per spectral element
  integer, parameter :: NGLLCUBE = NGLLX * NGLLY * NGLLZ

! for optimized routines by Deville et al. (2002)
  integer, parameter :: m1 = NGLLX, m2 = NGLLX * NGLLY

! ouput format of seismograms, ASCII or binary
  logical, parameter :: SEISMOGRAMS_BINARY = .false.
! output format of seismograms, Seismic Unix (binary with 240-byte-headers)
  logical, parameter :: SU_FORMAT=.false.

! input, output and main MPI I/O files
  integer, parameter :: ISTANDARD_OUTPUT = 6
  integer, parameter :: IIN = 40,IOUT = 41
! uncomment this to write messages to a text file
  integer, parameter :: IMAIN = 42
! uncomment this to write messages to the screen
! integer, parameter :: IMAIN = ISTANDARD_OUTPUT
! I/O unit for source and receiver vtk file
  integer, parameter :: IOVTK = 98
! I/O unit for absorbing boundary snapshots
!  integer, parameter :: IOABS = 31
!  integer, parameter :: IOABS_AC = 32
! I/O unit for plotting source time function
  integer, Parameter :: IOSTF = 71
! I/O unit for interface file
  integer, parameter :: IIN_INTERFACES = 43
! I/O unit for noise simulations
  integer, parameter :: IIN_NOISE = 44,IOUT_NOISE = 45
! I/O unit for SU formatted seismograms
  integer, parameter :: IOUT_SU = 46
! I/O unit for SU formatted adjoint sources
  integer, parameter :: IIN_SU1 = 47, IIN_SU2 = 48, IIN_SU3 = 49

! ignore variable name field (junk) at the beginning of each input line
  logical, parameter :: IGNORE_JUNK = .true.,DONT_IGNORE_JUNK = .false.

! flag to print the details of source location
  logical, parameter :: SHOW_DETAILS_LOCATE_SOURCE = .false.

! maximum length of station and network name for receivers
  integer, parameter :: MAX_LENGTH_STATION_NAME = 32
  integer, parameter :: MAX_LENGTH_NETWORK_NAME = 8

! number of sources to be gathered by MPI_Gather
  integer, parameter :: NGATHER_SOURCES = 100

! we mimic a triangle of half duration equal to half_duration_triangle
! using a Gaussian having a very close shape, as explained in Figure 4.2
! of the manual. This source decay rate to mimic an equivalent triangle
! was found by trial and error
  double precision, parameter :: SOURCE_DECAY_MIMIC_TRIANGLE = 1.628d0

! decide if master process writes all the seismograms or if all processes do it in parallel
  logical, parameter :: WRITE_SEISMOGRAMS_BY_MASTER = .false.

! use directory OUTPUT_FILES/ for seismogram output
  logical,parameter :: USE_OUTPUT_FILES_PATH = .true.

! absorb top surface
! (defined in mesh as 'free_surface_file')
  logical,parameter :: ABSORB_FREE_SURFACE = .false.

! absorb boundaries using a PML region
! (EXPERIMENTAL feature)
! (only acoustic domains supported...)
! (user parameters can be specified in PML_init.f90)
  logical,parameter :: ABSORB_USE_PML = .false.

! paths for inputs and outputs files
  character(len=256), parameter :: IN_DATA_FILES_PATH = './DATA/'
  character(len=256), parameter :: MF_IN_DATA_FILES_PATH = './DATA/meshfem3D_files/'
  character(len=256), parameter :: OUTPUT_FILES_PATH = './OUTPUT_FILES/'


! ---------------------------------------------------------------------------------------
! LQY -- Following 3 variables stays here temporarily,
!        we need to move them to Par_file at a proper time
! ---------------------------------------------------------------------------------------
! save moho mesh and compute Moho boundary kernels
  logical, parameter :: SAVE_MOHO_MESH = .false.

! outputs approximate hessian for preconditioning
  logical, parameter :: APPROXIMATE_HESS_KL = .false.

! number of steps to save the state variables in the forward simulation,
! to be used in the backward reconstruction in the presence of attenuation
  integer, parameter :: NSTEP_Q_SAVE = 50 ! depending on stability of reconstruction, up to 200

! the scratch disk to save the state variables saved in the forward
! simulation, this can be a global scratch disk in case you run out of
! space on the local scratch disk, e.g. '/ibrixfs1/scratch/lqy/DATABASES_MPI_Q/'
  character(len=256), parameter :: LOCAL_PATH_Q = './OUTPUT_FILES/DATABASES_MPI/'

!------------------------------------------------------
! nlegoff -- Variables that should be read/computed elsewhere.
!            Temporarily declared here.
!------------------------------------------------------

! no lagrange interpolation on seismograms (we take the value on one NGLL point)
  logical, parameter :: FASTER_RECEIVERS_POINTS_ONLY = .false.

! use a force source located exactly at a grid point instead of a CMTSOLUTION source
! this can be useful e.g. for oil industry foothills simulations or asteroid simulations
! in which the source is a vertical force, normal force, impact etc.
  logical, parameter :: USE_FORCE_POINT_SOURCE = .false.
  double precision, parameter :: FACTOR_FORCE_SOURCE = 1.d15
  integer, parameter :: COMPONENT_FORCE_SOURCE = 3  ! takes direction in comp E/N/Z = 1/2/3

! use this t0 as earliest starting time rather than the automatically calculated one
! (must be positive and bigger than the automatically one to be effective;
!  simulation will start at t = - t0)
  double precision, parameter :: USER_T0 = 0.0d0

! the receivers can be located inside the model
  logical, parameter :: RECVS_CAN_BE_BURIED_EXT_MESH = .true.
  logical, parameter :: SOURCES_CAN_BE_BURIED_EXT_MESH = .true.

! sources and receivers Z coordinates given directly instead of with depth
  logical, parameter :: USE_SOURCES_RECVS_Z = .false.

! the seismograms are normal to surface
! Z record corresponds to the normal, while E and N are two tangent vectors
! that completes an orthonormal.
  logical, parameter :: EXT_MESH_RECV_NORMAL = .false.

! shakemaps and movies can not be generated during the same run. Mutually exclusive.
  logical, parameter :: EXTERNAL_MESH_MOVIE_SURFACE = .false.
  logical, parameter :: EXTERNAL_MESH_CREATE_SHAKEMAP = .false.

! plots VTK cross-section planes instead of model surface
! (EXPERIMENTAL feature)
! (requires EXTERNAL_MESH_MOVIE_SURFACE set to true)
  logical, parameter :: PLOT_CROSS_SECTIONS = .false.
  real(kind=CUSTOM_REAL),parameter :: CROSS_SECTION_X = 67000.0
  real(kind=CUSTOM_REAL),parameter :: CROSS_SECTION_Y = 65500.0
  real(kind=CUSTOM_REAL),parameter :: CROSS_SECTION_Z = -30000.0

! plots GIF cross-section image
! (EXPERIMENTAL feature)
! (cross-section plane parameters can be specified in create_color_image.f90)
  logical, parameter :: PNM_GIF_IMAGE = .false.

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

! for vectorization of loops
  integer, parameter :: NGLLSQUARE_NDIM = NGLLSQUARE * NDIM
  integer, parameter :: NGLLCUBE_NDIM = NGLLCUBE * NDIM

! number of nodes for 2D and 3D shape functions for hexahedra
! we use 8-node mesh bricks, which are more stable than 27-node elements
  integer, parameter :: NGNOD = 8, NGNOD2D = 4

! a few useful constants
  double precision, parameter :: ZERO = 0.d0,ONE = 1.d0,TWO = 2.d0 !,HALF = 0.5d0

  real(kind=CUSTOM_REAL), parameter :: &
    ONE_THIRD   = 1._CUSTOM_REAL/3._CUSTOM_REAL, &
    FOUR_THIRDS = 4._CUSTOM_REAL/3._CUSTOM_REAL

! very large and very small values
  double precision, parameter :: HUGEVAL = 1.d+30,TINYVAL = 1.d-9

! tiny real value declared independently of the machine
  real(kind=CUSTOM_REAL), parameter :: TINYVAL_SNGL = 1.e-25_CUSTOM_REAL

! Olsen's constant for Q_mu = constant * v_s attenuation rule
  real, parameter :: OLSEN_ATTENUATION_RATIO = 0.05

! number of standard linear solids in parallel for attenuation
  integer, parameter :: N_SLS = 3

! computation of standard linear solids
! ATTENUATION_COMP_RESOLUTION: Number of Digits after decimal
! ATTENUATION_COMP_MAXIMUM:    Maximum Q Value
  integer, parameter :: ATTENUATION_COMP_RESOLUTION = 1
  integer, parameter :: ATTENUATION_COMP_MAXIMUM    = 9000

! reference frequency for target velocity values in velocity model
! arbitrarily set to typical resolution of model (3 sec)
  double precision, parameter :: ATTENUATION_f0_REFERENCE = 0.3d0

! define flag for regions for anisotropy
  integer, parameter :: IANISOTROPY_MODEL1 = 1
  integer, parameter :: IANISOTROPY_MODEL2 = 2

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

! ! size of topography and bathymetry file for Piero Basini's model
!   integer, parameter :: NX_TOPO = 787, NY_TOPO = 793
!   double precision, parameter :: ORIG_LAT_TOPO = -102352.d0
!   double precision, parameter :: ORIG_LONG_TOPO = 729806.d0
! ! for Piero Basini's model this is the resolution in meters of the topo file
!   double precision, parameter :: DEGREES_PER_CELL_TOPO = 250.d0
!   character(len=256), parameter :: TOPO_FILE = 'DATA/piero_model/dem_EV_UTM_regular_250_reordered.dat'

! flag for projection from latitude/longitude to UTM, and back
  integer, parameter :: ILONGLAT2UTM = 0, IUTM2LONGLAT = 1

! minimum thickness in meters to include the effect of the oceans
! to avoid taking into account spurious oscillations in topography model
  double precision, parameter :: MINIMUM_THICKNESS_3D_OCEANS = 10.d0
! density of sea water
  real(kind=CUSTOM_REAL), parameter :: RHO_OCEANS = 1020.0

! material domain ids
  integer, parameter :: IDOMAIN_ACOUSTIC    = 1
  integer, parameter :: IDOMAIN_ELASTIC     = 2
  integer, parameter :: IDOMAIN_POROELASTIC = 3

! number of points in each AVS or OpenDX quadrangular cell for movies
  integer, parameter :: NGNOD2D_AVS_DX = 4
