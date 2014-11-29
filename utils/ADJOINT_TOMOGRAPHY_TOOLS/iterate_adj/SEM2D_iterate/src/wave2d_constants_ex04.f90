module wave2d_constants

  ! This file is copied into the output directory when wave2d.f90 is run
  ! There are two basic options: membrane (surface) waves and body waves.
  ! Several lines need to be commented/uncommented to go from one to the other.

  ! index of the reference directory for the simulation output
  integer, parameter :: IRUNZ = 400

  !========================
  ! GRID, TIME-STEP, AND SOURCE PARAMETERS

  ! NFRAME : number of frames to save: membrane, 10; body, 20
  ! NSAVE  : timestep increment to save the wavefield: membrane, 400; body, 400
  ! NSTEP  : number of timesteps
  integer, parameter :: NFRAME = 10     ! membrane
  !integer, parameter :: NFRAME = 11     ! body
  integer, parameter :: NSAVE  = 400   ! 200,400
  integer, parameter :: NSTEP  = NFRAME*NSAVE

  ! time step in seconds
  !double precision, parameter :: DT = 2.0d-2 ! body waves
  double precision, parameter :: DT = 6.0d-2 ! membrane surface waves

  ! temporal properties of source (source time function)
  integer, parameter :: ISRC_TIME = 1                   ! type (1)
  double precision, parameter :: hdur = 10.0            ! HALF-duration (s), membrane
  !double precision, parameter :: hdur = 2.0            ! HALF-duration (s), body
  double precision, parameter :: tshift = 2.0*DT*dble(NSAVE) ! time shift (s)
  !double precision, parameter :: tshift = 8.0*hdur
  logical, parameter :: SRC_TAPER  = .true.             ! taper the endpoints of the time series

  ! normalization factor of point source force
  double precision, parameter :: FNORM = 1.0d10

  ! forward wavefield source-time function
  ! (For membrane waves, FOR_Y = 1 is all that matters.)
  integer, parameter :: FOR_X = 1       ! boolean
  integer, parameter :: FOR_Y = 1       ! boolean
  integer, parameter :: FOR_Z = 1       ! boolean

  ! adjoint wavefield source-time function
  ! (For membrane waves, REV_Y = 1 is all that matters.)
  !  integer, parameter :: IPSV  = 0       ! boolean: plot PSV or SH kernels
  integer, parameter :: REV_X = 1       ! boolean
  integer, parameter :: REV_Y = 1       ! boolean
  integer, parameter :: REV_Z = 1       ! boolean

  ! spatial properties of sources
  ! BODY WAVE OPTIONS : (6) a point source EVENT (GJI paper)
  ! SURFACE WAVE OPTIONS : (1) point source, (2) finite segment, (3) CA shelf boundary
  !                        (4) CA coast, (5) finite circle, (6) a point source EVENT (GJI paper)
  integer, parameter :: ISRC_SPACE = 6 ! see wave2d.f90

  ! spatial properties of receivers
  ! IREC_SPACE
  !   (1) individual station(s)
  !   (2) SoCal (used for GJI paper)
  !   (3) regular mesh on land
  !   (4) regular mesh
  ! NMESH_REC : determines the number of receivers in a regular mesh (IREC_SPACE > 3)
  ! STATION_GRID_BUFFER : exclude stations within this distance from edge of grid
  ! STATION_COAST_BUFFER : exclude stations within this distance from edge of coast
  integer, parameter          :: IREC_SPACE = 2 ! see wave2d.f90
  integer, parameter          :: NMESH_REC = 17
  double precision, parameter :: SOURCE_GRID_BUFFER   =  4.0d3  ! m
  double precision, parameter :: STATION_GRID_BUFFER  = 15.0d3  ! m
  double precision, parameter :: STATION_COAST_BUFFER =  0.0d3  ! m

  ! lower right corner for membrane surface waves plotting grid
  double precision, parameter :: LAT_MIN = 32.0d0
  double precision, parameter :: LON_MIN = -120.0d0
  integer, parameter :: UTM_PROJECTION_ZONE = 11     ! southern California

  ! mesh specifications: membrane surface waves
  double precision, parameter :: LENGTH = 480.0d3 ! m
  double precision, parameter :: HEIGHT = 480.0d3 ! m
  double precision, parameter :: AREA = LENGTH*HEIGHT  
  integer, parameter :: NEX = 40
  integer, parameter :: NEZ = 40
!!$
!!$! mesh specifications: body waves
!!$  double precision, parameter :: LENGTH = 200.0d3 ! m      ! 400 for 1D body waves
!!$  double precision, parameter :: HEIGHT = 80.0d3 ! m
!!$  integer, parameter :: NEX = 80   ! 160
!!$  integer, parameter :: NEZ = 32   ! 32

  !========================
  ! MODEL SPECIFICATIONS (REFERENCE MODEL AND TARGET MODEL)

  ! model perturbations for HOMOGENEOUS model (or perturbation)
  ! scaling from beta to alpha
  ! value is from Masters et al. (2000), "The relative behavior of shear velocity..."
  double precision, parameter :: R_BETA_OVER_ALPHA = 1.3d0
  double precision, parameter :: PBETA  = 10.0d0
  double precision, parameter :: PALPHA = PBETA / R_BETA_OVER_ALPHA
  double precision, parameter :: PRHO   = 0.0d0

  ! reference model and target model choice
  integer, parameter :: IMODEL_SYN = 0
  integer, parameter :: IMODEL_DAT = 2
  !-----------------------------------------------------------------------------------------
  !                               0           1                2              3          
  !-----------------------------------------------------------------------------------------
  ! ISURFACE=1, IMODEL_SYN : homo                         checker          het
  ! ISURFACE=1, IMODEL_DAT : homo pert                    checker pert     het
  !-----------------------------------------------------------------------------------------
  ! ISURFACE=0, IMODEL_SYN : homo         1D model        checker          NA       
  ! ISURFACE=0, IMODEL_DAT : homo pert    1D model        checker pert     NA  
  !-----------------------------------------------------------------------------------------

  ! smooth the structure model or kernels
  ! For the CG algorithm, ISMOOTH_EVENT_KERNEL --> use smoothed event kernels
  integer, parameter :: ISMOOTH_EVENT_KERNEL  = 0
  integer, parameter :: ISMOOTH_MISFIT_KERNEL = 1
  integer, parameter :: ISMOOTH_INITIAL_MODEL = 0
  integer, parameter :: ISMOOTH_MODEL_UPDATE  = 0

  ! scalelength of smoothing Gaussian
  ! GJI paper : 30,60,90
  ! body waves : 10
  double precision, parameter :: SIGMA_SMOOTH_KERNEL = 2.121320d4
  double precision, parameter :: SIGMA_SMOOTH_MODEL  = 1.060660d4
  double precision, parameter :: GAMMA_SMOOTH_KERNEL = sqrt(8.0)*SIGMA_SMOOTH_KERNEL 
  double precision, parameter :: GAMMA_SMOOTH_MODEL  = sqrt(8.0)*SIGMA_SMOOTH_MODEL

  ! parameters for smoothing
  logical, parameter :: HIGH_RES_SMOOTHING = .true.  ! smooth at high resolution
  logical, parameter :: EXAMPLE_GAUSSIAN = .false.   ! plot an example Gaussian

  !========================

  ! boolean parameters
  ! IKER: (0) waveform
  !       (1) traveltime, cross-correlation, misfit
  !       (2) amplitude, cross-correlation, misfit
  !       (3) traveltime, multitaper
  !       (4) amplitude, multitaper
  !       (5) traveltime, cross-correlation, sampling -- NOT AN OPTION
  !       (6) amplitude, cross-correlation, sampling -- NOT AN OPTION
  integer, parameter :: IKER = 1
  integer, parameter :: IAMP_VEL = 0  ! measure amplitudes between velocity traces

  ! KEY: USE ONE LINE OR THE OTHER
  integer, parameter :: ISURFACE = 1, NCOMP = 1, NABSORB = 4   ! surface waves
  !  integer, parameter :: ISURFACE = 0, NCOMP = 3, NABSORB = 3   ! body waves

  ! parameters controlling what to write to file
  ! NOTE: for the tomography simulations, ALL of these can be .false.

  logical, parameter :: WRITE_STF_F           = .false.
  logical, parameter :: WRITE_SEISMO_F        = .false.
  !  logical, parameter :: WRITE_SPECTRA_F       = .false.
  !  logical, parameter :: WRITE_SPECTRAL_MAP_F  = .false.

  logical, parameter :: WRITE_SEISMO_RECONSTRUCT = .false.  ! multitaper
  !  logical, parameter :: WRITE_MTM_FILES = .false.           ! multitaper

  logical, parameter :: WRITE_STF_A           = .false.
  logical, parameter :: WRITE_SEISMO_A        = .false.     ! source inversions
  !  logical, parameter :: WRITE_SPECTRA_A       = .false. 
  !  logical, parameter :: WRITE_SPECTRAL_MAP_A  = .false.

  logical, parameter :: WRITE_KERNELS = .true.             ! write all nine kernels
  logical, parameter :: WRITE_KERNEL_SNAPSHOTS = .false.    ! kernel snapshots
  logical, parameter :: WRITE_WAVFIELD_SNAPSHOTS = .false.  ! wavefield snapshots

  !--------------------------------------
  ! INVERSION PARAMETERS

  ! whether you want to compute kernels or simply the misfit function
  logical, parameter :: COMPUTE_KERNELS = .true.

  ! whether to read in models generated from outside wave2d.f90
  logical, parameter :: READ_IN = .false.
  logical, parameter :: READ_SINGLE = .false.  ! read single or multiple files

  ! stopping criteria
  ! NITERATION  : number of iterations
  ! VAR_RED_MIN : minimum variance reduction
  ! SIGMA_FAC   : stop if a model value exceeds SIGMA_FAC * sigma_m 
  ! CONV_STOP   : stop when the misfit value is this fraction of the INITIAL misfit value
  integer, parameter :: NITERATION = 16
  double precision, parameter :: VAR_RED_MIN = 0.05d0
  !double precision, parameter :: SIGMA_FAC = 2.0d0 
  !double precision, parameter :: CONV_STOP = 1.0d-4

  ! Gaussian errors containted in input file
  ! see wave2d_sigmas.m and INPUT/sigma_0p1_pert.dat
  double precision, parameter :: SIGMA_DT = 0.10d0
  double precision, parameter :: SIGMA_DLNA = 1.0d0
  double precision, parameter :: SIGMA_WAVEFORM = 1.0d0
  logical, parameter :: ADD_DATA_ERRORS = .true.

  ! order of interpolating polynomial in conjugate gradient algorithm
  ! using POLY_ORDER = 3 required computing the gradient of the test model
  integer, parameter :: POLY_ORDER = 2      ! 2 (preferred) or 3

  ! use inversion and scaling used for the GJI paper
  !integer, parameter :: GJI_PAPER = 1

  ! what to perturb, what to invert
  ! (For the inverse tests, we only allow perturbations in beta.)
  integer, parameter :: PERT_STRUCT_BETA = 1
  integer, parameter :: PERT_SOURCE_T = 1
  integer, parameter :: PERT_SOURCE_X = 1

  integer, parameter ::  INV_STRUCT_BETA = 1
  integer, parameter ::  INV_SOURCE_T = 1
  integer, parameter ::  INV_SOURCE_X = 1

  ! whether to include the model norm term in the misfit function, which acts like damping
  logical, parameter :: INCLUDE_MODEL_NORM = .true.  

  ! log file showing source loactions
  logical, parameter :: ISOURCE_LOG = .true.

  ! DO NOT CHANGE THESE
  !  integer, parameter :: NVAR_STRUCT = 2   ! alpha, beta (seismic velocities only)
  integer, parameter :: NVAR_STRUCT = 1   ! beta only (24-Jan-2010)
  integer, parameter :: NVAR_SOURCE = 3   ! x position, z position, origin time
  integer, parameter :: NVAR = NVAR_STRUCT + NVAR_SOURCE

  ! parameterization for structure in the inversion
  ! 1 : kappa-mu-rho
  ! 2 : alpha-beta-rho
  ! 3 : c-beta-rho
  integer, parameter :: STRUCTURE_PARAMETER_TYPE = 2

  ! homogeneous background model (S.I. units)
  double precision, parameter :: DENSITY           = 2.60d3 ! kg/m^3
  double precision, parameter :: INCOMPRESSIBILITY = 4.50d10 ! Pa
  double precision, parameter :: RIGIDITY          = 3.185d10 ! Pa
  !  double precision, parameter :: INCOMPRESSIBILITY = 5.20d10 ! Pa
  !  double precision, parameter :: RIGIDITY          = 2.66d10 ! Pa

  ! compute additional parameters
  double precision, parameter :: PWAVESPEED = &
       sqrt( (INCOMPRESSIBILITY + (4.0d0/3.0d0)*RIGIDITY)/DENSITY )
  double precision, parameter :: SWAVESPEED = sqrt( RIGIDITY/DENSITY )
  double precision, parameter :: BWAVESPEED = sqrt( INCOMPRESSIBILITY/DENSITY )

  !---------------------------------------------------------------
  ! measurement windows

  double precision, parameter :: HWIN1 = 2.5*hdur   ! half-window width for pulse (2.5)
  double precision, parameter :: HWIN2 = 5.0*hdur   ! half-window with for MTM measurement (zero-padding)

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

  !========================
  ! GRID AND GLL POINTS

  integer, parameter :: NELE = MAX(NEX,NEZ)
  integer, parameter :: NSPEC = NEX*NEZ

  ! number of GLL points (polynomial degree plus one)
  integer, parameter :: NGLLX = 5
  integer, parameter :: NGLLZ = 5
  integer, parameter :: NGLL = MAX(NGLLX,NGLLZ)

  ! number of global points
  integer, parameter :: NGLLSQUARE = NGLLX * NGLLZ                       ! element GLL points
  integer, parameter :: NGLOB = ((NGLLX-1)*NEX + 1)*((NGLLZ-1)*NEZ +1)   ! global GLL points
  integer, parameter :: NLOCAL = NGLLX * NGLLZ * NSPEC                   ! local GLL points
  integer, parameter :: NSPEC_CORNER = (NEX+1) * (NEZ+1)                 ! element corner points

  ! number of nodes for 2D and 3D shape functions for hexahedra
  ! we use 8-node mesh bricks, which are more stable than 27-node elements
  integer, parameter :: NGNOD = 8, NGNOD2D = 4

  ! number of iterations to solve the system for xi and eta
  integer, parameter :: NUM_ITER = 1

  ! very large and very small values
  double precision, parameter :: HUGEVAL = 1.0d30, TINYVAL = 1.0d-9

  ! for the Gauss-Lobatto-Legendre points and weights
  double precision, parameter :: GAUSSALPHA = 0.0d0,GAUSSBETA = 0.0d0

  !
  ! CONSTANTS
  !
  double precision, parameter :: PI = 3.141592653589793d0
  double precision, parameter :: FOUR_THIRDS = 4.0d0 / 3.0d0
  double precision, parameter :: ONE_THIRD = 1.0d0 / 3.0d0
  double precision, parameter :: ONEOVERTWO = 0.5d0
  double precision, parameter :: DEG = 180.0d0/PI
  !  double precision, parameter :: EPS = 1.0d-35

!!$! parameter for FFTW
!!$  integer, parameter :: NOUT = NSTEP/2 + 1
!!$
!!$! bounds for bandpass filter
!!$  double precision, parameter :: hwid = 3.0  ! HALF-width of window (s)
!!$  double precision, parameter :: tmin = 2.0*hdur-hwid
!!$  double precision, parameter :: tmax = 2.0*hdur+hwid
!!$  double precision, parameter :: fmin = 1.0/tmax, fmax = 1.0/tmin
!!$  double precision, parameter :: trbdndw = 0.3, a = 30.0
!!$  integer, parameter :: passes = 2, iord = 4
!!$
!!$! 
!!$! MULTI-TAPER PARAMETERS
!!$!
!!$! Ying Zhou: The fit between the recovered data and the data can be improved
!!$! by either increasing the window width (HWIN above) or by decreasing NPI.
!!$! In her experience, NPI = 2.5 is good for noisy data.
!!$! For synthetic data, we can use a lower NPI.
!!$! number of tapers should be fixed as twice NPI -- see Latex notes
!!$  double precision, parameter :: WTR = 0.02
!!$  double precision, parameter :: NPI = 2.5
!!$  integer, parameter :: NTAPER = int(2.0*NPI)
!!$
!!$! FFT parameters
!!$! FOR SOME REASON, I CANNOT SET NDIM = npt when lnpt <= 11 -- WHY NOT?
!!$  integer, parameter :: lnpt = 13, npt = 2**lnpt, NDIM = 20000
!!$  !integer, parameter :: lnpt = 14, npt = 2**lnpt, NDIM = npt
!!$
!!$  !double precision, parameter :: ZZIGN = -1.0   ! sign convention -- WATCH OUT
!!$  double precision, parameter :: FORWARD_FFT = 1.0
!!$  double precision, parameter :: REVERSE_FFT = -1.0

  integer, parameter :: NDIM = 20000

end module wave2d_constants
