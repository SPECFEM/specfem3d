!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
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


module constants

  implicit none

  include "constants.h"

  ! proc number for MPI process
  integer :: myrank

  ! a negative initial value is a convention that indicates that groups (i.e. sub-communicators, one per run) are off by default
  integer :: mygroup = -1

  ! MPI status size (will be set in init_mpi()
  integer :: my_status_size   = 1
  integer :: my_status_source = 1
  integer :: my_status_tag    = 1

  ! create a copy of the original output file path, to which we may add a "run0001/", "run0002/", "run0003/" prefix later
  ! if NUMBER_OF_SIMULTANEOUS_RUNS > 1
  character(len=MAX_STRING_LEN) :: OUTPUT_FILES = OUTPUT_FILES_BASE

  ! if doing simultaneous runs for the same mesh and model, see who should read the mesh and the model and broadcast it to others
  ! we put a default value here
  logical :: I_should_read_the_database = .true.

end module constants

!
!-------------------------------------------------------------------------------------------------
!

  module shared_input_parameters

! holds input parameters given in DATA/Par_file

  use constants, only: MAX_STRING_LEN,STANDARD_GRAVITY,NDIM, &
                       DEFAULT_ROTATION_OMEGA,DEFAULT_ROTATION_ORIGIN

  implicit none

  ! parameters read from parameter file
  integer :: NPROC = 0

  ! simulation parameters
  integer :: SIMULATION_TYPE = 1
  integer :: NOISE_TOMOGRAPHY = 0
  logical :: SAVE_FORWARD = .false.
  logical :: INVERSE_FWI_FULL_PROBLEM = .false.

  integer :: UTM_PROJECTION_ZONE = 0
  logical :: SUPPRESS_UTM_PROJECTION = .true.

  logical :: UNDO_ATTENUATION_AND_OR_PML = .false.
  integer :: NT_DUMP_ATTENUATION = 0

  ! number of time steps
  integer :: NSTEP = 0
  double precision :: DT = 0.d0

  ! number of time step for external source time function
  integer :: NSTEP_STF = 0

  ! Local Time Stepping (LTS)
  logical :: LTS_MODE = .false.

  ! partitioning scheme
  integer :: PARTITIONING_TYPE = 1

  ! LDD Runge-Kutta time scheme
  logical :: USE_LDDRK = .false.
  logical :: INCREASE_CFL_FOR_LDDRK = .false.
  double precision :: RATIO_BY_WHICH_TO_INCREASE_IT = 1.4d0

  ! mesh
  integer :: NGNOD = 8

  character(len=MAX_STRING_LEN) :: MODEL = 'default'
  character(len=MAX_STRING_LEN) :: SEP_MODEL_DIRECTORY = './DATA/sep_model'

  ! physical parameters
  logical :: APPROXIMATE_OCEAN_LOAD = .false.
  logical :: TOPOGRAPHY = .false.
  logical :: ATTENUATION = .false.
  logical :: ANISOTROPY = .false.
  logical :: GRAVITY = .false.
  logical :: ROTATION = .false.

  ! (optional) rotation angular velocity and origin of center of rotation (relative to mesh coordinates)
  double precision, dimension(NDIM) :: ROTATION_OMEGA = DEFAULT_ROTATION_OMEGA
  double precision, dimension(NDIM) :: ROTATION_ORIGIN = DEFAULT_ROTATION_ORIGIN

  ! tomography model file path
  character(len=MAX_STRING_LEN) :: TOMOGRAPHY_PATH = './DATA/tomo_files'

  ! attenuation
  ! reference frequency of seismic model
  double precision :: ATTENUATION_f0_REFERENCE = 1.d0
  ! Olsen attenuation (scaling from Vs)
  logical :: USE_OLSEN_ATTENUATION = .false.
  double precision :: OLSEN_ATTENUATION_RATIO = 0.05d0
  ! automatic frequency band selection
  logical :: COMPUTE_FREQ_BAND_AUTOMATIC = .true.
  ! attenuation period range over which we try to mimic a constant Q factor
  double precision :: MIN_ATTENUATION_PERIOD = 999999.d0
  double precision :: MAX_ATTENUATION_PERIOD = 999999.d0
  ! logarithmic center frequency (center of attenuation band)
  double precision :: ATT_F_C_SOURCE = 1.d0

  ! absorbing boundaries
  ! PML
  logical :: PML_CONDITIONS = .false.
  logical :: PML_INSTEAD_OF_FREE_SURFACE = .false.
  double precision :: f0_FOR_PML = 1.d0
  ! Stacey
  logical :: STACEY_ABSORBING_CONDITIONS = .false.
  logical :: STACEY_INSTEAD_OF_FREE_SURFACE = .false.
  ! To use a bottom free surface instead of absorbing Stacey or PML condition
  logical :: BOTTOM_FREE_SURFACE = .false.

  ! sources and receivers Z coordinates given directly instead of with depth
  logical :: USE_SOURCES_RECEIVERS_Z = .false.

  ! for simultaneous runs from the same batch job
  integer :: NUMBER_OF_SIMULTANEOUS_RUNS = 1
  logical :: BROADCAST_SAME_MESH_AND_MODEL = .false.

  ! movies
  logical :: CREATE_SHAKEMAP = .false.
  logical :: MOVIE_SURFACE = .false.
  logical :: MOVIE_VOLUME = .false.
  logical :: SAVE_DISPLACEMENT = .false.
  logical :: USE_HIGHRES_FOR_MOVIES = .false.
  logical :: MOVIE_VOLUME_STRESS = .false.
  integer :: MOVIE_TYPE = 1

  integer :: NTSTEP_BETWEEN_FRAMES = 100
  double precision :: HDUR_MOVIE = 0.d0

  ! mesh
  logical :: SAVE_MESH_FILES = .false.
  character(len=MAX_STRING_LEN) :: LOCAL_PATH = './DATABASES_MPI'

  ! seismograms
  integer :: NTSTEP_BETWEEN_OUTPUT_INFO = 500
  integer :: NTSTEP_BETWEEN_OUTPUT_SEISMOS = 100000
  integer :: NTSTEP_BETWEEN_READ_ADJSRC = 0
  integer :: NTSTEP_BETWEEN_OUTPUT_SAMPLE = 1 ! subsamp_seismos is deprecated and renamed to NTSTEP_BETWEEN_OUTPUT_SAMPLE
  logical :: SAVE_SEISMOGRAMS_DISPLACEMENT = .true.
  logical :: SAVE_SEISMOGRAMS_VELOCITY = .false.
  logical :: SAVE_SEISMOGRAMS_ACCELERATION = .false.
  logical :: SAVE_SEISMOGRAMS_PRESSURE = .false.
  logical :: SAVE_SEISMOGRAMS_STRAIN = .false.
  logical :: SAVE_SEISMOGRAMS_IN_ADJOINT_RUN = .false.
  logical :: WRITE_SEISMOGRAMS_BY_MAIN = .true.
  logical :: SAVE_ALL_SEISMOS_IN_ONE_FILE = .false.
  logical :: USE_BINARY_FOR_SEISMOGRAMS = .false.
  logical :: SU_FORMAT = .false.
  logical :: ASDF_FORMAT = .false.
  logical :: READ_ADJSRC_ASDF = .false.

  ! sources
  logical :: USE_FORCE_POINT_SOURCE = .false.
  logical :: USE_RICKER_TIME_FUNCTION = .false.
  logical :: PRINT_SOURCE_TIME_FUNCTION = .false.
  character(len=MAX_STRING_LEN) :: CMT_CONVENTION_FORMAT = 'USE'  ! Up-South-East (Harvard) convention

  ! cmt + point force simulation
  logical :: USE_CMT_AND_FORCE_SOURCE = .false.
  logical :: USE_BINARY_SOURCE_FILE = .false.

  ! external source time function
  logical :: USE_EXTERNAL_SOURCE_FILE = .false.

  logical :: USE_TRICK_FOR_BETTER_PRESSURE = .false.
  logical :: USE_SOURCE_ENCODING = .false.
  logical :: OUTPUT_ENERGY = .false.
  logical :: ANISOTROPIC_KL = .false.
  logical :: SAVE_TRANSVERSE_KL = .false.
  logical :: APPROXIMATE_HESS_KL = .false.
  logical :: SAVE_MOHO_MESH = .false.
  logical :: ANISOTROPIC_VELOCITY_KL = .false.
  integer :: NTSTEP_BETWEEN_OUTPUT_ENERGY = 10

  ! GPU simulations
  logical :: GPU_MODE = .false.

  ! adios file output
  logical :: ADIOS_ENABLED = .false.
  logical :: ADIOS_FOR_DATABASES = .false.
  logical :: ADIOS_FOR_MESH = .false.
  logical :: ADIOS_FOR_FORWARD_ARRAYS = .false.
  logical :: ADIOS_FOR_KERNELS = .false.
  logical :: ADIOS_FOR_UNDO_ATTENUATION = .false.

  ! HDF5 file i/o
  logical :: HDF5_ENABLED = .false.              ! for all databases i/o in hdf5
  logical :: HDF5_FOR_MOVIES = .false.           ! for movies (shakemap, surface movies, volume movies)

  ! HDF5 seismogram output
  logical :: HDF5_FORMAT  = .false.           ! for seismograms output in hdf5

  ! HDF5 IO server
  ! number of io dedicated nodes
  integer :: HDF5_IO_NODES = 0

  ! HDF5 IO writing mode (collective or independent)
  logical :: H5_COL = .true.

  ! flag for io-dedicated/compute node.
  logical :: IO_storage_task = .false.
  logical :: IO_compute_task = .true.

  ! external code coupling (DSM, AxiSEM)
  logical :: COUPLE_WITH_INJECTION_TECHNIQUE = .false.
  integer :: INJECTION_TECHNIQUE_TYPE = 0
  character(len=MAX_STRING_LEN) :: TRACTION_PATH = 'DATA/tractions'
  character(len=MAX_STRING_LEN) :: FKMODEL_FILE = 'DATA/FKMODEL'
  logical :: MESH_A_CHUNK_OF_THE_EARTH = .false.
  logical :: RECIPROCITY_AND_KH_INTEGRAL = .false.
  double precision :: INJECTION_START_TIME = -999999.d0

  ! prescribed wavefield discontinuity on an interface
  logical :: IS_WAVEFIELD_DISCONTINUITY = .false. ! if .true. then wavefield discontinuity is turned on (default is false)

  ! (optional) scattering perturbations
  logical :: ADD_SCATTERING_PERTURBATIONS = .false.
  double precision :: SCATTERING_STRENGTH = 0.d0
  double precision :: SCATTERING_CORRELATION = 1.d0
  character(len=MAX_STRING_LEN) :: SCATTERING_MATERIAL_IDS = ""

  ! Moon's Lunar Projections (LTM/LPS) instead of UTM
  logical :: USE_LUNAR_PROJECTIONS = .false.

  ! (optional) gravity min/max values
  logical :: USE_GRAVITY_MINMAX = .false.
  double precision :: GRAVITY_MINMAX_TOP = STANDARD_GRAVITY      ! default in m/s^2
  double precision :: GRAVITY_MINMAX_BOTTOM = STANDARD_GRAVITY

  end module shared_input_parameters

!
!-------------------------------------------------------------------------------------------------
!

  module shared_compute_parameters

  ! parameters to be computed based upon parameters above read from file
  use constants, only: CUSTOM_REAL

  implicit none

  ! number of sources given in CMTSOLUTION file
  integer :: NSOURCES = 0
  logical :: HAS_FINITE_FAULT_SOURCE = .false.

  !number of sources in CMTSOLUTION/FORCESOLUTION
  integer :: NSOURCES_CMT = 0
  integer :: NSOURCES_FORCE = 0

  ! anchor points
  integer :: NGNOD2D = 0

  ! model
  integer :: IMODEL = 0

  !! VM VM number of source for external source time function
  integer :: NSOURCES_STF = 0

  ! simulation type
  logical :: ACOUSTIC_SIMULATION = .false.
  logical :: ELASTIC_SIMULATION = .false.
  logical :: POROELASTIC_SIMULATION = .false.

  ! fault rupture simulation
  logical :: FAULT_SIMULATION = .false.

  ! free surface
  ! for elevation search: x/y coordinates of free surface element midpoints
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: free_surface_xy_midpoints
  real(kind=CUSTOM_REAL) :: free_surface_typical_size
  ! flag to calculate typical size only once
  logical :: free_surface_has_typical_size = .false.

  end module shared_compute_parameters

!
!-------------------------------------------------------------------------------------------------
!

  module shared_parameters

  use shared_input_parameters
  use shared_compute_parameters

  implicit none

  end module shared_parameters

