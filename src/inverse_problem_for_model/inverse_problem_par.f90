!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
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


module inverse_problem_par

  !! IMPORT VARIABLES
  use specfem_par, only: CUSTOM_REAL, MAX_STRING_LEN, MAX_LENGTH_STATION_NAME, MAX_LENGTH_NETWORK_NAME

  implicit none

  !! ------------------------------ compilation config parameters ---------------------------------

  !! maximum line length allowed in input files
  integer,                       public, parameter  :: MAX_LEN_STRING=256
  !! log file for inversion
  integer,                       public, parameter  :: INVERSE_LOG_FILE=6666, OUTPUT_ITERATION_FILE=6667
  integer,                       public, parameter  :: OUTPUT_FWI_LOG=6668
  !! useful files id
  integer,                       public, parameter  :: IINN=667, IIDD=668
  !!! verbose debug mode
  logical,                       public, parameter  :: DEBUG_MODE=.true.
  !!! verbose mode for with outputs for checking the FWI
  logical,                       public, parameter  :: VERBOSE_MODE=.true.
  !!! write kernels on disk
  logical,                       public, parameter  :: SAVE_KERNEL=.false.
  !!! use fast code and undoing_attenuation for adjoints (under development ... does not work yet)
  logical,                       public, parameter  :: USE_UNDO_ATTENUATION_AND_OR_PML=.false.
  !!! projection on FD grid for outputs (in developmement if useful need to move elsewhere)
  logical,                       public             :: PROJ_ON_FD=.false.
  !!! test for some preconditionners (in developmement if useful need to move elsewhere)
  !!logical,                       public             :: USE_PRECOND_OIL_INDUSTRY=.false.
  !!! if needed use steepest descent instead of l-bfgs (if useful need to move elsewhere)
  logical,                       public             :: USE_GRADIENT_OPTIM=.false.
  !! use simplified station location instead of specfem subroutine which can be problematic with a
  !! big number of stations
  logical,                       public             :: USE_LIGHT_STATIONS=.true.
  ! ------------------------------  global parameters for fwi ---------------------------------------
  !! name for outputs files
  character(len=MAX_STRING_LEN), public             :: prname_specfem
  character(len=8),              public             :: prefix_to_path='./'
  character(len=MAX_STRING_LEN), public             :: type_input='exploration'

!################################################# STRUCTURES ###############################################

  ! PROJECTION IN FD GRID STRUCTURE
  type, public :: profd

     integer                                                                  :: nx, ny, nz
     real(kind=CUSTOM_REAL)                                                   :: hx, hy, hz
     real(kind=CUSTOM_REAL)                                                   :: ox, oy, oz

     integer                                                                  :: nb_fd_point
     integer,                               dimension(:),         allocatable :: ispec_selected
     integer,                               dimension(:,:),       allocatable :: index_on_fd_grid
     double precision,                      dimension(:,:),       allocatable :: hxi, heta, hgamma

  end type profd

  ! INVERSION PARAMETERS STRUCTURE
  type, public :: inver

     !! inputs files to read -------------------------------------------------------------------------------
     character(len= MAX_LEN_STRING)                                           :: input_acqui_file
     character(len= MAX_LEN_STRING)                                           :: input_inver_file

     !! managing parameters family -------------------------------------------------------------------------
     !! choice of family parameters  :
     !! ISO : rho vp vs
     !! VTI : rho vp vs ep gm de
     !! ... todo add more ...
     character(len=MAX_LEN_STRING)                                            :: parameter_family_name="ISO"
     character(len=MAX_LEN_STRING), dimension(50)                             :: param_inv_name
     character(len=MAX_LEN_STRING), dimension(50)                             :: param_ref_name
     integer                                                                  :: NfamilyPar = 3
     integer                                                                  :: NinvPar = 3
     integer, dimension(:), allocatable                                       :: Index_Invert
     !! this is not useful todo remove it. -------------------------------------------------------
     logical                                                                  :: use_log=.true.
     !! ---------------------------------------------------------------------------------------------
     integer                                                                  :: parameter_metric=2 ! see below :
     !! choice of metric in parameter :
     !! 0 : directly use the parameter P
     !! 1 : use P / Pref
     !! 2 : use log(P)
     !! 3 : use log(P/Pref)

     !! stopping criteria ----------------------------------------------------------------------------------
     integer                                                                  :: Niter = 100
     integer                                                                  :: Niter_wolfe = 10
     real(kind=CUSTOM_REAL)                                                   :: relat_grad = 1.e-3
     real(kind=CUSTOM_REAL)                                                   :: relat_cost = 1.e-1

     !! use filter band pass for inversion
     integer                                                                  :: Nifrq=1
     logical                                                                  :: use_band_pass_filter=.false.

     !!---- TUNING PARAMETERS FOR INVERSION -------------
     !! --- for L-BFGS
     integer                                                                  :: max_history_bfgs = 10

     !!--- for Wolfe line search
     real(kind=CUSTOM_REAL)                                                   :: m1 = 1.e-4, m2 = 0.9
     real(kind=CUSTOM_REAL)                                                   :: current_step_length = 1.

     !! weight for regularization on family parameters
     real(kind=CUSTOM_REAL), dimension(:), allocatable                        :: reg_family, damp_family

     !! max perturbation allowed for initial guest step length relative to gradient max value
     real(kind=CUSTOM_REAL)                                                   :: max_relative_pert = 0.05 ! (5%)

     !! we can import model from disk
     logical                                                                  :: input_sem_model = .false.
     logical                                                                  :: input_sem_prior = .false.
     logical                                                                  :: input_fd_model = .false.

     !! write FWI solution to disk (or not)
     logical                                                                  :: output_model = .true.

     !! what to do: FWI, or just forward modeling
     logical                                                                  :: only_forward = .false.

     !! for preconditioning of FWI
     logical                                                                  :: use_taper=.false.
     logical                                                                  :: shin_precond=.false.
     logical                                                                  :: z2_precond=.false.
     logical                                                                  :: z_precond=.false.
     logical                                                                  :: energy_precond=.false.

     !! flags to choose if the user wants to dump some wave fields
     logical                                                                  :: dump_model_at_each_iteration=.true.
     logical                                                                  :: dump_gradient_at_each_iteration=.true.
     logical                                                                  :: dump_descent_direction_at_each_iteration=.true.

     !! user-defined taper
     real(kind=CUSTOM_REAL)                                                   :: xmin_taper, xmax_taper
     real(kind=CUSTOM_REAL)                                                   :: ymin_taper, ymax_taper
     real(kind=CUSTOM_REAL)                                                   :: zmin_taper, zmax_taper

     !! parameters for z_precond
     real(kind=CUSTOM_REAL)                                                   :: zPrc1=-200., zPrc2=-400., aPrc=3.

     !! domain boundary
     real(kind=CUSTOM_REAL)                                                   :: xmin, xmax
     real(kind=CUSTOM_REAL)                                                   :: ymin, ymax
     real(kind=CUSTOM_REAL)                                                   :: zmin, zmax

     !! cost
     real(kind=CUSTOM_REAL)                                                   :: total_current_cost, total_previous_cost
     real(kind=CUSTOM_REAL)                                                   :: total_current_prim, total_previous_prim
     real(kind=CUSTOM_REAL)                                                   :: penalty_term, damping_term
     real(kind=CUSTOM_REAL)                                                   :: adjust, penalty
     real(kind=CUSTOM_REAL)                                                   :: Cost_init, Norm_grad_init
     real(kind=CUSTOM_REAL)                                                   :: data_std, nb_data_std

     !!  cost function for each event
     real(kind=CUSTOM_REAL), dimension(:), allocatable                        :: current_cost, previous_cost
     real(kind=CUSTOM_REAL), dimension(:), allocatable                        :: current_cost_prime, previous_cost_prime
     real(kind=CUSTOM_REAL)                                                   :: prior_data_std=1.
     integer                                                                  :: current_iteration = 0
     integer                                                                  :: current_ifrq = 0
     real(kind=CUSTOM_REAL)                                                   :: nb_traces_tot, window_lenght

     !! projection in fd grid
     type(profd)                                                              :: projection_fd

     !! regularization ----------------------------------------------------------------------------------------------------
     logical                                                                  :: use_regularization_FD_Tikonov=.false.
     logical                                                                  :: use_regularization_SEM_Tikonov=.false.
     logical                                                                  :: use_damping_SEM_Tikonov=.false.
     logical                                                                  :: use_variable_SEM_damping=.false.
     real(kind=CUSTOM_REAL)                                                   :: weight_Tikonov=0.1
     real(kind=CUSTOM_REAL)                                                   :: cost_penalty
     real(kind=CUSTOM_REAL)                                                   :: volume_domain
     real(kind=CUSTOM_REAL)                                                   :: min_damp=1., max_damp=10.
     real(kind=CUSTOM_REAL)                                                   :: distance_from_source=100.
     real(kind=CUSTOM_REAL),  dimension(:), allocatable                       :: smooth_weight, damp_weight
     !! prior model
     real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable                :: prior_model
     !!-------------------------------------------------------------------------------------------------------------------

     !! inverted components
     character(len=2),                       dimension(3)                     :: component
     character(len=3)                                                         :: inverted_data_sys
     logical, dimension(3)                                                    :: inverted_data_comp
     character                                                                :: inverted_data_type
     logical                                                                  :: is_src_weigh_gradient=.false.
     logical                                                                  :: convolution_by_wavelet

     !! --- here add parameters for other methods (further developments)
     logical                                                                  :: get_synthetic_pressure=.true.
     logical                                                                  :: get_synthetic_displacement=.true.
     logical                                                                  :: get_synthetic_velocity=.true.
     logical                                                                  :: get_synthetic_acceleration=.true.

  end type inver

  ! ACQUISITION STRUCTURE
  type, public :: acqui

     !!------------------  event general parameters ----------------------
     !! number total of events
     integer                                                                   :: nevent_tot
     !! id for the event
     character(len= MAX_LEN_STRING)                                            :: event_name
     !! name of pif event repository
     character(len= MAX_LEN_STRING)                                            :: event_rep
     !! name for outputs files
     character(len= MAX_LEN_STRING)                                            :: prname_inversion
     !! file contains source parameter for 'moment' or 'fk' or axisem traction
     character(len= MAX_LEN_STRING)                                            :: source_file
     !! kind of source to be used ('moment', 'force', 'axisem', 'dsm', 'fk')
     character(len=MAX_LEN_STRING)                                             :: source_type
     !!
     !! SB SB add source_type_physical and source_type_modeling to distinguish
     !!       between the method used for modeling (local point source(s), injection)
     !!       and the physical type of source ('moment' or 'force') describe by cmt
     !!       or force files. I will use them instead of the above source_type that mix
     !!       these informations.
     !!
     !! kind of source to be used ('moment', 'force')
     character(len=MAX_LEN_STRING)                                             :: source_type_physical
     !! kind of source to be used ('pointsource', 'finitefault','axisem', 'dsm', 'fk')
     character(len=MAX_LEN_STRING)                                             :: source_type_modeling
     !! position of source in case of internal point source
     double precision, dimension(:), allocatable                               :: Xs,Ys,Zs
     !! source time function
     logical                                                                   :: band_pass_filter=.false.
     integer                                                                   :: Nfrq=1
     real(kind=CUSTOM_REAL), dimension(:), allocatable                         :: fl_event, fh_event
     real(kind=CUSTOM_REAL), dimension(:,:), allocatable                       :: user_source_time_function
     !! use external source time function
     character(len= MAX_LEN_STRING)                                            :: source_wavelet_file
     logical                                                                   :: external_source_wavelet=.false.
     !! in case of exploration geophysics,
     !! saving temporary shot point to be able to read it directly in acqui_file
     real(kind=CUSTOM_REAL)                                                    :: xshot, yshot, zshot, shot_ampl
     !! --------------------- source parameter specific for Specfem ---------------------
     !! time parameters needed for specfem
     double precision, dimension(:), allocatable                               :: tshift, hdur, hdur_Gaussian
     !! total number of sources
     integer                                                                   :: nsources_tot
     !! number of sources in my slice
     integer                                                                   :: nsources_local
     !! MPI slice contains source
     integer, dimension(:), allocatable                                        :: islice_selected_source
     !! ispec element contains source
     integer, dimension(:), allocatable                                        :: ispec_selected_source
     real(kind=CUSTOM_REAL),dimension(:,:,:,:,:), allocatable                  :: sourcearrays
     double precision                                                          :: t0
     !! ------------------- station general parameters ----------------------------------
     !! stations (in specfem format)
     character(len= MAX_LEN_STRING)                                            :: station_file
     !! id of stations
     character(len= MAX_LENGTH_STATION_NAME), dimension(:),       allocatable  :: station_name
     !! id of stations
     character(len= MAX_LENGTH_NETWORK_NAME), dimension(:),       allocatable  :: network_name
     !! total number of stations
     integer                                                                   :: nsta_tot
     !! coordinates of stations
     real(kind=CUSTOM_REAL),                  dimension(:,:),     allocatable  :: position_station

     !! ------------------ station parameters specific for specfem ----------------------
     !! general coordinate of receiver
     double precision,                        dimension(:),       allocatable  :: xi_rec,eta_rec,gamma_rec
     !! MPI slice contains receiver
     integer,                                 dimension(:),       allocatable  :: islice_selected_rec
     !! ispec element contains receiver
     integer,                                 dimension(:),       allocatable  :: ispec_selected_rec

     !! ----------------- local in MPI slice ----------------
     !! number of station in my MPI slice
     integer                                                                   :: nsta_slice
     !! global numbering of stations located in my slice
     integer,                                 dimension(:),       allocatable  :: number_receiver_global
     !! rotation matrix for the seismograms
     double precision,                        dimension(:,:,:),   allocatable  :: nu_rec
     !! interpolation along x direction
     double precision,                        dimension(:,:),     allocatable  :: hxi, hpxi
     !! interpolation along y direction
     double precision,                        dimension(:,:),     allocatable  :: heta, hpeta
     !! interpolation along z direction
     double precision,                        dimension(:,:),     allocatable  :: hgamma, hpgamma

     !! ----------------- waveform data local in MPI slice -----------------
     !! file contains the waveform
     character(len= MAX_LEN_STRING)                                            :: data_file_gather
     !! traces stored in memory (NCOMP,NSTA, NT)
     real(kind=CUSTOM_REAL),                  dimension(:,:,:),   allocatable  :: data_traces
     !! synthetics stored in memory (NCOMP,NSTA, NT)
     real(kind=CUSTOM_REAL),                  dimension(:,:,:),   allocatable  :: synt_traces
     !! stored adjoint sources (NCOMP,NSTA, NT)
     real(kind=CUSTOM_REAL),                  dimension(:,:,:),   allocatable  :: adjoint_sources
     !! window used for FWI (NCOMP,2,NSTA)
     real(kind=CUSTOM_REAL),                  dimension(:,:,:),   allocatable  :: window_to_invert
     !! low-high frequency used for FWI (NCOMP,2,NSTA)
     real(kind=CUSTOM_REAL),                  dimension(:,:,:),   allocatable  :: freqcy_to_invert
     !! weigth on each trace used for FWI (NCOMP,NSTA,NT)
     real(kind=CUSTOM_REAL),                  dimension(:,:,:),   allocatable  :: weight_trace

     !! gather time sampling
     real(kind=CUSTOM_REAL)                                                    :: dt_data
     !! number of samples
     integer                                                                   :: Nt_data
     !! components used
     character(len=2),                        dimension(3)                     :: component
     character(len=3)                                                          :: read_data_sys
     logical, dimension(3)                                                     :: read_data_comp
     character                                                                 :: read_data_type
     !! adjoint source to use
     character(len= MAX_LEN_STRING)                                            :: adjoint_source_type

     !! ---------------------- information needed for teleseismic fwi -----------------------------
     !! for rotation matrices
     real(kind=CUSTOM_REAL)                                                    :: Origin_chunk_lat, &
          Origin_chunk_lon, Origin_chunk_azi
     real(kind=CUSTOM_REAL)                                                    :: event_lat,  event_lon, event_depth
     !! window for inversion
     logical                                                                   :: is_time_pick
     logical                                                                   :: time_window
     real(kind=CUSTOM_REAL)                                                    :: time_before_pick, time_after_pick

     !! stations network
     character(len= MAX_LEN_STRING)                                            :: station_coord_system
     !! station list
     real(kind=CUSTOM_REAL), dimension(:,:), allocatable                       :: read_station_position
     real(kind=CUSTOM_REAL), dimension(:), allocatable                         :: time_pick, baz, inc
     real(kind=CUSTOM_REAL), dimension(:), allocatable                         :: dist, gcarc

     !! traction directory in case of AxiSem or DSM coupling
     character(len= MAX_LEN_STRING)                                            :: traction_dir


  end type acqui

  ! SET OF POINT SOURCES  (in development test not working yet)
  type, public :: source_points

     !! number of sources
     integer                                                                   :: NSOURCES
     !! dimension of computed field (eg acoustic, elastic)
     integer                                                                   :: ND
     !! number of time steps
     integer                                                                   :: NtStep
     !! size(NSOURCES)
     integer,                               dimension(:),         allocatable  :: ispec
     !! size(ND,NGLLX,NGLLY,NGLLZ,NSOURCES)
     real(kind=CUSTOM_REAL),                dimension(:,:,:,:,:), allocatable  :: array
     !! size(NSOURCES,ND,NtStep)
     real(kind=CUSTOM_REAL),                dimension(:,:,:),     allocatable  :: stf
     !! size(NDIM,NSOURCES)
     !real(kind=CUSTOM_REAL),                dimension(:,:)                     :: source_position

  end type source_points

  ! STATIONS POINTS (in development test not working yet)
  type, public :: station_network

     integer                                                                   :: NREC
     integer,                               dimension(:),      allocatable     :: ispec
     real(kind=CUSTOM_REAL),                dimension(:,:),    allocatable     :: hxir, hetar, hgammar

  end type station_network

  ! REGULARIZATION STRUCTURE
  type, public :: regul

     integer                                                                  :: iglob
     integer,                               dimension(:),        allocatable  :: iglob_regular_point_to_use
     integer,                               dimension(:),        allocatable  :: iglob_neighbo_point_to_use
     integer,                               dimension(:),        allocatable  :: iglob_neighbo_rank_slice
     real(kind=CUSTOM_REAL),                dimension(:,:),      allocatable  :: Deriv_FD_Matrix
     integer                                                                  :: nReg, nNei, MaxOrder

  end type regul

contains
  ! -------------------------------------------------------------------------------------------------
  !
  ! I/O wrapper function
  !
  !-------------------------------------------------------------------------------------------------

  subroutine flush_iunit(iunit)

    implicit none

    integer, intent(in) :: iunit

    ! note: Fortran2003 includes a FLUSH statement
    !          which is implemented by most compilers by now
    !
    ! otherwise:
    !   a) comment out the line below
    !   b) try to use instead: call flush(iunit)

    flush(iunit)

  end subroutine flush_iunit

end module inverse_problem_par

