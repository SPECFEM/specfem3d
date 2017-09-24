!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
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
!
! United States and French Government Sponsorship Acknowledged.

module inverse_problem_par

  !! IMPORT VARIABLES ------------------------------------------------------------------------------------------------
  use specfem_par, only: CUSTOM_REAL, MAX_STRING_LEN, MAX_LENGTH_STATION_NAME, MAX_LENGTH_NETWORK_NAME
  !-------------------------------------------------------------------------------------------------------------------

  implicit none

  !! ------------------------------ compialtion config parameters -----------------------------------------------------------------

  !! maximum line length allowed in input files
  integer,                       public, parameter  :: MAX_LEN_STRING=256
  !! log file for inversion
  integer,                       public, parameter  :: INVERSE_LOG_FILE=6666, OUTPUT_ITERATION_FILE=6667
  integer,                       public, parameter  :: OUTPUT_FWI_LOG=6668
  !! useful files id
  integer,                       public, parameter  :: IINN=667, IIDD=668
  !!! verbose debug mode
  logical,                       public, parameter  :: DEBUG_MODE=.false.
  !!! verbose mode for with outputs for checking the FWI
  logical,                       public, parameter  :: VERBOSE_MODE=.true.
  !!! write kernels on disk
  logical,                       public, parameter  :: SAVE_KERNEL=.false.
  !!! use fast code and undoing_attenuation for adjoints (in developmement ... not working yet)
  logical,                       public, parameter  :: USE_UNDO_ATT=.false.
  !!! projection on FD grid for outputs (in developmement if useful need to move elsewhere)
  logical,                       public             :: PROJ_ON_FD=.false.
  !!! test for some preconditionners (in developmement if useful need to move elsewhere)
  !!logical,                       public             :: USE_PRECOND_OIL_INDUSTRY=.false.
  !!! if needed use steepest descent instead of l-bfgs (if useful need to move elsewhere)
  logical,                       public             :: USE_GRADIENT_OPTIM=.false.

  ! ------------------------------  global parameters for fwi ---------------------------------------------------------------------
  !! name for outputs files
  character(len=MAX_STRING_LEN), public             :: prname_specfem
  character(len=8),              public             :: prefix_to_path='./'

!################################################# STRUCTURES ######################################################################

  ! INVERSION PARAMETERS STRUCTURE -------------------------------------------------------------------------------------------------
  type, public :: inver

     character(len= MAX_LEN_STRING)                                           :: input_acqui_file
     character(len= MAX_LEN_STRING)                                           :: input_inver_file
     character(len= MAX_LEN_STRING)                                           :: param_family='rho_vp_vs'
     integer                                                                  :: NfamilyPar=3
     integer                                                                  :: NinvPar=3
     integer, dimension(:), allocatable                                       :: Index_Invert

     !! stopping criteria
     integer                                                                  :: Niter = 100
     integer                                                                  :: Niter_wolfe = 10
     real(kind=CUSTOM_REAL)                                                   :: relat_grad=1.e-3
     real(kind=CUSTOM_REAL)                                                   :: relat_cost=1.e-1

     !!---- TUNING PARAMETERS FOR INVERSION -------------
     !! --- for L-BFGS
     integer                                                                  :: max_history_bfgs = 10

     !!--- for Wolfe line search
     real(kind=CUSTOM_REAL)                                                   :: m1 = 0.01, m2 = 0.7
     real(kind=CUSTOM_REAL)                                                   :: current_step_length = 1.

     !! weight for regularization on family parameters
     real(kind=CUSTOM_REAL), dimension(:), allocatable                        :: reg_family, damp_family

     !! max perturbation allowed for initial guest step length relative to gradient max value
     real(kind=CUSTOM_REAL)                                                   :: max_relative_pert = 0.1 ! (1%)

     !! we can import model from disk
     logical                                                                  :: input_sem_model=.false.
     logical                                                                  :: input_fd_model=.false.

     !! write FWI solution to disk (or not)
     logical                                                                  :: output_model=.true.

     !! what to do: FWI, or just forward modeling
     logical                                                                  :: only_forward=.false.

     !! for preconditioning of FWI
     logical                                                                  :: use_taper=.false.
     logical                                                                  :: shin_precond=.false.
     logical                                                                  :: z2_precond=.false.
     logical                                                                  :: energy_precond=.false.

     !! flags to choose if the user wants to dump some wave fields
     logical                                                                  :: dump_model_at_each_iteration=.true.
     logical                                                                  :: dump_gradient_at_each_iteration=.true.
     logical                                                                  :: dump_descent_direction_at_each_iteration=.true.

     !! user-defined taper
     real(kind=CUSTOM_REAL)                                                   :: xmin_taper, xmax_taper
     real(kind=CUSTOM_REAL)                                                   :: ymin_taper, ymax_taper
     real(kind=CUSTOM_REAL)                                                   :: zmin_taper, zmax_taper

     !! domain boundary
     real(kind=CUSTOM_REAL)                                                   :: xmin, xmax
     real(kind=CUSTOM_REAL)                                                   :: ymin, ymax
     real(kind=CUSTOM_REAL)                                                   :: zmin, zmax

     !! cost (need to move it in type inver one day)
     real(kind=CUSTOM_REAL)                                                   :: total_current_cost, total_previous_cost
     real(kind=CUSTOM_REAL)                                                   :: total_current_prim, total_previous_prim
     real(kind=CUSTOM_REAL)                                                   :: penalty_term, damping_term
     real(kind=CUSTOM_REAL)                                                   :: adjust, penalty
     real(kind=CUSTOM_REAL)                                                   :: Cost_init, Norm_grad_init

     !!  cost function for each event
     real(kind=CUSTOM_REAL), dimension(:), allocatable                        :: current_cost, previous_cost
     real(kind=CUSTOM_REAL), dimension(:), allocatable                        :: current_cost_prime, previous_cost_prime
     integer                                                                  :: current_iteration=0

     !! --- here add parameters for other methods (further developments)

  end type inver !------------------------------------------------------------------------------------------------------------------

  ! ACQUISITION STRUCTURE  ---------------------------------------------------------------------------------------------------------
  type, public :: acqui

     !!------------------  event general parameters ----------------------
     !! number total of events
     integer                                                                   :: nevent_tot
     !! id for the event
     character(len= MAX_LEN_STRING)                                            :: event_name
     !! name for outputs files
     character(len= MAX_LEN_STRING)                                            :: prname_inversion
     !! file contains source parameter for 'moment' or 'fk' or axisem traction
     character(len= MAX_LEN_STRING)                                            :: source_file
     !! kind of source to be used ('moment', 'force', 'axisem', 'dsm', 'fk')
     character(len=10)                                                         :: source_type
     !! traction directory in case of AxiSem or DSM coupling
     character(len= MAX_LEN_STRING)                                            :: traction_dir
     !! source time function
     real(kind=CUSTOM_REAL)                                                    :: fl_event, fh_event
     real(kind=CUSTOM_REAL), dimension(:,:), allocatable                       :: user_source_time_function
     !! --------------------- source parameter specific for Specfem ---------------------
     !! time parameters needed for specfem
     double precision, dimension(:), allocatable                               :: tshift, hdur, hdur_gaussian
     !! total number of sources
     integer                                                                   :: nsources_tot
     !! number of sources in my slice
     integer                                                                   :: nsources_local
     !! MPI slice contains source
     integer, dimension(:), allocatable                                        :: islice_selected_source
     !! ispec element contains source
     integer, dimension(:), allocatable                                        :: ispec_selected_source
     real(kind=CUSTOM_REAL),dimension(:,:,:,:,:), allocatable                   :: sourcearrays
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
     double precision,                        dimension(:,:,:),   allocatable  :: nu
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
     !! weigth on each trace used for FWI (NCOMP,NSTA)
     real(kind=CUSTOM_REAL),                  dimension(:,:),     allocatable  :: weight_trace
     !! gather time sampling
     real(kind=CUSTOM_REAL)                                                    :: dt_data
     !! number of samples
     integer                                                                   :: Nt_data
     !! components used
     character(len=2),                        dimension(3)                     :: component

     !! adjoint source to use
     character(len= MAX_LEN_STRING)                                            :: adjoint_source_type

  end type acqui !-----------------------------------------------------------------------------------------------------------------

  ! SET OF POINT SOURCES  (in development test not working yet)--------------------------------------------------------------------
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

  end type source_points  !--------------------------------------------------------------------------------------------------------

  ! STATIONS POINTS (in development test not working yet)--------------------------------------------------------------------------
  type, public :: station_network

     integer                                                                   :: NREC
     integer,                               dimension(:),      allocatable     :: ispec
     real(kind=CUSTOM_REAL),                dimension(:,:),    allocatable     :: hxir, hetar, hgammar

  end type station_network !-------------------------------------------------------------------------------------------------------

  ! REGULARIZATION STRUCTURE ------------------------------------------------------------------------------------------------------
  type, public :: regul

     integer                                                                  :: iglob
     integer,                               dimension(:),        allocatable  :: iglob_regular_point_to_use
     integer,                               dimension(:),        allocatable  :: iglob_neighbo_point_to_use
     integer,                               dimension(:),        allocatable  :: iglob_neighbo_rank_slice
     real(kind=CUSTOM_REAL),                dimension(:,:),      allocatable  :: Deriv_FD_Matrix
     integer                                                                  :: nReg, nNei, MaxOrder

  end type regul !-----------------------------------------------------------------------------------------------------------------

  ! PROJECTION IN FD GRID STRUCTURE -----------------------------------------------------------------------------------------------
  type, public :: profd

     integer                                                                  :: nx, ny, nz
     real(kind=CUSTOM_REAL)                                                   :: hx, hy, hz
     real(kind=CUSTOM_REAL)                                                   :: ox, oy, oz

     integer                                                                  :: nb_fd_point
     integer,                               dimension(:),         allocatable :: ispec_selected
     integer,                               dimension(:,:),       allocatable :: index_on_fd_grid
     double precision,                      dimension(:,:),       allocatable :: hxi, heta, hgamma

  end type profd !-----------------------------------------------------------------------------------------------------------------

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

