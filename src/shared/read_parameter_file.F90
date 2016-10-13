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

  subroutine read_parameter_file(myrank,BROADCAST_AFTER_READ)

  use constants

  use shared_parameters

  implicit none

  integer, intent(in) :: myrank
  logical, intent(in) :: BROADCAST_AFTER_READ

  ! local variables
  integer :: ier
  integer :: nproc_eta_old,nproc_xi_old

  character(len=MAX_STRING_LEN) :: MODEL

  double precision :: minval_hdur
  character(len=MAX_STRING_LEN) :: CMTSOLUTION,FORCESOLUTION

  character(len=MAX_STRING_LEN) :: path_to_add

  !logical :: sep_dir_exists
  integer :: i,irange

  !LDDRK
  logical :: INCREASE_CFL_FOR_LDDRK
  double precision :: RATIO_BY_WHICH_TO_INCREASE_IT

  ! read from a single processor (the master) and then use MPI to broadcast to others
  ! to avoid an I/O bottleneck in the case of very large runs
  if (myrank == 0) then

    ! opens file Par_file
    call open_parameter_file(ier)

    ! reads in mandatory simulation parameters

    !-------------------------------------------------------
    ! Simulation input parameters
    !-------------------------------------------------------
    call read_value_integer(SIMULATION_TYPE, 'SIMULATION_TYPE', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter SIMULATION_TYPE'
    call read_value_integer(NOISE_TOMOGRAPHY, 'NOISE_TOMOGRAPHY', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter NOISE_TOMOGRAPHY'
    call read_value_logical(SAVE_FORWARD, 'SAVE_FORWARD', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter SAVE_FORWARD'

    call read_value_integer(UTM_PROJECTION_ZONE, 'UTM_PROJECTION_ZONE', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter UTM_PROJECTION_ZONE'
    call read_value_logical(SUPPRESS_UTM_PROJECTION, 'SUPPRESS_UTM_PROJECTION', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter SUPPRESS_UTM_PROJECTION'

    ! total number of processors
    call read_value_integer(NPROC, 'NPROC', ier)
    if (ier /= 0) then
      ! checks if it uses an old Par_file format
      call read_value_integer(nproc_eta_old, 'NPROC_ETA', ier)
      if (ier /= 0) then
        print *,'please specify the number of processes in Par_file as:'
        print *,'NPROC           = my_number_of_desired_processes'
        stop 'Error reading Par_file parameter NPROC'
      endif
      ! checks if it uses an old Par_file format
      call read_value_integer(nproc_xi_old, 'NPROC_XI', ier)
      if (ier /= 0) then
        print *,'please specify the number of processes in Par_file as:'
        print *,'NPROC           = my_number_of_desired_processes'
        stop 'Error reading Par_file parameter NPROC'
      endif
      NPROC = nproc_eta_old * nproc_xi_old
    endif

    call read_value_integer(NSTEP, 'NSTEP', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter NSTEP'
    call read_value_double_precision(DT, 'DT', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter DT'

    !-------------------------------------------------------
    ! LDDRK time scheme
    !-------------------------------------------------------
    call read_value_logical(USE_LDDRK, 'USE_LDDRK', ier)   ! low-memory Runge-Kutta time scheme
    if (ier /= 0) stop 'an error occurred while reading the parameter file: USE_LDDRK'
    call read_value_logical(INCREASE_CFL_FOR_LDDRK, 'INCREASE_CFL_FOR_LDDRK', ier)
    if (ier /= 0) stop 'an error occurred while reading the parameter file: INCREASE_CFL_FOR_LDDRK'
    call read_value_double_precision(RATIO_BY_WHICH_TO_INCREASE_IT, 'RATIO_BY_WHICH_TO_INCREASE_IT', ier)
    if (ier /= 0) stop 'an error occurred while reading the parameter file: RATIO_BY_WHICH_TO_INCREASE_IT'

    !-------------------------------------------------------
    ! Mesh
    !-------------------------------------------------------
    call read_value_integer(NGNOD, 'NGNOD', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter NGNOD'
    call read_value_string(MODEL, 'MODEL', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter MODEL'

    call read_value_string(TOMOGRAPHY_PATH, 'TOMOGRAPHY_PATH', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter TOMOGRAPHY_PATH'
    call read_value_string(SEP_MODEL_DIRECTORY, 'SEP_MODEL_DIRECTORY', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter SEP_MODEL_DIRECTORY'

    !-------------------------------------------------------

    call read_value_logical(APPROXIMATE_OCEAN_LOAD, 'APPROXIMATE_OCEAN_LOAD', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter APPROXIMATE_OCEAN_LOAD'
    call read_value_logical(TOPOGRAPHY, 'TOPOGRAPHY', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter TOPOGRAPHY'
    call read_value_logical(ATTENUATION, 'ATTENUATION', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter ATTENUATION'
    call read_value_logical(ANISOTROPY, 'ANISOTROPY', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter ANISOTROPY'
    call read_value_logical(GRAVITY, 'GRAVITY', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter GRAVITY'

    call read_value_double_precision(ATTENUATION_f0_REFERENCE, 'ATTENUATION_f0_REFERENCE', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter ATTENUATION_f0_REFERENCE'

    call read_value_logical(USE_OLSEN_ATTENUATION, 'USE_OLSEN_ATTENUATION', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter USE_OLSEN_ATTENUATION'
    call read_value_double_precision(OLSEN_ATTENUATION_RATIO, 'OLSEN_ATTENUATION_RATIO', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter OLSEN_ATTENUATION_RATIO'

    !-------------------------------------------------------
    ! Absorbing boundary conditions
    !-------------------------------------------------------
    call read_value_logical(PML_CONDITIONS, 'PML_CONDITIONS', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter PML_CONDITIONS'
    call read_value_logical(PML_INSTEAD_OF_FREE_SURFACE, 'PML_INSTEAD_OF_FREE_SURFACE', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter PML_INSTEAD_OF_FREE_SURFACE'
    call read_value_double_precision(f0_FOR_PML, 'f0_FOR_PML', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter f0_FOR_PML'

    call read_value_logical(STACEY_ABSORBING_CONDITIONS, 'STACEY_ABSORBING_CONDITIONS', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter STACEY_ABSORBING_CONDITIONS'
    call read_value_logical(STACEY_INSTEAD_OF_FREE_SURFACE, 'STACEY_INSTEAD_OF_FREE_SURFACE', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter STACEY_INSTEAD_OF_FREE_SURFACE'
    call read_value_logical(BOTTOM_FREE_SURFACE, 'BOTTOM_FREE_SURFACE', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter BOTTOM_FREE_SURFACE'

    !-------------------------------------------------------
    ! Visualization
    !-------------------------------------------------------
    call read_value_logical(CREATE_SHAKEMAP, 'CREATE_SHAKEMAP', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter CREATE_SHAKEMAP'
    call read_value_logical(MOVIE_SURFACE, 'MOVIE_SURFACE', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter MOVIE_SURFACE'
    call read_value_integer(MOVIE_TYPE, 'MOVIE_TYPE', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter MOVIE_TYPE'
    call read_value_logical(MOVIE_VOLUME, 'MOVIE_VOLUME', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter MOVIE_VOLUME'
    call read_value_logical(SAVE_DISPLACEMENT, 'SAVE_DISPLACEMENT', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter SAVE_DISPLACEMENT'
    call read_value_logical(USE_HIGHRES_FOR_MOVIES, 'USE_HIGHRES_FOR_MOVIES', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter USE_HIGHRES_FOR_MOVIES'
    call read_value_integer(NTSTEP_BETWEEN_FRAMES, 'NTSTEP_BETWEEN_FRAMES', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter NTSTEP_BETWEEN_FRAMES'
    call read_value_double_precision(HDUR_MOVIE, 'HDUR_MOVIE', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter HDUR_MOVIE'

    call read_value_logical(SAVE_MESH_FILES, 'SAVE_MESH_FILES', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter SAVE_MESH_FILES'
    call read_value_string(LOCAL_PATH, 'LOCAL_PATH', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter LOCAL_PATH'
    call read_value_integer(NTSTEP_BETWEEN_OUTPUT_INFO, 'NTSTEP_BETWEEN_OUTPUT_INFO', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter NTSTEP_BETWEEN_OUTPUT_INFO'

    !-------------------------------------------------------
    ! Sources
    !-------------------------------------------------------
    call read_value_logical(USE_FORCE_POINT_SOURCE, 'USE_FORCE_POINT_SOURCE', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter USE_FORCE_POINT_SOURCE'
    call read_value_logical(USE_RICKER_TIME_FUNCTION, 'USE_RICKER_TIME_FUNCTION', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter USE_RICKER_TIME_FUNCTION'

    call read_value_logical(PRINT_SOURCE_TIME_FUNCTION, 'PRINT_SOURCE_TIME_FUNCTION', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter PRINT_SOURCE_TIME_FUNCTION'

    !-------------------------------------------------------
    ! Seismograms
    !-------------------------------------------------------
    call read_value_integer(NTSTEP_BETWEEN_OUTPUT_SEISMOS, 'NTSTEP_BETWEEN_OUTPUT_SEISMOS', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter NTSTEP_BETWEEN_OUTPUT_SEISMOS'

    call read_value_logical(SAVE_SEISMOGRAMS_DISPLACEMENT, 'SAVE_SEISMOGRAMS_DISPLACEMENT', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter SAVE_SEISMOGRAMS_DISPLACEMENT'
    call read_value_logical(SAVE_SEISMOGRAMS_VELOCITY, 'SAVE_SEISMOGRAMS_VELOCITY', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter SAVE_SEISMOGRAMS_VELOCITY'
    call read_value_logical(SAVE_SEISMOGRAMS_ACCELERATION, 'SAVE_SEISMOGRAMS_ACCELERATION', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter SAVE_SEISMOGRAMS_ACCELERATION'
    call read_value_logical(SAVE_SEISMOGRAMS_PRESSURE, 'SAVE_SEISMOGRAMS_PRESSURE', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter SAVE_SEISMOGRAMS_PRESSURE'

    call read_value_logical(USE_BINARY_FOR_SEISMOGRAMS, 'USE_BINARY_FOR_SEISMOGRAMS', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter USE_BINARY_FOR_SEISMOGRAMS'
    call read_value_logical(SU_FORMAT, 'SU_FORMAT', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter SU_FORMAT'
    call read_value_logical(WRITE_SEISMOGRAMS_BY_MASTER, 'WRITE_SEISMOGRAMS_BY_MASTER', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter WRITE_SEISMOGRAMS_BY_MASTER'
    call read_value_logical(SAVE_ALL_SEISMOS_IN_ONE_FILE, 'SAVE_ALL_SEISMOS_IN_ONE_FILE', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter SAVE_ALL_SEISMOS_IN_ONE_FILE'
    call read_value_logical(USE_TRICK_FOR_BETTER_PRESSURE, 'USE_TRICK_FOR_BETTER_PRESSURE', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter USE_TRICK_FOR_BETTER_PRESSURE'

    !-------------------------------------------------------
    ! Source encoding
    !-------------------------------------------------------
    call read_value_logical(USE_SOURCE_ENCODING, 'USE_SOURCE_ENCODING', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter USE_SOURCE_ENCODING'

    !-------------------------------------------------------
    ! Total energy calculation
    !-------------------------------------------------------
    call read_value_logical(OUTPUT_ENERGY, 'OUTPUT_ENERGY', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter OUTPUT_ENERGY'

    call read_value_integer(NTSTEP_BETWEEN_OUTPUT_ENERGY, 'NTSTEP_BETWEEN_OUTPUT_ENERGY', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter NTSTEP_BETWEEN_OUTPUT_ENERGY'

    !-------------------------------------------------------
    ! Adjoint kernel outputs
    !-------------------------------------------------------
    call read_value_integer(NTSTEP_BETWEEN_READ_ADJSRC, 'NTSTEP_BETWEEN_READ_ADJSRC', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter NTSTEP_BETWEEN_READ_ADJSRC'


    call read_value_logical(ANISOTROPIC_KL, 'ANISOTROPIC_KL', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter ANISOTROPIC_KL'
    call read_value_logical(SAVE_TRANSVERSE_KL, 'SAVE_TRANSVERSE_KL', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter SAVE_TRANSVERSE_KL'
    call read_value_logical(APPROXIMATE_HESS_KL, 'APPROXIMATE_HESS_KL', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter APPROXIMATE_HESS_KL'
    call read_value_logical(SAVE_MOHO_MESH, 'SAVE_MOHO_MESH', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter SAVE_MOHO_MESH'

    !-------------------------------------------------------

    ! for simultaneous runs from the same batch job
    call read_value_integer(NUMBER_OF_SIMULTANEOUS_RUNS, 'NUMBER_OF_SIMULTANEOUS_RUNS', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter NUMBER_OF_SIMULTANEOUS_RUNS'
    call read_value_logical(BROADCAST_SAME_MESH_AND_MODEL, 'BROADCAST_SAME_MESH_AND_MODEL', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter BROADCAST_SAME_MESH_AND_MODEL'

    !-------------------------------------------------------

    call read_value_logical(GPU_MODE, 'GPU_MODE', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter GPU_MODE'

  !> Read ADIOS related flags from the Par_file
  !! \param ADIOS_ENABLED Main flag to decide if ADIOS is used. If setted to
  !!                      false no other parameter is taken into account.
  !! \param ADIOS_FOR_DATABASES Flag to indicate if the databases are written
  !!                            and read with the help of ADIOS.
  !! \param ADIOS_FOR_MESH flag to indicate if the mesh (generate database) is
  !!                       written using ADIOS.
  !! \param ADIOS_FOR_FORWARD_ARRAYS flag to indicate if the solver forward arrays
  !!                                 are written using ADIOS.
  !! \param ADIOS_FOR_KERNELS flag to indicate if the kernels are saved using
  !!                          adios
  !! \author MPBL
    call read_value_logical(ADIOS_ENABLED, 'ADIOS_ENABLED', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter ADIOS_ENABLED'
    call read_value_logical(ADIOS_FOR_DATABASES, 'ADIOS_FOR_DATABASES', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter ADIOS_FOR_DATABASES'
    call read_value_logical(ADIOS_FOR_MESH, 'ADIOS_FOR_MESH', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter ADIOS_FOR_MESH'
    call read_value_logical(ADIOS_FOR_FORWARD_ARRAYS, 'ADIOS_FOR_FORWARD_ARRAYS', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter ADIOS_FOR_FORWARD_ARRAYS'
    call read_value_logical(ADIOS_FOR_KERNELS, 'ADIOS_FOR_KERNELS', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter ADIOS_FOR_KERNELS'

    call read_value_logical(EXTERNAL_STF, 'EXTERNAL_SOURCE_FILE', ier)
    if (ier /= 0) stop 'Error reading Par_file parameter EXTERNAL_SOURCE_FILE'

#ifdef DEBUG_COUPLED
      include "../../../add_to_read_parameter_file_1.F90"
#endif

    ! closes parameter file
    call close_parameter_file()

    ! checks parameter settings
    call check_simulation_parameters()

    ! updates values for simulation settings
    ! LDDRK scheme
    if (USE_LDDRK .and. INCREASE_CFL_FOR_LDDRK) DT = DT * RATIO_BY_WHICH_TO_INCREASE_IT

    ! to be consistent with external source file option turned to on
    if (EXTERNAL_STF) USE_RICKER_TIME_FUNCTION = .false.

    ! see if we are running several independent runs in parallel
    ! if so, add the right directory for that run
    ! (group numbers start at zero, but directory names start at run0001, thus we add one)
    ! a negative value for "mygroup" is a convention that indicates that groups (i.e. sub-communicators, one per run) are off
    if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
      write(path_to_add,"('run',i4.4,'/')") mygroup + 1
      LOCAL_PATH = path_to_add(1:len_trim(path_to_add))//LOCAL_PATH(1:len_trim(LOCAL_PATH))
      TOMOGRAPHY_PATH = path_to_add(1:len_trim(path_to_add))//TOMOGRAPHY_PATH(1:len_trim(TOMOGRAPHY_PATH))
#ifdef DEBUG_COUPLED
      include "../../../add_to_read_parameter_file_2.F90"
#endif
    endif

    ! noise simulations:
    if (NOISE_TOMOGRAPHY /= 0) then
      ! double the number of time steps, if running noise simulations (+/- branches)
      NSTEP = 2 * NSTEP - 1

      ! for noise simulations, we need to save movies at the surface (where the noise is generated)
      ! and thus we force MOVIE_SURFACE to be .true., in order to use variables defined for surface movies later
      MOVIE_TYPE = 1
      MOVIE_SURFACE = .true.
      USE_HIGHRES_FOR_MOVIES = .true.     ! we need to save surface movie everywhere, i.e. at all GLL points on the surface
    endif

    ! the default value of NTSTEP_BETWEEN_READ_ADJSRC (0) is to read the whole trace at the same time
    if (NTSTEP_BETWEEN_READ_ADJSRC == 0)  NTSTEP_BETWEEN_READ_ADJSRC = NSTEP

    ! total times steps must be dividable by adjoint source chunks/blocks
    if (mod(NSTEP,NTSTEP_BETWEEN_READ_ADJSRC) /= 0) then
      print *,'When NOISE_TOMOGRAPHY is not equal to zero, ACTUAL_NSTEP=2*NSTEP-1'
      stop 'Error: mod(NSTEP,NTSTEP_BETWEEN_READ_ADJSRC) must be zero! Please modify Par_file and recompile solver'
    endif

    ! checks number of nodes for 2D and 3D shape functions for quadrilaterals and hexahedra
    ! curvature (i.e. HEX27 elements) is not handled by our internal mesher, for that use Gmsh (CUBIT does not handle it either)
    if (NGNOD == 8) then
      NGNOD2D = 4
    else if (NGNOD == 27) then
      NGNOD2D = 9
    else if (NGNOD /= 8 .and. NGNOD /= 27) then
      stop 'Error elements should have 8 or 27 control nodes, please modify NGNOD in Par_file and recompile solver'
    endif

    ! determines number of sources depending on number of lines in source file
    call get_number_of_sources(minval_hdur,FORCESOLUTION,CMTSOLUTION)

    ! converts all string characters to lowercase
    irange = iachar('a') - iachar('A')
    do i = 1,len_trim(MODEL)
      if (lge(MODEL(i:i),'A') .and. lle(MODEL(i:i),'Z')) then
        MODEL(i:i) = achar(iachar(MODEL(i:i)) + irange)
      endif
    enddo

    ! determines velocity model
    select case (trim(MODEL))

    ! default mesh model
    case ('default')
      IMODEL = IMODEL_DEFAULT

    ! 1-D models
    case ('1d_prem')
      IMODEL = IMODEL_1D_PREM
    case ('1d_socal')
      IMODEL = IMODEL_1D_SOCAL
    case ('1d_cascadia')
      IMODEL = IMODEL_1D_CASCADIA

    ! user models
    case ('1d_prem_pb')
      IMODEL = IMODEL_1D_PREM_PB
    case ('aniso')
      IMODEL = IMODEL_DEFAULT
      ANISOTROPY = .true.
    case ('external')
      IMODEL = IMODEL_USER_EXTERNAL
    case ('ipati')
      IMODEL = IMODEL_IPATI
    case ('ipati_water')
      IMODEL = IMODEL_IPATI_WATER
    case ('gll')
      IMODEL = IMODEL_GLL
    case ('salton_trough')
      IMODEL = IMODEL_SALTON_TROUGH
    case ('tomo')
      IMODEL = IMODEL_TOMO
    case ('sep')
      IMODEL = IMODEL_SEP
      if (trim(SEP_MODEL_DIRECTORY) == '') then
        stop 'Error using sep model requires defining a SEP_MODEL_DIRECTORY.'
      endif
      !inquire(directory=trim(SEP_MODEL_DIRECTORY), exists=sep_dir_exists)
      !if (.not. sep_dir_exists) then
      !  stop 'Error: SEP_MODEL_DIRECTORY should exist.'
      !endif
    case default
      print *
      print *,'********** model not recognized: ',trim(MODEL),' **************'
      print *,'********** using model: default',' **************'
      print *
      IMODEL = IMODEL_DEFAULT
    end select

    ! check
    if (IMODEL == IMODEL_IPATI .or. IMODEL == IMODEL_IPATI_WATER) then
      if (USE_RICKER_TIME_FUNCTION .eqv. .false.) &
        stop 'Error for IPATI model, please set USE_RICKER_TIME_FUNCTION to .true. in Par_file and recompile solver'
    endif

  endif ! of if (myrank == 0) then

! read from a single processor (the master) and then use MPI to broadcast to others
! to avoid an I/O bottleneck in the case of very large runs
  if (BROADCAST_AFTER_READ) then

    call bcast_all_singlei_world(NPROC)
    call bcast_all_singlei_world(SIMULATION_TYPE)
    call bcast_all_singlei_world(NOISE_TOMOGRAPHY)
    call bcast_all_singlel_world(SAVE_FORWARD)
    call bcast_all_singlei_world(UTM_PROJECTION_ZONE)
    call bcast_all_singlel_world(SUPPRESS_UTM_PROJECTION)
    call bcast_all_singlei_world(NSTEP)
    call bcast_all_singledp_world(DT)
    call bcast_all_singlel_world(USE_LDDRK)
    call bcast_all_singlei_world(NGNOD)
    call bcast_all_string_world(MODEL)
    call bcast_all_string_world(SEP_MODEL_DIRECTORY)
    call bcast_all_singlel_world(APPROXIMATE_OCEAN_LOAD)
    call bcast_all_singlel_world(TOPOGRAPHY)
    call bcast_all_singlel_world(ATTENUATION)
    call bcast_all_singlel_world(ANISOTROPY)
    call bcast_all_singlel_world(GRAVITY)
    call bcast_all_singledp_world(ATTENUATION_f0_REFERENCE)
    call bcast_all_singlel_world(USE_OLSEN_ATTENUATION)
    call bcast_all_singledp_world(OLSEN_ATTENUATION_RATIO)
    call bcast_all_string_world(TOMOGRAPHY_PATH)
    call bcast_all_singlel_world(PML_CONDITIONS)
    call bcast_all_singlel_world(PML_INSTEAD_OF_FREE_SURFACE)
    call bcast_all_singledp_world(f0_FOR_PML)
    call bcast_all_singlel_world(STACEY_ABSORBING_CONDITIONS)
    call bcast_all_singlel_world(STACEY_INSTEAD_OF_FREE_SURFACE)
    call bcast_all_singlel_world(BOTTOM_FREE_SURFACE)
    call bcast_all_singlel_world(CREATE_SHAKEMAP)
    call bcast_all_singlel_world(MOVIE_SURFACE)
    call bcast_all_singlei_world(MOVIE_TYPE)
    call bcast_all_singlel_world(MOVIE_VOLUME)
    call bcast_all_singlel_world(SAVE_DISPLACEMENT)
    call bcast_all_singlel_world(USE_HIGHRES_FOR_MOVIES)
    call bcast_all_singlei_world(NTSTEP_BETWEEN_FRAMES)
    call bcast_all_singledp_world(HDUR_MOVIE)
    call bcast_all_singlel_world(SAVE_MESH_FILES)
    call bcast_all_string_world(LOCAL_PATH)
    call bcast_all_singlei_world(NTSTEP_BETWEEN_OUTPUT_INFO)
    call bcast_all_singlei_world(NTSTEP_BETWEEN_OUTPUT_SEISMOS)
    call bcast_all_singlei_world(NTSTEP_BETWEEN_READ_ADJSRC)
    call bcast_all_singlel_world(USE_FORCE_POINT_SOURCE)
    call bcast_all_singlel_world(USE_RICKER_TIME_FUNCTION)
    call bcast_all_singlel_world(SAVE_SEISMOGRAMS_DISPLACEMENT)
    call bcast_all_singlel_world(SAVE_SEISMOGRAMS_VELOCITY)
    call bcast_all_singlel_world(SAVE_SEISMOGRAMS_ACCELERATION)
    call bcast_all_singlel_world(SAVE_SEISMOGRAMS_PRESSURE)
    call bcast_all_singlel_world(USE_BINARY_FOR_SEISMOGRAMS)
    call bcast_all_singlel_world(SU_FORMAT)
    call bcast_all_singlel_world(WRITE_SEISMOGRAMS_BY_MASTER)
    call bcast_all_singlel_world(SAVE_ALL_SEISMOS_IN_ONE_FILE)
    call bcast_all_singlel_world(USE_TRICK_FOR_BETTER_PRESSURE)
    call bcast_all_singlel_world(USE_SOURCE_ENCODING)
    call bcast_all_singlel_world(OUTPUT_ENERGY)
    call bcast_all_singlei_world(NTSTEP_BETWEEN_OUTPUT_ENERGY)
    call bcast_all_singlel_world(ANISOTROPIC_KL)
    call bcast_all_singlel_world(SAVE_TRANSVERSE_KL)
    call bcast_all_singlel_world(APPROXIMATE_HESS_KL)
    call bcast_all_singlel_world(SAVE_MOHO_MESH)
    call bcast_all_singlel_world(PRINT_SOURCE_TIME_FUNCTION)
    call bcast_all_singlei_world(NUMBER_OF_SIMULTANEOUS_RUNS)
    call bcast_all_singlel_world(BROADCAST_SAME_MESH_AND_MODEL)
    call bcast_all_singlel_world(GPU_MODE)
    call bcast_all_singlel_world(ADIOS_ENABLED)
    call bcast_all_singlel_world(ADIOS_FOR_DATABASES)
    call bcast_all_singlel_world(ADIOS_FOR_MESH)
    call bcast_all_singlel_world(ADIOS_FOR_FORWARD_ARRAYS)
    call bcast_all_singlel_world(ADIOS_FOR_KERNELS)
    call bcast_all_singlel_world(EXTERNAL_STF)

! broadcast all parameters computed from others
    call bcast_all_singlei_world(IMODEL)
    call bcast_all_singlei_world(NGNOD2D)

    call bcast_all_singlei_world(NSOURCES)
    call bcast_all_singledp_world(minval_hdur)
    call bcast_all_string_world(FORCESOLUTION)
    call bcast_all_string_world(CMTSOLUTION)

  endif ! of if (BROADCAST_AFTER_READ) then

  end subroutine read_parameter_file

!
!-------------------------------------------------------------------------------------------------
!

  subroutine check_simulation_parameters()

! safety checks
! some features might not be implemented depending on setup

!  use constants
  use shared_parameters

  implicit none

  ! time steps
  if (NSTEP <= 0) &
    stop 'NSTEP must be > 0 for any simulation'

  ! seismogram output
  if (.not. SAVE_SEISMOGRAMS_DISPLACEMENT .and. .not. SAVE_SEISMOGRAMS_VELOCITY .and. &
     .not. SAVE_SEISMOGRAMS_ACCELERATION .and. .not. SAVE_SEISMOGRAMS_PRESSURE) &
   stop 'Error: at least one of SAVE_SEISMOGRAMS_DISPLACEMENT SAVE_SEISMOGRAMS_VELOCITY SAVE_SEISMOGRAMS_ACCELERATION &
             &SAVE_SEISMOGRAMS_PRESSURE must be true'

  ! this could be implemented in the future if needed,
  ! see comments in the source code around the USE_TRICK_FOR_BETTER_PRESSURE
  ! option (use a "grep" command to find them) to see how this could/should be done
  if (USE_TRICK_FOR_BETTER_PRESSURE .and. (SAVE_SEISMOGRAMS_DISPLACEMENT .or. SAVE_SEISMOGRAMS_VELOCITY .or. &
        SAVE_SEISMOGRAMS_ACCELERATION)) stop 'USE_TRICK_FOR_BETTER_PRESSURE is currently incompatible with &
        &SAVE_SEISMOGRAMS_DISPLACEMENT .or. SAVE_SEISMOGRAMS_VELOCITY .or. SAVE_SEISMOGRAMS_ACCELERATION, &
        &only SAVE_SEISMOGRAMS_PRESSURE can be used'

  ! LDDRK
  if (USE_LDDRK) then
    if (SIMULATION_TYPE == 3 ) &
      stop 'USE_LDDRK support not implemented yet for SIMULATION_TYPE == 3'

    ! GPU mode
    if (GPU_MODE ) &
      stop 'USE_LDDRK support not implemented yet for GPU simulations'
  endif

  ! ADIOS safety check
  if (ADIOS_ENABLED) then
    if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. BROADCAST_SAME_MESH_AND_MODEL) &
      stop 'Error ADIOS not yet supported by option BROADCAST_SAME_MESH_AND_MODEL'
  endif

  ! stacey absorbing free surface must have also stacey turned on
  if (STACEY_INSTEAD_OF_FREE_SURFACE) then
    if (.not. STACEY_ABSORBING_CONDITIONS) &
      stop 'STACEY_INSTEAD_OF_FREE_SURFACE must have also STACEY_ABSORBING_CONDITIONS turned on'
  endif

  ! PML
  if (PML_CONDITIONS) then
!! DK DK added this for now (March 2013)
!! DK DK we will soon add it
    if (SAVE_FORWARD .or. SIMULATION_TYPE == 3) &
      stop 'PML_CONDITIONS is still under test for adjoint simulation'

    ! only one absorbing condition is possible
    if (STACEY_ABSORBING_CONDITIONS) &
      stop 'Error for PML, please set STACEY_ABSORBING_CONDITIONS and STACEY_INSTEAD_OF_FREE_SURFACE to .false. in Par_file'
  endif

  ! external STF
  if (EXTERNAL_STF .and. GPU_MODE) &
    stop 'EXTERNAL_SOURCE_FILE in GPU_MODE simulation not supported yet'

  end subroutine check_simulation_parameters

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_number_of_sources(minval_hdur,FORCESOLUTION,CMTSOLUTION)

! determines number of sources depending on number of lines in source file

  use constants, only: mygroup,IN_DATA_FILES,HUGEVAL,TINYVAL, &
    NLINES_PER_CMTSOLUTION_SOURCE,NLINES_PER_FORCESOLUTION_SOURCE
  use shared_parameters

  implicit none

  double precision,intent(out) :: minval_hdur
  character(len=MAX_STRING_LEN),intent(out) :: CMTSOLUTION,FORCESOLUTION

  ! local variables
  integer :: icounter,isource,idummy,ier
  double precision :: hdur
  character(len=MAX_STRING_LEN) :: dummystring
  character(len=MAX_STRING_LEN) :: path_to_add

  if (USE_FORCE_POINT_SOURCE) then
    ! compute the total number of sources in the FORCESOLUTION file
    ! there are NLINES_PER_FORCESOLUTION_SOURCE lines per source in that file
    FORCESOLUTION = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'FORCESOLUTION'
! see if we are running several independent runs in parallel
! if so, add the right directory for that run
! (group numbers start at zero, but directory names start at run0001, thus we add one)
! a negative value for "mygroup" is a convention that indicates that groups (i.e. sub-communicators, one per run) are off
    if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
      write(path_to_add,"('run',i4.4,'/')") mygroup + 1
      FORCESOLUTION = path_to_add(1:len_trim(path_to_add))//FORCESOLUTION(1:len_trim(FORCESOLUTION))
    endif

    open(unit=21,file=trim(FORCESOLUTION),status='old',action='read',iostat=ier)
    if (ier /= 0) stop 'Error opening FORCESOLUTION file'

    icounter = 0
    do while (ier == 0)
      read(21,"(a)",iostat=ier) dummystring
      if (ier == 0) icounter = icounter + 1
    enddo
    close(21)

    if (.not. EXTERNAL_STF) then
       if (mod(icounter,NLINES_PER_FORCESOLUTION_SOURCE) /= 0) &
            stop 'Error total number of lines in FORCESOLUTION file should be a multiple of NLINES_PER_FORCESOLUTION_SOURCE'
       NSOURCES = icounter / NLINES_PER_FORCESOLUTION_SOURCE
    else !! VM VM in case of EXTERNAL_STF we have to read one additional line per source (the name of external source file)
       NSOURCES = icounter / (NLINES_PER_FORCESOLUTION_SOURCE+1)
    endif

    if (NSOURCES < 1) stop 'Error need at least one source in FORCESOLUTION file'

  else
    ! compute the total number of sources in the CMTSOLUTION file
    ! there are NLINES_PER_CMTSOLUTION_SOURCE lines per source in that file
    CMTSOLUTION = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'CMTSOLUTION'
! see if we are running several independent runs in parallel
! if so, add the right directory for that run
! (group numbers start at zero, but directory names start at run0001, thus we add one)
! a negative value for "mygroup" is a convention that indicates that groups (i.e. sub-communicators, one per run) are off
    if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
      write(path_to_add,"('run',i4.4,'/')") mygroup + 1
      CMTSOLUTION = path_to_add(1:len_trim(path_to_add))//CMTSOLUTION(1:len_trim(CMTSOLUTION))
    endif

    open(unit=21,file=trim(CMTSOLUTION),status='old',action='read',iostat=ier)
    if (ier /= 0) stop 'Error opening CMTSOLUTION file'

    icounter = 0
    do while (ier == 0)
      read(21,"(a)",iostat=ier) dummystring
      if (ier == 0) icounter = icounter + 1
    enddo
    close(21)
    if (.not. EXTERNAL_STF) then
       if (mod(icounter,NLINES_PER_CMTSOLUTION_SOURCE) /= 0) &
            stop 'Error total number of lines in CMTSOLUTION file should be a multiple of NLINES_PER_CMTSOLUTION_SOURCE'
    else
       !! VM VM in case of EXTERNAL_STF we have to read one additional line per source (the name of external source file)
       NSOURCES = icounter / (NLINES_PER_FORCESOLUTION_SOURCE+1)
    endif

    NSOURCES = icounter / NLINES_PER_CMTSOLUTION_SOURCE
    if (NSOURCES < 1) stop 'Error need at least one source in CMTSOLUTION file'

    ! compute the minimum value of hdur in CMTSOLUTION file
    open(unit=21,file=trim(CMTSOLUTION),status='old',action='read')
    minval_hdur = HUGEVAL
    do isource = 1,NSOURCES

      ! skip other information
      do idummy = 1,3
        read(21,"(a)") dummystring
      enddo

      ! read half duration and compute minimum
      read(21,"(a)") dummystring
      read(dummystring(15:len_trim(dummystring)),*) hdur
      minval_hdur = min(minval_hdur,hdur)

      ! skip other information
      do idummy = 1,9
        read(21,"(a)") dummystring
      enddo

    enddo
    close(21)

    ! one cannot use a Heaviside source for the movies
    if ((MOVIE_SURFACE .or. MOVIE_VOLUME) .and. sqrt(minval_hdur**2 + HDUR_MOVIE**2) < TINYVAL) &
      stop 'Error hdur too small for movie creation, movies do not make sense for Heaviside source'
  endif

  end subroutine get_number_of_sources


