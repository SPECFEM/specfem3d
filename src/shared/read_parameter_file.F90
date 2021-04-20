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

  subroutine read_parameter_file(BROADCAST_AFTER_READ)

  use constants, only: &
    myrank, &
    INJECTION_TECHNIQUE_IS_AXISEM,INJECTION_TECHNIQUE_IS_DSM,INJECTION_TECHNIQUE_IS_FK

  use shared_parameters

  implicit none

  logical, intent(in) :: BROADCAST_AFTER_READ

  ! local variables
  integer :: ier
  logical :: some_parameters_missing_from_Par_file

  ! read from a single processor (the main) and then use MPI to broadcast to others
  ! to avoid an I/O bottleneck in the case of very large runs
  if (myrank == 0) then

    ! opens file Par_file
    call open_parameter_file(ier)

    ! reads in mandatory simulation parameters

    some_parameters_missing_from_Par_file = .false.

    !-------------------------------------------------------
    ! Simulation input parameters
    !-------------------------------------------------------
    call read_value_integer(SIMULATION_TYPE, 'SIMULATION_TYPE', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'SIMULATION_TYPE                 = 1'
      write(*,*)
    endif

    call read_value_integer(NOISE_TOMOGRAPHY, 'NOISE_TOMOGRAPHY', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'NOISE_TOMOGRAPHY                = 0'
      write(*,*)
    endif

    call read_value_logical(SAVE_FORWARD, 'SAVE_FORWARD', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'SAVE_FORWARD                    = .false.'
      write(*,*)
    endif

    call read_value_logical(INVERSE_FWI_FULL_PROBLEM, 'INVERSE_FWI_FULL_PROBLEM', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'INVERSE_FWI_FULL_PROBLEM        = .false.'
      write(*,*)
    endif

    call read_value_integer(UTM_PROJECTION_ZONE, 'UTM_PROJECTION_ZONE', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'UTM_PROJECTION_ZONE             = 11'
      write(*,*)
    endif

    call read_value_logical(SUPPRESS_UTM_PROJECTION, 'SUPPRESS_UTM_PROJECTION', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'SUPPRESS_UTM_PROJECTION         = .true.'
      write(*,*)
    endif

    ! total number of processors
    call read_value_integer(NPROC, 'NPROC', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'NPROC                           = 4'
      write(*,*)
    endif

    call read_value_integer(NSTEP, 'NSTEP', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'NSTEP                           = 5000'
      write(*,*)
    endif

    call read_value_double_precision(DT, 'DT', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'DT                              = 0.05'
      write(*,*)
    endif

    call read_value_logical(LTS_MODE, 'LTS_MODE', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'LTS_MODE                        = .false.'
      write(*,*)
    endif

    call read_value_integer(PARTITIONING_TYPE, 'PARTITIONING_TYPE', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'PARTITIONING_TYPE               = 1'
      write(*,*)
    endif

    !-------------------------------------------------------
    ! LDDRK time scheme
    !-------------------------------------------------------
    call read_value_logical(USE_LDDRK, 'USE_LDDRK', ier)   ! low-memory Runge-Kutta time scheme
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'USE_LDDRK                       = .false.'
      write(*,*)
    endif

    call read_value_logical(INCREASE_CFL_FOR_LDDRK, 'INCREASE_CFL_FOR_LDDRK', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'INCREASE_CFL_FOR_LDDRK          = .false.'
      write(*,*)
    endif

    call read_value_double_precision(RATIO_BY_WHICH_TO_INCREASE_IT, 'RATIO_BY_WHICH_TO_INCREASE_IT', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'RATIO_BY_WHICH_TO_INCREASE_IT   = 1.4'
      write(*,*)
    endif

    !-------------------------------------------------------
    ! Mesh
    !-------------------------------------------------------
    call read_value_integer(NGNOD, 'NGNOD', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'NGNOD                           = 8'
      write(*,*)
    endif

    call read_value_string(MODEL, 'MODEL', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'MODEL                           = default'
      write(*,*)
    endif

    call read_value_string(TOMOGRAPHY_PATH, 'TOMOGRAPHY_PATH', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'TOMOGRAPHY_PATH                 = ./DATA/tomo_files/'
      write(*,*)
    endif

    call read_value_string(SEP_MODEL_DIRECTORY, 'SEP_MODEL_DIRECTORY', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'SEP_MODEL_DIRECTORY             = ./DATA/my_SEP_model/'
      write(*,*)
    endif

    !-------------------------------------------------------

    call read_value_logical(APPROXIMATE_OCEAN_LOAD, 'APPROXIMATE_OCEAN_LOAD', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'APPROXIMATE_OCEAN_LOAD          = .false.'
      write(*,*)
    endif

    call read_value_logical(TOPOGRAPHY, 'TOPOGRAPHY', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'TOPOGRAPHY                      = .false.'
      write(*,*)
    endif

    call read_value_logical(ATTENUATION, 'ATTENUATION', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'ATTENUATION                     = .false.'
      write(*,*)
    endif

    call read_value_logical(ANISOTROPY, 'ANISOTROPY', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'ANISOTROPY                      = .false.'
      write(*,*)
    endif

    call read_value_logical(GRAVITY, 'GRAVITY', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'GRAVITY                         = .false.'
      write(*,*)
    endif

    call read_value_double_precision(ATTENUATION_f0_REFERENCE, 'ATTENUATION_f0_REFERENCE', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'ATTENUATION_f0_REFERENCE        = 0.33333d0'
      write(*,*)
    endif

    call read_value_double_precision(MIN_ATTENUATION_PERIOD, 'MIN_ATTENUATION_PERIOD', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'MIN_ATTENUATION_PERIOD          = 999999998.d0'
      write(*,*)
    endif

    call read_value_double_precision(MAX_ATTENUATION_PERIOD, 'MAX_ATTENUATION_PERIOD', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'MAX_ATTENUATION_PERIOD          = 999999999.d0'
      write(*,*)
    endif

    call read_value_logical(COMPUTE_FREQ_BAND_AUTOMATIC, 'COMPUTE_FREQ_BAND_AUTOMATIC', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'COMPUTE_FREQ_BAND_AUTOMATIC     = .true.'
      write(*,*)
    endif

    call read_value_logical(USE_OLSEN_ATTENUATION, 'USE_OLSEN_ATTENUATION', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'USE_OLSEN_ATTENUATION           = .false.'
      write(*,*)
    endif

    call read_value_double_precision(OLSEN_ATTENUATION_RATIO, 'OLSEN_ATTENUATION_RATIO', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'OLSEN_ATTENUATION_RATIO         = 0.05'
      write(*,*)
    endif

    !-------------------------------------------------------
    ! Absorbing boundary conditions
    !-------------------------------------------------------
    call read_value_logical(PML_CONDITIONS, 'PML_CONDITIONS', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'PML_CONDITIONS                  = .false.'
      write(*,*)
    endif

    call read_value_logical(PML_INSTEAD_OF_FREE_SURFACE, 'PML_INSTEAD_OF_FREE_SURFACE', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'PML_INSTEAD_OF_FREE_SURFACE     = .false.'
      write(*,*)
    endif

    call read_value_double_precision(f0_FOR_PML, 'f0_FOR_PML', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'f0_FOR_PML                      = 0.05555'
      write(*,*)
    endif

    call read_value_logical(STACEY_ABSORBING_CONDITIONS, 'STACEY_ABSORBING_CONDITIONS', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'STACEY_ABSORBING_CONDITIONS     = .false.'
      write(*,*)
    endif

    call read_value_logical(STACEY_INSTEAD_OF_FREE_SURFACE, 'STACEY_INSTEAD_OF_FREE_SURFACE', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'STACEY_INSTEAD_OF_FREE_SURFACE  = .false.'
      write(*,*)
    endif

    call read_value_logical(BOTTOM_FREE_SURFACE, 'BOTTOM_FREE_SURFACE', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'BOTTOM_FREE_SURFACE             = .false.'
      write(*,*)
    endif

    call read_value_logical(UNDO_ATTENUATION_AND_OR_PML, 'UNDO_ATTENUATION_AND_OR_PML', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'UNDO_ATTENUATION_AND_OR_PML     = .false.'
      write(*,*)
    endif

    call read_value_integer(NT_DUMP_ATTENUATION, 'NT_DUMP_ATTENUATION', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'NT_DUMP_ATTENUATION             = 500'
      write(*,*)
    endif

    !-------------------------------------------------------
    ! Visualization
    !-------------------------------------------------------
    call read_value_logical(CREATE_SHAKEMAP, 'CREATE_SHAKEMAP', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'CREATE_SHAKEMAP                 = .false.'
      write(*,*)
    endif

    call read_value_logical(MOVIE_SURFACE, 'MOVIE_SURFACE', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'MOVIE_SURFACE                   = .false.'
      write(*,*)
    endif

    call read_value_integer(MOVIE_TYPE, 'MOVIE_TYPE', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'MOVIE_TYPE                      = 1'
      write(*,*)
    endif

    call read_value_logical(MOVIE_VOLUME, 'MOVIE_VOLUME', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'MOVIE_VOLUME                    = .false.'
      write(*,*)
    endif

    call read_value_logical(SAVE_DISPLACEMENT, 'SAVE_DISPLACEMENT', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'SAVE_DISPLACEMENT               = .false.'
      write(*,*)
    endif

    call read_value_logical(USE_HIGHRES_FOR_MOVIES, 'USE_HIGHRES_FOR_MOVIES', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'USE_HIGHRES_FOR_MOVIES          = .false.'
      write(*,*)
    endif

    call read_value_integer(NTSTEP_BETWEEN_FRAMES, 'NTSTEP_BETWEEN_FRAMES', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'NTSTEP_BETWEEN_FRAMES           = 200'
      write(*,*)
    endif

    call read_value_double_precision(HDUR_MOVIE, 'HDUR_MOVIE', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'HDUR_MOVIE                      = 0.0'
      write(*,*)
    endif

    call read_value_logical(SAVE_MESH_FILES, 'SAVE_MESH_FILES', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'SAVE_MESH_FILES                 = .false.'
      write(*,*)
    endif

    call read_value_string(LOCAL_PATH, 'LOCAL_PATH', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'LOCAL_PATH                      = OUTPUT_FILES/DATABASES_MPI'
      write(*,*)
    endif

    call read_value_integer(NTSTEP_BETWEEN_OUTPUT_INFO, 'NTSTEP_BETWEEN_OUTPUT_INFO', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'NTSTEP_BETWEEN_OUTPUT_INFO      = 500'
      write(*,*)
    endif

    !-------------------------------------------------------
    ! Sources
    !-------------------------------------------------------
    call read_value_logical(USE_SOURCES_RECEIVERS_Z, 'USE_SOURCES_RECEIVERS_Z', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'USE_SOURCES_RECEIVERS_Z         = .false.'
      write(*,*)
    endif

    call read_value_logical(USE_FORCE_POINT_SOURCE, 'USE_FORCE_POINT_SOURCE', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'USE_FORCE_POINT_SOURCE          = .false.'
      write(*,*)
    endif

    call read_value_logical(USE_RICKER_TIME_FUNCTION, 'USE_RICKER_TIME_FUNCTION', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'USE_RICKER_TIME_FUNCTION        = .false.'
      write(*,*)
    endif

    call read_value_logical(USE_EXTERNAL_SOURCE_FILE, 'USE_EXTERNAL_SOURCE_FILE', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'USE_EXTERNAL_SOURCE_FILE        = .false.'
      write(*,*)
    endif

    call read_value_logical(PRINT_SOURCE_TIME_FUNCTION, 'PRINT_SOURCE_TIME_FUNCTION', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'PRINT_SOURCE_TIME_FUNCTION      = .false.'
      write(*,*)
    endif

    !-------------------------------------------------------
    ! Seismograms
    !-------------------------------------------------------
    call read_value_integer(NTSTEP_BETWEEN_OUTPUT_SEISMOS, 'NTSTEP_BETWEEN_OUTPUT_SEISMOS', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'NTSTEP_BETWEEN_OUTPUT_SEISMOS   = 10000'
      write(*,*)
    endif

    call read_value_logical(SAVE_SEISMOGRAMS_DISPLACEMENT, 'SAVE_SEISMOGRAMS_DISPLACEMENT', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'SAVE_SEISMOGRAMS_DISPLACEMENT   = .true.'
      write(*,*)
    endif

    call read_value_logical(SAVE_SEISMOGRAMS_VELOCITY, 'SAVE_SEISMOGRAMS_VELOCITY', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'SAVE_SEISMOGRAMS_VELOCITY       = .false.'
      write(*,*)
    endif

    call read_value_logical(SAVE_SEISMOGRAMS_ACCELERATION, 'SAVE_SEISMOGRAMS_ACCELERATION', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'SAVE_SEISMOGRAMS_ACCELERATION   = .false.'
      write(*,*)
    endif

    call read_value_logical(SAVE_SEISMOGRAMS_PRESSURE, 'SAVE_SEISMOGRAMS_PRESSURE', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'SAVE_SEISMOGRAMS_PRESSURE       = .false.'
      write(*,*)
    endif

    call read_value_logical(SAVE_SEISMOGRAMS_IN_ADJOINT_RUN, 'SAVE_SEISMOGRAMS_IN_ADJOINT_RUN', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'SAVE_SEISMOGRAMS_IN_ADJOINT_RUN = .false.'
      write(*,*)
    endif

    call read_value_integer(subsamp_seismos, 'subsamp_seismos', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'subsamp_seismos                 = 1'
      write(*,*)
    endif

    call read_value_logical(USE_BINARY_FOR_SEISMOGRAMS, 'USE_BINARY_FOR_SEISMOGRAMS', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'USE_BINARY_FOR_SEISMOGRAMS      = .false.'
      write(*,*)
    endif

    call read_value_logical(SU_FORMAT, 'SU_FORMAT', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'SU_FORMAT                       = .false.'
      write(*,*)
    endif

    call read_value_logical(ASDF_FORMAT, 'ASDF_FORMAT', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'ASDF_FORMAT                       = .false.'
      write(*,*)
    endif

    call read_value_logical(WRITE_SEISMOGRAMS_BY_MAIN, 'WRITE_SEISMOGRAMS_BY_MAIN', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'WRITE_SEISMOGRAMS_BY_MAIN     = .false.'
      write(*,*)
    endif

    call read_value_logical(SAVE_ALL_SEISMOS_IN_ONE_FILE, 'SAVE_ALL_SEISMOS_IN_ONE_FILE', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'SAVE_ALL_SEISMOS_IN_ONE_FILE    = .false.'
      write(*,*)
    endif

    call read_value_logical(USE_TRICK_FOR_BETTER_PRESSURE, 'USE_TRICK_FOR_BETTER_PRESSURE', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'USE_TRICK_FOR_BETTER_PRESSURE   = .false.'
      write(*,*)
    endif

    !-------------------------------------------------------
    ! Source encoding
    !-------------------------------------------------------
    call read_value_logical(USE_SOURCE_ENCODING, 'USE_SOURCE_ENCODING', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'USE_SOURCE_ENCODING             = .false.'
      write(*,*)
    endif

    !-------------------------------------------------------
    ! Energy calculation
    !-------------------------------------------------------
    call read_value_logical(OUTPUT_ENERGY, 'OUTPUT_ENERGY', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'OUTPUT_ENERGY                   = .false.'
      write(*,*)
    endif

    call read_value_integer(NTSTEP_BETWEEN_OUTPUT_ENERGY, 'NTSTEP_BETWEEN_OUTPUT_ENERGY', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'NTSTEP_BETWEEN_OUTPUT_ENERGY    = 10'
      write(*,*)
    endif

    !-------------------------------------------------------
    ! Adjoint kernel outputs
    !-------------------------------------------------------
    call read_value_integer(NTSTEP_BETWEEN_READ_ADJSRC, 'NTSTEP_BETWEEN_READ_ADJSRC', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'NTSTEP_BETWEEN_READ_ADJSRC      = 0'
      write(*,*)
    endif

    call read_value_logical(READ_ADJSRC_ASDF, 'READ_ADJSRC_ASDF', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'READ_ADJSRC_ASDF               = .false.'
      write(*,*)
    endif

    call read_value_logical(ANISOTROPIC_KL, 'ANISOTROPIC_KL', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'ANISOTROPIC_KL                  = .false.'
      write(*,*)
    endif

    call read_value_logical(SAVE_TRANSVERSE_KL, 'SAVE_TRANSVERSE_KL', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'SAVE_TRANSVERSE_KL              = .false.'
      write(*,*)
    endif

    call read_value_logical(ANISOTROPIC_VELOCITY_KL, 'ANISOTROPIC_VELOCITY_KL', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'ANISOTROPIC_VELOCITY_KL         = .false.'
      write(*,*)
    endif

    call read_value_logical(APPROXIMATE_HESS_KL, 'APPROXIMATE_HESS_KL', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'APPROXIMATE_HESS_KL             = .false.'
      write(*,*)
    endif

    call read_value_logical(SAVE_MOHO_MESH, 'SAVE_MOHO_MESH', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'SAVE_MOHO_MESH                  = .false.'
      write(*,*)
    endif

    !-------------------------------------------------------

    call read_value_logical(COUPLE_WITH_INJECTION_TECHNIQUE, 'COUPLE_WITH_INJECTION_TECHNIQUE', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'COUPLE_WITH_INJECTION_TECHNIQUE = .false.'
      write(*,*)
    endif

    call read_value_integer(INJECTION_TECHNIQUE_TYPE,'INJECTION_TECHNIQUE_TYPE',ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'INJECTION_TECHNIQUE_TYPE        = 3'
      write(*,*)
    endif

    call read_value_logical(MESH_A_CHUNK_OF_THE_EARTH,'MESH_A_CHUNK_OF_THE_EARTH',ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'MESH_A_CHUNK_OF_THE_EARTH       = .false.'
      write(*,*)
    endif

    call read_value_string(TRACTION_PATH, 'TRACTION_PATH', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'TRACTION_PATH                   = ./DATA/AxiSEM_tractions/3/'
      write(*,*)
    endif

    ! this one is only used when INJECTION_TECHNIQUE_IS_FK, but we read it anyway because we later broadcast it to the other nodes
    call read_value_string(FKMODEL_FILE,'FKMODEL_FILE',ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'FKMODEL_FILE                    = FKmodel'
      write(*,*)
    endif

    call read_value_logical(RECIPROCITY_AND_KH_INTEGRAL,'RECIPROCITY_AND_KH_INTEGRAL',ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'RECIPROCITY_AND_KH_INTEGRAL     = .false.'
      write(*,*)
    endif

    ! check the type of external code to couple with, if any
    !if (MESH_A_CHUNK_OF_THE_EARTH .and. .not. COUPLE_WITH_INJECTION_TECHNIQUE) &
    !  stop 'MESH_A_CHUNK_OF_THE_EARTH only available with COUPLE_WITH_INJECTION_TECHNIQUE for now, easy to change but not done yet'

    if (COUPLE_WITH_INJECTION_TECHNIQUE) then
      if (INJECTION_TECHNIQUE_TYPE /= INJECTION_TECHNIQUE_IS_DSM .and. &
         INJECTION_TECHNIQUE_TYPE /= INJECTION_TECHNIQUE_IS_AXISEM .and. &
         INJECTION_TECHNIQUE_TYPE /= INJECTION_TECHNIQUE_IS_FK) stop 'Error incorrect value of INJECTION_TECHNIQUE_TYPE read'

      if ( (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_DSM ) .and. &
           (.not. MESH_A_CHUNK_OF_THE_EARTH) ) stop 'Error, coupling with DSM only works with a Earth chunk mesh'

      if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_FK .and. MESH_A_CHUNK_OF_THE_EARTH) &
           stop 'Error: coupling with F-K is for models with a flat surface (Earth flattening), &
                       &thus turn MESH_A_CHUNK_OF_THE_EARTH off'

      if ((INJECTION_TECHNIQUE_TYPE /= INJECTION_TECHNIQUE_IS_AXISEM) .and. RECIPROCITY_AND_KH_INTEGRAL) &
           stop 'Error: the use of RECIPROCITY_AND_KH_INTEGRAL is only available for coupling with AxiSEM for now'
    endif

    !-------------------------------------------------------

    ! for simultaneous runs from the same batch job
    call read_value_integer(NUMBER_OF_SIMULTANEOUS_RUNS, 'NUMBER_OF_SIMULTANEOUS_RUNS', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'NUMBER_OF_SIMULTANEOUS_RUNS     = 1'
      write(*,*)
    endif

    call read_value_logical(BROADCAST_SAME_MESH_AND_MODEL, 'BROADCAST_SAME_MESH_AND_MODEL', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'BROADCAST_SAME_MESH_AND_MODEL   = .true.'
      write(*,*)
    endif

    !-------------------------------------------------------

    call read_value_logical(GPU_MODE, 'GPU_MODE', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'GPU_MODE                        = .false.'
      write(*,*)
    endif

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
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'ADIOS_ENABLED                   = .false.'
      write(*,*)
    endif

    call read_value_logical(ADIOS_FOR_DATABASES, 'ADIOS_FOR_DATABASES', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'ADIOS_FOR_DATABASES             = .false.'
      write(*,*)
    endif

    call read_value_logical(ADIOS_FOR_MESH, 'ADIOS_FOR_MESH', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'ADIOS_FOR_MESH                  = .false.'
      write(*,*)
    endif

    call read_value_logical(ADIOS_FOR_FORWARD_ARRAYS, 'ADIOS_FOR_FORWARD_ARRAYS', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'ADIOS_FOR_FORWARD_ARRAYS        = .false.'
      write(*,*)
    endif

    call read_value_logical(ADIOS_FOR_KERNELS, 'ADIOS_FOR_KERNELS', ier)
    if (ier /= 0) then
      some_parameters_missing_from_Par_file = .true.
      write(*,'(a)') 'ADIOS_FOR_KERNELS               = .false.'
      write(*,*)
    endif

    ! closes parameter file
    call close_parameter_file()

    if (some_parameters_missing_from_Par_file) then
      write(*,*)
      write(*,*) 'All the above parameters are missing from your Par_file.'
      write(*,*) 'Please cut and paste them somewhere in your Par_file (any place is fine), change their values if needed'
      write(*,*) '(the above values are just default values), and restart your run.'
      write(*,*)
      stop 'Error: some parameters are missing in your Par_file, it is incomplete or in an older format, &
         &see at the end of the standard output file of the run for detailed and easy instructions about how to fix that'
    endif

    ! re-sets attenuation flags
    !if (.not. ATTENUATION) then
    !  ! turns off UNDO_ATTENUATION when ATTENUATION is off in the Par_file
    !  UNDO_ATTENUATION_AND_OR_PML = .false.
    !endif
    ! for pure forward simulation, no need to store undo_attenuation arrays; uses default iteration routine
    !if (SIMULATION_TYPE == 1 .and. .not. SAVE_FORWARD) then
    !  UNDO_ATTENUATION_AND_OR_PML = .false.
    !endif

    ! re-sets ADIOS flags
    if (.not. ADIOS_ENABLED) then
      ADIOS_FOR_DATABASES = .false.
      ADIOS_FOR_MESH = .false.
      ADIOS_FOR_FORWARD_ARRAYS = .false.
      ADIOS_FOR_KERNELS = .false.
      ! ADIOS_FOR_UNDO_ATTENUATION = .false. ! not implemented yet
    endif

    ! re-sets PML free surface flag
    ! PML absorbing free surface must have also PML turned on
    if (.not. PML_CONDITIONS) then
      PML_INSTEAD_OF_FREE_SURFACE = .false.
    endif

    ! re-sets stacey free surface flag
    ! stacey absorbing free surface must have also stacey turned on
    if (.not. STACEY_ABSORBING_CONDITIONS) then
      STACEY_INSTEAD_OF_FREE_SURFACE = .false.
    endif

    ! checks parameter settings
    call check_simulation_parameters()

    ! computes additional parameters depending on setting
    call read_compute_parameters()

  endif ! of if (myrank == 0) then

  ! read from a single processor (the main) and then use MPI to broadcast to others
  ! to avoid an I/O bottleneck in the case of very large runs
  if (BROADCAST_AFTER_READ) call broadcast_computed_parameters()

  ! Cray compilers
#if _CRAYFTN
#warning "Warning: using Cray compiler assign function for un-compressed file output"
  ! Cray uses compressed formats by default for list-directed output, for example:
  !     write(*,*) 10,10            leads to output -> 2*10           compressed, instead of:       10        10
  !     write(*,*) 1,1.78e-5                        -> 1,   1.78e-5   comma delimiter, instead of:   1         1.78e-5
  ! this leads to problems when writing seismograms in ASCII-format.
  ! to circumvent this behaviour, one can use cray's assign environment:
  !   $ setenv FILENV ASGTMP
  !   $ assign -U on g:all        (g:all  - all file open requests)
  ! see: https://pubs.cray.com/bundle/Cray_Fortran_Reference_Manual_100_S-3901_Fortran_ditaval.xml/..
  !             ..page/Cray_Fortran_Implementation_Specifics.html
  ! or use the Fortran statement here below:

  !debug
  !if (myrank == 0) print *,'...compiled by Cray compilers'

  ! assigns -U (uncompressed format) for all subsequent file opens
  ! includes seismograms, but not the already opened IMAIN file output
  call assign('assign -U on g:all',ier)
#endif

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

  ! period band over which we mimic a constant Q factor for attenuation
  if (.not. COMPUTE_FREQ_BAND_AUTOMATIC .and. MIN_ATTENUATION_PERIOD >= MAX_ATTENUATION_PERIOD) &
    stop 'must have MIN_ATTENUATION_PERIOD < MAX_ATTENUATION_PERIOD'

  ! seismogram output
  if (.not. SAVE_SEISMOGRAMS_DISPLACEMENT .and. .not. SAVE_SEISMOGRAMS_VELOCITY .and. &
     .not. SAVE_SEISMOGRAMS_ACCELERATION .and. .not. SAVE_SEISMOGRAMS_PRESSURE) &
   stop 'Error: at least one of SAVE_SEISMOGRAMS_DISPLACEMENT SAVE_SEISMOGRAMS_VELOCITY SAVE_SEISMOGRAMS_ACCELERATION &
             &SAVE_SEISMOGRAMS_PRESSURE must be true'

  if (subsamp_seismos < 1) &
    stop 'Error: subsamp_seismos must be >= 1'

  ! this could be implemented in the future if needed,
  ! see comments in the source code around the USE_TRICK_FOR_BETTER_PRESSURE
  ! option (use a "grep" command to find them) to see how this could/should be done
  if (USE_TRICK_FOR_BETTER_PRESSURE .and. (SAVE_SEISMOGRAMS_DISPLACEMENT .or. SAVE_SEISMOGRAMS_VELOCITY .or. &
        SAVE_SEISMOGRAMS_ACCELERATION)) &
    stop 'USE_TRICK_FOR_BETTER_PRESSURE is currently incompatible with &
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

  ! UNDO_ATT
  !if (UNDO_ATTENUATION_AND_OR_PML .and. GPU_MODE) &
  !  stop 'For GPU_MODE, UNDO_ATTENUATION_AND_OR_PML is not implemented in this code yet'

  ! attenuation for backward simulation
  if (SIMULATION_TYPE == 3 .and. ATTENUATION .and. .not. UNDO_ATTENUATION_AND_OR_PML) &
    stop 'For SIMULATION_TYPE == 3 and attenuation, simulations need flag UNDO_ATTENUATION_AND_OR_PML set to .true.'

  ! local time stepping (still in experimental implementation phase...)
  if (LTS_MODE) then
    if (USE_LDDRK) &
      stop 'USE_LDDRK in LTS_MODE not supported yet'
    if (PML_CONDITIONS) &
      stop 'PML_CONDITIONS in LTS_MODE not supported yet'
    if (GPU_MODE) &
      stop 'GPU_MODE and LTS_MODE together not supported yet'
  endif

  ! PARTITIONING_TYPE
  if (PARTITIONING_TYPE < 1 .or. PARTITIONING_TYPE > 4) &
    stop 'PARTITIONING_TYPE must be 1,2,3 or 4 (for SCOTCH, METIS, PATOH or ROW_PARTS partitioner)'

  ! Warnings

  ! ADIOS is very useful for very large simulations (say using 2000 MPI tasks or more)
  ! but slows down the code if used for simulations that are small or medium size, because of the overhead any library has.
  if (ADIOS_ENABLED .and. NPROC < 2000) then
    print *
    print *,'**************'
    print *,'**************'
    print *,'ADIOS significantly slows down small or medium-size runs, which is the case here, please consider turning it off'
    print *,'**************'
    print *,'**************'
    print *
  endif

  end subroutine check_simulation_parameters

!
!-------------------------------------------------------------------------------------------------
!


  subroutine read_compute_parameters()

! computes additional parameters
! (only executed by main process)

  use constants
  use shared_parameters

  implicit none

  ! local parameters
  integer :: i,irange
  character(len=MAX_STRING_LEN) :: sources_filename
  character(len=MAX_STRING_LEN) :: path_to_add
  character(len=MAX_STRING_LEN) :: tmp_TOMOGRAPHY_PATH,tmp_LOCAL_PATH

  ! updates values for simulation settings
  ! LDDRK scheme
  if (USE_LDDRK .and. INCREASE_CFL_FOR_LDDRK) DT = DT * RATIO_BY_WHICH_TO_INCREASE_IT

  ! to be consistent with external source file option turned to on
  if (USE_EXTERNAL_SOURCE_FILE) USE_RICKER_TIME_FUNCTION = .false.

  ! see if we are running several independent runs in parallel
  ! if so, add the right directory for that run
  ! (group numbers start at zero, but directory names start at run0001, thus we add one)
  ! a negative value for "mygroup" is a convention that indicates that groups (i.e. sub-communicators, one per run) are off
  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then

!! DK DK remove leading ./ if any, Paul Cristini said it could lead to problems when NUMBER_OF_SIMULTANEOUS_RUNS > 1
    tmp_LOCAL_PATH = adjustl(LOCAL_PATH)
    if (index (tmp_LOCAL_PATH, './') == 1) then
      LOCAL_PATH = tmp_LOCAL_PATH(3:)
    endif

    tmp_TOMOGRAPHY_PATH = adjustl(TOMOGRAPHY_PATH)
    if (index (tmp_TOMOGRAPHY_PATH, './') == 1) then
      TOMOGRAPHY_PATH = tmp_TOMOGRAPHY_PATH(3:)
    endif

    TRACTION_PATH_new = adjustl(TRACTION_PATH)
    if (index (TRACTION_PATH_new, './') == 1) then
      TRACTION_PATH = TRACTION_PATH_new(3:)
    endif

    write(path_to_add,"('run',i4.4,'/')") mygroup + 1

    LOCAL_PATH = path_to_add(1:len_trim(path_to_add))//LOCAL_PATH(1:len_trim(LOCAL_PATH))
    TOMOGRAPHY_PATH = path_to_add(1:len_trim(path_to_add))//TOMOGRAPHY_PATH(1:len_trim(TOMOGRAPHY_PATH))
    TRACTION_PATH = path_to_add(1:len_trim(path_to_add))//TRACTION_PATH(1:len_trim(TRACTION_PATH))
  endif

  ! noise simulations:
  if (NOISE_TOMOGRAPHY /= 0) then
    ! double the number of time steps, if running noise simulations (+/- branches)
    NSTEP = 2 * NSTEP - 1

    ! for noise simulations, we need to save movies at the surface (where the noise is generated)
    ! and thus we force MOVIE_SURFACE in the mesher to be .true., in order to use variables defined for surface movies later
    ! the noise surface wavefield will be stored/load in files DATABASES_MPI/proc***_surface_movie
    !
    ! for the solver, setting MOVIE_SURFACE can be used to plot/visualize the wavefield, but is not needed for noise simulations.
    ! however, noise simulations require movie type 1 and highres
    ! defaults
    MOVIE_TYPE = 1                      ! 1 == only top surface (no side/bottom faces)
    USE_HIGHRES_FOR_MOVIES = .true.     ! we need to save surface movie everywhere, i.e. at all GLL points on the surface
    ! let user decide
    !MOVIE_SURFACE = .true.             ! (not necessary) to store/load generating wavefield for plotting
    !SAVE_DISPLACEMENT = .true.         ! (not necessary) stores displacement (flag not necessary, but to avoid confusion)
  endif

  ! make sure NSTEP is a multiple of subsamp_seismos
  ! if not, increase it a little bit, to the next multiple
  if (mod(NSTEP,subsamp_seismos) /= 0) then
    if (NOISE_TOMOGRAPHY /= 0) then
      if (myrank == 0) then
        print *,'Noise simulation: Invalid number of NSTEP = ',NSTEP
        print *,'Must be a multiple of subsamp_seismos = ',subsamp_seismos
      endif
      stop 'Error: NSTEP must be a multiple of subsamp_seismos'
    else
      NSTEP = (NSTEP/subsamp_seismos + 1)*subsamp_seismos
      ! user output
      if (myrank == 0) then
        print *
        print *,'NSTEP is not a multiple of subsamp_seismos'
        print *,'thus increasing it automatically to the next multiple, which is ',NSTEP
        print *
      endif
    endif
  endif

  ! output seismograms at least once at the end of the simulation
  NTSTEP_BETWEEN_OUTPUT_SEISMOS = min(NSTEP,NTSTEP_BETWEEN_OUTPUT_SEISMOS)

  ! make sure NSTEP_BETWEEN_OUTPUT_SEISMOS is a multiple of subsamp_seismos
  if (mod(NTSTEP_BETWEEN_OUTPUT_SEISMOS,subsamp_seismos) /= 0) then
    if (myrank == 0) then
      print *,'Invalid number of NTSTEP_BETWEEN_OUTPUT_SEISMOS = ',NTSTEP_BETWEEN_OUTPUT_SEISMOS
      print *,'Must be a multiple of subsamp_seismos = ',subsamp_seismos
    endif
    stop 'Error: NTSTEP_BETWEEN_OUTPUT_SEISMOS must be a multiple of subsamp_seismos'
  endif

  ! the default value of NTSTEP_BETWEEN_READ_ADJSRC (0) is to read the whole trace at the same time
  if (NTSTEP_BETWEEN_READ_ADJSRC == 0)  NTSTEP_BETWEEN_READ_ADJSRC = NSTEP

  ! total times steps must be dividable by adjoint source chunks/blocks
  if (mod(NSTEP,NTSTEP_BETWEEN_READ_ADJSRC) /= 0) then
    print *,'Error: NSTEP ',NSTEP,' not a multiple of NTSTEP_BETWEEN_READ_ADJSRC ',NTSTEP_BETWEEN_READ_ADJSRC
    print *,'       Please change NTSTEP_BETWEEN_READ_ADJSRC in the Par_file!'
    print *,'       (in case NOISE_TOMOGRAPHY is not equal to zero, NSTEP from Par_file becomes 2*NSTEP-1)'
    stop 'Error: mod(NSTEP,NTSTEP_BETWEEN_READ_ADJSRC) must be zero! Please modify Par_file and rerun solver'
  endif

  ! checks number of nodes for 2D and 3D shape functions for quadrilaterals and hexahedra
  ! curvature (i.e. HEX27 elements) is not handled by our internal mesher, for that use CUBIT/Trelis or Gmsh for instance
  if (NGNOD == 8) then
    NGNOD2D = 4
  else if (NGNOD == 27) then
    NGNOD2D = 9
  else if (NGNOD /= 8 .and. NGNOD /= 27) then
    stop 'Error elements should have 8 or 27 control nodes, please modify NGNOD in Par_file and recompile solver'
  endif

  ! get the name of the file describing the sources
  if (USE_FORCE_POINT_SOURCE) then
    sources_filename = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'FORCESOLUTION'
  else
    sources_filename = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'CMTSOLUTION'
  endif
  ! see if we are running several independent runs in parallel
  ! if so, add the right directory for that run
  ! (group numbers start at zero, but directory names start at run0001, thus we add one)
  ! a negative value for "mygroup" is a convention that indicates that groups (i.e. sub-communicators, one per run) are off
  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
    write(path_to_add,"('run',i4.4,'/')") mygroup + 1
    sources_filename = path_to_add(1:len_trim(path_to_add))//sources_filename(1:len_trim(sources_filename))
  endif

  ! determines number of sources depending on number of lines in sources file
  call get_number_of_sources(sources_filename)

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
  case ('coupled')
    IMODEL = IMODEL_COUPLED
  case default
    print *
    print *,'********** model not recognized: ',trim(MODEL),' **************'
    ! requires a valid name
    stop 'Invalid MODEL name'
    ! allows to continue with default
    !print *,'********** using model: default',' **************'
    !print *
    !IMODEL = IMODEL_DEFAULT
  end select

  ! check
  if (IMODEL == IMODEL_IPATI .or. IMODEL == IMODEL_IPATI_WATER) then
    if (USE_RICKER_TIME_FUNCTION .eqv. .false.) &
      stop 'Error for IPATI model, please set USE_RICKER_TIME_FUNCTION to .true. in Par_file and recompile solver'
  endif

  end subroutine read_compute_parameters

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_number_of_sources(sources_filename)

! determines number of sources depending on number of lines in source file
! (only executed by main process)

  use constants, only: IIN,IIN_PAR,IN_DATA_FILES,HUGEVAL,TINYVAL, &
    NLINES_PER_CMTSOLUTION_SOURCE,NLINES_PER_FORCESOLUTION_SOURCE

  use shared_parameters

  implicit none

  character(len=MAX_STRING_LEN),intent(in) :: sources_filename

  ! local variables
  integer :: icounter,isource,idummy,ier,nlines_per_source
  double precision :: hdur, minval_hdur
  character(len=MAX_STRING_LEN) :: dummystring

  ! initializes
  NSOURCES = 0

  ! checks if finite fault source
  open(unit=IIN_PAR,file=IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'Par_file_faults',status='old',iostat=ier)
  if (ier == 0) then
    HAS_FINITE_FAULT_SOURCE = .true.
    !write(IMAIN,*) 'provides finite faults'
    close(IIN_PAR)
  else
    HAS_FINITE_FAULT_SOURCE = .false.
  endif

  ! gets number of point sources
  if (USE_FORCE_POINT_SOURCE) then
    ! compute the total number of sources in the FORCESOLUTION file
    ! there are NLINES_PER_FORCESOLUTION_SOURCE lines per source in that file
    open(unit=IIN,file=trim(sources_filename),status='old',action='read',iostat=ier)
    if (ier /= 0) then
      if (HAS_FINITE_FAULT_SOURCE) then
        ! no need for FORCESOLUTION file
        return
      else
        stop 'Error opening FORCESOLUTION file'
      endif
    endif
    !write(IMAIN,*) 'provides force solution'

    icounter = 0
    do while (ier == 0)
      read(IIN,"(a)",iostat=ier) dummystring
      if (ier == 0) icounter = icounter + 1
    enddo
    close(IIN)

    ! number of lines for source description
    if (USE_EXTERNAL_SOURCE_FILE) then
      !! VM VM in case of USE_EXTERNAL_SOURCE_FILE we have to read one additional line per source (the name of external source file)
      nlines_per_source = NLINES_PER_FORCESOLUTION_SOURCE + 1
    else
      nlines_per_source = NLINES_PER_FORCESOLUTION_SOURCE
    endif

    ! checks lines are a multiple
    if (mod(icounter,nlines_per_source) /= 0) then
      print *,'Error: total number of lines in FORCESOLUTION file should be a multiple of ',nlines_per_source
      stop 'Error total number of lines in FORCESOLUTION file should be a multiple of NLINES_PER_FORCESOLUTION_SOURCE'
    endif

    ! number of sources in file
    NSOURCES = icounter / nlines_per_source

    ! checks if any
    if (NSOURCES < 1) stop 'Error need at least one source in FORCESOLUTION file'

  else
    ! compute the total number of sources in the CMTSOLUTION file
    ! there are NLINES_PER_CMTSOLUTION_SOURCE lines per source in that file
    open(unit=IIN,file=trim(sources_filename),status='old',action='read',iostat=ier)
    if (ier /= 0) then
      if (HAS_FINITE_FAULT_SOURCE) then
        ! no need for CMTSOLUTION file
        return
      else
        stop 'Error opening CMTSOLUTION file'
      endif
    endif
    !write(IMAIN,*) 'provides CMT solution'

    icounter = 0
    do while (ier == 0)
      read(IIN,"(a)",iostat=ier) dummystring
      if (ier == 0) icounter = icounter + 1
    enddo
    close(IIN)

    ! number of lines for source description
    if (USE_EXTERNAL_SOURCE_FILE) then
      !! VM VM in case of USE_EXTERNAL_SOURCE_FILE we have to read one additional line per source (the name of external source file)
      nlines_per_source = NLINES_PER_CMTSOLUTION_SOURCE + 1
    else
      nlines_per_source = NLINES_PER_CMTSOLUTION_SOURCE
    endif

    ! checks number of lines
    if (mod(icounter,nlines_per_source) /= 0) then
      print *,'Error: total number of lines in CMTSOLUTION file should be a multiple of ',nlines_per_source
      stop 'Error total number of lines in CMTSOLUTION file should be a multiple of NLINES_PER_CMTSOLUTION_SOURCE'
    endif

    ! number of sources in file
    NSOURCES = icounter / nlines_per_source

    ! checks if any
    if (NSOURCES < 1) stop 'Error need at least one source in CMTSOLUTION file'

    ! compute the minimum value of hdur in CMTSOLUTION file
    open(unit=IIN,file=trim(sources_filename),status='old',action='read')
    minval_hdur = HUGEVAL
    do isource = 1,NSOURCES

      ! skip other information
      do idummy = 1,3
        read(IIN,"(a)") dummystring
      enddo

      ! read half duration and compute minimum
      read(IIN,"(a)") dummystring
      read(dummystring(15:len_trim(dummystring)),*) hdur
      minval_hdur = min(minval_hdur,hdur)

      ! reads till the end of this source
      do idummy = 5,nlines_per_source
        read(IIN,"(a)") dummystring
      enddo

    enddo
    close(IIN)

    ! one cannot use a Heaviside source for the movies
    if ((MOVIE_SURFACE .or. MOVIE_VOLUME) .and. sqrt(minval_hdur**2 + HDUR_MOVIE**2) < TINYVAL) &
      stop 'Error hdur too small for movie creation, movies do not make sense for Heaviside source'
  endif

  end subroutine get_number_of_sources

!
!-------------------------------------------------------------------------------------------------
!


  subroutine broadcast_computed_parameters()

! broadcasts parameters to all processes

  use constants
  use shared_parameters

  implicit none

  ! broadcasts setting
  ! simulation parameters
  call bcast_all_singlei(SIMULATION_TYPE)
  call bcast_all_singlei(NOISE_TOMOGRAPHY)
  call bcast_all_singlel(SAVE_FORWARD)

  call bcast_all_singlel(INVERSE_FWI_FULL_PROBLEM)
  call bcast_all_singlei(UTM_PROJECTION_ZONE)
  call bcast_all_singlel(SUPPRESS_UTM_PROJECTION)

  call bcast_all_singlei(NPROC)
  call bcast_all_singlei(NSTEP)
  call bcast_all_singledp(DT)

  call bcast_all_singlel(LTS_MODE)
  call bcast_all_singlei(PARTITIONING_TYPE)

  ! LDDRK
  call bcast_all_singlel(USE_LDDRK)
  ! call bcast_all_singlel(INCREASE_CFL_FOR_LDDRK) ! not needed any further
  ! call bcast_all_singledp(RATIO_BY_WHICH_TO_INCREASE_IT) ! not needed any further

  ! mesh
  call bcast_all_singlei(NGNOD)
  call bcast_all_string(MODEL)

  call bcast_all_string(TOMOGRAPHY_PATH)
  call bcast_all_string(SEP_MODEL_DIRECTORY)

  call bcast_all_singlel(APPROXIMATE_OCEAN_LOAD)
  call bcast_all_singlel(TOPOGRAPHY)
  call bcast_all_singlel(ATTENUATION)
  call bcast_all_singlel(ANISOTROPY)
  call bcast_all_singlel(GRAVITY)

  call bcast_all_singledp(ATTENUATION_f0_REFERENCE)
  call bcast_all_singledp(MIN_ATTENUATION_PERIOD)
  call bcast_all_singledp(MAX_ATTENUATION_PERIOD)
  call bcast_all_singlel(COMPUTE_FREQ_BAND_AUTOMATIC)

  call bcast_all_singlel(USE_OLSEN_ATTENUATION)
  call bcast_all_singledp(OLSEN_ATTENUATION_RATIO)

  ! absorbing boundaries
  call bcast_all_singlel(PML_CONDITIONS)
  call bcast_all_singlel(PML_INSTEAD_OF_FREE_SURFACE)
  call bcast_all_singledp(f0_FOR_PML)

  call bcast_all_singlel(STACEY_ABSORBING_CONDITIONS)
  call bcast_all_singlel(STACEY_INSTEAD_OF_FREE_SURFACE)
  call bcast_all_singlel(BOTTOM_FREE_SURFACE)

  ! undoing att
  call bcast_all_singlel(UNDO_ATTENUATION_AND_OR_PML)
  call bcast_all_singlei(NT_DUMP_ATTENUATION)

  ! visualization
  call bcast_all_singlel(CREATE_SHAKEMAP)
  call bcast_all_singlel(MOVIE_SURFACE)
  call bcast_all_singlei(MOVIE_TYPE)
  call bcast_all_singlel(MOVIE_VOLUME)
  call bcast_all_singlel(SAVE_DISPLACEMENT)
  call bcast_all_singlel(USE_HIGHRES_FOR_MOVIES)
  call bcast_all_singlei(NTSTEP_BETWEEN_FRAMES)
  call bcast_all_singledp(HDUR_MOVIE)

  call bcast_all_singlel(SAVE_MESH_FILES)
  call bcast_all_string(LOCAL_PATH)
  call bcast_all_singlei(NTSTEP_BETWEEN_OUTPUT_INFO)

  ! sources
  call bcast_all_singlel(USE_SOURCES_RECEIVERS_Z)
  call bcast_all_singlel(USE_FORCE_POINT_SOURCE)
  call bcast_all_singlel(USE_RICKER_TIME_FUNCTION)
  call bcast_all_singlel(USE_EXTERNAL_SOURCE_FILE)
  call bcast_all_singlel(PRINT_SOURCE_TIME_FUNCTION)

  ! seismograms
  call bcast_all_singlei(NTSTEP_BETWEEN_OUTPUT_SEISMOS)
  call bcast_all_singlel(SAVE_SEISMOGRAMS_DISPLACEMENT)
  call bcast_all_singlel(SAVE_SEISMOGRAMS_VELOCITY)
  call bcast_all_singlel(SAVE_SEISMOGRAMS_ACCELERATION)
  call bcast_all_singlel(SAVE_SEISMOGRAMS_PRESSURE)
  call bcast_all_singlel(SAVE_SEISMOGRAMS_IN_ADJOINT_RUN)
  call bcast_all_singlei(subsamp_seismos)
  call bcast_all_singlel(USE_BINARY_FOR_SEISMOGRAMS)
  call bcast_all_singlel(SU_FORMAT)
  call bcast_all_singlel(ASDF_FORMAT)
  call bcast_all_singlel(WRITE_SEISMOGRAMS_BY_MAIN)
  call bcast_all_singlel(SAVE_ALL_SEISMOS_IN_ONE_FILE)
  call bcast_all_singlel(USE_TRICK_FOR_BETTER_PRESSURE)

  ! source encoding
  call bcast_all_singlel(USE_SOURCE_ENCODING)

  ! energy
  call bcast_all_singlel(OUTPUT_ENERGY)
  call bcast_all_singlei(NTSTEP_BETWEEN_OUTPUT_ENERGY)

  ! adjoint kernels
  call bcast_all_singlei(NTSTEP_BETWEEN_READ_ADJSRC)
  call bcast_all_singlel(READ_ADJSRC_ASDF)
  call bcast_all_singlel(ANISOTROPIC_KL)
  call bcast_all_singlel(SAVE_TRANSVERSE_KL)
  call bcast_all_singlel(ANISOTROPIC_VELOCITY_KL)
  call bcast_all_singlel(APPROXIMATE_HESS_KL)
  call bcast_all_singlel(SAVE_MOHO_MESH)

  ! coupling
  call bcast_all_singlel(COUPLE_WITH_INJECTION_TECHNIQUE)
  call bcast_all_singlei(INJECTION_TECHNIQUE_TYPE)
  call bcast_all_singlel(MESH_A_CHUNK_OF_THE_EARTH)
  call bcast_all_string(TRACTION_PATH)
  call bcast_all_string(FKMODEL_FILE)
  call bcast_all_singlel(RECIPROCITY_AND_KH_INTEGRAL)

  ! simultaneous runs
  call bcast_all_singlei(NUMBER_OF_SIMULTANEOUS_RUNS)
  call bcast_all_singlel(BROADCAST_SAME_MESH_AND_MODEL)

  ! GPU
  call bcast_all_singlel(GPU_MODE)

  ! ADIOS
  call bcast_all_singlel(ADIOS_ENABLED)
  call bcast_all_singlel(ADIOS_FOR_DATABASES)
  call bcast_all_singlel(ADIOS_FOR_MESH)
  call bcast_all_singlel(ADIOS_FOR_FORWARD_ARRAYS)
  call bcast_all_singlel(ADIOS_FOR_KERNELS)

  ! broadcast all parameters computed from others
  call bcast_all_singlei(IMODEL)
  call bcast_all_singlei(NGNOD2D)

  call bcast_all_singlei(NSOURCES)
  call bcast_all_singlel(HAS_FINITE_FAULT_SOURCE)

  end subroutine broadcast_computed_parameters

