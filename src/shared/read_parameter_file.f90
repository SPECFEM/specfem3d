!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
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

  subroutine read_parameter_file()

  use constants

  use shared_parameters

  implicit none

! local variables
  integer :: icounter,isource,idummy,ier
  integer :: nproc_eta_old,nproc_xi_old
  double precision :: hdur,minval_hdur
  character(len=MAX_STRING_LEN) :: dummystring

  character(len=MAX_STRING_LEN) :: MODEL
  character(len=MAX_STRING_LEN) :: CMTSOLUTION,FORCESOLUTION
  character(len=MAX_STRING_LEN) :: path_to_add

  !logical :: sep_dir_exists
  integer :: i,irange

  ! opens file Par_file
  call open_parameter_file(ier)

  ! total number of processors
  call read_value_integer(NPROC, 'NPROC', ier)
  if (ier /= 0) then
    ! checks if it's using an old Par_file format
    call read_value_integer(nproc_eta_old, 'NPROC_ETA', ier)
    if (ier /= 0) then
      print*,'please specify the number of processes in Par_file as:'
      print*,'NPROC           =    <my_number_of_desired_processes> '
      return
    endif
    ! checks if it's using an old Par_file format
    call read_value_integer(nproc_xi_old, 'NPROC_XI', ier)
    if (ier /= 0) then
      print*,'please specify the number of processes in Par_file as:'
      print*,'NPROC           =    <my_number_of_desired_processes> '
      return
    endif
    NPROC = nproc_eta_old * nproc_xi_old
  endif

  ! reads in mandatory simulation parameters
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
  call read_value_integer(NSTEP, 'NSTEP', ier)
  if (ier /= 0) stop 'Error reading Par_file parameter NSTEP'
  call read_value_double_precision(DT, 'DT', ier)
  if (ier /= 0) stop 'Error reading Par_file parameter DT'
  call read_value_integer(NGNOD, 'NGNOD', ier)
  if (ier /= 0) stop 'Error reading Par_file parameter NGNOD'
  call read_value_string(MODEL, 'MODEL', ier)
  if (ier /= 0) stop 'Error reading Par_file parameter MODEL'

  write(SEP_MODEL_DIRECTORY, '(a)') ''
  call read_value_string(SEP_MODEL_DIRECTORY, 'SEP_MODEL_DIRECTORY', ier)
  if (ier /= 0) write (0, '(a)') 'No SEP_MODEL_DIRECTORY defined in Par_file.'
  !if (ier /= 0) stop 'Error reading Par_file parameter SEP_MODEL_DIRECTORY'

  call read_value_logical(APPROXIMATE_OCEAN_LOAD, 'APPROXIMATE_OCEAN_LOAD', ier)
  if (ier /= 0) stop 'Error reading Par_file parameter APPROXIMATE_OCEAN_LOAD'
  call read_value_logical(TOPOGRAPHY, 'TOPOGRAPHY', ier)
  if (ier /= 0) stop 'Error reading Par_file parameter TOPOGRAPHY'
  call read_value_logical(ATTENUATION, 'ATTENUATION', ier)
  if (ier /= 0) stop 'Error reading Par_file parameter ATTENUATION'
  call read_value_logical(FULL_ATTENUATION_SOLID, 'FULL_ATTENUATION_SOLID', ier)
  if (ier /= 0) stop 'Error reading Par_file parameter FULL_ATTENUATION_SOLID'
  call read_value_logical(ANISOTROPY, 'ANISOTROPY', ier)
  if (ier /= 0) stop 'Error reading Par_file parameter ANISOTROPY'
  call read_value_string(TOMOGRAPHY_PATH, 'TOMOGRAPHY_PATH', ier)
  if (ier /= 0) stop 'Error reading Par_file parameter TOMOGRAPHY_PATH'
  call read_value_logical(USE_OLSEN_ATTENUATION, 'USE_OLSEN_ATTENUATION', ier)
  if (ier /= 0) stop 'Error reading Par_file parameter USE_OLSEN_ATTENUATION'
  call read_value_double_precision(OLSEN_ATTENUATION_RATIO, 'OLSEN_ATTENUATION_RATIO', ier)
  if (ier /= 0) stop 'Error reading Par_file parameter OLSEN_ATTENUATION_RATIO'
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
  call read_value_integer(NTSTEP_BETWEEN_OUTPUT_SEISMOS, 'NTSTEP_BETWEEN_OUTPUT_SEISMOS', ier)
  if (ier /= 0) stop 'Error reading Par_file parameter NTSTEP_BETWEEN_OUTPUT_SEISMOS'
  call read_value_integer(NTSTEP_BETWEEN_READ_ADJSRC, 'NTSTEP_BETWEEN_READ_ADJSRC', ier)
  if (ier /= 0) stop 'Error reading Par_file parameter NTSTEP_BETWEEN_READ_ADJSRC'
  call read_value_logical(USE_FORCE_POINT_SOURCE, 'USE_FORCE_POINT_SOURCE', ier)
  if (ier /= 0) stop 'Error reading Par_file parameter USE_FORCE_POINT_SOURCE'
  call read_value_logical(USE_RICKER_TIME_FUNCTION, 'USE_RICKER_TIME_FUNCTION', ier)
  if (ier /= 0) stop 'Error reading Par_file parameter USE_RICKER_TIME_FUNCTION'
  call read_value_logical(PRINT_SOURCE_TIME_FUNCTION, 'PRINT_SOURCE_TIME_FUNCTION', ier)
  if (ier /= 0) stop 'Error reading Par_file parameter PRINT_SOURCE_TIME_FUNCTION'
  call read_value_logical(COUPLE_WITH_EXTERNAL_CODE, 'COUPLE_WITH_EXTERNAL_CODE', ier)
  if (ier /= 0) stop 'Error reading Par_file parameter COUPLE_WITH_EXTERNAL_CODE'
  call read_value_integer(EXTERNAL_CODE_TYPE, 'EXTERNAL_CODE_TYPE', ier)
  if (ier /= 0) stop 'Error reading Par_file parameter EXTERNAL_CODE_TYPE'
  call read_value_string(TRACTION_PATH, 'TRACTION_PATH', ier)
  if (ier /= 0) stop 'Error reading Par_file parameter TRACTION_PATH'
  call read_value_logical(MESH_A_CHUNK_OF_THE_EARTH, 'MESH_A_CHUNK_OF_THE_EARTH', ier)
  if (ier /= 0) stop 'Error reading Par_file parameter MESH_A_CHUNK_OF_THE_EARTH'

  ! close parameter file
  call close_parameter_file()

! check the type of external code to couple with, if any
  if (COUPLE_WITH_EXTERNAL_CODE) then
    if (EXTERNAL_CODE_TYPE /= EXTERNAL_CODE_IS_DSM .and. &
       EXTERNAL_CODE_TYPE /= EXTERNAL_CODE_IS_AXISEM .and. &
       EXTERNAL_CODE_TYPE /= EXTERNAL_CODE_IS_FK) stop 'incorrect value of EXTERNAL_CODE_TYPE read'

    if (EXTERNAL_CODE_TYPE == EXTERNAL_CODE_IS_AXISEM) &
         stop 'coupling with AxiSEM not implemented yet, but will soon be'

    if (EXTERNAL_CODE_TYPE == EXTERNAL_CODE_IS_FK) &
         stop 'coupling with F-K not implemented yet, but see work by Ping et al. (GJI 2014, GRL 2015)'
  endif

! see if we are running several independent runs in parallel
! if so, add the right directory for that run (group numbers start at zero, but directory names start at run0001, thus we add one)
! a negative value for "mygroup" is a convention that indicates that groups (i.e. sub-communicators, one per run) are off
  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
    write(path_to_add,"('run',i4.4,'/')") mygroup + 1
    LOCAL_PATH = path_to_add(1:len_trim(path_to_add))//LOCAL_PATH(1:len_trim(LOCAL_PATH))
    TOMOGRAPHY_PATH = path_to_add(1:len_trim(path_to_add))//TOMOGRAPHY_PATH(1:len_trim(TOMOGRAPHY_PATH))
    TRACTION_PATH = path_to_add(1:len_trim(path_to_add))//TRACTION_PATH(1:len_trim(TRACTION_PATH))
  endif

  ! noise simulations:
  ! double the number of time steps, if running noise simulations (+/- branches)
  if (NOISE_TOMOGRAPHY /= 0)   NSTEP = 2*NSTEP-1

  ! for noise simulations, we need to save movies at the surface (where the noise is generated)
  ! and thus we force MOVIE_SURFACE to be .true., in order to use variables defined for surface movies later
  if (NOISE_TOMOGRAPHY /= 0) then
    MOVIE_TYPE = 1
    MOVIE_SURFACE = .true.
    USE_HIGHRES_FOR_MOVIES = .true.     ! we need to save surface movie everywhere, i.e. at all GLL points on the surface
  endif

  ! the default value of NTSTEP_BETWEEN_READ_ADJSRC (0) is to read the whole trace at the same time
  if (NTSTEP_BETWEEN_READ_ADJSRC == 0)  NTSTEP_BETWEEN_READ_ADJSRC = NSTEP

  ! total times steps must be dividable by adjoint source chunks/blocks
  if (mod(NSTEP,NTSTEP_BETWEEN_READ_ADJSRC) /= 0) then
    print*,'When NOISE_TOMOGRAPHY is not equal to zero, ACTUAL_NSTEP=2*NSTEP-1'
    stop 'error: mod(NSTEP,NTSTEP_BETWEEN_READ_ADJSRC) must be zero! Please modify Par_file and recompile solver'
  endif

  ! checks number of nodes for 2D and 3D shape functions for quadrilaterals and hexahedra
  ! curvature (i.e. HEX27 elements) is not handled by our internal mesher, for that use Gmsh (CUBIT does not handle it either)
  if (NGNOD == 8) then
    NGNOD2D = 4
  else if (NGNOD == 27) then
    NGNOD2D = 9
  else if (NGNOD /= 8 .and. NGNOD /= 27) then
    stop 'elements should have 8 or 27 control nodes, please modify NGNOD in Par_file and recompile solver'
  endif

  if (USE_FORCE_POINT_SOURCE) then
    ! compute the total number of sources in the FORCESOLUTION file
    ! there are NLINES_PER_FORCESOLUTION_SOURCE lines per source in that file
    FORCESOLUTION = IN_DATA_FILES_PATH(1:len_trim(IN_DATA_FILES_PATH))//'FORCESOLUTION'
! see if we are running several independent runs in parallel
! if so, add the right directory for that run (group numbers start at zero, but directory names start at run0001, thus we add one)
! a negative value for "mygroup" is a convention that indicates that groups (i.e. sub-communicators, one per run) are off
    if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
      write(path_to_add,"('run',i4.4,'/')") mygroup + 1
      FORCESOLUTION = path_to_add(1:len_trim(path_to_add))//FORCESOLUTION(1:len_trim(FORCESOLUTION))
    endif

    open(unit=21,file=trim(FORCESOLUTION),status='old',action='read',iostat=ier)
    if (ier /= 0) stop 'error opening FORCESOLUTION file'

    icounter = 0
    do while (ier == 0)
      read(21,"(a)",iostat=ier) dummystring
      if (ier == 0) icounter = icounter + 1
    enddo
    close(21)

    if (mod(icounter,NLINES_PER_FORCESOLUTION_SOURCE) /= 0) &
      stop 'error: total number of lines in FORCESOLUTION file should be a multiple of NLINES_PER_FORCESOLUTION_SOURCE'

    NSOURCES = icounter / NLINES_PER_FORCESOLUTION_SOURCE
    if (NSOURCES < 1) stop 'error: need at least one source in FORCESOLUTION file'

  else
    ! compute the total number of sources in the CMTSOLUTION file
    ! there are NLINES_PER_CMTSOLUTION_SOURCE lines per source in that file
    CMTSOLUTION = IN_DATA_FILES_PATH(1:len_trim(IN_DATA_FILES_PATH))//'CMTSOLUTION'
! see if we are running several independent runs in parallel
! if so, add the right directory for that run (group numbers start at zero, but directory names start at run0001, thus we add one)
! a negative value for "mygroup" is a convention that indicates that groups (i.e. sub-communicators, one per run) are off
    if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
      write(path_to_add,"('run',i4.4,'/')") mygroup + 1
      CMTSOLUTION = path_to_add(1:len_trim(path_to_add))//CMTSOLUTION(1:len_trim(CMTSOLUTION))
    endif

    open(unit=21,file=trim(CMTSOLUTION),status='old',action='read',iostat=ier)
    if (ier /= 0) stop 'error opening CMTSOLUTION file'

    icounter = 0
    do while (ier == 0)
      read(21,"(a)",iostat=ier) dummystring
      if (ier == 0) icounter = icounter + 1
    enddo
    close(21)

    if (mod(icounter,NLINES_PER_CMTSOLUTION_SOURCE) /= 0) &
      stop 'error: total number of lines in CMTSOLUTION file should be a multiple of NLINES_PER_CMTSOLUTION_SOURCE'

    NSOURCES = icounter / NLINES_PER_CMTSOLUTION_SOURCE
    if (NSOURCES < 1) stop 'error: need at least one source in CMTSOLUTION file'

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
      stop 'error: hdur too small for movie creation, movies do not make sense for Heaviside source'
  endif

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
      stop 'Error: Using sep model requires defining a SEP_MODEL_DIRECTORY.'
    endif
    !inquire(directory=trim(SEP_MODEL_DIRECTORY), exists=sep_dir_exists)
    !if (.not. sep_dir_exists) then
    !  stop 'Error: SEP_MODEL_DIRECTORY should exist.'
    !endif


  case default
    print*
    print*,'********** model not recognized: ',trim(MODEL),' **************'
    print*,'********** using model: default',' **************'
    print*
    IMODEL = IMODEL_DEFAULT
  end select

  ! check
  if (IMODEL == IMODEL_IPATI .or. IMODEL == IMODEL_IPATI_WATER) then
    if (USE_RICKER_TIME_FUNCTION .eqv. .false.) &
      stop 'error: please set USE_RICKER_TIME_FUNCTION to .true. in Par_file and recompile solver'
  endif


  ! absorbing conditions
  ! stacey absorbing free surface must have also stacey turned on
  if (STACEY_INSTEAD_OF_FREE_SURFACE) STACEY_ABSORBING_CONDITIONS = .true.

  ! only one absorbing condition is possible
  if (PML_CONDITIONS) then
    if (STACEY_ABSORBING_CONDITIONS) &
      stop 'error: please set STACEY_ABSORBING_CONDITIONS and STACEY_INSTEAD_OF_FREE_SURFACE to .false. in Par_file for PML'
  endif

  end subroutine read_parameter_file

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_gpu_mode(GPU_MODE,GRAVITY)

  implicit none

  logical,intent(out) :: GPU_MODE
  logical,intent(out) :: GRAVITY

  ! local parameters
  integer :: ier

  ! initializes flags
  GPU_MODE = .false.
  GRAVITY = .false.

  ! opens file Par_file
  call open_parameter_file(ier)

  call read_value_logical(GPU_MODE, 'GPU_MODE', ier)
  call read_value_logical(GRAVITY, 'GRAVITY', ier)

  ! close parameter file
  call close_parameter_file()

  end subroutine read_gpu_mode

!===============================================================================
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

subroutine read_adios_parameters()

  use constants,only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

  use shared_input_parameters

  implicit none

  ! local parameters
  integer :: ier

  ! initialize flags to false
  ADIOS_ENABLED            = .false.
  ADIOS_FOR_DATABASES      = .false.
  ADIOS_FOR_MESH           = .false.
  ADIOS_FOR_FORWARD_ARRAYS = .false.
  ADIOS_FOR_KERNELS        = .false.

  ! opens file Par_file
  call open_parameter_file(ier)
  call read_value_logical(ADIOS_ENABLED, 'ADIOS_ENABLED', ier)

  if (ier == 0 .and. ADIOS_ENABLED) then
    call read_value_logical(ADIOS_FOR_DATABASES, 'ADIOS_FOR_DATABASES', ier)
    call read_value_logical(ADIOS_FOR_MESH, 'ADIOS_FOR_MESH', ier)
    call read_value_logical(ADIOS_FOR_FORWARD_ARRAYS, 'ADIOS_FOR_FORWARD_ARRAYS', ier)
    call read_value_logical(ADIOS_FOR_KERNELS, 'ADIOS_FOR_KERNELS', ier)
  endif

  call close_parameter_file()

  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. BROADCAST_SAME_MESH_AND_MODEL .and. ADIOS_ENABLED) &
    stop 'ADIOS not yet supported by option BROADCAST_SAME_MESH_AND_MODEL'

end subroutine read_adios_parameters

