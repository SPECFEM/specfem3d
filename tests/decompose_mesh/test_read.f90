program test_read

  use decompose_mesh

  implicit none

  ! some test values from the default Par_file
  integer,parameter :: PAR_FILE_NPROC = 4
  integer,parameter :: PAR_FILE_NSTEP = 5000
  integer,parameter :: PAR_FILE_NT = 5000
  integer,parameter :: PAR_FILE_NTSTEP_BETWEEN_OUTPUT_INFO = 100
  double precision,parameter :: PAR_FILE_DT = 0.05
  logical,parameter :: PAR_FILE_USE_RICKER_TIME_FUNCTION = .false.

  ! reads ../DATA/Par_file
  call read_parameter_file(NPROC,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,DT,NGNOD,NGNOD2D, &
                           UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION,TOMOGRAPHY_PATH, &
                           ATTENUATION,USE_OLSEN_ATTENUATION,LOCAL_PATH,NSOURCES, &
                           APPROXIMATE_OCEAN_LOAD,TOPOGRAPHY,ANISOTROPY,STACEY_ABSORBING_CONDITIONS,MOVIE_TYPE, &
                           MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
                           NTSTEP_BETWEEN_FRAMES,USE_HIGHRES_FOR_MOVIES,HDUR_MOVIE, &
                           SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION, &
                           NTSTEP_BETWEEN_OUTPUT_INFO,SIMULATION_TYPE,SAVE_FORWARD, &
                           NTSTEP_BETWEEN_READ_ADJSRC,NOISE_TOMOGRAPHY, &
                           USE_FORCE_POINT_SOURCE,STACEY_INSTEAD_OF_FREE_SURFACE, &
                           USE_RICKER_TIME_FUNCTION,OLSEN_ATTENUATION_RATIO,PML_CONDITIONS, &
                           PML_INSTEAD_OF_FREE_SURFACE,f0_FOR_PML,IMODEL,SEP_MODEL_DIRECTORY, &
                           FULL_ATTENUATION_SOLID,TRACTION_PATH,COUPLE_WITH_EXTERNAL_CODE,EXTERNAL_CODE_TYPE, &
                           MESH_A_CHUNK_OF_THE_EARTH)

  ! punctual check of values for given default Par_file in SPECFEM3D/DATA/ directory
  print*,'NPROC = ',NPROC
  if (NPROC /= PAR_FILE_NPROC) then
    print*,'ERROR: NPROC value invalid'
    stop 1
  else
    print*,'  result is correct'
  endif

  print*,'NSTEP = ',NSTEP
  if (NSTEP /= PAR_FILE_NSTEP) then
    print*,'ERROR: NSTEP value invalid'
    stop 1
  else
    print*,'  result is correct'
  endif

  print*,'DT = ',DT
  if (abs(DT - PAR_FILE_DT) > 1.e-9) then
    print*,'ERROR: DT value invalid'
    stop 1
  else
    print*,'  result is correct'
  endif

  print*,'NTSTEP_BETWEEN_OUTPUT_INFO = ',NTSTEP_BETWEEN_OUTPUT_INFO
  if (NTSTEP_BETWEEN_OUTPUT_INFO /= PAR_FILE_NTSTEP_BETWEEN_OUTPUT_INFO) then
    print*,'ERROR: NPROC value invalid'
    stop 1
  else
    print*,'  result is correct'
  endif

  print*,'USE_RICKER_TIME_FUNCTION = ',USE_RICKER_TIME_FUNCTION
  if (USE_RICKER_TIME_FUNCTION .neqv. PAR_FILE_USE_RICKER_TIME_FUNCTION) then
    print*,'ERROR: NPROC value invalid'
    stop 1
  else
    print*,'  result is correct'
  endif

  ! done
  print*,'test_read done successfully'

end program test_read

