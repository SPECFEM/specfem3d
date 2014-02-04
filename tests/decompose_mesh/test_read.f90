program test_read

  use decompose_mesh

  implicit none

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
                          PML_INSTEAD_OF_FREE_SURFACE,f0_FOR_PML,IMODEL,FULL_ATTENUATION_SOLID,TRAC_PATH)

  ! punctual check of values for given default Par_file in SPECFEM3D/DATA/ directory
  print*,'NPROC = ',NPROC
  if( NPROC /= 4 ) then
    print*,'ERROR: NPROC value invalid'
    stop 1
  else
    print*,'  result is correct'
  endif

  print*,'NSTEP = ',NSTEP
  if( NSTEP /= 200 ) then
    print*,'ERROR: NSTEP value invalid'
    stop 1
  else
    print*,'  result is correct'
  endif

  print*,'DT = ',DT
  if( abs(DT - 0.03) > 1.e-9 ) then
    print*,'ERROR: DT value invalid'
    stop 1
  else
    print*,'  result is correct'
  endif

  print*,'NTSTEP_BETWEEN_OUTPUT_INFO = ',NTSTEP_BETWEEN_OUTPUT_INFO
  if( NTSTEP_BETWEEN_OUTPUT_INFO /= 100 ) then
    print*,'ERROR: NPROC value invalid'
    stop 1
  else
    print*,'  result is correct'
  endif

  print*,'USE_RICKER_TIME_FUNCTION = ',USE_RICKER_TIME_FUNCTION
  if( USE_RICKER_TIME_FUNCTION .neqv. .false. ) then
    print*,'ERROR: NPROC value invalid'
    stop 1
  else
    print*,'  result is correct'
  endif

  ! done
  print*,'test_read done successfully'

end program test_read

