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
  call read_parameter_file()

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

