! file "adios_helpers_readers_1d_generic.inc"
!
!
! common include section for adios read schedule helper functions

#if defined(USE_ADIOS)
  ! adios
  integer(kind=8), intent(in) :: adios_handle
  integer(kind=8), intent(in) :: adios_group
  ! for step
  integer :: istep
#elif defined(USE_ADIOS2)
  ! adios2
  type(adios2_engine), intent(in) :: adios_handle
  type(adios2_io), intent(in) :: adios_group
  ! variable
  type(adios2_variable) :: v
  integer(kind=8) :: step_start
#endif
  integer(kind=8), intent(in) :: sel
  integer(kind=8), dimension(1), intent(in) :: start,count

  character(len=*), intent(in) :: array_name
  integer(kind=8), intent(in), optional :: step

  ! local parameters
  integer :: ier
  integer(kind=8) :: idummy

  TRACE_ADIOS_L2_ARG('read_adios_schedule_array_global_1d adios_helpers_readers_1d_generic: ',trim(array_name))

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in adios_helpers_readers_1d_generic'

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! only schedules read
  if (present(step)) then
    ! reads single step at specified step position
    istep = int(step)
    call adios_schedule_read(adios_handle, sel, trim(array_name), istep, 1, array, ier)
    if (ier /= 0 ) then
      print *,'Error adios: scheduling read of array ',trim(array_name),' at step ',istep,step,' failed'
      stop 'Error adios helper schedule read array at step'
    endif
  else
    ! default read
    call adios_schedule_read(adios_handle, sel, trim(array_name), 0, 1, array, ier)
    if (ier /= 0 ) then
      print *,'Error adios: scheduling read of array ',trim(array_name),' failed'
      stop 'Error adios helper schedule read array'
    endif
  endif

  ! to avoid compiler warnings
  idummy = start(1)
  idummy = count(1)
  idummy = adios_group

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! gets associated variable for array
  call adios2_inquire_variable(v, adios_group, trim(array_name), ier)
  call check_adios_err(ier,"Error adios2 readers_1d_generic: inquire variable "//trim(array_name)//" failed")

  if (.not. v%valid) stop 'Error adios2 variable invalid'

  call adios2_set_selection(v, 1, start, count, ier)
  call check_adios_err(ier,"Error adios2 set selection for "//trim(array_name)//" failed")

  ! selection
  if (present(step)) then
    step_start = step
  else
    step_start = 0
  endif
  call adios2_set_step_selection(v, step_start, int(1,kind=8), ier)
  call check_adios_err(ier, "Error adios2 set step variable for "//trim(array_name)//" failed")

  ! reads array data
  ! default: deferred mode adios2_mode_deferred, otherwise sync mode adios2_mode_sync
  call adios2_get(adios_handle, v, array, myadios_reader_mode, ier)
  call check_adios_err(ier,"Error adios2 get for array "//trim(array_name)//" failed")

  ! to avoid compiler warnings
  idummy = sel

#endif
