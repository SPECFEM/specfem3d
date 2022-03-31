! file "adios_helpers_writers_1d_generic_offset.inc"
!
!
! common include section for adios write 1d helper functions with offset

  ! common parameters
#if defined(USE_ADIOS)
  integer(kind=8), intent(in) :: adios_handle
  integer(kind=8), intent(in) :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_engine), intent(in) :: adios_handle
  type(adios2_io), intent(in) :: adios_group
  type(adios2_variable) :: v
#endif
  integer, intent(in) :: myrank, sizeprocs
  integer(kind=8), intent(in) :: local_dim
  integer(kind=8), intent(in) :: global_dim, offset
  character(len=*) :: array_name
  ! Variables
  integer :: adios_err
  integer(kind=8) :: array_size

  TRACE_ADIOS_L2_ARG('write_adios_global_1d adios_helpers_writers_1d_generic_offset: ',trim(array_name))

  ! checks name
  if (len_trim(array_name) == 0) stop 'Error: array_name has zero length in adios_helpers_writers_1d_generic_offset'

  ! gets array size
  array_size = size(array,kind=8)

  ! adds local_dim/global_dim/offset infos
  call write_1D_global_array_adios_dims_offset(adios_handle, adios_group, myrank, &
                                               local_dim, global_dim, offset, sizeprocs, array_name, array_size)

#if defined(USE_ADIOS)
  ! ADIOS 1
  call adios_write(adios_handle, trim(array_name)// "/array", array, adios_err)
  call check_adios_err(adios_err,"Error writing "//trim(array_name) // "/array")

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! gets associated variable for array from current group
  call adios2_inquire_variable(v, adios_group, trim(array_name)// "/array", adios_err)
  call check_adios_err(adios_err,"Error adios2 writers offset: inquire variable "//trim(array_name)// "/array" //" failed")

  if (.not. v%valid) stop 'Error adios2 variable invalid'

  ! no special offset case ... yet
  call adios2_put(adios_handle, v, array, myadios_writer_mode, adios_err)
  call check_adios_err(adios_err,"Error writing adios2 array "//trim(array_name) // "/array")
#endif
