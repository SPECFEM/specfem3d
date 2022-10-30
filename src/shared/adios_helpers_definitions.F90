!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
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


!-------------------------------------------------------------------------------
!> Helpers to set up adios features.
!! * Scalar definition
!! * Global arrays definition
!!
!! \author MPBL
!-------------------------------------------------------------------------------

#include "config.fh"


module adios_helpers_definitions_mod

#if defined(USE_ADIOS)
  use adios_write_mod, only: adios_define_var, &
                             adios_long, adios_integer, adios_double, adios_real, adios_byte, adios_string
#elif defined(USE_ADIOS2)
  use adios2
#endif

  use manager_adios, only: check_adios_err
#if defined(USE_ADIOS2)
  use manager_adios, only: sizeprocs_adios,myrank_adios
#endif

#if defined(USE_ADIOS)
  use adios_write_mod, only: adios_set_transform
#elif defined(USE_ADIOS2)
  use manager_adios, only: myadios2_obj
#endif

  ! compression
  use manager_adios, only: use_adios_compression
  use constants, only: ADIOS_COMPRESSION_ALGORITHM,ADIOS_COMPRESSION_MODE,ADIOS_COMPRESSION_MODE_VALUE

  implicit none

  private

  public :: define_adios_scalar
  public :: define_adios_global_real_1d_array
  public :: define_adios_global_double_1d_array
  public :: define_adios_global_integer_1d_array
  public :: define_adios_global_long_1d_array
  public :: define_adios_global_logical_1d_array
  public :: define_adios_global_string_1d_array
  public :: define_adios_local_string_1d_array
  public :: define_adios_global_array1D
  public :: define_adios_compression

  ! Generic interface to define scalar variables in ADIOS
  interface define_adios_scalar
    module procedure define_adios_scalar_double
    module procedure define_adios_scalar_float
    module procedure define_adios_scalar_integer
    module procedure define_adios_scalar_long
    module procedure define_adios_scalar_byte
  end interface define_adios_scalar

  interface define_adios_global_real_1d_array
    module procedure define_adios_global_1d_real_1d
    ! unused so far...
    !module procedure define_adios_global_1d_real_2d
    !module procedure define_adios_global_1d_real_3d
    !module procedure define_adios_global_1d_real_4d
    !module procedure define_adios_global_1d_real_5d
  end interface define_adios_global_real_1d_array

  interface define_adios_global_double_1d_array
    module procedure define_adios_global_1d_double_1d
    ! unused so far..
    !module procedure define_adios_global_1d_double_2d
    !module procedure define_adios_global_1d_double_3d
    !module procedure define_adios_global_1d_double_4d
    !module procedure define_adios_global_1d_double_5d
  end interface define_adios_global_double_1d_array

  interface define_adios_global_integer_1d_array
    module procedure define_adios_global_1d_int_1d
    ! unused so far..
    !module procedure define_adios_global_1d_int_2d
    !module procedure define_adios_global_1d_int_3d
    !module procedure define_adios_global_1d_int_4d
    !module procedure define_adios_global_1d_int_5d
  end interface define_adios_global_integer_1d_array

  interface define_adios_global_long_1d_array
    module procedure define_adios_global_1d_long_1d
    ! unused so far..
    !module procedure define_adios_global_1d_long_2d
    !module procedure define_adios_global_1d_long_3d
    !module procedure define_adios_global_1d_long_4d
    !module procedure define_adios_global_1d_long_5d
  end interface define_adios_global_long_1d_array

  interface define_adios_global_logical_1d_array
    module procedure define_adios_global_1d_logical_1d
    ! unused so far..
    !module procedure define_adios_global_1d_logical_2d
    !module procedure define_adios_global_1d_logical_3d
    !module procedure define_adios_global_1d_logical_4d
    !module procedure define_adios_global_1d_logical_5d
  end interface define_adios_global_logical_1d_array

  interface define_adios_global_string_1d_array
    module procedure define_adios_global_1d_string_1d
  end interface define_adios_global_string_1d_array

  interface define_adios_local_string_1d_array
    module procedure define_adios_local_1d_string_1d
  end interface define_adios_local_string_1d_array

  ! Cannot include an interface in another interface
  interface define_adios_global_array1D
    module procedure define_adios_global_1d_int_1d
    module procedure define_adios_global_1d_int_2d
    module procedure define_adios_global_1d_int_3d
    module procedure define_adios_global_1d_int_4d
    module procedure define_adios_global_1d_int_5d

    module procedure define_adios_global_1d_long_1d
    module procedure define_adios_global_1d_long_2d
    module procedure define_adios_global_1d_long_3d
    module procedure define_adios_global_1d_long_4d
    module procedure define_adios_global_1d_long_5d

    module procedure define_adios_global_1d_logical_1d
    module procedure define_adios_global_1d_logical_2d
    module procedure define_adios_global_1d_logical_3d
    module procedure define_adios_global_1d_logical_4d
    module procedure define_adios_global_1d_logical_5d

    module procedure define_adios_global_1d_real_1d
    module procedure define_adios_global_1d_real_2d
    module procedure define_adios_global_1d_real_3d
    module procedure define_adios_global_1d_real_4d
    module procedure define_adios_global_1d_real_5d

    module procedure define_adios_global_1d_double_1d
    module procedure define_adios_global_1d_double_2d
    module procedure define_adios_global_1d_double_3d
    module procedure define_adios_global_1d_double_4d
    module procedure define_adios_global_1d_double_5d

    module procedure define_adios_global_1d_string_1d

  end interface define_adios_global_array1D

  ! compression
#if defined(USE_ADIOS)
  ! for undo_att snapshots compression - compression transform string
  character(len=128) :: myadios_comp_operator
#elif defined(USE_ADIOS2)
  ! for undo_att snapshots compression - compression operator
  type(adios2_operator), save :: myadios_comp_operator        ! see note about save attribute in adios_manager.F90
#endif

contains


!===============================================================================
!
! scalars
!
!===============================================================================

!===============================================================================
!> Define an ADIOS scalar double precision variable and autoincrement
!! the adios group size by (8).
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param name The variable name in the ADIOS file.
!! \param var The variable to be defined. Used for type inference. Can be
!!            ignored.
!!
!! \note 'name' and 'var' are written as successive arguments on purpose.
!!       One should be able to define a macro such as:
!!       #define STRINGIFY_VAR(x) #x, x
!!       Calling define_adi os_double_scalar with such a macro will be done as:
!!       call define_adios_scalar_double(group, size, path, STRINGIFY_VAR(x))
!!       as STRINGIFY_VAR(x) expand as:
!!       "x", x
!!       x being the variable name inside the code.
subroutine define_adios_scalar_double (adios_group, group_size_inc, path, name, var)

  implicit none
  ! Arguments
#if defined(USE_ADIOS)
  integer(kind=8),  intent(in)    :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_io),  intent(in)     :: adios_group
#endif
  integer(kind=8),  intent(inout) :: group_size_inc
  character(len=*), intent(in)    :: name, path
  real(kind=8),     intent(in)    :: var

  ! Local Variables
#if defined(USE_ADIOS)
  integer(kind=8) :: varid ! dummy variable, adios use var name
#elif defined(USE_ADIOS2)
  type(adios2_variable) :: v
  integer :: ier
  !integer(kind=8), parameter :: count_dims(1) = (/ 1 /)
  integer(kind=8) :: ldim(1),gdim(1),offs(1)
#endif
  real(kind=8) :: idummy
  character(len=256) :: full_name

  TRACE_ADIOS_L2_ARG('define_adios_scalar_double: ',trim(name))

  ! check
  if (len_trim(name) == 0) stop 'Error adios: invalid name in define_adios_scalar_double()'

  ! sets full variable name
  if (len_trim(path) == 0) then
    full_name = trim(name)
  else
    full_name = trim(path) // '/' // trim(name)
  endif

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! adios: 6 == real(kind=8)
  call adios_define_var (adios_group, trim(name), trim(path), adios_double,  '', '', '', varid)

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! note: we won't store the variable object, but get it back by adios2_inquire_** for writing
  ! defines scalar as global variable (same for all processes, would be enough to be stored by main process myrank==0)
  ! > call adios2_define_variable(v, adios_group, trim(full_name), adios2_type_real8, ier)
  ! defines scalar as local variable (value can vary between processes)
  ! > call adios2_define_variable(v, adios_group, trim(full_name), adios2_type_real8, &
  !                               1, adios2_null_dims, adios2_null_dims, count_dims, adios2_constant_dims, ier)
  !
  ! defines scalar as 1-d array with single entry
  ldim(1) = 1
  gdim(:) = int(sizeprocs_adios,kind=8) * ldim(:)
  offs(:) = int(myrank_adios,kind=8) * ldim(:)

  call adios2_define_variable(v, adios_group, trim(full_name), adios2_type_real8, &
                              1, gdim, offs, ldim, adios2_constant_dims, ier)
  call check_adios_err(ier,"Error adios2 could not define parameter: "//trim(full_name))

  if (.not. v%valid) stop 'Error adios2 defined variable is invalid'

#endif

  group_size_inc = group_size_inc + 8

  ! to avoid compiler warnings
  idummy = var

end subroutine define_adios_scalar_double


!===============================================================================
!> Define an ADIOS scalar single precision variable and autoincrement
!! the adios group size by (4).
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param name The variable name in the ADIOS file.
!! \param var The variable to be defined. Used for type inference. Can be
!             ignored.
!!
!! \note See define_adios_scalar_double()
subroutine define_adios_scalar_float(adios_group, group_size_inc, path, name, var)

  implicit none
  ! Arguments
#if defined(USE_ADIOS)
  integer(kind=8),  intent(in)    :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_io),  intent(in)     :: adios_group
#endif
  integer(kind=8),  intent(inout) :: group_size_inc
  character(len=*), intent(in)    :: name, path
  real(kind=4),     intent(in)    :: var

  ! Local Variables
#if defined(USE_ADIOS)
  integer(kind=8) :: varid ! dummy variable, adios use var name
#elif defined(USE_ADIOS2)
  type(adios2_variable) :: v
  integer :: ier
  !integer(kind=8), parameter :: count_dims(1) = (/ 1 /)
  integer(kind=8) :: ldim(1),gdim(1),offs(1)
#endif
  real(kind=4) :: idummy
  character(len=256) :: full_name

  TRACE_ADIOS_L2_ARG('define_adios_scalar_float: ',trim(name))

  ! check
  if (len_trim(name) == 0) stop 'Error adios: invalid name in define_adios_scalar_float()'

  ! sets full variable name
  if (len_trim(path) == 0) then
    full_name = trim(name)
  else
    full_name = trim(path) // '/' // trim(name)
  endif

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! adios: 6 == real(kind=8), 5 == real(kind=4)
  call adios_define_var (adios_group, trim(name), trim(path), adios_real,  '', '', '', varid)

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! note: we won't store the variable object, but get it back by adios2_inquire_** for writing
  ! defines scalar as global variable (same for all processes, would be enough to be stored by main process myrank==0)
  ! > call adios2_define_variable(v, adios_group, trim(full_name), adios2_type_real4, ier)
  ! defines scalar as local variable (value can vary between processes)
  ! > call adios2_define_variable(v, adios_group, trim(full_name), adios2_type_real4, &
  !                               1, adios2_null_dims, adios2_null_dims, count_dims, adios2_constant_dims, ier)
  !
  ! defines scalar as 1-d array with single entry
  ldim(1) = 1
  gdim(:) = int(sizeprocs_adios,kind=8) * ldim(:)
  offs(:) = int(myrank_adios,kind=8) * ldim(:)

  call adios2_define_variable(v, adios_group, trim(full_name), adios2_type_real4, &
                              1, gdim, offs, ldim, adios2_constant_dims, ier)
  call check_adios_err(ier,"Error adios2 could not define parameter: "//trim(full_name))

  if (.not. v%valid) stop 'Error adios2 defined variable is invalid'

#endif

  group_size_inc = group_size_inc + 4

  ! to avoid compiler warnings
  idummy = var

end subroutine define_adios_scalar_float


!===============================================================================
!> Define an ADIOS scalar integer variable and autoincrement the adios
!! group size by (4).
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param name The variable name in the ADIOS file.
!! \param var The variable to be defined. Used for type inference. Can be
!             ignored.
!!
!! \note See define_adios_scalar_double()
subroutine define_adios_scalar_integer(adios_group, group_size_inc, path, name, var)

  implicit none
  ! Arguments
#if defined(USE_ADIOS)
  integer(kind=8),  intent(in)    :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_io),  intent(in)    :: adios_group
#endif
  character(len=*), intent(in)    :: name, path
  integer(kind=8),  intent(inout) :: group_size_inc
  integer(kind=4),  intent(in)    :: var

  ! Local Variables
#if defined(USE_ADIOS)
  integer(kind=8) :: varid ! dummy variable, adios use var name
#elif defined(USE_ADIOS2)
  type(adios2_variable) :: v
  integer :: ier
  !integer(kind=8), parameter :: count_dims(1) = (/ 1 /)
  integer(kind=8) :: ldim(1),gdim(1),offs(1)
#endif
  integer(kind=4) :: idummy
  character(len=256) :: full_name

  TRACE_ADIOS_L2_ARG('define_adios_scalar_integer: ',trim(name))

  ! check
  if (len_trim(name) == 0) stop 'Error adios: invalid name in define_adios_scalar_integer()'

  ! sets full variable name
  if (len_trim(path) == 0) then
    full_name = trim(name)
  else
    full_name = trim(path) // '/' // trim(name)
  endif

  !debug
  !print *,'debug adios: ',myrank_adios,' define integer scalar: ',trim(full_name),' path: ',trim(path),' name: ',trim(name)

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! adios: 2 ~ integer(kind=4)
  call adios_define_var (adios_group, trim(name), trim(path), adios_integer, '', '', '', varid)

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! note 1: we won't store the variable object, but get it back by adios2_inquire_** for writing
  !
  ! note 2: we will store a scalar as a 1-D array with a single entry instead.
  !         this is due to appending to a file will increase the step count for variables.
  !         retrieving local scalar variables would use adios2_set_block_selection() which then fails in such cases with an error:
  !          " ERROR: invalid blockID 0 from steps start 0 in variable reg2/nspec,
  !                   check argument to Variable < T>::SetBlockID, in call to Get "
  !
  !         however, using 1-D arrays will use adios2_set_selection() which succeeds also for appended variables.
  !         until adios2 fixes this, we will use the 1-D work-around.

  ! defines scalar as global variable (same for all processes, would be enough to be stored by main process myrank==0)
  ! > call adios2_define_variable(v, adios_group, trim(full_name), adios2_type_integer4, ier)
  !
  ! defines scalar as local variable (value can vary between processes)
  ! > call adios2_define_variable(v, adios_group, trim(full_name), adios2_type_integer4, &
  !                                1, adios2_null_dims, adios2_null_dims, count_dims, adios2_constant_dims, ier)
  !
  ! defines scalar as 1-d array with single entry
  ldim(1) = 1
  gdim(:) = int(sizeprocs_adios,kind=8) * ldim(:)
  offs(:) = int(myrank_adios,kind=8) * ldim(:)

  call adios2_define_variable(v, adios_group, trim(full_name), adios2_type_integer4, &
                              1, gdim, offs, ldim, adios2_constant_dims, ier)
  call check_adios_err(ier,"Error adios2 could not define parameter: "//trim(full_name))

  if (.not. v%valid) stop 'Error adios2 defined variable is invalid'

#endif

  group_size_inc = group_size_inc + 4

  ! to avoid compiler warnings
  idummy = var

end subroutine define_adios_scalar_integer


!===============================================================================
!> Define an ADIOS scalar long integer variable and autoincrement the adios
!! group size by (8).
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param name The variable name in the ADIOS file.
!! \param var The variable to be defined. Used for type inference. Can be
!             ignored.
!!
!! \note See define_adios_scalar_double()
subroutine define_adios_scalar_long(adios_group, group_size_inc, path, name, var)

  implicit none
  ! Arguments
#if defined(USE_ADIOS)
  integer(kind=8),  intent(in)    :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_io),  intent(in)     :: adios_group
#endif
  character(len=*), intent(in)     :: name, path
  integer(kind=8),  intent(inout)  :: group_size_inc
  integer(kind=8),  intent(in)  :: var

  ! Local Variables
#if defined(USE_ADIOS)
  integer(kind=8) :: varid ! dummy variable, adios use var name
#elif defined(USE_ADIOS2)
  type(adios2_variable) :: v
  integer :: ier
  !integer(kind=8), parameter :: count_dims(1) = (/ 1 /)
  integer(kind=8) :: ldim(1),gdim(1),offs(1)
#endif
  integer(kind=8) :: idummy
  character(len=256) :: full_name

  TRACE_ADIOS_L2_ARG('define_adios_scalar_long: ',trim(name))

  ! check
  if (len_trim(name) == 0) stop 'Error adios: invalid name in define_adios_scalar_long()'

  ! sets full variable name
  if (len_trim(path) == 0) then
    full_name = trim(name)
  else
    full_name = trim(path) // '/' // trim(name)
  endif

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! adios: 2 ~ integer(kind=8)
  call adios_define_var (adios_group, trim(name), trim(path), adios_long, '', '', '', varid)

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! note: we won't store the variable object, but get it back by adios2_inquire_** for writing
  ! defines global variable
  ! defines scalar as global variable (same for all processes, would be enough to be stored by main process myrank==0)
  ! > call adios2_define_variable(v, adios_group, trim(full_name), adios2_type_integer8, ier)
  ! defines scalar as local variable (value can vary between processes)
  ! > call adios2_define_variable(v, adios_group, trim(full_name), adios2_type_integer8, &
  !                               1, adios2_null_dims, adios2_null_dims, count_dims, adios2_constant_dims, ier)
  !
  ! defines scalar as 1-d array with single entry
  ldim(1) = 1
  gdim(:) = int(sizeprocs_adios,kind=8) * ldim(:)
  offs(:) = int(myrank_adios,kind=8) * ldim(:)

  call adios2_define_variable(v, adios_group, trim(full_name), adios2_type_integer8, &
                              1, gdim, offs, ldim, adios2_constant_dims, ier)
  call check_adios_err(ier,"Error adios2 could not define parameter: "//trim(full_name))

  if (.not. v%valid) stop 'Error adios2 defined variable is invalid'

#endif

  group_size_inc = group_size_inc + 8

  ! to avoid compiler warnings
  idummy = var

end subroutine define_adios_scalar_long

!===============================================================================
!> Define an ADIOS scalar byte variable and autoincrement the adios
!! group size by (1).
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param name The variable name in the ADIOS file.
!! \param var The variable to be defined. Used for type inference. Can be
!             ignored.
!!
!! \note See define_adios_scalar_double()
subroutine define_adios_scalar_byte (adios_group, group_size_inc, name, path, var)

  implicit none
  ! Arguments
#if defined(USE_ADIOS)
  integer(kind=8),  intent(in)    :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_io),  intent(in)     :: adios_group
#endif
  character(len=*), intent(in)     :: name, path
  integer(kind=8),  intent(inout)  :: group_size_inc
  ! note: byte is non-standard gnu Fortran
  !byte,     intent(in)             :: var
  integer(kind=1),  intent(in)     :: var

  ! Local Variables
#if defined(USE_ADIOS)
  integer(kind=8) :: varid ! dummy variable, adios use var name
#elif defined(USE_ADIOS2)
  type(adios2_variable) :: v
  integer :: ier
  !integer(kind=8), parameter :: count_dims(1) = (/ 1 /)
  integer(kind=8) :: ldim(1),gdim(1),offs(1)
#endif
  integer(kind=1) :: idummy
  character(len=256) :: full_name

  TRACE_ADIOS_L2_ARG('define_adios_scalar_byte: ',trim(name))

  ! check
  if (len_trim(name) == 0) stop 'Error adios: invalid name in define_adios_scalar_byte()'

  ! sets full variable name
  if (len_trim(path) == 0) then
    full_name = trim(name)
  else
    full_name = trim(path) // '/' // trim(name)
  endif

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! adios: 0 == byte == any_data_type(kind=1)
  call adios_define_var (adios_group, trim(name), trim(path), adios_byte,  '', '', '', varid)

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! note: we won't store the variable object, but get it back by adios2_inquire_** for writing
  ! defines global variable
  ! defines scalar as global variable (same for all processes, would be enough to be stored by main process myrank==0)
  ! > call adios2_define_variable(v, adios_group, trim(full_name), adios2_type_integer1, ier)
  ! defines scalar as local variable (value can vary between processes)
  ! > call adios2_define_variable(v, adios_group, trim(full_name), adios2_type_integer1, &
  !                               1, adios2_null_dims, adios2_null_dims, count_dims, adios2_constant_dims, ier)
  !
  ! defines scalar as 1-d array with single entry
  ldim(1) = 1
  gdim(:) = int(sizeprocs_adios,kind=8) * ldim(:)
  offs(:) = int(myrank_adios,kind=8) * ldim(:)

  call adios2_define_variable(v, adios_group, trim(full_name), adios2_type_integer1, &
                              1, gdim, offs, ldim, adios2_constant_dims, ier)
  call check_adios_err(ier,"Error adios2 could not define parameter: "//trim(full_name))

  if (.not. v%valid) stop 'Error adios2 defined variable is invalid'

#endif

  group_size_inc = group_size_inc + 1

  ! to avoid compiler warnings
  idummy = var

end subroutine define_adios_scalar_byte


!===============================================================================
!
! arrays
!
!===============================================================================

!===============================================================================
!> Define the dimensions that will be written along a global array in ADIOS.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param local_dim The local dimension of the array.
subroutine define_adios_global_dims_1d(adios_group, group_size_inc, array_name, local_dim)

  implicit none
  ! Parameters
#if defined(USE_ADIOS)
  integer(kind=8),  intent(in) :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_io),  intent(in) :: adios_group
#endif
  character(len=*), intent(in) :: array_name
  integer(kind=8),  intent(in) :: local_dim
  integer(kind=8),  intent(inout) :: group_size_inc

  TRACE_ADIOS_L2_ARG('define_adios_global_dims_1d: ',trim(array_name))

  ! array_name should be defined
  if (len_trim(array_name) == 0) stop 'Error adios: invalid array_name in define_adios_global_dims_1d()'

  ! uses local_dim as dummy variable
  call define_adios_scalar(adios_group, group_size_inc, trim(array_name), "local_dim", local_dim)  ! scalar long type
  call define_adios_scalar(adios_group, group_size_inc, trim(array_name), "global_dim", local_dim)
  call define_adios_scalar(adios_group, group_size_inc, trim(array_name), "offset", local_dim)

  ! additional scalar to specify actual array size, especially needed for ADIOS1 when arrays have variable sizes across ranks.
  ! note: in ADIOS2, the local_dim array is a separate info from the actual array, so we could only store local_dim there
  !       and use the actual array size dimension in the adios2_define_var() call below.
  !       in ADIOS1, the adios_define_var() call needs a local array size specifier, for which we now use array/size.
  call define_adios_scalar(adios_group, group_size_inc, trim(array_name), "size", local_dim)

end subroutine define_adios_global_dims_1d


!===============================================================================
!> Define a real global array in ADIOS regardless of the array shape
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param local_dim The local dimension of the array.
  subroutine define_adios_global_1d_generic_real(adios_group, group_size_inc, array_name, local_dim, array_size)

  implicit none
  ! Parameters
#if defined(USE_ADIOS)
  integer(kind=8),  intent(in)    :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_io),  intent(in)     :: adios_group
#endif
  character(len=*), intent(in) :: array_name
  integer(kind=8),  intent(in) :: local_dim
  integer(kind=8),  intent(inout) :: group_size_inc
  integer(kind=8),  intent(in) :: array_size
  ! Variables
#if defined(USE_ADIOS)
  integer(kind=8) :: var_id
#elif defined(USE_ADIOS2)
  integer(kind=8), dimension(1) :: gdim  ! Global shape of array
  integer(kind=8), dimension(1) :: ldim  ! Local size of array
  integer(kind=8), dimension(1) :: offs  ! Starting offset in global array
  type(adios2_variable) :: v
  character(len=256) :: full_name
  ! compression
  integer :: operation_id
#endif
  integer :: ier

  TRACE_ADIOS_L2_ARG('define_adios_global_1d_generic_real: ',trim(array_name))

  ! Define the dimensions of the array.
  ! local_dim used as a dummy variable to call the integer routine.
  call define_adios_global_dims_1d(adios_group, group_size_inc, trim(array_name), local_dim)

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! we specify the actual size of the array in array_name/size which differs to array_name/local_dim for cases
  ! where ranks have different slice sizes.
  ! we will still use local_dim to calculate memory offsets for different ranks when reading arrays.
  call adios_define_var(adios_group, "array", trim(array_name), adios_real, &
                        trim(array_name) // "/size", &
                        trim(array_name) // "/global_dim", &
                        trim(array_name) // "/offset", var_id)

  ! old: assumes all ranks will have same array sizes
  !call adios_define_var(adios_group, "array", trim(array_name), adios_real, &
  !                      trim(array_name) // "/local_dim", &
  !                      trim(array_name) // "/global_dim", &
  !                      trim(array_name) // "/offset", var_id)

  ! compression
  if (use_adios_compression) then
    ! adds compression transform
    call adios_set_transform(var_id,myadios_comp_operator,ier)
    call check_adios_err(ier,"Error adios compression: set transform for variable "//trim(array_name)//" failed")
  endif

  ! to avoid compiler warning
  ier = array_size

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! note: variable parameter v will not be stored globally, instead will be retrieved by _inquire** function again
  ldim(1) = array_size   ! instead of local_dim, we assign the actual size of the array to write out
  gdim(1) = int(sizeprocs_adios,kind=8) * local_dim
  offs(1) = int(myrank_adios,kind=8) * local_dim

  full_name = trim(array_name) // '/array'

  ! defines variable  (adios2_constant_dims = .true. argument for constant dimensions, won't change)
  call adios2_define_variable(v, adios_group, trim(full_name), &
                              adios2_type_real4, 1, gdim, offs, ldim, adios2_constant_dims, ier)
  call check_adios_err(ier,"Error adios2 define variable "//trim(full_name)//" failed")

  if (.not. v%valid) stop 'Error adios2 defined variable is invalid'

  ! compression
  if (use_adios_compression) then
    ! adds compression operation to array variable
    call adios2_add_operation(operation_id, v, myadios_comp_operator, &
                              ADIOS_COMPRESSION_MODE, ADIOS_COMPRESSION_MODE_VALUE, ier)
    call check_adios_err(ier,"Error adios2 compression: add operation for variable "//trim(full_name)//" failed")
    if (operation_id /= 0) stop 'Error adios2 operation_id not added for array'
  endif

#endif

  group_size_inc = group_size_inc + local_dim * 4

  end subroutine define_adios_global_1d_generic_real

!===============================================================================
!> Define a global ADIOS 1D real array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_real_1d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  real, dimension(:), intent(in) :: var

#define VAR_TYPE real
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_real_1d


!===============================================================================
!> Define a global ADIOS 1D real array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_real_2d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  real, dimension(:,:), intent(in) :: var

#define VAR_TYPE real
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_real_2d


!===============================================================================
!> Define a global ADIOS 1D real array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_real_3d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  real, dimension(:,:,:), intent(in) :: var

#define VAR_TYPE real
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_real_3d


!===============================================================================
!> Define a global ADIOS 1D real array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_real_4d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  real, dimension(:,:,:,:), intent(in) :: var

#define VAR_TYPE real
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_real_4d


!===============================================================================
!> Define a global ADIOS 1D real array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_real_5d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  real, dimension(:,:,:,:,:), intent(in) :: var

#define VAR_TYPE real
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_real_5d


!===============================================================================
!> Define a double global array in ADIOS regardless of the array shape
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param local_dim The local dimension of the array.
  subroutine define_adios_global_1d_generic_double(adios_group, group_size_inc, array_name, local_dim, array_size)

  implicit none
  ! Parameters
#if defined(USE_ADIOS)
  integer(kind=8),  intent(in)    :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_io),  intent(in)     :: adios_group
#endif
  character(len=*), intent(in) :: array_name
  integer(kind=8),  intent(in) :: local_dim
  integer(kind=8),  intent(inout) :: group_size_inc
  integer(kind=8),  intent(in) :: array_size
  ! Variables
#if defined(USE_ADIOS)
  integer(kind=8) :: var_id
#elif defined(USE_ADIOS2)
  integer(kind=8), dimension(1) :: gdim  ! Global shape of array
  integer(kind=8), dimension(1) :: ldim  ! Local size of array
  integer(kind=8), dimension(1) :: offs  ! Starting offset in global array
  type(adios2_variable) :: v
  character(len=256) :: full_name
  ! compression
  integer :: operation_id
#endif
  integer :: ier

  TRACE_ADIOS_L2_ARG('define_adios_global_1d_generic_double: ',trim(array_name))

  ! Define the dimensions of the array. local_dim used as a dummy
  ! variable to call the integer routine.
  call define_adios_global_dims_1d(adios_group, group_size_inc, trim(array_name), local_dim)

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! we specify the actual size of the array in array_name/size which differs to array_name/local_dim for cases
  ! where ranks have different slice sizes.
  ! we will still use local_dim to calculate memory offsets for different ranks when reading arrays.
  call adios_define_var(adios_group, "array", trim(array_name), adios_double, &
                        trim(array_name) // "/size", &
                        trim(array_name) // "/global_dim", &
                        trim(array_name) // "/offset", var_id)

  ! old: assumes all ranks will have same array sizes
  !call adios_define_var(adios_group, "array", trim(array_name), adios_double, &
  !                      trim(array_name) // "/local_dim", &
  !                      trim(array_name) // "/global_dim", &
  !                      trim(array_name) // "/offset", var_id)

  ! compression
  if (use_adios_compression) then
    ! adds compression transform
    call adios_set_transform(var_id,myadios_comp_operator,ier)
    call check_adios_err(ier,"Error adios compression: set transform for variable "//trim(array_name)//" failed")
  endif

  ! to avoid compiler warning
  ier = array_size

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! note: variable parameter v will not be stored globally, instead will be retrieved by _inquire** function again
  ldim(1) = array_size   ! instead of local_dim, we assign the actual size of the array to write out
  gdim(1) = int(sizeprocs_adios,kind=8) * local_dim
  offs(1) = int(myrank_adios,kind=8) * local_dim

  full_name = trim(array_name) // '/array'

  ! defines variable  (adios2_constant_dims = .true. argument for constant dimensions, won't change)
  call adios2_define_variable(v, adios_group, trim(full_name), &
                              adios2_type_real8, 1, gdim, offs, ldim, adios2_constant_dims, ier)
  call check_adios_err(ier,"Error adios2 define variable "//trim(full_name)//" failed")

  if (.not. v%valid) stop 'Error adios2 defined variable is invalid'

  ! compression
  if (use_adios_compression) then
    ! adds compression operation to array variable
    call adios2_add_operation(operation_id, v, myadios_comp_operator, &
                              ADIOS_COMPRESSION_MODE, ADIOS_COMPRESSION_MODE_VALUE, ier)
    call check_adios_err(ier,"Error adios2 compression: add operation for variable "//trim(full_name)//" failed")
    if (operation_id /= 0) stop 'Error adios2 operation_id not added for array'
  endif

#endif

  group_size_inc = group_size_inc + local_dim * 8

  end subroutine define_adios_global_1d_generic_double

!===============================================================================
!> Define a global ADIOS 1D double array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_double_1d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  real(kind=8), dimension(:), intent(in) :: var

#define VAR_TYPE double
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_double_1d


!===============================================================================
!> Define a global ADIOS 1D double array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_double_2d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  real(kind=8), dimension(:,:), intent(in) :: var

#define VAR_TYPE double
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_double_2d


!===============================================================================
!> Define a global ADIOS 1D double array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_double_3d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  real(kind=8), dimension(:,:,:), intent(in) :: var

#define VAR_TYPE double
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_double_3d


!===============================================================================
!> Define a global ADIOS 1D double array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
subroutine define_adios_global_1d_double_4d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  real(kind=8), dimension(:,:,:,:), intent(in) :: var

#define VAR_TYPE double
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

end subroutine define_adios_global_1d_double_4d


!===============================================================================
!> Define a global ADIOS 1D double array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_double_5d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  real(kind=8), dimension(:,:,:,:,:), intent(in) :: var

#define VAR_TYPE double
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_double_5d


!===============================================================================
!> Define a integer global array in ADIOS regardless of the array shape
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param local_dim The local dimension of the array.
  subroutine define_adios_global_1d_generic_int(adios_group, group_size_inc, array_name, local_dim, array_size)

  implicit none
  ! Parameters
#if defined(USE_ADIOS)
  integer(kind=8),  intent(in)    :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_io),  intent(in)     :: adios_group
#endif
  character(len=*), intent(in) :: array_name
  integer(kind=8),  intent(in) :: local_dim
  integer(kind=8),  intent(inout) :: group_size_inc
  integer(kind=8),  intent(in) :: array_size
  ! Variables
#if defined(USE_ADIOS)
  integer(kind=8) :: var_id
#elif defined(USE_ADIOS2)
  integer(kind=8), dimension(1) :: gdim  ! Global shape of array
  integer(kind=8), dimension(1) :: ldim  ! Local size of array
  integer(kind=8), dimension(1) :: offs  ! Starting offset in global array
  type(adios2_variable) :: v
  character(len=256) :: full_name
#endif
  integer :: ier

  TRACE_ADIOS_L2_ARG('define_adios_global_1d_generic_int: ',trim(array_name))

  ! Define the dimensions of the array. local_dim used as a dummy
  ! variable to call the integer routine.
  call define_adios_global_dims_1d(adios_group, group_size_inc, trim(array_name), local_dim)

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! we specify the actual size of the array in array_name/size which differs to array_name/local_dim for cases
  ! where ranks have different slice sizes.
  ! we will still use local_dim to calculate memory offsets for different ranks when reading arrays.
  call adios_define_var(adios_group, "array", trim(array_name), adios_integer, &
                        trim(array_name) // "/size", &
                        trim(array_name) // "/global_dim", &
                        trim(array_name) // "/offset", var_id)

  ! old: assumes all ranks will have same array sizes
  !call adios_define_var(adios_group, "array", trim(array_name), adios_integer, &
  !                      trim(array_name) // "/local_dim", &
  !                      trim(array_name) // "/global_dim", &
  !                      trim(array_name) // "/offset", var_id)

  ! to avoid compiler warning
  ier = array_size

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! note: variable will not be stored globally, instead will be retrieved by _inquire** function again
  ldim(1) = array_size   ! instead of local_dim, we assign the actual size of the array to write out
  gdim(1) = int(sizeprocs_adios,kind=8) * local_dim
  offs(1) = int(myrank_adios,kind=8) * local_dim

  full_name = trim(array_name) // "/array"

  ! defines variable  (adios2_constant_dims = .true. argument for constant dimensions, won't change)
  call adios2_define_variable(v, adios_group, trim(full_name), &
                              adios2_type_integer4, 1, gdim, offs, ldim, adios2_constant_dims, ier)
  call check_adios_err(ier,"Error adios2 define variable "//trim(full_name)//" failed")

  if (.not. v%valid) stop 'Error adios2 defined variable is invalid'

#endif

  group_size_inc = group_size_inc + local_dim * 4

  end subroutine define_adios_global_1d_generic_int

!===============================================================================
!> Define a global ADIOS 1D int array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_int_1d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  integer(kind=4), dimension(:), intent(in) :: var

#define VAR_TYPE int
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_int_1d


!===============================================================================
!> Define a global ADIOS 1D int array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_int_2d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  integer(kind=4), dimension(:,:), intent(in) :: var

#define VAR_TYPE int
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_int_2d


!===============================================================================
!> Define a global ADIOS 1D int array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_int_3d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  integer(kind=4), dimension(:,:,:), intent(in) :: var

#define VAR_TYPE int
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_int_3d


!===============================================================================
!> Define a global ADIOS 1D int array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_int_4d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  integer(kind=4), dimension(:,:,:,:), intent(in) :: var

#define VAR_TYPE int
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_int_4d


!===============================================================================
!> Define a global ADIOS 1D int array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_int_5d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  integer(kind=4), dimension(:,:,:,:,:), intent(in) :: var

#define VAR_TYPE int
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_int_5d


!===============================================================================
!> Define a long integer global array in ADIOS regardless of the array shape
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param local_dim The local dimension of the array.
  subroutine define_adios_global_1d_generic_long(adios_group, group_size_inc, array_name, local_dim, array_size)

  implicit none
  ! Parameters
#if defined(USE_ADIOS)
  integer(kind=8),  intent(in)    :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_io),  intent(in)     :: adios_group
#endif
  character(len=*), intent(in) :: array_name
  integer(kind=8),  intent(in) :: local_dim
  integer(kind=8),  intent(inout) :: group_size_inc
  integer(kind=8),  intent(in) :: array_size
  ! Variables
#if defined(USE_ADIOS)
  integer(kind=8) :: var_id
#elif defined(USE_ADIOS2)
  integer(kind=8), dimension(1) :: gdim  ! Global shape of array
  integer(kind=8), dimension(1) :: ldim  ! Local size of array
  integer(kind=8), dimension(1) :: offs  ! Starting offset in global array
  type(adios2_variable) :: v
  character(len=256) :: full_name
#endif
  integer :: ier

  TRACE_ADIOS_L2_ARG('define_adios_global_1d_generic_long: ',trim(array_name))

  ! Define the dimensions of the array. local_dim used as a dummy
  ! variable to call the integer routine.
  call define_adios_global_dims_1d(adios_group, group_size_inc, trim(array_name), local_dim)

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! we specify the actual size of the array in array_name/size which differs to array_name/local_dim for cases
  ! where ranks have different slice sizes.
  ! we will still use local_dim to calculate memory offsets for different ranks when reading arrays.
  call adios_define_var(adios_group, "array", trim(array_name), adios_long, &
                        trim(array_name) // "/size", &
                        trim(array_name) // "/global_dim", &
                        trim(array_name) // "/offset", var_id)

  ! old: assumes all ranks will have same array sizes
  !call adios_define_var(adios_group, "array", trim(array_name), adios_long, &
  !                      trim(array_name) // "/local_dim", &
  !                      trim(array_name) // "/global_dim", &
  !                      trim(array_name) // "/offset", var_id)

  ! to avoid compiler warning
  ier = array_size
#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! note: variable will not be stored globally, instead will be retrieved by _inquire** function again
  ldim(1) = array_size   ! instead of local_dim, we assign the actual size of the array to write out
  gdim(1) = int(sizeprocs_adios,kind=8) * local_dim
  offs(1) = int(myrank_adios,kind=8) * local_dim

  full_name = trim(array_name) // '/array'

  ! defines variable  (adios2_constant_dims = .true. argument for constant dimensions, won't change)
  call adios2_define_variable(v, adios_group, trim(full_name), &
                              adios2_type_integer8, 1, gdim, offs, ldim, adios2_constant_dims, ier)
  call check_adios_err(ier,"Error adios2 define variable "//trim(full_name)//" failed")

  if (.not. v%valid) stop 'Error adios2 defined variable is invalid'

#endif

  group_size_inc = group_size_inc + local_dim * 8

  end subroutine define_adios_global_1d_generic_long

!===============================================================================
!> Define a global ADIOS 1D long array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_long_1d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  integer(kind=8), dimension(:), intent(in) :: var

#define VAR_TYPE long
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_long_1d


!===============================================================================
!> Define a global ADIOS 1D long array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_long_2d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  integer(kind=8), dimension(:,:), intent(in) :: var

#define VAR_TYPE long
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_long_2d


!===============================================================================
!> Define a global ADIOS 1D int array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_long_3d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  integer(kind=8), dimension(:,:,:), intent(in) :: var

#define VAR_TYPE long
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_long_3d


!===============================================================================
!> Define a global ADIOS 1D int array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_long_4d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  integer(kind=8), dimension(:,:,:,:), intent(in) :: var

#define VAR_TYPE long
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_long_4d


!===============================================================================
!> Define a global ADIOS 1D long array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_long_5d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  integer(kind=8), dimension(:,:,:,:,:), intent(in) :: var

#define VAR_TYPE long
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_long_5d

!===============================================================================
!> Define a logical global array in ADIOS regardless of the array shape
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param local_dim The local dimension of the array.
subroutine define_adios_global_1d_generic_logical(adios_group, group_size_inc, array_name, local_dim, array_size)

  implicit none
  ! Parameters
#if defined(USE_ADIOS)
  integer(kind=8),  intent(in)    :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_io),  intent(in)     :: adios_group
#endif
  character(len=*), intent(in) :: array_name
  integer(kind=8),  intent(in) :: local_dim
  integer(kind=8),  intent(inout) :: group_size_inc
  integer(kind=8),  intent(in) :: array_size
  ! Variables
#if defined(USE_ADIOS)
  integer(kind=8) :: var_id
#elif defined(USE_ADIOS2)
  integer(kind=8), dimension(1) :: gdim  ! Global shape of array
  integer(kind=8), dimension(1) :: ldim  ! Local size of array
  integer(kind=8), dimension(1) :: offs  ! Starting offset in global array
  type(adios2_variable) :: v
  character(len=256) :: full_name
#endif
  integer :: ier

  TRACE_ADIOS_L2_ARG('define_adios_global_1d_generic_logical: ',trim(array_name))

  ! Define the dimensions of the array. local_dim used as a dummy
  ! variable to call the integer routine.
  call define_adios_global_dims_1d(adios_group, group_size_inc, trim(array_name), local_dim)

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! The Fortran standard does not specify how variables of LOGICAL type are
  ! represented, beyond requiring that LOGICAL variables of default kind
  ! have the same storage size as default INTEGER and REAL variables.
  ! Hence the 'adios_integer' (2) data type to store logical values

  ! we specify the actual size of the array in array_name/size which differs to array_name/local_dim for cases
  ! where ranks have different slice sizes.
  ! we will still use local_dim to calculate memory offsets for different ranks when reading arrays.
  call adios_define_var(adios_group, "array", trim(array_name), adios_integer, &
                        trim(array_name) // "/size", &
                        trim(array_name) // "/global_dim", &
                        trim(array_name) // "/offset", var_id)

  ! old: assumes all ranks will have same array sizes
  !call adios_define_var(adios_group, "array", trim(array_name), adios_integer, &
  !                      trim(array_name) // "/local_dim", &
  !                      trim(array_name) // "/global_dim", &
  !                      trim(array_name) // "/offset", var_id)

  ! to avoid compiler warning
  ier = array_size

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! note: variable parameter v will not be stored globally, instead will be retrieved by _inquire** function again
  !
  !       also, adios2 has no adios2_get()/adios2_put() routine for logicals.
  !       we need to use integer array instead for storing/reading.
  ldim(1) = array_size   ! instead of local_dim, we assign the actual size of the array to write out
  gdim(1) = int(sizeprocs_adios,kind=8) * local_dim
  offs(1) = int(myrank_adios,kind=8) * local_dim

  full_name = trim(array_name) // '/array'

  ! defines variable  (adios2_constant_dims = .true. argument for constant dimensions, won't change)
  call adios2_define_variable(v, adios_group, trim(full_name), &
                              adios2_type_integer4, 1, gdim, offs, ldim, adios2_constant_dims, ier)
  call check_adios_err(ier,"Error adios2 define variable "//trim(full_name)//" failed")

  if (.not. v%valid) stop 'Error adios2 defined variable is invalid'

#endif

  group_size_inc = group_size_inc + local_dim * 4

end subroutine define_adios_global_1d_generic_logical

!===============================================================================
!> Define a global ADIOS 1D logical array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_logical_1d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  logical, dimension(:), intent(in) :: var

#define VAR_TYPE logical
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_logical_1d


!===============================================================================
!> Define a global ADIOS 1D logical array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_logical_2d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  logical, dimension(:,:), intent(in) :: var

#define VAR_TYPE logical
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_logical_2d


!===============================================================================
!> Define a global ADIOS 1D logical array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_logical_3d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  logical, dimension(:,:,:), intent(in) :: var

#define VAR_TYPE logical
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_logical_3d


!===============================================================================
!> Define a global ADIOS 1D logical array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_logical_4d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  logical, dimension(:,:,:,:), intent(in) :: var

#define VAR_TYPE logical
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_logical_4d


!===============================================================================
!> Define a global ADIOS 1D logical array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_logical_5d(adios_group, group_size_inc,local_dim, path, array_name, var)

  implicit none
  ! Parameters
  logical, dimension(:,:,:,:,:), intent(in) :: var

#define VAR_TYPE logical
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_logical_5d

!===============================================================================

!string added
subroutine define_adios_global_1d_string_generic(adios_group, group_size_inc, array_name, local_dim, array_size)

  implicit none
  ! Parameters
#if defined(USE_ADIOS)
  integer(kind=8),  intent(in)    :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_io),  intent(in)     :: adios_group
#endif
  character(len=*), intent(in) :: array_name
  integer(kind=8),  intent(in) :: local_dim
  integer(kind=8),  intent(inout) :: group_size_inc
  integer(kind=8),  intent(in) :: array_size
  ! Variables
#if defined(USE_ADIOS)
  integer(kind=8) :: var_id
#elif defined(USE_ADIOS2)
  integer(kind=8), dimension(1) :: gdim  ! Global shape of array
  integer(kind=8), dimension(1) :: ldim  ! Local size of array
  integer(kind=8), dimension(1) :: offs  ! Starting offset in global array
  type(adios2_variable) :: v
  character(len=256) :: full_name
#endif
  integer :: ier

  TRACE_ADIOS_L2_ARG('define_adios_global_1d_string_generic: ',trim(array_name))

  ! Define the dimensions of the array. local_dim used as a dummy
  ! variable to call the integer routine.
  call define_adios_global_dims_1d(adios_group, group_size_inc, trim(array_name), local_dim)

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! we specify the actual size of the array in array_name/size which differs to array_name/local_dim for cases
  ! where ranks have different slice sizes.
  ! we will still use local_dim to calculate memory offsets for different ranks when reading arrays.
  call adios_define_var(adios_group, "array", trim(array_name), adios_string, &
                        trim(array_name) // "/size", &
                        trim(array_name) // "/global_dim", &
                        trim(array_name) // "/offset", var_id)

  ! old: assumes all ranks will have same array sizes
  !call adios_define_var(adios_group, "array", trim(array_name), adios_string, &
  !                      trim(array_name) // "/local_dim", &
  !                      trim(array_name) // "/global_dim", &
  !                      trim(array_name) // "/offset", var_id)

  ! to avoid compiler warning
  ier = array_size

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! note: variable parameter v will not be stored globally, instead will be retrieved by _inquire** function again
  !
  !       we define the array as a global array (adios2_constant_dims = .true.) over all processes,
  !       but specify the corresponding sub-array dimension for each process using the dims

  ! local array dimensions
  ldim(1) = array_size
  gdim(:) = int(sizeprocs_adios,kind=8) * local_dim
  offs(:) = int(myrank_adios,kind=8) * local_dim

  full_name = trim(array_name) // '/array'

  ! defines variable  (adios2_constant_dims = .true. argument for constant dimensions, won't change)
  call adios2_define_variable(v, adios_group, trim(full_name), &
                              adios2_type_string, 1, gdim, offs, ldim, adios2_constant_dims, ier)
  call check_adios_err(ier,"Error adios2 define variable "//trim(full_name)//" failed")

  if (.not. v%valid) stop 'Error adios2 defined variable is invalid'

#endif

  group_size_inc = group_size_inc + local_dim * 1

end subroutine define_adios_global_1d_string_generic

!===============================================================================

subroutine define_adios_global_1d_string_1d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
#if defined(USE_ADIOS)
  integer(kind=8),  intent(in)    :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_io),  intent(in)     :: adios_group
#endif
  character(len=*), intent(in) :: path, array_name
  integer(kind=8),  intent(in) :: local_dim
  integer(kind=8),  intent(inout) :: group_size_inc
  character(len=*), intent(in) :: var
  ! Local vars
  character(len=256) :: full_name
  integer(kind=8) :: array_size

  TRACE_ADIOS_L2_ARG('define_adios_global_1d_string_1d: ',trim(array_name))

  ! check
  if (len_trim(array_name) == 0) stop 'Error adios: invalid name in define_adios_global_1d_string_1d()'

  ! sets full variable name
  if (len_trim(path) == 0) then
    full_name = trim(array_name)
  else
    full_name = trim(path) // '/' // trim(array_name)
  endif

  ! debug
  !print *,"full name", trim(full_name),"local_dim:",local_dim

  ! size
  array_size = len(var)

  call define_adios_global_1d_string_generic(adios_group, group_size_inc, full_name, local_dim, array_size)

end subroutine define_adios_global_1d_string_1d

!
!------------
!

subroutine  define_adios_local_1d_string_1d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
#if defined(USE_ADIOS)
  integer(kind=8),  intent(in)    :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_io),  intent(in)     :: adios_group
#endif
  character(len=*), intent(in) :: path, array_name
  integer(kind=8),  intent(in) :: local_dim
  integer(kind=8),  intent(inout) :: group_size_inc
  character(len=*), intent(in) :: var
  ! Local
#if defined(USE_ADIOS)
  integer(kind=8) :: var_id
#elif defined(USE_ADIOS2)
  type(adios2_variable) :: v
  integer :: ier
#endif
  integer :: idummy
  character(len=256) :: full_name

  TRACE_ADIOS_L2_ARG('define_adios_local_1d_string_1d: ',trim(array_name))

  ! check
  if (len_trim(array_name) == 0) stop 'Error adios: invalid name in define_adios_local_1d_string_1d()'

  ! sets full variable name
  if (len_trim(path) == 0) then
    full_name = trim(array_name)
  else
    full_name = trim(path) // '/' // trim(array_name)
  endif

  !debug
  !print *,"in define local: full_name:", trim(full_name)

#if defined(USE_ADIOS)
  ! ADIOS 1
  call adios_define_var(adios_group, trim(array_name), trim(path), adios_string, '', '', '', var_id )

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! note: we won't store the variable object, but get it back by adios2_inquire_** for writing
  ! defines global variable
  call adios2_define_variable(v, adios_group, trim(full_name), adios2_type_string, ier)
  call check_adios_err(ier,"Error adios2 could not define parameter: "//trim(full_name))
  if (.not. v%valid) stop 'Error adios2 defined variable is invalid'

#endif

  group_size_inc = group_size_inc + local_dim * 1

  ! to avoid compiler warnings
  idummy = len(var)

end subroutine define_adios_local_1d_string_1d


!-------------------------------------------------------------------------------
!
! compression
!
!-------------------------------------------------------------------------------

  subroutine define_adios_compression()

  use constants, only: IMAIN
  use manager_adios, only: myrank_adios

  implicit none

  ! local parameters
  integer :: ier
#if defined(USE_ADIOS)
  character(len=128) :: transform
#endif

! note: We add compression operations before saving the undo att forward arrays to file disk.
!       this should help reducing file storage sizes and help mitigating I/O issues for larger runs.
!
!       SZ and ZFP are the two leading lossy compressors available to compress scientific data sets.
!       Their performance can vary however across different data sets:
!         https://www.osti.gov/servlets/purl/1657917
!
!       ZFP compression is lossy compression, and works well with high compression rates for floating-point arrays:
!         https://computing.llnl.gov/projects/zfp
!         https://github.com/LLNL/zfp
!         Please cite this reference if using ZFP:
!           Peter Lindstrom. Fixed-Rate Compressed Floating-Point Arrays.
!           IEEE Transactions on Visualization and Computer Graphics, 20(12):2674-2683, December 2014.
!           doi:10.1109/TVCG.2014.2346458
!
!       SZ compression would be another possibility for floating-point arrays:
!         https://szcompressor.org
!         https://github.com/szcompressor/SZ
!
!       LZ4 lossless compression (see ADIOS 1.13.1 manual)
!
!       ADIOS2 will assign the compression operation to arrays and use compression before writing/putting them to disk.
!       Reading back these arrays will not need explicit assignment of the compression operator anymore.
!       Thus, the only change is required in the saving routine.
!
! note: crust_mantle arrays will dominate file sizes, where in particular attenuation arrays will be largest
!       as a small example NEX=48:  displ/veloc/accel crust_mantle   ~ 11.21 MB
!                                                     inner_core     ~  0.19 MB
!                                                     outer_core     ~  0.31 MB
!
!                                   rotation          outer_core     ~  0.36 MB
!
!                                   attenuation       crust_mantle   ~ 34.76 MB
!                                                     inner_core     ~  0.51 MB
!
!       to avoid spending too much time for the compression operation on small arrays,
!       one could apply it only to the dominant ones.
!       at the moment, we apply compression to all arrays to shrink the snapshot file size as much as possible.
!       -> todo in future, re-evaluate and optimize...

  ! initializes flag
  use_adios_compression = .true.
  ier = 0

  ! defines compression operator
  select case(ADIOS_COMPRESSION_ALGORITHM)
  case (0)
    ! no compression
    use_adios_compression = .false.

  case (1)
    ! ZFP compression
    ! options: precision w/ value 12
#if defined(USE_ADIOS)
    ! ADIOS 1
    !transform = "zfp:precision=12"
    write(transform,'("zfp:",a,"=",a)') trim(ADIOS_COMPRESSION_MODE),trim(ADIOS_COMPRESSION_MODE_VALUE)
#elif defined(USE_ADIOS2)
    ! ADIOS 2
    ! options: name = 'CompressorZfp', type = 'zfp'
    call adios2_define_operator(myadios_comp_operator, myadios2_obj, 'CompressorZfp', 'zfp', ier)
#endif

  case (2)
    ! SZ compression
    ! options: relative error w/ value 1.e-4
#if defined(USE_ADIOS)
    ! ADIOS 1
    !transform = "sz:relative=0.0001"
    write(transform,'("sz:",a,"=",a)') trim(ADIOS_COMPRESSION_MODE),trim(ADIOS_COMPRESSION_MODE_VALUE)
#elif defined(USE_ADIOS2)
    ! ADIOS 2
    call adios2_define_operator(myadios_comp_operator, myadios2_obj, 'CompressorSZ', 'sz', ier)
#endif

  case (3)
    ! LZ4 compression (lossless)
#if defined(USE_ADIOS)
    ! ADIOS 1
    !transform = "lz4:threshold=4096,lvl=9"     ! see ADIOS-UsersManual-1.13.1, page 75
    write(transform,'("lz4:",a,"=",a)') trim(ADIOS_COMPRESSION_MODE),trim(ADIOS_COMPRESSION_MODE_VALUE)
#elif defined(USE_ADIOS2)
    ! ADIOS 2
    call adios2_define_operator(myadios_comp_operator, myadios2_obj, 'CompressorLZ4', 'lz4', ier)
#endif

  case default
    ! no compression
    stop 'Invalid ADIOS compression algorithm selection, please choose 0 == none / 1 == ZFP / 2 == SZ / 3 == LZ4 compression'
  end select

#if defined(USE_ADIOS2)
  ! checks returned result
  if (ier /= 0) then
    ! no compression supported, resets flag
    ! (might need to install ZFP/SZ/LZ4 library first and re-configure/compile ADIOS2)
    use_adios_compression = .false.
    myadios_comp_operator%valid = .false.
    ! user output
    if (myrank_adios == 0) then
      print *
      print *,'WARNING: Selected ADIOS2 compression is not supported by ADIOS2 library. Please check your installation.'
      print *,'         Run will proceed without compression...'
      print *
      ! file output
      write(IMAIN,*)
      if (ADIOS_COMPRESSION_ALGORITHM == 1) then
        write(IMAIN,*) 'ADIOS2 compression algorithm ZFP is not supported by ADIOS2 library. Please check your installation.'
      else if (ADIOS_COMPRESSION_ALGORITHM == 2) then
        write(IMAIN,*) 'ADIOS2 compression algorithm SZ  is not supported by ADIOS2 library. Please check your installation.'
      else if (ADIOS_COMPRESSION_ALGORITHM == 3) then
        write(IMAIN,*) 'ADIOS2 compression algorithm LZ4 is not supported by ADIOS2 library. Please check your installation.'
      endif
      write(IMAIN,*) 'Run will proceed without compression...'
      write(IMAIN,*)
      call flush_IMAIN()
    endif
    ! or safety stop
    !stop 'ADIOS2 compression not supported. Please set ADIOS_COMPRESSION_ALGORITHM == 0 in constants.h file'
  endif
#endif

  ! applies compression to undo_att defined array variables
  if (use_adios_compression) then
#if defined(USE_ADIOS)
    ! checks if transform valid
    if (len_trim(transform) == 0) &
      stop 'Invalid ADIOS transform in define_adios_compression()'

    ! sets compression string
    myadios_comp_operator = trim(transform)

#elif defined(USE_ADIOS2)
    ! check if operator valid
    if (.not. myadios_comp_operator%valid) &
      stop 'Invalid ADIOS2 compression operator handle returned in define_adios_compression()'
#endif

    ! user output
    if (myrank_adios == 0) then
      write(IMAIN,*)
      write(IMAIN,*) "undoing attenuation:"
#if defined(USE_ADIOS)
      write(IMAIN,*) "  adding ADIOS compression operation for snapshot wavefield storage"
#elif defined(USE_ADIOS2)
      write(IMAIN,*) "  adding ADIOS2 compression operation for snapshot wavefield storage"
#endif
      if (ADIOS_COMPRESSION_ALGORITHM == 1) then
#if defined(USE_ADIOS)
        write(IMAIN,*) "  ZFP compression mode: ",trim(myadios_comp_operator)
#elif defined(USE_ADIOS2)
        write(IMAIN,*) "  ZFP compression mode: {'",trim(ADIOS_COMPRESSION_MODE),"','",trim(ADIOS_COMPRESSION_MODE_VALUE),"'}"
#endif
      else if (ADIOS_COMPRESSION_ALGORITHM == 2) then
#if defined(USE_ADIOS)
        write(IMAIN,*) "  SZ compression  mode: ",trim(myadios_comp_operator)
#elif defined(USE_ADIOS2)
        write(IMAIN,*) "  SZ  compression mode: {'",trim(ADIOS_COMPRESSION_MODE),"','",trim(ADIOS_COMPRESSION_MODE_VALUE),"'}"
#endif
      else if (ADIOS_COMPRESSION_ALGORITHM == 3) then
#if defined(USE_ADIOS)
        write(IMAIN,*) "  LZ4 compression mode: ",trim(myadios_comp_operator)
#elif defined(USE_ADIOS2)
        write(IMAIN,*) "  LZ4 compression mode: {'",trim(ADIOS_COMPRESSION_MODE),"','",trim(ADIOS_COMPRESSION_MODE_VALUE),"'}"
#endif
      endif
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

  end subroutine define_adios_compression


end module adios_helpers_definitions_mod
