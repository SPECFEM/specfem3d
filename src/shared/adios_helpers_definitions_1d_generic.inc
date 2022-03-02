! file "adios_helpers_definitions_1d_generic.inc"
!
! generic type
!
! common include section for adios define 1d helper functions

#if _CRAYFTN
! Cray preprocessor has problems with the default macro
! work-around see: https://stackoverflow.com/questions/1489932/..
!                  ../how-to-concatenate-twice-with-the-c-preprocessor-and-expand-a-macro-as-in-arg/1489985#1489985
! to check: compiler with ftn -eP .. and look at processed source file.i
!!#warning "Compiling with Cray compiler"
#define CONCAT(a,b) a/**/b
#define EVAL(x,y) CONCAT(CONCAT(x,_),y)
#define FUNC_DEFINE_1D_GENERIC(a) EVAL(define_adios_global_1d_generic,a)
#else
! standard compiler
#define STR(a) a
#define CONCAT(a,b) STR(a)STR(b)
#define FUNC_DEFINE_1D_GENERIC(a) CONCAT(define_adios_global_1d_generic_,a)
#endif

  ! common parameters
#if defined(USE_ADIOS)
  integer(kind=8),  intent(in)    :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_io),  intent(in)     :: adios_group
#endif
  character(len=*), intent(in) :: path, array_name
  integer(kind=8),  intent(in) :: local_dim
  integer(kind=8),  intent(inout) :: group_size_inc
  ! Local vars
  integer :: idummy
  character(len=256) :: full_name

  TRACE_ADIOS_L2_ARG('define_adios_global_1d adios_helpers_definitions_1d_generic: ',trim(path) // trim(array_name))

  ! sets full variable name
  full_name = trim(path) // trim(array_name)

  ! checks name
  if (len_trim(full_name) == 0) stop 'Error adios: invalid full_name in adios_helpers_definitions()'
  if (len_trim(array_name) == 0) then
    print *,'Error adios: invalid path in adios_helpers_definitions, provides only: ',trim(full_name)
    stop 'Error adios: invalid path in adios_helpers_definitions()'
  endif

!! in case preprocessor macro fails with specific compiler...
!!#if _CRAYFTN
!!  ! Cray compiler has problems with the macro
!!#warning "Compiling with Cray compiler"
!!  ! expands call to individual cases
!!#ifdef VAR_TYPE == int
!!  call define_adios_global_1d_generic_int(adios_group, group_size_inc, trim(full_name), local_dim)
!!#endif
!!#ifdef VAR_TYPE == long
!!  call define_adios_global_1d_generic_long(adios_group, group_size_inc, trim(full_name), local_dim)
!!#endif
!!#ifdef VAR_TYPE == real
!!  call define_adios_global_1d_generic_real(adios_group, group_size_inc, trim(full_name), local_dim)
!!#endif
!!#ifdef VAR_TYPE == double
!!  call define_adios_global_1d_generic_double(adios_group, group_size_inc, trim(full_name), local_dim)
!!#endif
!!#ifdef VAR_TYPE == logical
!!  call define_adios_global_1d_generic_logical(adios_group, group_size_inc, trim(full_name), local_dim)
!!#endif
!!#else
!!  ! other compilers can use the macro definition
!!  call FUNC_DEFINE_1D_GENERIC(VAR_TYPE) (adios_group, group_size_inc, trim(full_name), local_dim)
!!#endif

  call FUNC_DEFINE_1D_GENERIC(VAR_TYPE) (adios_group, group_size_inc, trim(full_name), local_dim)

  ! to avoid compiler warnings
  idummy = size(var)

