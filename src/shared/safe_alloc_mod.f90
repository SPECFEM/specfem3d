

!==============================================================================
!> \file  safe_alloc_mod.f90
!! \brief Helper module for (de)allocation error check.
!!
!! \author MPBL
!==============================================================================

!==============================================================================
!> Helpers (de)allocate arrays of various types and dimensions. Check errors.
!!
!! \note The routine of this module only allow to allocate arrays starting at
!!       index 1. If you want custom bound, write your own allocation. You can
!!       still use the error checking routines, though.
!!
!! \author MPBL
!------------------------------------------------------------------------------
module safe_alloc_mod

  implicit none

  private

  public :: safe_alloc
  public :: safe_dealloc

  public :: check_alloc_err
  public :: check_dealloc_err

  interface safe_alloc
    module procedure safe_alloc_float_1d
    module procedure safe_alloc_float_2d
    module procedure safe_alloc_float_3d
    module procedure safe_alloc_float_4d
    module procedure safe_alloc_float_5d

    module procedure safe_alloc_double_1d
    module procedure safe_alloc_double_2d
    module procedure safe_alloc_double_3d
    module procedure safe_alloc_double_4d
    module procedure safe_alloc_double_5d

    module procedure safe_alloc_int_1d
    module procedure safe_alloc_int_2d
    module procedure safe_alloc_int_3d
    module procedure safe_alloc_int_4d
    module procedure safe_alloc_int_5d

    module procedure safe_alloc_long_1d
    module procedure safe_alloc_long_2d
    module procedure safe_alloc_long_3d
    module procedure safe_alloc_long_4d
    module procedure safe_alloc_long_5d

    !module procedure safe_alloc_byte_1d
    !module procedure safe_alloc_byte_2d
    !module procedure safe_alloc_byte_3d
    !module procedure safe_alloc_byte_4d
    !module procedure safe_alloc_byte_5d

    module procedure safe_alloc_logical_1d
    module procedure safe_alloc_logical_2d
    module procedure safe_alloc_logical_3d
    module procedure safe_alloc_logical_4d
    module procedure safe_alloc_logical_5d
  end interface safe_alloc

  interface safe_dealloc
    module procedure safe_dealloc_float_1d
    module procedure safe_dealloc_float_2d
    module procedure safe_dealloc_float_3d
    module procedure safe_dealloc_float_4d
    module procedure safe_dealloc_float_5d

    module procedure safe_dealloc_double_1d
    module procedure safe_dealloc_double_2d
    module procedure safe_dealloc_double_3d
    module procedure safe_dealloc_double_4d
    module procedure safe_dealloc_double_5d

    module procedure safe_dealloc_int_1d
    module procedure safe_dealloc_int_2d
    module procedure safe_dealloc_int_3d
    module procedure safe_dealloc_int_4d
    module procedure safe_dealloc_int_5d

    module procedure safe_dealloc_long_1d
    module procedure safe_dealloc_long_2d
    module procedure safe_dealloc_long_3d
    module procedure safe_dealloc_long_4d
    module procedure safe_dealloc_long_5d

    !module procedure safe_dealloc_byte_1d
    !module procedure safe_dealloc_byte_2d
    !module procedure safe_dealloc_byte_3d
    !module procedure safe_dealloc_byte_4d
    !module procedure safe_dealloc_byte_5d

    module procedure safe_dealloc_logical_1d
    module procedure safe_dealloc_logical_2d
    module procedure safe_dealloc_logical_3d
    module procedure safe_dealloc_logical_4d
    module procedure safe_dealloc_logical_5d
  end interface safe_dealloc

contains

!==============================================================================
!>
subroutine check_alloc_err(ier, usr_msg)
  use iso_fortran_env, only : error_unit
  use mpi

  integer, intent(in) :: ier
  character(len=*), intent(in), optional :: usr_msg

  integer :: myrank, mpi_er

  if(ier /= 0) then
    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, mpi_er)
    if (present(usr_msg)) then
      write(error_unit, "('Process ', i6.6, " // &
                        "': Allocation error. ', A)") myrank, usr_msg
    else
      write(error_unit, "('Process ', i6.6, ': Allocation error. " // &
                        " No user message specified.')") myrank
    endif
    !call exit(ier)
    call MPI_Abort(MPI_COMM_WORLD, ier, mpi_er)
  endif
end subroutine check_alloc_err

!==============================================================================
!>
subroutine check_dealloc_err(ier, usr_msg)
  use iso_fortran_env, only : error_unit
  use mpi

  integer, intent(in) :: ier
  character(len=*), intent(in), optional :: usr_msg

  integer :: myrank, mpi_er

  if(ier /= 0) then
    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, mpi_er)
    if (present(usr_msg)) then
      write(error_unit, "('Process ', i6.6, " // &
                        "': Deallocation error. ', A)") myrank, usr_msg
    else
      write(error_unit, "('Process ', i6.6, ': Deallocation error. " // &
                        " No user message specified.')") myrank
    endif
    !call exit(ier)
    call MPI_Abort(MPI_COMM_WORLD, ier, mpi_er)
  endif
end subroutine check_dealloc_err

!==============================================================================
!> Allocate a 1D float array and check for errors
!! \param array The array to allocate
!! \param dim1 The dimension of the array
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_alloc_float_1d(array, dim1, usr_msg)
  ! Arguments
  real(kind=4), dimension(:), allocatable, intent(inout) :: array
  integer, intent(in) :: dim1
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  allocate(array(dim1), stat=ier)
  call check_alloc_err(ier, usr_msg)
end subroutine safe_alloc_float_1d

!==============================================================================
!> Allocate a 2D float array and check for errors
!! \param array The array to allocate
!! \param dim1 The 1st dimension of the array
!! \param dim2 The 2nd dimension of the array
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_alloc_float_2d(array, dim1, dim2, usr_msg)
  ! Arguments
  real(kind=4), dimension(:,:), allocatable, intent(inout) :: array
  integer, intent(in) :: dim1, dim2
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  allocate(array(dim1, dim2), stat=ier)
  call check_alloc_err(ier, usr_msg)
end subroutine safe_alloc_float_2d

!==============================================================================
!> Allocate a 3D float array and check for errors
!! \param array The array to allocate
!! \param dim1 The 1st dimension of the array
!! \param dim2 The 2nd dimension of the array
!! \param dim3 The 3rd dimension of the array
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_alloc_float_3d(array, dim1, dim2, dim3, usr_msg)
  ! Arguments
  real(kind=4), dimension(:,:,:), allocatable, intent(inout) :: array
  integer, intent(in) :: dim1, dim2, dim3
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  allocate(array(dim1, dim2, dim3), stat=ier)
  call check_alloc_err(ier, usr_msg)
end subroutine safe_alloc_float_3d

!==============================================================================
!> Allocate a 4D float array and check for errors
!! \param array The array to allocate
!! \param dim1 The 1st dimension of the array
!! \param dim2 The 2nd dimension of the array
!! \param dim3 The 3rd dimension of the array
!! \param dim4 The 4th dimension of the array
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_alloc_float_4d(array, dim1, dim2, dim3, dim4, usr_msg)
  ! Arguments
  real(kind=4), dimension(:,:,:,:), allocatable, intent(inout) :: array
  integer, intent(in) :: dim1, dim2, dim3, dim4
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  allocate(array(dim1, dim2, dim3, dim4), stat=ier)
  call check_alloc_err(ier, usr_msg)
end subroutine safe_alloc_float_4d

!==============================================================================
!> Allocate a 1D float array and check for errors
!! \param array The array to allocate
!! \param dim1 The 1st dimension of the array
!! \param dim2 The 2nd dimension of the array
!! \param dim3 The 3rd dimension of the array
!! \param dim4 The 4th dimension of the array
!! \param dim5 The 5th dimension of the array
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_alloc_float_5d(array, dim1, dim2, dim3, dim4, dim5, usr_msg)
  ! Arguments
  real(kind=4), dimension(:,:,:,:,:), allocatable, intent(inout) :: array
  integer, intent(in) :: dim1, dim2, dim3, dim4, dim5
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  allocate(array(dim1, dim2, dim3, dim4, dim5), stat=ier)
  call check_alloc_err(ier, usr_msg)
end subroutine safe_alloc_float_5d

!==============================================================================
!> Allocate a 1D double array and check for errors
!! \param array The array to allocate
!! \param dim1 The dimension of the array
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_alloc_double_1d(array, dim1, usr_msg)
  ! Arguments
  real(kind=8), dimension(:), allocatable, intent(inout) :: array
  integer, intent(in) :: dim1
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  allocate(array(dim1), stat=ier)
  call check_alloc_err(ier, usr_msg)
end subroutine safe_alloc_double_1d

!==============================================================================
!> Allocate a 2D double array and check for errors
!! \param array The array to allocate
!! \param dim1 The 1st dimension of the array
!! \param dim2 The 2nd dimension of the array
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_alloc_double_2d(array, dim1, dim2, usr_msg)
  ! Arguments
  real(kind=8), dimension(:,:), allocatable, intent(inout) :: array
  integer, intent(in) :: dim1, dim2
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  allocate(array(dim1, dim2), stat=ier)
  call check_alloc_err(ier, usr_msg)
end subroutine safe_alloc_double_2d

!==============================================================================
!> Allocate a 3D double array and check for errors
!! \param array The array to allocate
!! \param dim1 The 1st dimension of the array
!! \param dim2 The 2nd dimension of the array
!! \param dim3 The 3rd dimension of the array
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_alloc_double_3d(array, dim1, dim2, dim3, usr_msg)
  ! Arguments
  real(kind=8), dimension(:,:,:), allocatable, intent(inout) :: array
  integer, intent(in) :: dim1, dim2, dim3
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  allocate(array(dim1, dim2, dim3), stat=ier)
  call check_alloc_err(ier, usr_msg)
end subroutine safe_alloc_double_3d

!==============================================================================
!> Allocate a 4D double array and check for errors
!! \param array The array to allocate
!! \param dim1 The 1st dimension of the array
!! \param dim2 The 2nd dimension of the array
!! \param dim3 The 3rd dimension of the array
!! \param dim4 The 4th dimension of the array
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_alloc_double_4d(array, dim1, dim2, dim3, dim4, usr_msg)
  ! Arguments
  real(kind=8), dimension(:,:,:,:), allocatable, intent(inout) :: array
  integer, intent(in) :: dim1, dim2, dim3, dim4
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  allocate(array(dim1, dim2, dim3, dim4), stat=ier)
  call check_alloc_err(ier, usr_msg)
end subroutine safe_alloc_double_4d

!==============================================================================
!> Allocate a 1D double array and check for errors
!! \param array The array to allocate
!! \param dim1 The 1st dimension of the array
!! \param dim2 The 2nd dimension of the array
!! \param dim3 The 3rd dimension of the array
!! \param dim4 The 4th dimension of the array
!! \param dim5 The 5th dimension of the array
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_alloc_double_5d(array, dim1, dim2, dim3, dim4, dim5, usr_msg)
  ! Arguments
  real(kind=8), dimension(:,:,:,:,:), allocatable, intent(inout) :: array
  integer, intent(in) :: dim1, dim2, dim3, dim4, dim5
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  allocate(array(dim1, dim2, dim3, dim4, dim5), stat=ier)
  call check_alloc_err(ier, usr_msg)
end subroutine safe_alloc_double_5d

!==============================================================================
!> Allocate a 1D int array and check for errors
!! \param array The array to allocate
!! \param dim1 The dimension of the array
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_alloc_int_1d(array, dim1, usr_msg)
  ! Arguments
  integer(kind=4), dimension(:), allocatable, intent(inout) :: array
  integer, intent(in) :: dim1
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  allocate(array(dim1), stat=ier)
  call check_alloc_err(ier, usr_msg)
end subroutine safe_alloc_int_1d

!==============================================================================
!> Allocate a 2D int array and check for errors
!! \param array The array to allocate
!! \param dim1 The 1st dimension of the array
!! \param dim2 The 2nd dimension of the array
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_alloc_int_2d(array, dim1, dim2, usr_msg)
  ! Arguments
  integer(kind=4), dimension(:,:), allocatable, intent(inout) :: array
  integer, intent(in) :: dim1, dim2
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  allocate(array(dim1, dim2), stat=ier)
  call check_alloc_err(ier, usr_msg)
end subroutine safe_alloc_int_2d

!==============================================================================
!> Allocate a 3D int array and check for errors
!! \param array The array to allocate
!! \param dim1 The 1st dimension of the array
!! \param dim2 The 2nd dimension of the array
!! \param dim3 The 3rd dimension of the array
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_alloc_int_3d(array, dim1, dim2, dim3, usr_msg)
  ! Arguments
  integer(kind=4), dimension(:,:,:), allocatable, intent(inout) :: array
  integer, intent(in) :: dim1, dim2, dim3
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  allocate(array(dim1, dim2, dim3), stat=ier)
  call check_alloc_err(ier, usr_msg)
end subroutine safe_alloc_int_3d

!==============================================================================
!> Allocate a 4D int array and check for errors
!! \param array The array to allocate
!! \param dim1 The 1st dimension of the array
!! \param dim2 The 2nd dimension of the array
!! \param dim3 The 3rd dimension of the array
!! \param dim4 The 4th dimension of the array
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_alloc_int_4d(array, dim1, dim2, dim3, dim4, usr_msg)
  ! Arguments
  integer(kind=4), dimension(:,:,:,:), allocatable, intent(inout) :: array
  integer, intent(in) :: dim1, dim2, dim3, dim4
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  allocate(array(dim1, dim2, dim3, dim4), stat=ier)
  call check_alloc_err(ier, usr_msg)
end subroutine safe_alloc_int_4d

!==============================================================================
!> Allocate a 1D int array and check for errors
!! \param array The array to allocate
!! \param dim1 The 1st dimension of the array
!! \param dim2 The 2nd dimension of the array
!! \param dim3 The 3rd dimension of the array
!! \param dim4 The 4th dimension of the array
!! \param dim5 The 5th dimension of the array
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_alloc_int_5d(array, dim1, dim2, dim3, dim4, dim5, usr_msg)
  ! Arguments
  integer(kind=4), dimension(:,:,:,:,:), allocatable, intent(inout) :: array
  integer, intent(in) :: dim1, dim2, dim3, dim4, dim5
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  allocate(array(dim1, dim2, dim3, dim4, dim5), stat=ier)
  call check_alloc_err(ier, usr_msg)
end subroutine safe_alloc_int_5d

!==============================================================================
!> Allocate a 1D long array and check for errors
!! \param array The array to allocate
!! \param dim1 The dimension of the array
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_alloc_long_1d(array, dim1, usr_msg)
  ! Arguments
  integer(kind=8), dimension(:), allocatable, intent(inout) :: array
  integer, intent(in) :: dim1
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  allocate(array(dim1), stat=ier)
  call check_alloc_err(ier, usr_msg)
end subroutine safe_alloc_long_1d

!==============================================================================
!> Allocate a 2D long array and check for errors
!! \param array The array to allocate
!! \param dim1 The 1st dimension of the array
!! \param dim2 The 2nd dimension of the array
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_alloc_long_2d(array, dim1, dim2, usr_msg)
  ! Arguments
  integer(kind=8), dimension(:,:), allocatable, intent(inout) :: array
  integer, intent(in) :: dim1, dim2
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  allocate(array(dim1, dim2), stat=ier)
  call check_alloc_err(ier, usr_msg)
end subroutine safe_alloc_long_2d

!==============================================================================
!> Allocate a 3D long array and check for errors
!! \param array The array to allocate
!! \param dim1 The 1st dimension of the array
!! \param dim2 The 2nd dimension of the array
!! \param dim3 The 3rd dimension of the array
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_alloc_long_3d(array, dim1, dim2, dim3, usr_msg)
  ! Arguments
  integer(kind=8), dimension(:,:,:), allocatable, intent(inout) :: array
  integer, intent(in) :: dim1, dim2, dim3
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  allocate(array(dim1, dim2, dim3), stat=ier)
  call check_alloc_err(ier, usr_msg)
end subroutine safe_alloc_long_3d

!==============================================================================
!> Allocate a 4D long array and check for errors
!! \param array The array to allocate
!! \param dim1 The 1st dimension of the array
!! \param dim2 The 2nd dimension of the array
!! \param dim3 The 3rd dimension of the array
!! \param dim4 The 4th dimension of the array
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_alloc_long_4d(array, dim1, dim2, dim3, dim4, usr_msg)
  ! Arguments
  integer(kind=8), dimension(:,:,:,:), allocatable, intent(inout) :: array
  integer, intent(in) :: dim1, dim2, dim3, dim4
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  allocate(array(dim1, dim2, dim3, dim4), stat=ier)
  call check_alloc_err(ier, usr_msg)
end subroutine safe_alloc_long_4d

!==============================================================================
!> Allocate a 1D long array and check for errors
!! \param array The array to allocate
!! \param dim1 The 1st dimension of the array
!! \param dim2 The 2nd dimension of the array
!! \param dim3 The 3rd dimension of the array
!! \param dim4 The 4th dimension of the array
!! \param dim5 The 5th dimension of the array
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_alloc_long_5d(array, dim1, dim2, dim3, dim4, dim5, usr_msg)
  ! Arguments
  integer(kind=8), dimension(:,:,:,:,:), allocatable, intent(inout) :: array
  integer, intent(in) :: dim1, dim2, dim3, dim4, dim5
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  allocate(array(dim1, dim2, dim3, dim4, dim5), stat=ier)
  call check_alloc_err(ier, usr_msg)
end subroutine safe_alloc_long_5d

!!==============================================================================
!!> Allocate a 1D byte array and check for errors
!!! \param array The array to allocate
!!! \param dim1 The dimension of the array
!!! \param usr_msg A custom error message to print if allocation fails.
!subroutine safe_alloc_byte_1d(array, dim1, usr_msg)
  !! Arguments
  !byte, dimension(:), allocatable, intent(inout) :: array
  !integer, intent(in) :: dim1
  !character(len=*), intent(in), optional :: usr_msg
  !! Local parameters
  !integer :: ier
!
  !allocate(array(dim1), stat=ier)
  !call check_alloc_err(ier, usr_msg)
!end subroutine safe_alloc_byte_1d
!
!!==============================================================================
!!> Allocate a 2D byte array and check for errors
!!! \param array The array to allocate
!!! \param dim1 The 1st dimension of the array
!!! \param dim2 The 2nd dimension of the array
!!! \param usr_msg A custom error message to print if allocation fails.
!subroutine safe_alloc_byte_2d(array, dim1, dim2, usr_msg)
  !! Arguments
  !byte, dimension(:,:), allocatable, intent(inout) :: array
  !integer, intent(in) :: dim1, dim2
  !character(len=*), intent(in), optional :: usr_msg
  !! Local parameters
  !integer :: ier
!
  !allocate(array(dim1, dim2), stat=ier)
  !call check_alloc_err(ier, usr_msg)
!end subroutine safe_alloc_byte_2d
!
!!==============================================================================
!!> Allocate a 3D byte array and check for errors
!!! \param array The array to allocate
!!! \param dim1 The 1st dimension of the array
!!! \param dim2 The 2nd dimension of the array
!!! \param dim3 The 3rd dimension of the array
!!! \param usr_msg A custom error message to print if allocation fails.
!subroutine safe_alloc_byte_3d(array, dim1, dim2, dim3, usr_msg)
  !! Arguments
  !byte, dimension(:,:,:), allocatable, intent(inout) :: array
  !integer, intent(in) :: dim1, dim2, dim3
  !character(len=*), intent(in), optional :: usr_msg
  !! Local parameters
  !integer :: ier
!
  !allocate(array(dim1, dim2, dim3), stat=ier)
  !call check_alloc_err(ier, usr_msg)
!end subroutine safe_alloc_byte_3d
!
!!==============================================================================
!!> Allocate a 4D byte array and check for errors
!!! \param array The array to allocate
!!! \param dim1 The 1st dimension of the array
!!! \param dim2 The 2nd dimension of the array
!!! \param dim3 The 3rd dimension of the array
!!! \param dim4 The 4th dimension of the array
!!! \param usr_msg A custom error message to print if allocation fails.
!subroutine safe_alloc_byte_4d(array, dim1, dim2, dim3, dim4, usr_msg)
  !! Arguments
  !byte, dimension(:,:,:,:), allocatable, intent(inout) :: array
  !integer, intent(in) :: dim1, dim2, dim3, dim4
  !character(len=*), intent(in), optional :: usr_msg
  !! Local parameters
  !integer :: ier
!
  !allocate(array(dim1, dim2, dim3, dim4), stat=ier)
  !call check_alloc_err(ier, usr_msg)
!end subroutine safe_alloc_byte_4d
!
!!==============================================================================
!!> Allocate a 1D byte array and check for errors
!!! \param array The array to allocate
!!! \param dim1 The 1st dimension of the array
!!! \param dim2 The 2nd dimension of the array
!!! \param dim3 The 3rd dimension of the array
!!! \param dim4 The 4th dimension of the array
!!! \param dim5 The 5th dimension of the array
!!! \param usr_msg A custom error message to print if allocation fails.
!subroutine safe_alloc_byte_5d(array, dim1, dim2, dim3, dim4, dim5, usr_msg)
  !! Arguments
  !byte, dimension(:,:,:,:,:), allocatable, intent(inout) :: array
  !integer, intent(in) :: dim1, dim2, dim3, dim4, dim5
  !character(len=*), intent(in), optional :: usr_msg
  !! Local parameters
  !integer :: ier
!
  !allocate(array(dim1, dim2, dim3, dim4, dim5), stat=ier)
  !call check_alloc_err(ier, usr_msg)
!end subroutine safe_alloc_byte_5d

!==============================================================================
!> Allocate a 1D logical array and check for errors
!! \param array The array to allocate
!! \param dim1 The dimension of the array
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_alloc_logical_1d(array, dim1, usr_msg)
  ! Arguments
  logical, dimension(:), allocatable, intent(inout) :: array
  integer, intent(in) :: dim1
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  allocate(array(dim1), stat=ier)
  call check_alloc_err(ier, usr_msg)
end subroutine safe_alloc_logical_1d

!==============================================================================
!> Allocate a 2D logical array and check for errors
!! \param array The array to allocate
!! \param dim1 The 1st dimension of the array
!! \param dim2 The 2nd dimension of the array
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_alloc_logical_2d(array, dim1, dim2, usr_msg)
  ! Arguments
  logical, dimension(:,:), allocatable, intent(inout) :: array
  integer, intent(in) :: dim1, dim2
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  allocate(array(dim1, dim2), stat=ier)
  call check_alloc_err(ier, usr_msg)
end subroutine safe_alloc_logical_2d

!==============================================================================
!> Allocate a 3D logical array and check for errors
!! \param array The array to allocate
!! \param dim1 The 1st dimension of the array
!! \param dim2 The 2nd dimension of the array
!! \param dim3 The 3rd dimension of the array
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_alloc_logical_3d(array, dim1, dim2, dim3, usr_msg)
  ! Arguments
  logical, dimension(:,:,:), allocatable, intent(inout) :: array
  integer, intent(in) :: dim1, dim2, dim3
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  allocate(array(dim1, dim2, dim3), stat=ier)
  call check_alloc_err(ier, usr_msg)
end subroutine safe_alloc_logical_3d

!==============================================================================
!> Allocate a 4D logical array and check for errors
!! \param array The array to allocate
!! \param dim1 The 1st dimension of the array
!! \param dim2 The 2nd dimension of the array
!! \param dim3 The 3rd dimension of the array
!! \param dim4 The 4th dimension of the array
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_alloc_logical_4d(array, dim1, dim2, dim3, dim4, usr_msg)
  ! Arguments
  logical, dimension(:,:,:,:), allocatable, intent(inout) :: array
  integer, intent(in) :: dim1, dim2, dim3, dim4
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  allocate(array(dim1, dim2, dim3, dim4), stat=ier)
  call check_alloc_err(ier, usr_msg)
end subroutine safe_alloc_logical_4d

!==============================================================================
!> Allocate a 1D logical array and check for errors
!! \param array The array to allocate
!! \param dim1 The 1st dimension of the array
!! \param dim2 The 2nd dimension of the array
!! \param dim3 The 3rd dimension of the array
!! \param dim4 The 4th dimension of the array
!! \param dim5 The 5th dimension of the array
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_alloc_logical_5d(array, dim1, dim2, dim3, dim4, dim5, usr_msg)
  ! Arguments
  logical, dimension(:,:,:,:,:), allocatable, intent(inout) :: array
  integer, intent(in) :: dim1, dim2, dim3, dim4, dim5
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  allocate(array(dim1, dim2, dim3, dim4, dim5), stat=ier)
  call check_alloc_err(ier, usr_msg)
end subroutine safe_alloc_logical_5d

!==============================================================================
!> Deallocate a 1D float array and check for errors
!! \param array The array to allocate
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_dealloc_float_1d(array, usr_msg)
  ! Arguments
  real(kind=4), dimension(:), allocatable, intent(inout) :: array
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  deallocate(array, stat=ier)
  call check_dealloc_err(ier, usr_msg)
end subroutine safe_dealloc_float_1d

!==============================================================================
!> Deallocate a 2D float array and check for errors
!! \param array The array to allocate
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_dealloc_float_2d(array, usr_msg)
  ! Arguments
  real(kind=4), dimension(:,:), allocatable, intent(inout) :: array
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  deallocate(array, stat=ier)
  call check_dealloc_err(ier, usr_msg)
end subroutine safe_dealloc_float_2d

!==============================================================================
!> Deallocate a 3D float array and check for errors
!! \param array The array to allocate
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_dealloc_float_3d(array, usr_msg)
  ! Arguments
  real(kind=4), dimension(:,:,:), allocatable, intent(inout) :: array
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  deallocate(array, stat=ier)
  call check_dealloc_err(ier, usr_msg)
end subroutine safe_dealloc_float_3d

!==============================================================================
!> Deallocate a 4D float array and check for errors
!! \param array The array to allocate
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_dealloc_float_4d(array, usr_msg)
  ! Arguments
  real(kind=4), dimension(:,:,:,:), allocatable, intent(inout) :: array
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  deallocate(array, stat=ier)
  call check_dealloc_err(ier, usr_msg)
end subroutine safe_dealloc_float_4d

!==============================================================================
!> Deallocate a 5D float array and check for errors
!! \param array The array to allocate
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_dealloc_float_5d(array, usr_msg)
  ! Arguments
  real(kind=4), dimension(:,:,:,:,:), allocatable, intent(inout) :: array
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  deallocate(array, stat=ier)
  call check_dealloc_err(ier, usr_msg)
end subroutine safe_dealloc_float_5d

!==============================================================================
!> Deallocate a 1D double array and check for errors
!! \param array The array to allocate
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_dealloc_double_1d(array, usr_msg)
  ! Arguments
  real(kind=8), dimension(:), allocatable, intent(inout) :: array
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  deallocate(array, stat=ier)
  call check_dealloc_err(ier, usr_msg)
end subroutine safe_dealloc_double_1d

!==============================================================================
!> Deallocate a 2D double array and check for errors
!! \param array The array to allocate
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_dealloc_double_2d(array, usr_msg)
  ! Arguments
  real(kind=8), dimension(:,:), allocatable, intent(inout) :: array
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  deallocate(array, stat=ier)
  call check_dealloc_err(ier, usr_msg)
end subroutine safe_dealloc_double_2d

!==============================================================================
!> Deallocate a 3D double array and check for errors
!! \param array The array to allocate
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_dealloc_double_3d(array, usr_msg)
  ! Arguments
  real(kind=8), dimension(:,:,:), allocatable, intent(inout) :: array
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  deallocate(array, stat=ier)
  call check_dealloc_err(ier, usr_msg)
end subroutine safe_dealloc_double_3d

!==============================================================================
!> Deallocate a 4D double array and check for errors
!! \param array The array to allocate
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_dealloc_double_4d(array, usr_msg)
  ! Arguments
  real(kind=8), dimension(:,:,:,:), allocatable, intent(inout) :: array
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  deallocate(array, stat=ier)
  call check_dealloc_err(ier, usr_msg)
end subroutine safe_dealloc_double_4d

!==============================================================================
!> Deallocate a 5D double array and check for errors
!! \param array The array to allocate
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_dealloc_double_5d(array, usr_msg)
  ! Arguments
  real(kind=8), dimension(:,:,:,:,:), allocatable, intent(inout) :: array
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  deallocate(array, stat=ier)
  call check_dealloc_err(ier, usr_msg)
end subroutine safe_dealloc_double_5d

!==============================================================================
!> Deallocate a 1D int array and check for errors
!! \param array The array to allocate
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_dealloc_int_1d(array, usr_msg)
  ! Arguments
  integer(kind=4), dimension(:), allocatable, intent(inout) :: array
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  deallocate(array, stat=ier)
  call check_dealloc_err(ier, usr_msg)
end subroutine safe_dealloc_int_1d

!==============================================================================
!> Deallocate a 2D int array and check for errors
!! \param array The array to allocate
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_dealloc_int_2d(array, usr_msg)
  ! Arguments
  integer(kind=4), dimension(:,:), allocatable, intent(inout) :: array
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  deallocate(array, stat=ier)
  call check_dealloc_err(ier, usr_msg)
end subroutine safe_dealloc_int_2d

!==============================================================================
!> Deallocate a 3D int array and check for errors
!! \param array The array to allocate
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_dealloc_int_3d(array, usr_msg)
  ! Arguments
  integer(kind=4), dimension(:,:,:), allocatable, intent(inout) :: array
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  deallocate(array, stat=ier)
  call check_dealloc_err(ier, usr_msg)
end subroutine safe_dealloc_int_3d

!==============================================================================
!> Deallocate a 4D int array and check for errors
!! \param array The array to allocate
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_dealloc_int_4d(array, usr_msg)
  ! Arguments
  integer(kind=4), dimension(:,:,:,:), allocatable, intent(inout) :: array
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  deallocate(array, stat=ier)
  call check_dealloc_err(ier, usr_msg)
end subroutine safe_dealloc_int_4d

!==============================================================================
!> Deallocate a 5D int array and check for errors
!! \param array The array to allocate
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_dealloc_int_5d(array, usr_msg)
  ! Arguments
  integer(kind=4), dimension(:,:,:,:,:), allocatable, intent(inout) :: array
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  deallocate(array, stat=ier)
  call check_dealloc_err(ier, usr_msg)
end subroutine safe_dealloc_int_5d

!==============================================================================
!> Deallocate a 1D long array and check for errors
!! \param array The array to allocate
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_dealloc_long_1d(array, usr_msg)
  ! Arguments
  integer(kind=8), dimension(:), allocatable, intent(inout) :: array
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  deallocate(array, stat=ier)
  call check_dealloc_err(ier, usr_msg)
end subroutine safe_dealloc_long_1d

!==============================================================================
!> Deallocate a 2D long array and check for errors
!! \param array The array to allocate
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_dealloc_long_2d(array, usr_msg)
  ! Arguments
  integer(kind=8), dimension(:,:), allocatable, intent(inout) :: array
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  deallocate(array, stat=ier)
  call check_dealloc_err(ier, usr_msg)
end subroutine safe_dealloc_long_2d

!==============================================================================
!> Deallocate a 3D long array and check for errors
!! \param array The array to allocate
!!dim1 The 1st dimension of the array
!! \param dim2 The 2nd dimension of the array
!! \param dim3 The 3rd dimension of the array
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_dealloc_long_3d(array, usr_msg)
  ! Arguments
  integer(kind=8), dimension(:,:,:), allocatable, intent(inout) :: array
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  deallocate(array, stat=ier)
  call check_dealloc_err(ier, usr_msg)
end subroutine safe_dealloc_long_3d

!==============================================================================
!> Deallocate a 4D long array and check for errors
!! \param array The array to allocate
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_dealloc_long_4d(array, usr_msg)
  ! Arguments
  integer(kind=8), dimension(:,:,:,:), allocatable, intent(inout) :: array
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  deallocate(array, stat=ier)
  call check_dealloc_err(ier, usr_msg)
end subroutine safe_dealloc_long_4d

!==============================================================================
!> Deallocate a 5D long array and check for errors
!! \param array The array to allocate
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_dealloc_long_5d(array, usr_msg)
  ! Arguments
  integer(kind=8), dimension(:,:,:,:,:), allocatable, intent(inout) :: array
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  deallocate(array, stat=ier)
  call check_dealloc_err(ier, usr_msg)
end subroutine safe_dealloc_long_5d

!!==============================================================================
!!> Deallocate a 1D byte array and check for errors
!!! \param array The array to allocate
!!! \param usr_msg A custom error message to print if allocation fails.
!subroutine safe_dealloc_byte_1d(array, usr_msg)
  !! Arguments
  !byte, dimension(:), allocatable, intent(inout) :: array
  !character(len=*), intent(in), optional :: usr_msg
  !! Local parameters
  !integer :: ier
!
  !deallocate(array, stat=ier)
  !call check_dealloc_err(ier, usr_msg)
!end subroutine safe_dealloc_byte_1d
!
!!==============================================================================
!!> Deallocate a 2D byte array and check for errors
!!! \param array The array to allocate
!!! \param usr_msg A custom error message to print if allocation fails.
!subroutine safe_dealloc_byte_2d(array, usr_msg)
  !! Arguments
  !byte, dimension(:,:), allocatable, intent(inout) :: array
  !character(len=*), intent(in), optional :: usr_msg
  !! Local parameters
  !integer :: ier
!
  !deallocate(array, stat=ier)
  !call check_dealloc_err(ier, usr_msg)
!end subroutine safe_dealloc_byte_2d
!
!!==============================================================================
!!> Deallocate a 3D byte array and check for errors
!!! \param array The array to allocate
!!! \param usr_msg A custom error message to print if allocation fails.
!subroutine safe_dealloc_byte_3d(array, usr_msg)
  !! Arguments
  !byte, dimension(:,:,:), allocatable, intent(inout) :: array
  !character(len=*), intent(in), optional :: usr_msg
  !! Local parameters
  !integer :: ier
!
  !deallocate(array, stat=ier)
  !call check_dealloc_err(ier, usr_msg)
!end subroutine safe_dealloc_byte_3d
!
!!==============================================================================
!!> Deallocate a 4D byte array and check for errors
!!! \param array The array to allocate
!!! \param usr_msg A custom error message to print if allocation fails.
!subroutine safe_dealloc_byte_4d(array, usr_msg)
  !! Arguments
  !byte, dimension(:,:,:,:), allocatable, intent(inout) :: array
  !character(len=*), intent(in), optional :: usr_msg
  !! Local parameters
  !integer :: ier
!
  !deallocate(array, stat=ier)
  !call check_dealloc_err(ier, usr_msg)
!end subroutine safe_dealloc_byte_4d
!
!!==============================================================================
!!> Deallocate a 5D byte array and check for errors
!!! \param array The array to allocate
!!! \param usr_msg A custom error message to print if allocation fails.
!subroutine safe_dealloc_byte_5d(array, usr_msg)
  !! Arguments
  !byte, dimension(:,:,:,:,:), allocatable, intent(inout) :: array
  !character(len=*), intent(in), optional :: usr_msg
  !! Local parameters
  !integer :: ier
!
  !deallocate(array, stat=ier)
  !call check_dealloc_err(ier, usr_msg)
!end subroutine safe_dealloc_byte_5d

!==============================================================================
!> Deallocate a 1D logical array and check for errors
!! \param array The array to allocate
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_dealloc_logical_1d(array, usr_msg)
  ! Arguments
  logical, dimension(:), allocatable, intent(inout) :: array
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  deallocate(array, stat=ier)
  call check_dealloc_err(ier, usr_msg)
end subroutine safe_dealloc_logical_1d

!==============================================================================
!> Deallocate a 2D logical array and check for errors
!! \param array The array to allocate
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_dealloc_logical_2d(array, usr_msg)
  ! Arguments
  logical, dimension(:,:), allocatable, intent(inout) :: array
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  deallocate(array, stat=ier)
  call check_dealloc_err(ier, usr_msg)
end subroutine safe_dealloc_logical_2d

!==============================================================================
!> Deallocate a 3D logical array and check for errors
!! \param array The array to allocate
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_dealloc_logical_3d(array, usr_msg)
  ! Arguments
  logical, dimension(:,:,:), allocatable, intent(inout) :: array
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  deallocate(array, stat=ier)
  call check_dealloc_err(ier, usr_msg)
end subroutine safe_dealloc_logical_3d

!==============================================================================
!> Deallocate a 4D logical array and check for errors
!! \param array The array to allocate
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_dealloc_logical_4d(array, usr_msg)
  ! Arguments
  logical, dimension(:,:,:,:), allocatable, intent(inout) :: array
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  deallocate(array, stat=ier)
  call check_dealloc_err(ier, usr_msg)
end subroutine safe_dealloc_logical_4d

!==============================================================================
!> Deallocate a 1D logical array and check for errors
!! \param array The array to allocate
!! \param usr_msg A custom error message to print if allocation fails.
subroutine safe_dealloc_logical_5d(array, usr_msg)
  ! Arguments
  logical, dimension(:,:,:,:,:), allocatable, intent(inout) :: array
  character(len=*), intent(in), optional :: usr_msg
  ! Local parameters
  integer :: ier

  deallocate(array, stat=ier)
  call check_dealloc_err(ier, usr_msg)
end subroutine safe_dealloc_logical_5d

end module safe_alloc_mod
