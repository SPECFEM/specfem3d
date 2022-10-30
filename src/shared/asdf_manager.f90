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

!==============================================================================
!> Initialize ASDF for reading the adjoint sources
subroutine asdf_setup(file_id, path_to_add, simul_run_flag)

  use iso_c_binding, only: C_NULL_CHAR
  implicit none

  ! asdf file handle
  !
  ! note: ASDF uses hdf5 file i/o routines. therefore, ASDF c routines define the file_id handle as hid_t.
  !       the datatype hid_t is defined by the hdf5 library, and for Fortran in file H5f90i_gen.h as:
  !          #define c_int_8 long long
  !          typedef c_int_8 hid_t_f
  !       in Fortran codes, one could use the hdf5 module for this
  !          use hdf5, only: HID_T
  !          integer(HID_T) :: file_id
  !       which will required the hdf5 library paths set for compilation and linking.
  !       instead here, the c_int_8 corresponds to long long, which in Fortran would be an 8-byte integer
  integer(kind=8),intent(inout) :: file_id

  character(len=*), intent(in) :: path_to_add
  logical, intent(in) :: simul_run_flag

  ! local parameters
  character(len=512) :: filename
  integer :: comm
  integer :: ier

  call world_duplicate(comm)
  call ASDF_initialize_hdf5_f(ier)

  if (simul_run_flag) then
    filename = trim(path_to_add) // 'SEM/adjoint.h5' // C_NULL_CHAR
  else
    filename = 'SEM/adjoint.h5' // C_NULL_CHAR
  endif

  call ASDF_open_read_only_f(filename, comm, file_id)

  ! checks file id
  if (file_id < 0) then
    ! negative file id indicates an opening error
    print *,'Error opening ASDF file   : ',trim(filename)
    print *,'       got invalid file id: ',file_id
    print *,'Please check if the file exists...'
    stop 'Error opening ASDF file'
  endif

end subroutine asdf_setup

!==============================================================================
!> Finalize ASDF. Called once all adjoint sources have been read from the file
subroutine asdf_cleanup()

  implicit none

  ! local parameters
  integer :: myrank
  integer :: ier

  call world_rank(myrank)
  call synchronize_all()

  call ASDF_finalize_hdf5_f(ier);

  if (ier /= 0 ) stop 'Error cleaning up ASDF: calling ASDF_finalize_hdf5_f() routine failed'

end subroutine asdf_cleanup
