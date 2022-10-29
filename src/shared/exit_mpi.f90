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

! end the simulation and exit MPI

  subroutine exit_MPI(myrank,error_msg)

  use constants, only: MAX_STRING_LEN,IMAIN,ISTANDARD_OUTPUT,OUTPUT_FILES

  implicit none

  ! identifier for error message file
  integer, parameter :: IERROR = 30

  integer, intent(in) :: myrank
  character(len=*),intent(in) :: error_msg

  character(len=MAX_STRING_LEN) :: outputname

  ! write error message to screen
  write(*,*) error_msg(1:len(error_msg))
  write(*,*) 'Error detected, aborting MPI... proc ',myrank

  ! write error message to file
  write(outputname,"('/error_message',i6.6,'.txt')") myrank
  open(unit=IERROR,file=trim(OUTPUT_FILES)//outputname,status='unknown')
  write(IERROR,*) error_msg(1:len(error_msg))
  write(IERROR,*) 'Error detected, aborting MPI... proc ',myrank
  close(IERROR)

  ! close output file
  if (myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) close(IMAIN)

  ! flushes possible left-overs from print-statements
  call flush_stdout()

  ! abort execution
  call abort_mpi()

  end subroutine exit_MPI

!
!-------------------------------------------------------------------------------------------------
!

! version without rank number printed in the error message

  subroutine exit_MPI_without_rank(error_msg)

  use constants

  implicit none

  character(len=*) :: error_msg

  ! write error message to screen
  write(*,*) error_msg(1:len(error_msg))
  write(*,*) 'Error detected, aborting MPI...'

  ! flushes possible left-overs from print-statements
  call flush_stdout()

  ! abort execution
  call abort_mpi()

  end subroutine exit_MPI_without_rank

!-------------------------------------------------------------------------------------------------
!
! I/O wrapper function
!
!-------------------------------------------------------------------------------------------------

  subroutine flush_IMAIN()

  use constants, only: IMAIN

  implicit none

  ! only main process writes out to main output file
  ! file I/O in Fortran is buffered by default
  !
  ! note: Fortran2003 includes a FLUSH statement
  !          which is implemented by most compilers by now
  !
  ! otherwise:
  !   a) comment out the line below
  !   b) try to use instead: call flush(IMAIN)

  flush(IMAIN)

  end subroutine flush_IMAIN

!
!-------------------------------------------------------------------------------------------------
!

  subroutine flush_stdout()

! flushes possible left-overs from print-statements

  implicit none

  logical :: is_connected

  ! note: Cray systems don't flush print statements before ending with an MPI abort,
  !       which often omits debugging statements with print before it.
  !
  !       to check which unit is used for standard output, one might also use a Fortran2003 module iso_Fortran_env:
  !         use, intrinsic :: iso_Fortran_env, only: output_unit

  ! checks default stdout unit 6
  inquire(unit=6,opened=is_connected)
  if (is_connected) &
    flush(6)

  ! checks Cray stdout unit 101
  inquire(unit=101,opened=is_connected)
  if (is_connected) &
    flush(101)

  end subroutine flush_stdout
