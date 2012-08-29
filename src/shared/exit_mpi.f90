!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
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

  implicit none

  include "constants.h"

! identifier for error message file
  integer, parameter :: IERROR = 30

  integer myrank
  character(len=*) error_msg

  character(len=80) outputname
  character(len=256) OUTPUT_FILES

! write error message to screen
  write(*,*) error_msg(1:len(error_msg))
  write(*,*) 'Error detected, aborting MPI... proc ',myrank

! write error message to file
  call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', OUTPUT_FILES_PATH(1:len_trim(OUTPUT_FILES_PATH)))
  write(outputname,"('/error_message',i6.6,'.txt')") myrank
  open(unit=IERROR,file=trim(OUTPUT_FILES)//outputname,status='unknown')
  write(IERROR,*) error_msg(1:len(error_msg))
  write(IERROR,*) 'Error detected, aborting MPI... proc ',myrank
  close(IERROR)

! close output file
  if(myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) close(IMAIN)

  call stop_all()

  end subroutine exit_MPI

!
!-------------------------------------------------------------------------------------------------
!

! version without rank number printed in the error message

  subroutine exit_MPI_without_rank(error_msg)

  implicit none

  include "constants.h"

  character(len=*) error_msg

! write error message to screen
  write(*,*) error_msg(1:len(error_msg))
  write(*,*) 'Error detected, aborting MPI...'

  call stop_all()

  end subroutine exit_MPI_without_rank


!-------------------------------------------------------------------------------------------------
!
! I/O wrapper function
!
!-------------------------------------------------------------------------------------------------

  subroutine flush_IMAIN()

  implicit none

  include "constants.h"

  ! only master process writes out to main output file
  ! file I/O in fortran is buffered by default
  !
  ! note: Fortran2003 includes a FLUSH statement
  !          which is implemented by most compilers by now
  !
  ! otherwise:
  !   a) comment out the line below
  !   b) try to use instead: call flush(IMAIN)

  flush(IMAIN)

  end subroutine flush_IMAIN

