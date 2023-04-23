!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
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

! version with rank number printed in the error message
  subroutine exit_MPI(myrank,error_msg)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "../../SHARE_FILES/HEADER_FILES/constants.h"

! identifier for error message file
  integer, parameter :: IERROR = 30

  integer myrank
  character(len=*) error_msg

  integer ier
  character(len=80) outputname
  character(len=150) OUTPUT_FILES

! write error message to screen
  write(*,*) error_msg(1:len(error_msg))
  write(*,*) 'Error detected, aborting MPI... proc ',myrank

! write error message to file
  call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', 'OUTPUT_FILES')
  write(outputname,"('/error_message',i6.6,'.txt')") myrank
  open(unit=IERROR,file=trim(OUTPUT_FILES)//outputname,status='unknown')
  write(IERROR,*) error_msg(1:len(error_msg))
  write(IERROR,*) 'Error detected, aborting MPI... proc ',myrank
  close(IERROR)

! close output file
  if (myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) close(IMAIN)

! stop all the MPI processes, and exit
! note: MPI_ABORT does not return, and does exit the
!          program with an error code of 30
  call MPI_ABORT(MPI_COMM_WORLD,30,ier)

! otherwise: there is no standard behaviour to exit with an error code in Fortran,
! however most compilers do recognize this as an error code stop statement;
! to check stop code in terminal: > echo $?
  stop 30

  ! or just exit with message:
  !stop 'error, program ended in exit_MPI'

  end subroutine exit_MPI

!
!----
!

! version without rank number printed in the error message
  subroutine exit_MPI_without_rank(error_msg)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "../../SHARE_FILES/HEADER_FILES/constants.h"

  character(len=*) error_msg

  integer ier

! write error message to screen
  write(*,*) error_msg(1:len(error_msg))
  write(*,*) 'Error detected, aborting MPI...'

! stop all the MPI processes, and exit
  call MPI_ABORT(MPI_COMM_WORLD,30,ier)
  stop 'error, program ended in exit_MPI'

  end subroutine exit_MPI_without_rank

