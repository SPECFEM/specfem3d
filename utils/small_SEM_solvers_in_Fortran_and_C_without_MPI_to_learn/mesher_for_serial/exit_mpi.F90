!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2008
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

! version with rank number printed in the error message
  subroutine exit_MPI(myrank,error_msg)

  implicit none

! standard include of the MPI library
#ifdef USE_MPI
  include 'mpif.h'
#endif

  include "constants.h"

! identifier for error message file
  integer, parameter :: IERROR = 30

  integer myrank
  character(len=*) error_msg

#ifdef USE_MPI
  integer ier
#endif
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
  if(myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) close(IMAIN)

! stop all the MPI processes, and exit
#ifdef USE_MPI
  call MPI_ABORT(MPI_COMM_WORLD,30,ier)
#endif
  stop 'error, program ended in exit_MPI'

  end subroutine exit_MPI

!
!----
!

! version without rank number printed in the error message
  subroutine exit_MPI_without_rank(error_msg)

  implicit none

! standard include of the MPI library
#ifdef USE_MPI
  include 'mpif.h'
#endif

  include "constants.h"

  character(len=*) error_msg

#ifdef USE_MPI
  integer ier
#endif

! write error message to screen
  write(*,*) error_msg(1:len(error_msg))
  write(*,*) 'Error detected, aborting MPI...'

! stop all the MPI processes, and exit
#ifdef USE_MPI
  call MPI_ABORT(MPI_COMM_WORLD,30,ier)
#endif
  stop 'error, program ended in exit_MPI'

  end subroutine exit_MPI_without_rank

