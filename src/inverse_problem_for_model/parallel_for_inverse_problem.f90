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


  subroutine sum_all_all_cr_for_simulatenous_runs(sendbuf, recvbuf, countval)

  use my_mpi
  use constants, only: CUSTOM_REAL
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS

  implicit none

  include "precision.h"

  integer,                intent(in)    :: countval
  real(kind=CUSTOM_REAL), intent(in)    :: sendbuf
  real(kind=CUSTOM_REAL), intent(inout) :: recvbuf
  integer                               :: ier

  if (NUMBER_OF_SIMULTANEOUS_RUNS <= 1) then
    recvbuf = sendbuf
    return
  endif

  call MPI_ALLREDUCE(sendbuf, recvbuf, countval, CUSTOM_MPI_TYPE, MPI_SUM, my_local_mpi_comm_for_bcast, ier)

  end subroutine sum_all_all_cr_for_simulatenous_runs

!
!-------------------------------------------------------------------------------------------------------------
!

  subroutine max_all_all_cr_for_simulatenous_runs(sendbuf, recvbuf, countval)

  use my_mpi
  use constants, only: CUSTOM_REAL
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS

  implicit none

  include "precision.h"

  integer,                intent(in)    :: countval
  real(kind=CUSTOM_REAL), intent(in)    :: sendbuf
  real(kind=CUSTOM_REAL), intent(inout) :: recvbuf
  integer                               :: ier

  if (NUMBER_OF_SIMULTANEOUS_RUNS <= 1) then
    recvbuf = sendbuf
    return
  endif

  call MPI_ALLREDUCE(sendbuf, recvbuf, countval, CUSTOM_MPI_TYPE, MPI_MAX, my_local_mpi_comm_for_bcast, ier)

  end subroutine max_all_all_cr_for_simulatenous_runs

!
!-------------------------------------------------------------------------------------------------------------
!

  subroutine min_all_all_cr_for_simulatenous_runs(sendbuf, recvbuf, countval)

  use my_mpi
  use constants, only: CUSTOM_REAL
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS

  implicit none

  include "precision.h"

  integer,                intent(in)    :: countval
  real(kind=CUSTOM_REAL), intent(in)    :: sendbuf
  real(kind=CUSTOM_REAL), intent(inout) :: recvbuf
  integer                               :: ier

  if (NUMBER_OF_SIMULTANEOUS_RUNS <= 1) then
    recvbuf = sendbuf
    return
  endif

  call MPI_ALLREDUCE(sendbuf, recvbuf, countval, CUSTOM_MPI_TYPE, MPI_MIN, my_local_mpi_comm_for_bcast, ier)

  end subroutine min_all_all_cr_for_simulatenous_runs

!
!-------------------------------------------------------------------------------------------------------------
!

  subroutine synchronize_all_world()

  use my_mpi

  implicit none

  integer :: ier

  ! synchronizes MPI processes
  call MPI_BARRIER(MPI_COMM_WORLD,ier)
  if (ier /= 0 ) stop 'Error synchronize MPI processes'

  end subroutine synchronize_all_world

!
!-------------------------------------------------------------------------------------------------------------
!

  subroutine synchronize_for_bcast()

  use my_mpi

  implicit none

  integer :: ier

  ! synchronizes MPI processes
  call MPI_BARRIER(my_local_mpi_comm_for_bcast,ier)
  if (ier /= 0 ) stop 'Error synchronize MPI processes'

  end subroutine synchronize_for_bcast

!
!-------------------------------------------------------------------------------------------------------------
!

  subroutine dummy_bcast(dummy_integer)

   use my_mpi

   implicit none
   integer :: dummy_integer, ier

   call mpi_bcast(dummy_integer, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)
   if (ier /= 0 ) stop 'Error synchronize MPI processes'

  end subroutine dummy_bcast
!
!-------------------------------------------------------------------------------------------------------------
!

  subroutine sum_all_all_cr_array(sendbuf, recvbuf, countval)

  use my_mpi
  use constants, only: CUSTOM_REAL

  implicit none

  include "precision.h"

  integer,                intent(in)    :: countval
  real(kind=CUSTOM_REAL), intent(in)    :: sendbuf
  real(kind=CUSTOM_REAL), intent(inout) :: recvbuf
  integer                               :: ier

  call MPI_ALLREDUCE(sendbuf, recvbuf, countval, CUSTOM_MPI_TYPE, MPI_SUM, my_local_mpi_comm_world, ier)

  end subroutine sum_all_all_cr_array
