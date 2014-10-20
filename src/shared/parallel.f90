!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
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

!! DK DK July 2014, CNRS Marseille, France:
!! DK DK added the ability to run several calculations (several earthquakes)
!! DK DK in an embarrassingly-parallel fashion from within the same run;
!! DK DK this can be useful when using a very large supercomputer to compute
!! DK DK many earthquakes in a catalog, in which case it can be better from
!! DK DK a batch job submission point of view to start fewer and much larger jobs,
!! DK DK each of them computing several earthquakes in parallel.
!! DK DK To turn that option on, set parameter NUMBER_OF_SIMULTANEOUS_RUNS
!! DK DK to a value greater than 1 in file setup/constants.h.in before
!! DK DK configuring and compiling the code.
!! DK DK To implement that, we create NUMBER_OF_SIMULTANEOUS_RUNS MPI sub-communicators,
!! DK DK each of them being labeled "my_local_mpi_comm_world", and we use them
!! DK DK in all the routines in "src/shared/parallel.f90", except in MPI_ABORT() because in that case
!! DK DK we need to kill the entire run.
!! DK DK When that option is on, of course the number of processor cores used to start
!! DK DK the code in the batch system must be a multiple of NUMBER_OF_SIMULTANEOUS_RUNS,
!! DK DK all the individual runs must use the same number of processor cores,
!! DK DK which as usual is NPROC in the input file DATA/Par_file,
!! DK DK and thus the total number of processor cores to request from the batch system
!! DK DK should be NUMBER_OF_SIMULTANEOUS_RUNS * NPROC.
!! DK DK All the runs to perform must be placed in directories called run0001, run0002, run0003 and so on
!! DK DK (with exactly four digits).

module my_mpi

! main parameter module for specfem simulations

  use mpi

  implicit none

  integer :: my_local_mpi_comm_world, my_local_mpi_comm_for_bcast

end module my_mpi

!----
!---- Parallel routines.  All MPI calls belong in this file!
!----

  subroutine stop_all()

  use mpi

  implicit none

  integer ier

! stop all the MPI processes, and exit
  call MPI_ABORT(MPI_COMM_WORLD,30,ier)
  stop 'error, program ended in exit_MPI'

  end subroutine stop_all

!
!----
!

  double precision function wtime()

  use my_mpi

  implicit none

  wtime = MPI_WTIME()

  end function wtime

!
!----
!

  subroutine synchronize_all()

  use my_mpi

  implicit none

  integer ier

  call MPI_BARRIER(my_local_mpi_comm_world,ier)

  end subroutine synchronize_all

!
!---- broadcast using the default communicator for the whole run
!

  subroutine bcast_all_i(buffer, countval)

  use my_mpi

  implicit none

  integer countval
  integer, dimension(countval) :: buffer

  integer ier

  call MPI_BCAST(buffer,countval,MPI_INTEGER,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_i

!
!----
!

  subroutine bcast_all_cr(buffer, countval)

  use my_mpi
  use constants,only: CUSTOM_REAL

  implicit none

  include "precision.h"

  integer countval
  real(kind=CUSTOM_REAL), dimension(countval) :: buffer

  integer ier

  call MPI_BCAST(buffer,countval,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_cr


!
!----
!

  subroutine bcast_all_singlecr(buffer)

  use my_mpi
  use constants,only: CUSTOM_REAL

  implicit none

  include "precision.h"

  real(kind=CUSTOM_REAL) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,1,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)

  end subroutine bcast_all_singlecr

!
!----
!

  subroutine bcast_all_dp(buffer, countval)

  use my_mpi

  implicit none

  integer countval
  double precision, dimension(countval) :: buffer

  integer ier

  call MPI_BCAST(buffer,countval,MPI_DOUBLE_PRECISION,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_dp


!
!----
!

  subroutine bcast_all_singledp(buffer)

  use my_mpi

  implicit none

  double precision :: buffer

  integer :: ier

  call MPI_BCAST(buffer,1,MPI_DOUBLE_PRECISION,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_singledp

!
!----
!

  subroutine bcast_all_r(buffer, countval)

  use my_mpi

  implicit none

  integer countval
  real, dimension(countval) :: buffer

  integer ier

  call MPI_BCAST(buffer,countval,MPI_REAL,0,my_local_mpi_comm_world,ier)

  end subroutine bcast_all_r

!
!---- broadcast using the communicator to send the mesh and model to other simultaneous runs
!

  subroutine bcast_all_i_for_database(buffer, countval)

  use my_mpi
  use constants,only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

  implicit none

  integer countval
  ! by not specifying any dimensions for the buffer here we can use this routine for arrays of any number
  ! of indices, provided we call the routine using the first memory cell of that multidimensional array,
  ! i.e. for instance buffer(1,1,1) if the array has three dimensions with indices that all start at 1.
  integer :: buffer

  integer ier

  if(.not. (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. BROADCAST_SAME_MESH_AND_MODEL)) return

  call MPI_BCAST(buffer,countval,MPI_INTEGER,0,my_local_mpi_comm_for_bcast,ier)

  end subroutine bcast_all_i_for_database

!
!----
!

  subroutine bcast_all_l_for_database(buffer, countval)

  use my_mpi
  use constants,only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

  implicit none

  integer countval
  ! by not specifying any dimensions for the buffer here we can use this routine for arrays of any number
  ! of indices, provided we call the routine using the first memory cell of that multidimensional array,
  ! i.e. for instance buffer(1,1,1) if the array has three dimensions with indices that all start at 1.
  logical :: buffer

  integer ier

  if(.not. (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. BROADCAST_SAME_MESH_AND_MODEL)) return

  call MPI_BCAST(buffer,countval,MPI_INTEGER,0,my_local_mpi_comm_for_bcast,ier)

  end subroutine bcast_all_l_for_database

!
!----
!

  subroutine bcast_all_cr_for_database(buffer, countval)

  use my_mpi
  use constants,only: CUSTOM_REAL,NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

  implicit none

  include "precision.h"

  integer countval
  ! by not specifying any dimensions for the buffer here we can use this routine for arrays of any number
  ! of indices, provided we call the routine using the first memory cell of that multidimensional array,
  ! i.e. for instance buffer(1,1,1) if the array has three dimensions with indices that all start at 1.
  real(kind=CUSTOM_REAL) :: buffer

  integer ier

  if(.not. (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. BROADCAST_SAME_MESH_AND_MODEL)) return

  call MPI_BCAST(buffer,countval,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_for_bcast,ier)

  end subroutine bcast_all_cr_for_database

!
!----
!

  subroutine bcast_all_dp_for_database(buffer, countval)

  use my_mpi
  use constants,only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

  implicit none

  integer countval
  ! by not specifying any dimensions for the buffer here we can use this routine for arrays of any number
  ! of indices, provided we call the routine using the first memory cell of that multidimensional array,
  ! i.e. for instance buffer(1,1,1) if the array has three dimensions with indices that all start at 1.
  double precision :: buffer

  integer ier

  if(.not. (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. BROADCAST_SAME_MESH_AND_MODEL)) return

  call MPI_BCAST(buffer,countval,MPI_DOUBLE_PRECISION,0,my_local_mpi_comm_for_bcast,ier)

  end subroutine bcast_all_dp_for_database

!
!----
!

  subroutine bcast_all_r_for_database(buffer, countval)

  use my_mpi
  use constants,only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

  implicit none

  integer countval
  ! by not specifying any dimensions for the buffer here we can use this routine for arrays of any number
  ! of indices, provided we call the routine using the first memory cell of that multidimensional array,
  ! i.e. for instance buffer(1,1,1) if the array has three dimensions with indices that all start at 1.
  real :: buffer

  integer ier

  if(.not. (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. BROADCAST_SAME_MESH_AND_MODEL)) return

  call MPI_BCAST(buffer,countval,MPI_REAL,0,my_local_mpi_comm_for_bcast,ier)

  end subroutine bcast_all_r_for_database

!
!----
!

  subroutine gather_all_i(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  use my_mpi

  implicit none

  integer sendcnt, recvcount, NPROC
  integer, dimension(sendcnt) :: sendbuf
  integer, dimension(recvcount,0:NPROC-1) :: recvbuf

  integer ier

  call MPI_GATHER(sendbuf,sendcnt,MPI_INTEGER, &
                  recvbuf,recvcount,MPI_INTEGER, &
                  0,my_local_mpi_comm_world,ier)

  end subroutine gather_all_i

!
!----
!

  subroutine gather_all_singlei(sendbuf, recvbuf, NPROC)

  use my_mpi

  implicit none

  integer NPROC
  integer :: sendbuf
  integer, dimension(0:NPROC-1) :: recvbuf

  integer ier

  call MPI_GATHER(sendbuf,1,MPI_INTEGER, &
                  recvbuf,1,MPI_INTEGER, &
                  0,my_local_mpi_comm_world,ier)

  end subroutine gather_all_singlei

!
!----
!

  subroutine gather_all_dp(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  use my_mpi

  implicit none

  integer sendcnt, recvcount, NPROC
  double precision, dimension(sendcnt) :: sendbuf
  double precision, dimension(recvcount,0:NPROC-1) :: recvbuf

  integer ier

  call MPI_GATHER(sendbuf,sendcnt,MPI_DOUBLE_PRECISION, &
                  recvbuf,recvcount,MPI_DOUBLE_PRECISION, &
                  0,my_local_mpi_comm_world,ier)

  end subroutine gather_all_dp

!
!----
!

  subroutine gather_all_cr(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  use my_mpi
  use constants,only: CUSTOM_REAL

  implicit none

  include "precision.h"

  integer sendcnt, recvcount, NPROC
  real(kind=CUSTOM_REAL), dimension(sendcnt) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(recvcount,0:NPROC-1) :: recvbuf

  integer ier

  call MPI_GATHER(sendbuf,sendcnt,CUSTOM_MPI_TYPE, &
                  recvbuf,recvcount,CUSTOM_MPI_TYPE, &
                  0,my_local_mpi_comm_world,ier)

  end subroutine gather_all_cr

!
!----
!

  subroutine gather_all_all_cr(sendbuf, recvbuf, counts, NPROC)

  use my_mpi
  use constants,only: CUSTOM_REAL

  implicit none

  include "precision.h"

  integer NPROC,counts
  real(kind=CUSTOM_REAL), dimension(counts) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(counts,0:NPROC-1) :: recvbuf

  integer ier

  call MPI_ALLGATHER(sendbuf,counts,CUSTOM_MPI_TYPE,recvbuf,counts,CUSTOM_MPI_TYPE,my_local_mpi_comm_world,ier)

  end subroutine gather_all_all_cr

!
!----
!

  subroutine gatherv_all_cr(sendbuf, sendcnt, recvbuf, recvcount, recvoffset,recvcounttot, NPROC)

  use my_mpi
  use constants,only: CUSTOM_REAL

  implicit none

  include "precision.h"

  integer sendcnt,recvcounttot,NPROC
  integer, dimension(NPROC) :: recvcount,recvoffset
  real(kind=CUSTOM_REAL), dimension(sendcnt) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(recvcounttot) :: recvbuf

  integer ier

  call MPI_GATHERV(sendbuf,sendcnt,CUSTOM_MPI_TYPE, &
                  recvbuf,recvcount,recvoffset,CUSTOM_MPI_TYPE, &
                  0,my_local_mpi_comm_world,ier)

  end subroutine gatherv_all_cr

!
!----
!

  subroutine init_mpi()

  use my_mpi

  implicit none

  integer ier

! initialize the MPI communicator and start the NPROCTOT MPI processes.
  call MPI_INIT(ier)

! create sub-communicators if needed, if running more than one earthquake from the same job
  call world_split()

  end subroutine init_mpi

!
!----
!

  subroutine finalize_mpi()

  use my_mpi

  implicit none

  integer ier

! close sub-communicators if needed, if running more than one earthquake from the same job
  call world_unsplit()

! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)

  end subroutine finalize_mpi

!
!----
!

  subroutine world_size(sizeval)

  use my_mpi

  implicit none

  integer sizeval
  integer ier

  call MPI_COMM_SIZE(my_local_mpi_comm_world,sizeval,ier)

  end subroutine world_size

!
!----
!

  subroutine world_rank(rank)

  use my_mpi

  implicit none

  integer rank
  integer ier

  call MPI_COMM_RANK(my_local_mpi_comm_world,rank,ier)

  end subroutine world_rank

!
!----
!

  subroutine min_all_dp(sendbuf, recvbuf)

  use my_mpi

  implicit none

  double precision sendbuf, recvbuf
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,my_local_mpi_comm_world,ier)

  end subroutine min_all_dp

!
!----
!

  subroutine max_all_dp(sendbuf, recvbuf)

  use my_mpi

  implicit none

  double precision sendbuf, recvbuf
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,my_local_mpi_comm_world,ier)

  end subroutine max_all_dp

!
!----
!

  subroutine max_all_cr(sendbuf, recvbuf)

  use my_mpi
  use constants,only: CUSTOM_REAL

  implicit none

  include "precision.h"

  real(kind=CUSTOM_REAL) sendbuf, recvbuf
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE,MPI_MAX,0,my_local_mpi_comm_world,ier)

  end subroutine max_all_cr

!
!----
!

  subroutine min_all_cr(sendbuf, recvbuf)

  use my_mpi
  use constants,only: CUSTOM_REAL

  implicit none

  include "precision.h"

  real(kind=CUSTOM_REAL) sendbuf, recvbuf
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE,MPI_MIN,0,my_local_mpi_comm_world,ier)

  end subroutine min_all_cr


!
!----
!

  subroutine min_all_all_cr(sendbuf, recvbuf)

  use my_mpi
  use constants,only: CUSTOM_REAL

  implicit none

  include "precision.h"

  real(kind=CUSTOM_REAL):: sendbuf, recvbuf
  integer ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE,MPI_MIN,my_local_mpi_comm_world,ier)

  end subroutine min_all_all_cr

!
!----
!
!
!
!  subroutine min_all_all_dp(sendbuf, recvbuf)
!
!  use my_mpi
!
!  implicit none
!
!  double precision :: sendbuf, recvbuf
!  integer ier
!
!  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION, &
!                  MPI_MIN,my_local_mpi_comm_world,ier)
!
!  end subroutine min_all_all_dp
!
!
!----
!

  subroutine max_all_i(sendbuf, recvbuf)

  use my_mpi

  implicit none

  integer :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_INTEGER,MPI_MAX,0,my_local_mpi_comm_world,ier)

  end subroutine max_all_i


!
!----
!

  subroutine max_allreduce_i(buffer,countval)

  use my_mpi

  implicit none

  integer :: countval
  integer,dimension(countval),intent(inout) :: buffer

  ! local parameters
  integer :: ier
  integer,dimension(countval) :: send

  ! seems not to be supported on all kind of MPI implementations...
  !call MPI_ALLREDUCE(MPI_IN_PLACE, buffer, countval, MPI_INTEGER, MPI_MAX, my_local_mpi_comm_world, ier)

  send(:) = buffer(:)

  call MPI_ALLREDUCE(send, buffer, countval, MPI_INTEGER, MPI_MAX, my_local_mpi_comm_world, ier)
  if( ier /= 0 ) stop 'Allreduce to get max values failed.'

  end subroutine max_allreduce_i

!
!----
!

  subroutine max_all_all_cr(sendbuf, recvbuf)

  use my_mpi
  use constants,only: CUSTOM_REAL

  implicit none

  include "precision.h"

  real(kind=CUSTOM_REAL):: sendbuf, recvbuf
  integer ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE,MPI_MAX,my_local_mpi_comm_world,ier)

  end subroutine max_all_all_cr


!
!----
!


  subroutine max_all_all_dp(sendbuf, recvbuf)

  use my_mpi

  implicit none

  double precision :: sendbuf, recvbuf
  integer ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_MAX,my_local_mpi_comm_world,ier)

  end subroutine max_all_all_dp


!
!----
!

  subroutine min_all_i(sendbuf, recvbuf)

  use my_mpi

  implicit none

  integer:: sendbuf, recvbuf
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_INTEGER,MPI_MIN,0,my_local_mpi_comm_world,ier)

  end subroutine min_all_i

!
!----
!

  subroutine maxloc_all_dp(sendbuf, recvbuf)

  use my_mpi

  implicit none

  double precision, dimension(2) :: sendbuf,recvbuf
  integer ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_2DOUBLE_PRECISION,MPI_MAXLOC,my_local_mpi_comm_world,ier)

  end subroutine maxloc_all_dp


!
!----
!


  subroutine sum_all_dp(sendbuf, recvbuf)

  use my_mpi

  implicit none

  double precision sendbuf, recvbuf
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,my_local_mpi_comm_world,ier)

  end subroutine sum_all_dp

!
!----
!

  subroutine sum_all_cr(sendbuf, recvbuf)

  use my_mpi
  use constants,only: CUSTOM_REAL

  implicit none

  include "precision.h"

  real(kind=CUSTOM_REAL) sendbuf, recvbuf
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE,MPI_SUM,0,my_local_mpi_comm_world,ier)

  end subroutine sum_all_cr

!
!----
!

  subroutine sum_all_1Darray_dp(sendbuf, recvbuf, nx)

  use my_mpi

  implicit none

  integer :: nx
  double precision, dimension(nx) :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,nx,MPI_DOUBLE_PRECISION,MPI_SUM,0,my_local_mpi_comm_world,ier)

  end subroutine sum_all_1Darray_dp

!
!----
!

  subroutine sum_all_all_cr(sendbuf, recvbuf)

  use my_mpi
  use constants,only: CUSTOM_REAL

  implicit none

  include "precision.h"

  real(kind=CUSTOM_REAL) sendbuf, recvbuf
  integer ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE,MPI_SUM,my_local_mpi_comm_world,ier)

  end subroutine sum_all_all_cr

!
!----
!

  subroutine sum_all_i(sendbuf, recvbuf)

  use my_mpi

  implicit none

  integer sendbuf, recvbuf
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_INTEGER,MPI_SUM,0,my_local_mpi_comm_world,ier)

  end subroutine sum_all_i

!
!----
!

  subroutine sum_all_all_i(sendbuf, recvbuf)

  use my_mpi

  implicit none

  integer sendbuf, recvbuf
  integer ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_INTEGER,MPI_SUM,my_local_mpi_comm_world,ier)

  end subroutine sum_all_all_i

!
!----
!

  subroutine any_all_l(sendbuf, recvbuf)

  use my_mpi

  implicit none

  logical sendbuf, recvbuf
  integer ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_LOGICAL,MPI_LOR,my_local_mpi_comm_world,ier)

  end subroutine any_all_l

!
!----
!

  subroutine sendrecv_all_cr(sendbuf, sendcount, dest, sendtag, &
                             recvbuf, recvcount, source, recvtag)

  use my_mpi
  use constants,only: CUSTOM_REAL

  implicit none

  include "precision.h"

  integer sendcount, recvcount, dest, sendtag, source, recvtag
  real(kind=CUSTOM_REAL), dimension(sendcount) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(recvcount) :: recvbuf

  integer ier

  call MPI_SENDRECV(sendbuf,sendcount,CUSTOM_MPI_TYPE,dest,sendtag, &
                    recvbuf,recvcount,CUSTOM_MPI_TYPE,source,recvtag, &
                    my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)

  end subroutine sendrecv_all_cr

!
!----
!

  integer function proc_null()

  use my_mpi

  implicit none

  proc_null = MPI_PROC_NULL

  end function proc_null

!
!----
!

  subroutine isend_cr(sendbuf, sendcount, dest, sendtag, req)

  use my_mpi
  use constants,only: CUSTOM_REAL

  implicit none

  include "precision.h"

  integer sendcount, dest, sendtag, req
  real(kind=CUSTOM_REAL), dimension(sendcount) :: sendbuf

  integer ier

  call MPI_ISEND(sendbuf,sendcount,CUSTOM_MPI_TYPE,dest,sendtag,my_local_mpi_comm_world,req,ier)

  end subroutine isend_cr

!
!----
!

  subroutine irecv_cr(recvbuf, recvcount, dest, recvtag, req)

  use my_mpi
  use constants,only: CUSTOM_REAL

  implicit none

  include "precision.h"

  integer recvcount, dest, recvtag, req
  real(kind=CUSTOM_REAL), dimension(recvcount) :: recvbuf

  integer ier

  call MPI_IRECV(recvbuf,recvcount,CUSTOM_MPI_TYPE,dest,recvtag,my_local_mpi_comm_world,req,ier)

  end subroutine irecv_cr

!
!----
!

  subroutine isend_i(sendbuf, sendcount, dest, sendtag, req)

  use my_mpi

  implicit none

  integer sendcount, dest, sendtag, req
  integer, dimension(sendcount) :: sendbuf

  integer ier

  call MPI_ISEND(sendbuf,sendcount,MPI_INTEGER,dest,sendtag,my_local_mpi_comm_world,req,ier)

  end subroutine isend_i

!
!----
!

  subroutine irecv_i(recvbuf, recvcount, dest, recvtag, req)

  use my_mpi

  implicit none

  integer recvcount, dest, recvtag, req
  integer, dimension(recvcount) :: recvbuf
  integer ier

  call MPI_IRECV(recvbuf,recvcount,MPI_INTEGER,dest,recvtag,my_local_mpi_comm_world,req,ier)

  end subroutine irecv_i


!
!----
!

  subroutine recv_i(recvbuf, recvcount, dest, recvtag )

  use my_mpi

  implicit none

  integer dest,recvtag
  integer recvcount
  !integer recvbuf
  integer,dimension(recvcount):: recvbuf
  integer ier

  call MPI_RECV(recvbuf,recvcount,MPI_INTEGER,dest,recvtag, &
                my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)

  end subroutine recv_i

!
!----
!

  subroutine recvv_cr(recvbuf, recvcount, dest, recvtag )

  use my_mpi
  use constants,only: CUSTOM_REAL

  implicit none

  include "precision.h"

  integer recvcount,dest,recvtag
  real(kind=CUSTOM_REAL),dimension(recvcount) :: recvbuf
  integer ier

  call MPI_RECV(recvbuf,recvcount,CUSTOM_MPI_TYPE,dest,recvtag, &
                my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)

  end subroutine recvv_cr


!
!----
!

  subroutine send_i(sendbuf, sendcount, dest, sendtag)

  use my_mpi

  implicit none

  !integer sendbuf,sendcount,dest,sendtag
  integer dest,sendtag
  integer sendcount
  integer,dimension(sendcount):: sendbuf
  integer ier

  call MPI_SEND(sendbuf,sendcount,MPI_INTEGER,dest,sendtag,my_local_mpi_comm_world,ier)

  end subroutine send_i


!
!----
!

  subroutine send_i_t(sendbuf,sendcount,dest)

  use my_mpi

  implicit none

  integer :: dest,sendcount,ier
  integer :: tag = 100
  integer, dimension(sendcount) :: sendbuf

  call MPI_SEND(sendbuf,sendcount,MPI_INTEGER,dest,tag,my_local_mpi_comm_world,ier)

  end subroutine send_i_t

!
!----
!


  subroutine recv_i_t(recvbuf,recvcount,source)

  use my_mpi

  implicit none

  integer :: source,recvcount,ier
  integer :: tag = 100
  integer, dimension(recvcount) :: recvbuf

  call MPI_RECV(recvbuf,recvcount,MPI_INTEGER,source,tag, &
                my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)

  end subroutine recv_i_t


!
!----
!
!
!  subroutine send_dp_t(sendbuf,sendcount,dest)
!
!  use my_mpi
!
!  implicit none
!
!  integer :: dest,sendcount,ier
!  integer :: tag = 100
!  double precision, dimension(sendcount) :: sendbuf
!
!  call MPI_SEND(sendbuf,sendcount,MPI_DOUBLE_PRECISION,dest,tag, &
!       my_local_mpi_comm_world,ier)
!
!  end subroutine send_dp_t
!
!
!----
!
!
!  subroutine recv_dp_t(recvbuf,recvcount,source)
!
!  use my_mpi
!
!  implicit none
!
!  integer :: recvcount,source,ier
!  integer :: tag = 100
!  double precision, dimension(recvcount) :: recvbuf
!
!  call MPI_RECV(recvbuf,recvcount,MPI_DOUBLE_PRECISION,source,tag, &
!                my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)
!
!  end subroutine recv_dp_t
!
!
!
!----
!

  subroutine send_dp(sendbuf, sendcount, dest, sendtag)

  use my_mpi

  implicit none

  integer dest,sendtag
  integer sendcount
  double precision,dimension(sendcount):: sendbuf
  integer ier

  call MPI_SEND(sendbuf,sendcount,MPI_DOUBLE_PRECISION,dest,sendtag,my_local_mpi_comm_world,ier)

  end subroutine send_dp

!
!----
!

  subroutine recv_dp(recvbuf, recvcount, dest, recvtag)

  use my_mpi

  implicit none

  integer dest,recvtag
  integer recvcount
  double precision,dimension(recvcount):: recvbuf
  integer ier

  call MPI_RECV(recvbuf,recvcount,MPI_DOUBLE_PRECISION,dest,recvtag, &
                my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)

  end subroutine recv_dp

!
!----
!

  subroutine sendv_cr(sendbuf, sendcount, dest, sendtag)

  use my_mpi
  use constants,only: CUSTOM_REAL

  implicit none

  include "precision.h"

  integer sendcount,dest,sendtag
  real(kind=CUSTOM_REAL),dimension(sendcount) :: sendbuf
  integer ier

  call MPI_SEND(sendbuf,sendcount,CUSTOM_MPI_TYPE,dest,sendtag,my_local_mpi_comm_world,ier)

  end subroutine sendv_cr

!
!----
!

  subroutine wait_req(req)

  use my_mpi

  implicit none

  integer :: req

  integer :: ier

  call mpi_wait(req,MPI_STATUS_IGNORE,ier)

  end subroutine wait_req

!
!----
!

  subroutine world_get_comm(comm)

  use my_mpi

  implicit none

  integer,intent(out) :: comm

  comm = my_local_mpi_comm_world

  end subroutine world_get_comm

!
!----
!

  subroutine world_duplicate(comm)

  use my_mpi

  implicit none

  integer,intent(out) :: comm
  integer :: ier

  call MPI_COMM_DUP(my_local_mpi_comm_world,comm,ier)
  if( ier /= 0 ) stop 'error duplicating my_local_mpi_comm_world communicator'

  end subroutine world_duplicate

!
!----
!

! create sub-communicators if needed, if running more than one earthquake from the same job.
!! DK DK create a sub-communicator for each independent run;
!! DK DK if there is a single run to do, then just copy the default communicator to the new one
  subroutine world_split()

  use my_mpi
  use constants,only: MAX_STRING_LEN,NUMBER_OF_SIMULTANEOUS_RUNS,OUTPUT_FILES_PATH, &
    IMAIN,ISTANDARD_OUTPUT,mygroup,BROADCAST_SAME_MESH_AND_MODEL,I_should_read_the_database

  implicit none

  integer :: sizeval,myrank,ier,key,my_group_for_bcast,my_local_rank_for_bcast,NPROC

  character(len=MAX_STRING_LEN) :: path_to_add

  if(NUMBER_OF_SIMULTANEOUS_RUNS <= 0) stop 'NUMBER_OF_SIMULTANEOUS_RUNS <= 0 makes no sense'

  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeval,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)

  if(NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mod(sizeval,NUMBER_OF_SIMULTANEOUS_RUNS) /= 0) &
    stop 'the number of MPI processes is not a multiple of NUMBER_OF_SIMULTANEOUS_RUNS'

  if(NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. IMAIN == ISTANDARD_OUTPUT) &
    stop 'must not have IMAIN == ISTANDARD_OUTPUT when NUMBER_OF_SIMULTANEOUS_RUNS > 1 otherwise output to screen is mingled'

  if(NUMBER_OF_SIMULTANEOUS_RUNS == 1) then

    my_local_mpi_comm_world = MPI_COMM_WORLD

! no broadcast of the mesh and model databases to other runs in that case
    my_group_for_bcast = 0
    my_local_mpi_comm_for_bcast = MPI_COMM_NULL

  else

!--- create a subcommunicator for each independent run

    NPROC = sizeval / NUMBER_OF_SIMULTANEOUS_RUNS

!   create the different groups of processes, one for each independent run
    mygroup = myrank / NPROC
    key = myrank
    if(mygroup < 0 .or. mygroup > NUMBER_OF_SIMULTANEOUS_RUNS-1) stop 'invalid value of mygroup'

!   build the sub-communicators
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, mygroup, key, my_local_mpi_comm_world, ier)
    if(ier /= 0) stop 'error while trying to create the sub-communicators'

!   add the right directory for that run (group numbers start at zero, but directory names start at run0001, thus we add one)
    write(path_to_add,"('run',i4.4,'/')") mygroup + 1
    OUTPUT_FILES_PATH = path_to_add(1:len_trim(path_to_add))//OUTPUT_FILES_PATH(1:len_trim(OUTPUT_FILES_PATH))

!--- create a subcommunicator to broadcast the identical mesh and model databases if needed
    if(BROADCAST_SAME_MESH_AND_MODEL) then

      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)
!     to broadcast the model, split along similar ranks per run instead
      my_group_for_bcast = mod(myrank,NPROC)
      key = myrank
      if(my_group_for_bcast < 0 .or. my_group_for_bcast > NPROC-1) stop 'invalid value of my_group_for_bcast'

!     build the sub-communicators
      call MPI_COMM_SPLIT(MPI_COMM_WORLD, my_group_for_bcast, key, my_local_mpi_comm_for_bcast, ier)
      if(ier /= 0) stop 'error while trying to create the sub-communicators'

!     see if that process will need to read the mesh and model database and then broadcast it to others
      call MPI_COMM_RANK(my_local_mpi_comm_for_bcast,my_local_rank_for_bcast,ier)
      if(my_local_rank_for_bcast > 0) I_should_read_the_database = .false.

    else

! no broadcast of the mesh and model databases to other runs in that case
      my_group_for_bcast = 0
      my_local_mpi_comm_for_bcast = MPI_COMM_NULL

    endif

  endif

  end subroutine world_split

!
!----
!

! close sub-communicators if needed, if running more than one earthquake from the same job.
  subroutine world_unsplit()

  use my_mpi
  use constants,only: NUMBER_OF_SIMULTANEOUS_RUNS,BROADCAST_SAME_MESH_AND_MODEL

  implicit none

  integer :: ier

  if(NUMBER_OF_SIMULTANEOUS_RUNS > 1) then
    call MPI_COMM_FREE(my_local_mpi_comm_world,ier)
    if(BROADCAST_SAME_MESH_AND_MODEL) call MPI_COMM_FREE(my_local_mpi_comm_for_bcast,ier)
  endif

  end subroutine world_unsplit

