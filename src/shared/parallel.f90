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

  use mpi

  implicit none

  wtime = MPI_WTIME()

  end function wtime

!
!----
!

  subroutine synchronize_all()

  use mpi

  implicit none

  integer ier

  call MPI_BARRIER(MPI_COMM_WORLD,ier)

  end subroutine synchronize_all

!
!----
!

  subroutine bcast_all_i(buffer, countval)

  use mpi

  implicit none

  integer countval
  integer, dimension(countval) :: buffer

  integer ier

  call MPI_BCAST(buffer,countval,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  end subroutine bcast_all_i

!
!----
!

  subroutine bcast_all_cr(buffer, countval)

  use constants
  use mpi

  implicit none

  include "precision.h"

  integer countval
  real(kind=CUSTOM_REAL), dimension(countval) :: buffer

  integer ier

  call MPI_BCAST(buffer,countval,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)

  end subroutine bcast_all_cr

!
!----
!

  subroutine bcast_all_dp(buffer, countval)

  use mpi

  implicit none

  integer countval
  double precision, dimension(countval) :: buffer

  integer ier

  call MPI_BCAST(buffer,countval,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

  end subroutine bcast_all_dp

!
!----
!

  subroutine bcast_all_r(buffer, countval)

  use constants
  use mpi

  implicit none

  include "precision.h"

  integer countval
  real, dimension(countval) :: buffer

  integer ier

  call MPI_BCAST(buffer,countval,MPI_REAL,0,MPI_COMM_WORLD,ier)

  end subroutine bcast_all_r


!
!----
!

  subroutine gather_all_i(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  use mpi

  implicit none

  integer sendcnt, recvcount, NPROC
  integer, dimension(sendcnt) :: sendbuf
  integer, dimension(recvcount,0:NPROC-1) :: recvbuf

  integer ier

  call MPI_GATHER(sendbuf,sendcnt,MPI_INTEGER, &
                  recvbuf,recvcount,MPI_INTEGER, &
                  0,MPI_COMM_WORLD,ier)

  end subroutine gather_all_i


!
!----
!

  subroutine gather_all_singlei(sendbuf, recvbuf, NPROC)

  use mpi

  implicit none

  integer NPROC
  integer :: sendbuf
  integer, dimension(0:NPROC-1) :: recvbuf

  integer ier

  call MPI_GATHER(sendbuf,1,MPI_INTEGER, &
                  recvbuf,1,MPI_INTEGER, &
                  0,MPI_COMM_WORLD,ier)

  end subroutine gather_all_singlei


!
!----
!

  subroutine gather_all_dp(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  use mpi

  implicit none

  integer sendcnt, recvcount, NPROC
  double precision, dimension(sendcnt) :: sendbuf
  double precision, dimension(recvcount,0:NPROC-1) :: recvbuf

  integer ier

  call MPI_GATHER(sendbuf,sendcnt,MPI_DOUBLE_PRECISION, &
                  recvbuf,recvcount,MPI_DOUBLE_PRECISION, &
                  0,MPI_COMM_WORLD,ier)

  end subroutine gather_all_dp

!
!----
!

  subroutine gather_all_cr(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  use constants
  use mpi

  implicit none

  include "precision.h"

  integer sendcnt, recvcount, NPROC
  real(kind=CUSTOM_REAL), dimension(sendcnt) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(recvcount,0:NPROC-1) :: recvbuf

  integer ier

  call MPI_GATHER(sendbuf,sendcnt,CUSTOM_MPI_TYPE, &
                  recvbuf,recvcount,CUSTOM_MPI_TYPE, &
                  0,MPI_COMM_WORLD,ier)

  end subroutine gather_all_cr

!
!----
!

  subroutine gather_all_all_cr(sendbuf, recvbuf, counts, NPROC)

  use constants
  use mpi

  implicit none

  include "precision.h"

  integer NPROC,counts
  real(kind=CUSTOM_REAL), dimension(counts) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(counts,0:NPROC-1) :: recvbuf

  integer ier

  call MPI_ALLGATHER(sendbuf,counts,CUSTOM_MPI_TYPE,recvbuf,counts,CUSTOM_MPI_TYPE, &
                 MPI_COMM_WORLD,ier)

  end subroutine gather_all_all_cr

!
!----
!

  subroutine gatherv_all_cr(sendbuf, sendcnt, recvbuf, recvcount, recvoffset,recvcounttot, NPROC)

  use constants
  use mpi

  implicit none

  include "precision.h"

  integer sendcnt,recvcounttot,NPROC
  integer, dimension(NPROC) :: recvcount,recvoffset
  real(kind=CUSTOM_REAL), dimension(sendcnt) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(recvcounttot) :: recvbuf

  integer ier

  call MPI_GATHERV(sendbuf,sendcnt,CUSTOM_MPI_TYPE, &
                  recvbuf,recvcount,recvoffset,CUSTOM_MPI_TYPE, &
                  0,MPI_COMM_WORLD,ier)

  end subroutine gatherv_all_cr

!
!----
!

  subroutine init()

  use mpi

  implicit none

  integer ier

! initialize the MPI communicator and start the NPROCTOT MPI processes.
  call MPI_INIT(ier)

  end subroutine init

!
!----
!

  subroutine finalize()

  use mpi

  implicit none

  integer ier

! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)

  end subroutine finalize

!
!----
!

  subroutine world_size(sizeval)

  use mpi

  implicit none

  integer sizeval
  integer ier

  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeval,ier)

  end subroutine world_size

!
!----
!

  subroutine world_rank(rank)

  use mpi

  implicit none

  integer rank
  integer ier

  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ier)

  end subroutine world_rank

!
!----
!

  subroutine min_all_dp(sendbuf, recvbuf)

  use mpi

  implicit none

  double precision sendbuf, recvbuf
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION, &
                  MPI_MIN,0,MPI_COMM_WORLD,ier)

  end subroutine min_all_dp

!
!----
!

  subroutine max_all_dp(sendbuf, recvbuf)

  use mpi

  implicit none

  double precision sendbuf, recvbuf
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION, &
                  MPI_MAX,0,MPI_COMM_WORLD,ier)

  end subroutine max_all_dp

!
!----
!

  subroutine max_all_cr(sendbuf, recvbuf)

  use constants
  use mpi

  implicit none

  include "precision.h"

  real(kind=CUSTOM_REAL) sendbuf, recvbuf
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE, &
                  MPI_MAX,0,MPI_COMM_WORLD,ier)

  end subroutine max_all_cr

!
!----
!

  subroutine min_all_cr(sendbuf, recvbuf)

  use constants
  use mpi

  implicit none

  include "precision.h"

  real(kind=CUSTOM_REAL) sendbuf, recvbuf
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE, &
                  MPI_MIN,0,MPI_COMM_WORLD,ier)

  end subroutine min_all_cr


!
!----
!

  subroutine min_all_all_cr(sendbuf, recvbuf)

  use constants
  use mpi

  implicit none

  include "precision.h"

  real(kind=CUSTOM_REAL):: sendbuf, recvbuf
  integer ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE, &
                  MPI_MIN,MPI_COMM_WORLD,ier)

  end subroutine min_all_all_cr

!
!----
!
!
!
!  subroutine min_all_all_dp(sendbuf, recvbuf)
!
!  use constants
!  use mpi
!
!  implicit none
!
!  include "precision.h"
!
!  double precision :: sendbuf, recvbuf
!  integer ier
!
!  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION, &
!                  MPI_MIN,MPI_COMM_WORLD,ier)
!
!  end subroutine min_all_all_dp
!
!
!----
!

  subroutine max_all_i(sendbuf, recvbuf)

  use constants
  use mpi

  implicit none

  include "precision.h"

  integer :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_INTEGER, &
                  MPI_MAX,0,MPI_COMM_WORLD,ier)

  end subroutine max_all_i


!
!----
!

  subroutine max_allreduce_i(buffer,countval)

  use mpi

  implicit none

  integer :: countval
  integer,dimension(countval),intent(inout) :: buffer

  ! local parameters
  integer :: ier
  integer,dimension(countval) :: send

  ! seems not to be supported on all kind of MPI implementations...
  !call MPI_ALLREDUCE(MPI_IN_PLACE, buffer, countval, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ier)

  send(:) = buffer(:)

  call MPI_ALLREDUCE(send, buffer, countval, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ier)
  if( ier /= 0 ) stop 'Allreduce to get max values failed.'

  end subroutine max_allreduce_i

!
!----
!

  subroutine max_all_all_cr(sendbuf, recvbuf)

  use constants
  use mpi

  implicit none

  include "precision.h"

  real(kind=CUSTOM_REAL):: sendbuf, recvbuf
  integer ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE, &
                  MPI_MAX,MPI_COMM_WORLD,ier)

  end subroutine max_all_all_cr


!
!----
!


  subroutine max_all_all_dp(sendbuf, recvbuf)

  use constants
  use mpi

  implicit none

  include "precision.h"

  double precision :: sendbuf, recvbuf
  integer ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION, &
                  MPI_MAX,MPI_COMM_WORLD,ier)

  end subroutine max_all_all_dp


!
!----
!

  subroutine min_all_i(sendbuf, recvbuf)

  use constants
  use mpi

  implicit none

  include "precision.h"

  integer:: sendbuf, recvbuf
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_INTEGER, &
                  MPI_MIN,0,MPI_COMM_WORLD,ier)

  end subroutine min_all_i

!
!----
!

  subroutine maxloc_all_dp(sendbuf, recvbuf)

  use mpi

  implicit none

  double precision, dimension(2) :: sendbuf,recvbuf
  integer ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_2DOUBLE_PRECISION, &
                  MPI_MAXLOC,MPI_COMM_WORLD,ier)

  end subroutine maxloc_all_dp


!
!----
!


  subroutine sum_all_dp(sendbuf, recvbuf)

  use mpi

  implicit none

  double precision sendbuf, recvbuf
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION, &
                  MPI_SUM,0,MPI_COMM_WORLD,ier)

  end subroutine sum_all_dp

!
!----
!

  subroutine sum_all_cr(sendbuf, recvbuf)

  use constants
  use mpi

  implicit none

  include "precision.h"

  real(kind=CUSTOM_REAL) sendbuf, recvbuf
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE, &
                  MPI_SUM,0,MPI_COMM_WORLD,ier)

  end subroutine sum_all_cr

!
!----
!

  subroutine sum_all_all_cr(sendbuf, recvbuf)

  use constants
  use mpi

  implicit none

  include "precision.h"

  real(kind=CUSTOM_REAL) sendbuf, recvbuf
  integer ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE, &
                  MPI_SUM,MPI_COMM_WORLD,ier)

  end subroutine sum_all_all_cr

!
!----
!

  subroutine sum_all_i(sendbuf, recvbuf)

  use mpi

  implicit none

  integer sendbuf, recvbuf
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_INTEGER, &
                  MPI_SUM,0,MPI_COMM_WORLD,ier)

  end subroutine sum_all_i

!
!----
!

  subroutine sum_all_all_i(sendbuf, recvbuf)

  use mpi

  implicit none

  integer sendbuf, recvbuf
  integer ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_INTEGER, &
                  MPI_SUM,MPI_COMM_WORLD,ier)

  end subroutine sum_all_all_i

!
!----
!

  subroutine any_all_l(sendbuf, recvbuf)

  use mpi

  implicit none

  logical sendbuf, recvbuf
  integer ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_LOGICAL, &
                  MPI_LOR,MPI_COMM_WORLD,ier)

  end subroutine any_all_l

!
!----
!

  subroutine sendrecv_all_cr(sendbuf, sendcount, dest, sendtag, &
                             recvbuf, recvcount, source, recvtag)

  use constants
  use mpi

  implicit none

  include "precision.h"

  integer sendcount, recvcount, dest, sendtag, source, recvtag
  real(kind=CUSTOM_REAL), dimension(sendcount) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(recvcount) :: recvbuf

  integer ier

  call MPI_SENDRECV(sendbuf,sendcount,CUSTOM_MPI_TYPE,dest,sendtag, &
                    recvbuf,recvcount,CUSTOM_MPI_TYPE,source,recvtag, &
                    MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)

  end subroutine sendrecv_all_cr

!
!----
!

  integer function proc_null()

  use mpi

  implicit none

  proc_null = MPI_PROC_NULL

  end function proc_null

!
!----
!

  subroutine isend_cr(sendbuf, sendcount, dest, sendtag, req)

  use constants
  use mpi

  implicit none

  include "precision.h"

  integer sendcount, dest, sendtag, req
  real(kind=CUSTOM_REAL), dimension(sendcount) :: sendbuf

  integer ier

  call MPI_ISEND(sendbuf,sendcount,CUSTOM_MPI_TYPE,dest,sendtag, &
                  MPI_COMM_WORLD,req,ier)

  end subroutine isend_cr

!
!----
!

  subroutine irecv_cr(recvbuf, recvcount, dest, recvtag, req)

  use constants
  use mpi

  implicit none

  include "precision.h"

  integer recvcount, dest, recvtag, req
  real(kind=CUSTOM_REAL), dimension(recvcount) :: recvbuf

  integer ier

  call MPI_IRECV(recvbuf,recvcount,CUSTOM_MPI_TYPE,dest,recvtag, &
                  MPI_COMM_WORLD,req,ier)

  end subroutine irecv_cr

!
!----
!

  subroutine isend_i(sendbuf, sendcount, dest, sendtag, req)

  use constants
  use mpi

  implicit none

  include "precision.h"

  integer sendcount, dest, sendtag, req
  integer, dimension(sendcount) :: sendbuf

  integer ier

  call MPI_ISEND(sendbuf,sendcount,MPI_INTEGER,dest,sendtag, &
                  MPI_COMM_WORLD,req,ier)

  end subroutine isend_i

!
!----
!

  subroutine irecv_i(recvbuf, recvcount, dest, recvtag, req)

  use constants
  use mpi

  implicit none

  include "precision.h"

  integer recvcount, dest, recvtag, req
  integer, dimension(recvcount) :: recvbuf
  integer ier

  call MPI_IRECV(recvbuf,recvcount,MPI_INTEGER,dest,recvtag, &
                  MPI_COMM_WORLD,req,ier)

  end subroutine irecv_i


!
!----
!

  subroutine recv_i(recvbuf, recvcount, dest, recvtag )

  use mpi

  implicit none

  integer dest,recvtag
  integer recvcount
  !integer recvbuf
  integer,dimension(recvcount):: recvbuf
  integer ier

  call MPI_RECV(recvbuf,recvcount,MPI_INTEGER,dest,recvtag, &
                MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)

  end subroutine recv_i

!
!----
!

  subroutine recvv_cr(recvbuf, recvcount, dest, recvtag )

  use constants
  use mpi

  implicit none

  include "precision.h"

  integer recvcount,dest,recvtag
  real(kind=CUSTOM_REAL),dimension(recvcount) :: recvbuf
  integer ier

  call MPI_RECV(recvbuf,recvcount,CUSTOM_MPI_TYPE,dest,recvtag, &
                MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)

  end subroutine recvv_cr


!
!----
!

  subroutine send_i(sendbuf, sendcount, dest, sendtag)

  use mpi

  implicit none

  !integer sendbuf,sendcount,dest,sendtag
  integer dest,sendtag
  integer sendcount
  integer,dimension(sendcount):: sendbuf
  integer ier

  call MPI_SEND(sendbuf,sendcount,MPI_INTEGER,dest,sendtag,MPI_COMM_WORLD,ier)

  end subroutine send_i


!
!----
!

  subroutine send_i_t(sendbuf,sendcount,dest)

  use mpi

  implicit none

  integer :: dest,sendcount,ier
  integer :: tag = 100
  integer, dimension(sendcount) :: sendbuf

  call MPI_SEND(sendbuf,sendcount,MPI_INTEGER,dest,tag, &
       MPI_COMM_WORLD,ier)

  end subroutine send_i_t

!
!----
!


  subroutine recv_i_t(recvbuf,recvcount,source)

  use mpi

  implicit none

  integer :: source,recvcount,ier
  integer :: tag = 100
  integer, dimension(recvcount) :: recvbuf

  call MPI_RECV(recvbuf,recvcount,MPI_INTEGER,source,tag, &
                MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)

  end subroutine recv_i_t


!
!----
!
!
!  subroutine send_dp_t(sendbuf,sendcount,dest)
!
!  use mpi
!
!  implicit none
!
!  integer :: dest,sendcount,ier
!  integer :: tag = 100
!  double precision, dimension(sendcount) :: sendbuf
!
!  call MPI_SEND(sendbuf,sendcount,MPI_DOUBLE_PRECISION,dest,tag, &
!       MPI_COMM_WORLD,ier)
!
!  end subroutine send_dp_t
!
!
!----
!
!
!  subroutine recv_dp_t(recvbuf,recvcount,source)
!
!  use mpi
!
!  implicit none
!
!  integer :: recvcount,source,ier
!  integer :: tag = 100
!  double precision, dimension(recvcount) :: recvbuf
!
!  call MPI_RECV(recvbuf,recvcount,MPI_DOUBLE_PRECISION,source,tag, &
!                MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)
!
!  end subroutine recv_dp_t
!
!
!
!----
!

  subroutine send_dp(sendbuf, sendcount, dest, sendtag)

  use mpi

  implicit none

  integer dest,sendtag
  integer sendcount
  double precision,dimension(sendcount):: sendbuf
  integer ier

  call MPI_SEND(sendbuf,sendcount,MPI_DOUBLE_PRECISION,dest,sendtag,MPI_COMM_WORLD,ier)

  end subroutine send_dp

!
!----
!

  subroutine recv_dp(recvbuf, recvcount, dest, recvtag)

  use mpi

  implicit none

  integer dest,recvtag
  integer recvcount
  double precision,dimension(recvcount):: recvbuf
  integer ier

  call MPI_RECV(recvbuf,recvcount,MPI_DOUBLE_PRECISION,dest,recvtag, &
                MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)

  end subroutine recv_dp

!
!----
!

  subroutine sendv_cr(sendbuf, sendcount, dest, sendtag)

  use constants
  use mpi

  implicit none

  include "precision.h"

  integer sendcount,dest,sendtag
  real(kind=CUSTOM_REAL),dimension(sendcount) :: sendbuf
  integer ier

  call MPI_SEND(sendbuf,sendcount,CUSTOM_MPI_TYPE,dest,sendtag,MPI_COMM_WORLD,ier)

  end subroutine sendv_cr

!
!----
!

  subroutine wait_req(req)

  use mpi

  implicit none

  integer :: req

  integer :: ier

  call mpi_wait(req,MPI_STATUS_IGNORE,ier)

  end subroutine wait_req

!
!----
!

  subroutine world_get_comm(comm)

  use mpi

  implicit none

  integer,intent(out) :: comm

  comm = MPI_COMM_WORLD

  end subroutine world_get_comm

!
!----
!

  subroutine world_duplicate(comm)

  use mpi

  implicit none

  integer,intent(out) :: comm
  integer :: ier

  call MPI_COMM_DUP(MPI_COMM_WORLD,comm,ier)
  if( ier /= 0 ) stop 'error duplicating MPI_COMM_WORLD communicator'

  end subroutine world_duplicate
